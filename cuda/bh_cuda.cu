/**
 * @file bh_serial.cpp
 * @author Martin Baer
 * @brief 
 * Particle-particle simulator serial implimentation
 * This code simulates a given particle system for a given number of steps, 
 * saving the particle positions of each step to a file.
 * It also prints the time it took to run the simulation.
 * 
 * Barnes-Hut simulator serial implimentation
 * @version 0.1
 * @date 2022-11-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <limits>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "bh_tree.h"

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2
#define BLOCK_SIZE 256

#define QUEUE_EXPANSION_FACTOR 1.5

void checkError(cudaError_t e)
{
   if (e != cudaSuccess)
   {
      std::cerr << "CUDA error: " << int(e) << " : " << cudaGetErrorString(e) << '\n';
      abort();
   }
}


typedef struct QueuedNode
{
   unsigned int tree_index;
   double s;
} QueuedNode;

typedef struct CircularArrayQueue
{
   QueuedNode *data;
   int num_elements;
   int max_elements;
   int head;
   int tail;
} CircularArrayQueue;

__device__ void enqueue(CircularArrayQueue &queue, unsigned int tree_index, double s)
{
	int new_tail = (queue.tail + 1) % queue.max_elements;
	// Expand the queue if needed (this should only occur a few times)
	if (new_tail == queue.head)
	{
		// Copy the queue data to a temporary array
		QueuedNode *temp_data = new QueuedNode[queue.max_elements];
		for (int i = 0; i < queue.max_elements; i++)
		{
			temp_data[i] = queue.data[i];
		}
		// Delete the queue data
		delete[] queue.data;
		// Reallocate the queue
		queue.max_elements = queue.max_elements * QUEUE_EXPANSION_FACTOR;
		queue.data = new QueuedNode[queue.max_elements];
		// Copy the data back from the temporary array
		for (int i = 0; i < queue.max_elements; i++)
		{
			queue.data[i] = temp_data[i];
		}
		// Delete the temporary array
		delete[] temp_data;
	}
	// Add the value to the queue
	queue.data[new_tail].tree_index = tree_index;
	queue.data[new_tail].s = s;
	queue.tail = new_tail;
	queue.num_elements++;
}

__device__ QueuedNode dequeue(CircularArrayQueue &queue)
{
	// Get the value at the head of the queue
	QueuedNode head = queue.data[queue.head];
	// Move the head of the queue
	queue.head = (queue.head - 1) % queue.max_elements;
	queue.num_elements--;
	return head;
}

// Calculate the node accleration and then multiply it by gravity 
__global__ void calculate_acceleration_kernel(ArrayVector2D pos, ArrayVector2D acc, BHTree bh_tree, double root_half_width, Constants constants)
{
	int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
	if (thread_id < constants.num_particles)
		return;
	// Create queue
	CircularArrayQueue queue;
	int expected_size = (int)__log2f(bh_tree.num_nodes);
	queue.data = new QueuedNode[expected_size];
	queue.num_elements = 0;
	queue.max_elements = expected_size;
	queue.head = 0;
	queue.tail = 0;
	// Add the root node to the queue
	enqueue(queue, ROOT_INDEX, root_half_width);
	// Add node acceleration to the acceleration of the particle
	while (queue.num_elements > 0)
	{
		// Get the next node from the queue
		QueuedNode node = dequeue(queue);
		// Calculate the distance between the particle and the node
		double dx = bh_tree.nodes[node.tree_index].centre_of_mass_x - pos.x[thread_id];
		double dy = bh_tree.nodes[node.tree_index].centre_of_mass_y - pos.y[thread_id];
		double d = sqrt(dx * dx + dy * dy);
		// If the node is a leaf, add the acceleration
		if (bh_tree.nodes[node.tree_index].mass == 1)
		{
			// Calculate and add the acceleration (mass is 1)
			acc.x[thread_id] += constants.gravity * dx / (d * d * d + constants.softening);
			acc.y[thread_id] += constants.gravity * dy / (d * d * d + constants.softening);
		}
		else // Check if the node is far enough to take its centre of mass
		{
			// Check the s/d ratio for the node
			if (node.s / d < constants.theta)
			{
				// Calculate and add the acceleration
				acc.x[thread_id] += constants.gravity * bh_tree.nodes[node.tree_index].mass * dx / (d * d * d + constants.softening);
				acc.y[thread_id] += constants.gravity * bh_tree.nodes[node.tree_index].mass * dy / (d * d * d + constants.softening);
			}
			else // Add the children to the queue
			{
				// Add the children to the queue
				double child_s = node.s / 2;
				if (bh_tree.nodes[node.tree_index].bottom_left)
					enqueue(queue, bh_tree.nodes[node.tree_index].bottom_left, child_s);
				if (bh_tree.nodes[node.tree_index].bottom_right)
					enqueue(queue, bh_tree.nodes[node.tree_index].bottom_right, child_s);
				if (bh_tree.nodes[node.tree_index].top_left)
					enqueue(queue, bh_tree.nodes[node.tree_index].top_left, child_s);
				if (bh_tree.nodes[node.tree_index].top_right)
					enqueue(queue, bh_tree.nodes[node.tree_index].top_right, child_s);
			}
		}
	}
	// Multiply the acceleration by gravity
	acc.x[thread_id] *= constants.gravity;
	acc.y[thread_id] *= constants.gravity;
}


/**
 * @brief 
 * Simulates a particle system under the given constants 
 * and saves the particle positions to a binary file.
 */
int main(int argc, char *argv[])
{
	// start timer
	auto start = std::chrono::high_resolution_clock::now();

	// Check if the correct number of arguments were given
	if (argc != NUM_ARGS)
	{
		std::cerr << USAGE << std::endl;
		return 1;
	}
	// Parse the constants file
	Constants constants;
	parse_constants(argv[1], constants);
	// Initialise the output file
	std::ofstream output_file(constants.output_filename);
	// Check if the file opened
	if (!output_file.is_open())
	{
		std::cerr << "Error opening file: " << argv[2] << std::endl;
		return 1;
	}

	// Initialise physical vectors
	ArrayVector2D pos;
	ArrayVector2D vel;
	ArrayVector2D acc;
	pos.x = new double[constants.num_particles];
	pos.y = new double[constants.num_particles];
	vel.x = new double[constants.num_particles];
	vel.y = new double[constants.num_particles];
	acc.x = new double[constants.num_particles];
	acc.y = new double[constants.num_particles];
	// Initialise physical vectors on the GPU
	ArrayVector2D pos_device;
	ArrayVector2D vel_device;
	ArrayVector2D acc_device;

	checkError(cudaMalloc(&pos_device.x, constants.num_particles * sizeof(double)));
	checkError(cudaMalloc(&pos_device.y, constants.num_particles * sizeof(double)));
	checkError(cudaMalloc(&vel_device.x, constants.num_particles * sizeof(double)));
	checkError(cudaMalloc(&vel_device.y, constants.num_particles * sizeof(double)));
	checkError(cudaMalloc(&acc_device.x, constants.num_particles * sizeof(double)));
	checkError(cudaMalloc(&acc_device.y, constants.num_particles * sizeof(double)));

	// Initialise the BH tree
	BHTree bh_tree;
	bh_tree.max_nodes = constants.num_particles * 2;
	bh_tree.nodes = new Node[bh_tree.max_nodes];
	// Initialise the BH tree on the GPU
	BHTree bh_tree_device;
	checkError(cudaMalloc(&bh_tree_device.nodes, bh_tree.max_nodes * sizeof(Node)));

	// Set the initial positions, velocities and accelerations
	for (int i = 0; i < constants.num_particles; i++)
	{
		pos.x[i] = constants.init_pos.x[i];
		pos.y[i] = constants.init_pos.y[i];
		vel.x[i] = 0;
		vel.y[i] = 0;
		acc.x[i] = 0;
		acc.y[i] = 0;
	}
	// Delete the initial positions in the constants struct
	delete[] constants.init_pos.x;
	delete[] constants.init_pos.y;

	// Loop over the number of steps
	for (int step = 0; step < constants.num_steps; step++)
	{
		// Write the positions to the binary output file
		if (step % constants.write_interval == 0)
			write_positions(output_file, pos, constants);
		// Check the energy conservation
		if (constants.log_energy_conservation) 
			log_energy_conservation(pos, vel, constants);

		// Calculate the bounds of the simulation : TODO AVX using _mm256_cmp_pd (and MPI)
		double min_x = pos.x[0];
		double max_x = pos.x[0];
		double min_y = pos.y[0];
		double max_y = pos.y[0];
		for (int i = 1; i < constants.num_particles; i++)
		{
			if (pos.x[i] < min_x)
				min_x = pos.x[i];
			else if (pos.x[i] > max_x)
				max_x = pos.x[i];
			if (pos.y[i] < min_y)
				min_y = pos.y[i];
			else if (pos.y[i] > max_y)
				max_y = pos.y[i];
		}

		// Create tree : TODO MPI by distrubiting particles, creating trees and merging them

		// Create root node info
		NodeDescriber root;
		root.index = ROOT_INDEX;
		root.centre_x = (max_x + min_x) / 2;
		root.centre_y = (max_y + min_y) / 2;
		root.half_width = max_x - min_x < max_y - min_y ? (max_y - min_y) / 2 : (max_x - min_x) / 2;
		// Reset the tree
		zero_node(bh_tree, ROOT_INDEX);
		bh_tree.num_nodes = 1;
		// Add each particle to the barnes-hut tree
		for (int i = 0; i < constants.num_particles; i++)
		{
			bh_tree_insert(pos.x[i], pos.y[i], bh_tree, root);
		}
		// Print tree
		// print_tree(0, bh_tree, TREE_ROOT);
		// exit(0);

		if (constants.log_tree_size)
			log_tree_size(bh_tree, constants);



		/* Phase 2: Calculate accleration */

		// Reallocate the tree on the GPU if it is too small (this will only happen a couple of times)
		if (bh_tree_device.max_nodes < bh_tree.num_nodes)
		{
			checkError(cudaFree(bh_tree_device.nodes));
			checkError(cudaMalloc(&bh_tree_device.nodes, bh_tree.num_nodes * sizeof(Node)));
			bh_tree_device.max_nodes = bh_tree.num_nodes;
		}
		// Copy the tree to the GPU
		checkError(cudaMemcpy(bh_tree_device.nodes, bh_tree.nodes, bh_tree.max_nodes * sizeof(Node), cudaMemcpyHostToDevice));

		// Copy the positions to the GPU
		checkError(cudaMemcpy(pos_device.x, pos.x, constants.num_particles * sizeof(double), cudaMemcpyHostToDevice));
		checkError(cudaMemcpy(pos_device.y, pos.y, constants.num_particles * sizeof(double), cudaMemcpyHostToDevice));

		// Call the CUDA kernel to calculate the acceleration pos_device.x, pos_device.y, acc_device.x, acc_device.y, bh_tree_device.nodes, root.half_width, 
		calculate_acceleration_kernel<<<constants.num_particles / BLOCK_SIZE + 1, BLOCK_SIZE>>>(pos_device, acc_device, bh_tree_device, root.half_width, constants);


		// Loop over the particles to update their velocities and positions
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Update the velocity
			vel.x[i] += acc.x[i] * constants.delta_t;
			vel.y[i] += acc.y[i] * constants.delta_t;
			// Update the position
			pos.x[i] += vel.x[i] * constants.delta_t;
			pos.y[i] += vel.y[i] * constants.delta_t;
			// Reset the acceleration
			acc.x[i] = 0;
			acc.y[i] = 0;
		}
	}

	// Close output file
	output_file.close();
	// Free memory
	delete[] pos.x;
	delete[] pos.y;
	delete[] vel.x;
	delete[] vel.y;
	delete[] acc.x;
	delete[] acc.y;
	delete[] bh_tree.nodes;
	// Free memory on the GPU
	checkError(cudaFree(pos_device.x));
	checkError(cudaFree(pos_device.y));
	checkError(cudaFree(vel_device.x));
	checkError(cudaFree(vel_device.y));
	checkError(cudaFree(acc_device.x));
	checkError(cudaFree(acc_device.y));
	checkError(cudaFree(bh_tree_device.nodes));

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the time taken
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	// Print the time taken
	std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

	return 0;
}