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

#include "bh_tree.h"

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2
#define BLOCK_SIZE 256

void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, unsigned int node_index, double s, BHTree bh_tree, Constants constants)
{
	Node node = bh_tree.nodes[node_index];
	// Calculate the distance between the particle and the node
	double dx = node.centre_of_mass_x - x;
	double dy = node.centre_of_mass_y - y;
	double d = sqrt(dx * dx + dy * dy);
	// If the node is a leaf, add the acceleration
	if (node.mass == 1)
	{
		// Calculate and add the acceleration (mass is 1)
		acc_x += dx / (d * d * d + constants.softening);
		acc_y += dy / (d * d * d + constants.softening);
	}
	// If the node is not a leaf, check if the node is far enough to take its centre of mass
	else
	{
		// Check the s/d ratio for the node
		if (s / d < constants.theta)
		{
			// Calculate and add the acceleration (mass is >1)
			acc_x += node.mass * dx / (d * d * d + constants.softening);
			acc_y += node.mass * dy / (d * d * d + constants.softening);
		}
		else
		{
			// Recursively calculate the acceleration
			double new_s = s / 2;
			if (node.bottom_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_left, new_s, bh_tree, constants);
			if (node.bottom_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_right, new_s, bh_tree, constants);
			if (node.top_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_left, new_s, bh_tree, constants);
			if (node.top_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_right, new_s, bh_tree, constants);
		}
	}
}


__device__ add_node_acceleration_kernel(double &acc_x, double &acc_y, double x, double y, unsigned int node_index, double s, BHTree bh_tree, Constants constants)
{
	Node node = bh_tree.nodes[node_index];
	// Calculate the distance between the particle and the node
	double dx = node.centre_of_mass_x - x;
	double dy = node.centre_of_mass_y - y;
	double d = sqrt(dx * dx + dy * dy);
	// If the node is a leaf, add the acceleration
	if (node.mass == 1)
	{
		// Calculate and add the acceleration (mass is 1)
		acc_x += dx / (d * d * d + constants.softening);
		acc_y += dy / (d * d * d + constants.softening);
	}
	// If the node is not a leaf, check if the node is far enough to take its centre of mass
	else
	{
		// Check the s/d ratio for the node
		if (s / d < constants.theta)
		{
			// Calculate and add the acceleration (mass is >1)
			acc_x += node.mass * dx / (d * d * d + constants.softening);
			acc_y += node.mass * dy / (d * d * d + constants.softening);
		}
		else
		{
			// Recursively calculate the acceleration
			double new_s = s / 2;
			if (node.bottom_left)
				add_node_acceleration_kernel(acc_x, acc_y, x, y, node.bottom_left, new_s, bh_tree, constants);
			if (node.bottom_right)
				add_node_acceleration_kernel(acc_x, acc_y, x, y, node.bottom_right, new_s, bh_tree, constants);
			if (node.top_left)
				add_node_acceleration_kernel(acc_x, acc_y, x, y, node.top_left, new_s, bh_tree, constants);
			if (node.top_right)
				add_node_acceleration_kernel(acc_x, acc_y, x, y, node.top_right, new_s, bh_tree, constants);
		}
	}
}

// Calculate the node accleration and then multiply it by gravity
__device__ calculate_acceleration_kernel(ArrayVector2D pos, ArrayVector2D acc, BHTree bh_tree, Constants constants)
{
	// starting index for the thread's row
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// Zero the acceleration
	acc[i].x = 0;
	acc[i].y = 0;
	// Calculate the acceleration for the particle using iteration instead of recursion
	
	// Multiply by gravity
	acc[i].x *= constants.gravity;
	acc[i].y *= constants.gravity;
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
		std::cerr << "Error opening file " << argv[2] << std::endl;
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


		// Loop over each particle to calculate th_ acceleration: TODO CUDA

		// Copy the tree to the GPU

		// Reallocate the tree on the GPU if it is too small (this will only happen a couple of times)
		if (bh_tree_device.max_nodes < bh_tree.num_nodes)
		{
			checkError(cudaFree(bh_tree_device.nodes));
			checkError(cudaMalloc(&bh_tree_device.nodes, bh_tree.num_nodes * sizeof(Node)));
			bh_tree_device.max_nodes = bh_tree.num_nodes;
		}

		// Copy the tree to the GPU
		checkError(cudaMemcpy(bh_tree_device.nodes, bh_tree.nodes, bh_tree.max_nodes * sizeof(Node), cudaMemcpyHostToDevice));

		// Call the CUDA kernel to calculate the acceleration
		calculate_acceleration_kernel<<<constants.num_particles / BLOCK_SIZE + 1, BLOCK_SIZE>>>(pos_device, acc_device, bh_tree_device, constants);

		// for (int i = 0; i < constants.num_particles; i++)
		// {
		// 	// Get the acceleration for the particle
		// 	calc_node_acceleration(acc.x[i], acc.y[i], pos.x[i], pos.y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
		// 	acc.x[i] *= constants.gravity;
		// 	acc.y[i] *= constants.gravity;
		// }

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