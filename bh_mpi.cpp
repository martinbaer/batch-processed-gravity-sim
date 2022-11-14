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
#include <chrono>

#include <mpi.h>

#include "bh_tree.h"
#include "helpers.h"

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2

#define ROOT 0


/**
 * MPI plan 1:
 * each process is exclusively responsible for a certain number of particles' positions, velocities and accelerations
 * the bh tree is shared between all processes
 * 
 * Initilise:
 * i. root reads constants and initialises the particle positions
 * ii. constants are broadcast to all processes
 * iii. particles' positions are scattered to all processes (each has n/p particles)
 * iv. collective output file is opened
 * 
 * Main loop:
 * 1. each process contructs a sub-tree with a subset of the particles
 * 2. repetatively until root is the only process left:
 * 		2.1 split the processes in half
 * 		2.2 the second half sends its sub-tree to the first half
 * 		2.3 the first half merges the sub-trees
 * 3. the root process now has the full tree, broadcast it to all processes
 * 4. each process uses the tree to calculate the acceleration on its subset of particles
 * 5. each process updates its subset of particles' velocities and positions
 * 6. each process writes its subset of particles' positions to the output file
 * 
 * ----------------
 * 
 * MPI plan 2:
 * each process is exclusively responsible for a sub-tree of the bh tree
 * the particles' positions, velocities and accelerations are shared between all processes
 * 
 * Initilise:
 * i. root reads constants and initialises the particle positions
 * ii. constants are broadcast to all processes
 * iii. particles' positions are broadcast to all processes
 * iv. collective output file is opened
 * 
 * Main loop:
 * forces
 * 1. each process contructs a sub-tree with a subset of the particles
 * 2. each process uses its partial tree to calculate the acceleration on all particles
 * 3. repetatively until root is the only process left:
 * 		3.1 split the processes in half
 * 		3.2 the second half sends its particle accelerations to the first half
 * 		3.3 the first half sums the accelerations
 * 4. the root node now has the total acceleration on all particles, scatter it to all processes
 * 5. each process updates its subset of the particles' velocities and positions
 * 6. each process writes its particles' positions to the output file
 * 7. the processes allgather their particles' positions (they now have all particles' positions for the next iteration)
 * 
 */
int main(int argc, char *argv[])
{
	// start timer
	auto start = std::chrono::high_resolution_clock::now();

	// Initialise MPI
	int world_size, my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int last_rank = world_size - 1;

	// Create MPI datatype for Node
	MPI_Datatype MPI_BH_NODE;
	MPI_Datatype type[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
	int blocklen[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Aint offsets[7];
	offsets[0] = offsetof(Node, centre_of_mass_x);
	offsets[1] = offsetof(Node, centre_of_mass_y);
	offsets[2] = offsetof(Node, mass);
	offsets[3] = offsetof(Node, top_left);
	offsets[4] = offsetof(Node, top_right);
	offsets[5] = offsetof(Node, bottom_left);
	offsets[6] = offsetof(Node, bottom_right);
	MPI_Type_create_struct(7, blocklen, offsets, type, &MPI_BH_NODE);
	MPI_Type_commit(&MPI_BH_NODE);


	// Check if the correct number of arguments were given
	if (argc != NUM_ARGS)
	{
		std::cerr << USAGE << std::endl;
		return 1;
	}

	// All processes parse the constants file
	Constants constants;
	parse_constants(argv[1], constants, my_rank == ROOT, 0);

	// Open the output file
	MPI_File output_file;
	MPI_File_open(MPI_COMM_WORLD, constants.output_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
	// Check if the file was opened successfully
	if (output_file == MPI_FILE_NULL)
	{
		std::cerr <<  "Error opening output file" << std::endl;
		return 1;
	}

	bool is_partial_tail = my_rank == world_size - 1 && constants.num_particles % world_size != 0;

	// Initialise physical vectors
	int num_particles_per_proc = constants.num_particles / world_size + (constants.num_particles % world_size != 0); // round up
	ArrayVector2D pos;
	ArrayVector2D vel;
	ArrayVector2D acc;
	pos.x = new double[num_particles_per_proc];
	pos.y = new double[num_particles_per_proc];
	vel.x = new double[num_particles_per_proc];
	vel.y = new double[num_particles_per_proc];
	acc.x = new double[num_particles_per_proc];
	acc.y = new double[num_particles_per_proc];

	// std::cout << "Process " << my_rank << " has " << num_particles_per_proc << " of " << constants.num_particles << " particles" << std::endl;
	// if (my_rank == ROOT)
	// {
	// 	// print the initial positions
	// 	std::cout << "Initial positions:" << std::endl;
	// 	for (int i = 0; i < constants.num_particles; i++)
	// 	{
	// 		std::cout << "(" << constants.init_pos.x[i] << ", " << constants.init_pos.y[i] << ")" << std::endl;
	// 	}
	// }

	// Scatter the particles' positions to all processes
	MPI_Scatter(constants.init_pos.x, num_particles_per_proc, MPI_DOUBLE, pos.x, num_particles_per_proc, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Scatter(constants.init_pos.y, num_particles_per_proc, MPI_DOUBLE, pos.y, num_particles_per_proc, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	std::cout << "Process " << my_rank << " has finished scattering" << std::endl;

	// Set the initial velocities and accelerations
	for (int i = 0; i < num_particles_per_proc; i++)
	{
		vel.x[i] = 0;
		vel.y[i] = 0;
		acc.x[i] = 0;
		acc.y[i] = 0;
	}
	// Delete the initial positions in the constants struct
	if (my_rank == ROOT)
	{
		delete[] constants.init_pos.x;
		delete[] constants.init_pos.y;
	}

	// Initialise the BH tree for this process - this data structure is used for both the sub-tree (during construction) and the full tree (during force calculation)
	BHTree bh_tree;
	bh_tree.max_nodes = constants.num_particles * 2;
	bh_tree.nodes = new Node[bh_tree.max_nodes];

	unsigned int output_file_offset = 0; // how many doubles are offset
	// Loop over the number of steps
	for (int step = 0; step < constants.num_steps; step++)
	{
		// print position of first particle
		std::cout << "rank: " << my_rank << ", step " << step << ", positions: (" << pos.x[0] << ", " << pos.y[0] << ")" << std::endl;
		// Write the positions to the binary output file
		if (step % constants.write_interval == 0)
		{
			MPI_File_write_at_all(output_file, (output_file_offset + my_rank * num_particles_per_proc) * sizeof(double), pos.x, num_particles_per_proc, MPI_DOUBLE, MPI_STATUS_IGNORE);
			MPI_File_write_at_all(output_file, (output_file_offset + constants.num_particles + my_rank * num_particles_per_proc) * sizeof(double), pos.y, num_particles_per_proc, MPI_DOUBLE, MPI_STATUS_IGNORE);
			output_file_offset += 2 * constants.num_particles; // how many doubles are offset
		}
		// Check the energy conservation
		if (constants.log_energy_conservation && my_rank == ROOT) 
			log_energy_conservation(pos, vel, constants);

		// Calculate the bounds of the simulation : TODO AVX using _mm256_cmp_pd (and MPI)
		double min_x = pos.x[0];
		double max_x = pos.x[0];
		double min_y = pos.y[0];
		double max_y = pos.y[0];
		for (int i = 1; i < num_particles_per_proc; i++)
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
		double lower_bounds[2] = {min_x, min_y};
		double upper_bounds[2] = {max_x, max_y};
		MPI_Allreduce(MPI_IN_PLACE, lower_bounds, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, upper_bounds, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		min_x = lower_bounds[0];
		min_y = lower_bounds[1];
		max_x = upper_bounds[0];
		max_y = upper_bounds[1];

		// Create tree : TODO MPI creating trees and merging them

		// Create root node info
		NodeDescriber root;
		root.index = ROOT_INDEX;
		root.centre_x = (max_x + min_x) / 2;
		root.centre_y = (max_y + min_y) / 2;
		root.half_width = max_x - min_x < max_y - min_y ? (max_y - min_y) / 2 : (max_x - min_x) / 2;
		// Reset the tree
		zero_node(bh_tree, ROOT_INDEX);
		bh_tree.num_nodes = 1;
		// Add each particle to the subset barnes-hut tree
		for (int i = 0; i < num_particles_per_proc; i++)
		{
			if (my_rank == last_rank && i > constants.num_particles % world_size && constants.num_particles % world_size != 0)
			 	break;
			bh_tree_insert(pos.x[i], pos.y[i], bh_tree, root);
		}
		// Print tree
		// print_tree(0, bh_tree, TREE_ROOT);
		// exit(0);

		if (constants.log_tree_size && my_rank == ROOT)
			log_tree_size(bh_tree, constants);
			
		// Merge the trees
		int active_procs = world_size;
		while (active_procs > 1)
		{
			// Halve the number of active procs (rounding up)
			active_procs = (active_procs + 1) / 2;
			// If this process is in the first half of the active processes
			if (my_rank < active_procs)
			{
				// Receive the tree from the other process
				int other_rank = my_rank + active_procs;
				// Receive the number of nodes in the other tree
				unsigned int num_nodes;
				MPI_Recv(&num_nodes, 1, MPI_UNSIGNED, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Receive the nodes in the other tree
				Node* other_nodes = new Node[num_nodes];
				MPI_Recv(other_nodes, num_nodes, MPI_BH_NODE, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Merge the trees
				merge_trees(bh_tree, other_nodes, ROOT_INDEX, ROOT_INDEX);
				// Delete the other nodes
				delete[] other_nodes;
			}
			// If this process is in the second half of the active processes
			else
			{
				// Send the tree to the other process
				int other_rank = my_rank - active_procs;
				// Send the number of nodes
				MPI_Send(&bh_tree.num_nodes, 1, MPI_UNSIGNED, other_rank, 0, MPI_COMM_WORLD);
				// Send the nodes
				MPI_Send(bh_tree.nodes, bh_tree.num_nodes, MPI_BH_NODE, other_rank, 0, MPI_COMM_WORLD);
				// Exit the loop
				break;
			}
		}

		// Now the root process has the full tree, broadcast it to all processes
		// Broadcast the number of nodes
		MPI_Bcast(&bh_tree.num_nodes, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
		// Make sure the nodes array is big enough
		while (bh_tree.num_nodes > bh_tree.max_nodes)
			reallocate_tree(bh_tree);
		// Broadcast the nodes
		MPI_Bcast(bh_tree.nodes, bh_tree.num_nodes, MPI_BH_NODE, ROOT, MPI_COMM_WORLD);


		// print_tree(0, bh_tree, ROOT_INDEX);

		// Loop over each particle in subset of particles to calculate the acceleration, velocity and position, and then zero the acceleration
		for (int i = 0; i < num_particles_per_proc; i++)
		{
			if (my_rank == last_rank && i > constants.num_particles % world_size && constants.num_particles % world_size != 0)
			 	break;
			// Get the acceleration for the particle
			add_node_acceleration(acc.x[i], acc.y[i], pos.x[i], pos.y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
			acc.x[i] *= constants.gravity;
			acc.y[i] *= constants.gravity;
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

	// Close MPI output file
	MPI_File_close(&output_file);
	
	// Free memory
	delete[] pos.x;
	delete[] pos.y;
	delete[] vel.x;
	delete[] vel.y;
	delete[] acc.x;
	delete[] acc.y;
	delete[] bh_tree.nodes;

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the time taken
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	// Print the time taken
	std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

	// Finalise MPI
	MPI_Finalize();

	return 0;
}