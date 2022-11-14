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

#define ROOT 0

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2

#define STARTING_PROBLEM_SIZE 1<<7
#define D_PROBLEM_SIZE 1<<7
#define N_TRIALS 1
#define TIME_LIMIT_MS 600000 //10 mins
#define BENCHMARKING_OUTPUT_FILE "benchmarking_bh_mpi2.csv"

template <class T> 
T avg(std::vector<T> arr) {
    T sum = 0;
    for (int i = 0; i < arr.size(); i++) {
        sum += arr[i];
    }
    return sum / arr.size();
}

template <class T>
T sum(std::vector<T> arr) {
	T sum = 0;
	for (int i = 0; i < arr.size(); i++) {
		sum += arr[i];
	}
	return sum;
}

/**
 * @brief 
 * Simulates a particle system under the given constants 
 * and saves the particle positions to a binary file.
 */
int main(int argc, char *argv[])
{

	// Initialise MPI
	int world_size, my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	std::ofstream benchmarking_file;

	// start timer
	auto start = std::chrono::high_resolution_clock::now();
	if (my_rank == ROOT)
	{
		// open benchmarking file
		benchmarking_file.open(BENCHMARKING_OUTPUT_FILE);
		benchmarking_file << "problem_size,total_time,write_time,tree_construction_time,acceleration_calc_time,integration_time" << std::endl;
	}

	// Check if the correct number of arguments were given
	if (argc != NUM_ARGS)
	{
		std::cerr << USAGE << std::endl;
		return 1;
	}

	


	/* BENCHMARKING LOOP */
	for (int problem_size = STARTING_PROBLEM_SIZE; std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() < TIME_LIMIT_MS; problem_size += D_PROBLEM_SIZE)
	{
		/* TRIAL LOOP */
		std::vector<int> total_times_trial;
		std::vector<int> write_times_trial;
		std::vector<int> tree_construction_times_trial;
		std::vector<int> acceleration_calc_times_trial;
		std::vector<int> integration_times_trial;

		
		for (int trial = 0; trial < N_TRIALS; trial++)
		{
			auto start_trial = std::chrono::high_resolution_clock::now();

			
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
			parse_constants(argv[1], constants, my_rank == ROOT, problem_size);

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


			// Scatter the particles' positions to all processes
			MPI_Scatter(constants.init_pos.x, num_particles_per_proc, MPI_DOUBLE, pos.x, num_particles_per_proc, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
			MPI_Scatter(constants.init_pos.y, num_particles_per_proc, MPI_DOUBLE, pos.y, num_particles_per_proc, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


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


			std::vector<int> write_times_inner;
			std::vector<int> tree_construction_times_inner;
			std::vector<int> acceleration_calc_times_inner;
			std::vector<int> integration_times_inner;


			unsigned int output_file_offset = 0; // how many doubles are offset
			// Loop over the number of steps
			for (int step = 0; step < constants.num_steps; step++)
			{
				// print position of first particle
				// Write the positions to the binary output file

				auto start_inner = std::chrono::high_resolution_clock::now();

				if (step % constants.write_interval == 0)
				{
					MPI_File_write_at_all(output_file, (output_file_offset + my_rank * num_particles_per_proc) * sizeof(double), pos.x, num_particles_per_proc, MPI_DOUBLE, MPI_STATUS_IGNORE);
					MPI_File_write_at_all(output_file, (output_file_offset + constants.num_particles + my_rank * num_particles_per_proc) * sizeof(double), pos.y, num_particles_per_proc, MPI_DOUBLE, MPI_STATUS_IGNORE);
					output_file_offset += 2 * constants.num_particles; // how many doubles are offset
				}

				auto finish_write = std::chrono::high_resolution_clock::now();
				write_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_write - start_inner).count());


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


				auto finish_tree_construction = std::chrono::high_resolution_clock::now();
				tree_construction_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_tree_construction - finish_write).count());

				std::vector<int> acceleration_calc_times_inner_inner;
				std::vector<int> integration_times_inner_inner;

				// print_tree(0, bh_tree, ROOT_INDEX);

				// Loop over each particle in subset of particles to calculate the acceleration, velocity and position, and then zero the acceleration
				for (int i = 0; i < num_particles_per_proc; i++)
				{
					if (my_rank == last_rank && i > constants.num_particles % world_size && constants.num_particles % world_size != 0)
						break;
					// Get the acceleration for the particle
					auto start_acceleration_calc = std::chrono::high_resolution_clock::now();
					add_node_acceleration(acc.x[i], acc.y[i], pos.x[i], pos.y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
					auto finish_acceleration_calc = std::chrono::high_resolution_clock::now();
					acceleration_calc_times_inner_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_acceleration_calc - start_acceleration_calc).count());
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

					integration_times_inner_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - finish_acceleration_calc).count());
				}

				acceleration_calc_times_inner.push_back(sum(acceleration_calc_times_inner_inner));
				integration_times_inner.push_back(sum(integration_times_inner_inner));
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


			

			auto finish_trial = std::chrono::high_resolution_clock::now();

			total_times_trial.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_trial - start_trial).count());
			write_times_trial.push_back(sum(write_times_inner));
			tree_construction_times_trial.push_back(sum(tree_construction_times_inner));
			acceleration_calc_times_trial.push_back(sum(acceleration_calc_times_inner));
			integration_times_trial.push_back(sum(integration_times_inner));
		}

		// save the time results to file
		// problem_size,write_time,tree_construction_time,acceleration_calc_time,integration_time
		if (my_rank == ROOT)
			benchmarking_file << problem_size << "," << avg(total_times_trial) << "," << avg(write_times_trial) << "," << avg(tree_construction_times_trial) << "," << avg(acceleration_calc_times_trial) << "," << avg(integration_times_trial) << "," << avg(total_times_trial) << std::endl;
	}
	// close the benchmarking file
	if (my_rank == ROOT)
		benchmarking_file.close();
	// Finalise MPI
	MPI_Finalize();
	return 0;
}