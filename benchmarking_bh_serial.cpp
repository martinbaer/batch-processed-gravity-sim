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

#include "bh_tree.h"

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2

#define STARTING_PROBLEM_SIZE 1<<7
#define D_PROBLEM_SIZE 1<<7
#define N_TRIALS 1 // this can be 1 as number of steps is high
#define TIME_LIMIT_MS 600000 //10 mins
#define BENCHMARKING_OUTPUT_FILE "benchmarking_bh_serial.csv"

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
	// start timer
	auto start = std::chrono::high_resolution_clock::now();
	// open benchmarking file
	std::ofstream benchmarking_file(BENCHMARKING_OUTPUT_FILE);
	benchmarking_file << "problem_size,total_timewrite_time,tree_construction_time,acceleration_calc_time,integration_time" << std::endl;

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
		std::vector<int> ingeration_times_trial;

		
		for (int trial = 0; trial < N_TRIALS; trial++)
		{
			auto start_trial = std::chrono::high_resolution_clock::now();
			// Parse the constants file
			Constants constants;
			parse_constants(argv[1], constants, true, problem_size);
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

			// Initialise the BH tree
			BHTree bh_tree;
			bh_tree.max_nodes = constants.num_particles * 2;
			bh_tree.nodes = new Node[bh_tree.max_nodes];

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


			std::vector<int> write_times_inner;
			std::vector<int> tree_construction_times_inner;
			std::vector<int> acceleration_calc_times_inner;
			std::vector<int> ingeration_times_inner;

			// Loop over the number of steps
			for (int step = 0; step < constants.num_steps; step++)
			{
				auto start_inner = std::chrono::high_resolution_clock::now();
				// Write the positions to the binary output file
				if (step % constants.write_interval == 0)
					write_positions(output_file, pos, constants);

				auto finish_write = std::chrono::high_resolution_clock::now();
				write_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_write - start_inner).count());


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

				auto finish_tree_construction = std::chrono::high_resolution_clock::now();
				tree_construction_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_tree_construction - finish_write).count());

				// Loop over each particle to calculate th_ acceleration: TODO CUDA
				for (int i = 0; i < constants.num_particles; i++)
				{
					// Get the acceleration for the particle
					add_node_acceleration(acc.x[i], acc.y[i], pos.x[i], pos.y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
					acc.x[i] *= constants.gravity;
					acc.y[i] *= constants.gravity;
				}

				auto finish_acceleration_calc = std::chrono::high_resolution_clock::now();
				acceleration_calc_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_acceleration_calc - finish_tree_construction).count());

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

				auto finish_integration = std::chrono::high_resolution_clock::now();
				ingeration_times_inner.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_integration - finish_acceleration_calc).count());
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

			auto finish_trial = std::chrono::high_resolution_clock::now();

			total_times_trial.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(finish_trial - start_trial).count());
			write_times_trial.push_back(sum(write_times_inner));
			tree_construction_times_trial.push_back(sum(tree_construction_times_inner));
			acceleration_calc_times_trial.push_back(sum(acceleration_calc_times_inner));
			ingeration_times_trial.push_back(sum(ingeration_times_inner));
		}

		// save the time results to file
		// problem_size,write_time,tree_construction_time,acceleration_calc_time,integration_time
		benchmarking_file << problem_size << "," << avg(total_times_trial) << "," << avg(write_times_trial) << "," << avg(tree_construction_times_trial) << "," << avg(acceleration_calc_times_trial) << "," << avg(ingeration_times_trial) << "," << avg(total_times_trial) << std::endl;
	}
	// close the benchmarking file
	benchmarking_file.close();
	return 0;
}