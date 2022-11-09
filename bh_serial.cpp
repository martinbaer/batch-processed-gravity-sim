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
	// Check that the number of particles does not exceed the limitations of data types
	// (the BH tree must store pointers to the children of each node)
	if (2 * constants.num_particles > USHRT_MAX)
	{
		std::cerr << "Error: Too many particles for BH tree: " << constants.num_particles << " * 2 > " << USHRT_MAX << std::endl;
		return 1;
	}
	// Initialise the output file
	std::ofstream output_file(constants.output_filename);
	// Check if the file opened
	if (!output_file.is_open())
	{
		std::cerr << "Error opening file " << argv[2] << std::endl;
		return 1;
	}
	// Initialise physical vectors
	std::vector<double> pos_x = std::vector<double>(constants.num_particles);
	std::vector<double> pos_y = std::vector<double>(constants.num_particles);
	std::vector<double> vel_x = std::vector<double>(constants.num_particles);
	std::vector<double> vel_y = std::vector<double>(constants.num_particles);
	std::vector<double> acc_x = std::vector<double>(constants.num_particles);
	std::vector<double> acc_y = std::vector<double>(constants.num_particles);

	// Set the initial positions, velocities and accelerations
	for (int i = 0; i < constants.num_particles; i++)
	{
		pos_x[i] = constants.init_pos_x[i];
		pos_y[i] = constants.init_pos_y[i];
		vel_x[i] = 0;
		vel_y[i] = 0;
		acc_x[i] = 0;
		acc_y[i] = 0;
	}
	constants.init_pos_x.clear();
	constants.init_pos_y.clear();

	// Loop over the number of steps
	for (int step = 0; step < constants.num_steps; step++)
	{
		// Write the positions to the binary output file
		if (step % constants.write_interval == 0)
			write_positions(output_file, pos_x, pos_y, constants);
		// Check the energy conservation
		if (constants.log_energy_conservation) 
			log_energy_conservation(pos_x, pos_y, vel_x, vel_y, constants);

		// Calculate the bounds of the simulation : TODO AVX using _mm256_cmp_pd (and MPI)
		double min_x = pos_x[0];
		double max_x = pos_x[0];
		double min_y = pos_y[0];
		double max_y = pos_y[0];
		for (int i = 1; i < constants.num_particles; i++)
		{
			if (pos_x[i] < min_x)
				min_x = pos_x[i];
			else if (pos_x[i] > max_x)
				max_x = pos_x[i];
			if (pos_y[i] < min_y)
				min_y = pos_y[i];
			else if (pos_y[i] > max_y)
				max_y = pos_y[i];
		}

		// Create tree : TODO MPI by distrubiting particles, creating trees and merging them
		std::vector<Node> bh_tree;
		bh_tree.push_back(Node());

		// Create root node info
		NodeDescriber root;
		root.centre_x = (max_x + min_x) / 2;
		root.centre_y = (max_y + min_y) / 2;
		root.half_width = max_x - min_x < max_y - min_y ? (max_y - min_y) / 2 : (max_x - min_x) / 2;
		// Add each particle to the barnes-hut tree
		for (int i = 0; i < constants.num_particles; i++)
		{
			bh_tree_insert(pos_x[i], pos_y[i], bh_tree, root);
		}
		// Print tree
		// print_tree(0, bh_tree, TREE_ROOT);
		// exit(0);

		if (constants.log_tree_size)
			log_tree_size(bh_tree, constants);

		// Loop over each particle to calculate th_ acceleration: TODO CUDA
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Get the acceleration for the particle
			add_node_acceleration(acc_x[i], acc_y[i], pos_x[i], pos_y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
			acc_x[i] *= constants.gravity;
			acc_y[i] *= constants.gravity;
		}

		// Loop over the particles to update their velocities and positions
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Update the velocity
			vel_x[i] += acc_x[i] * constants.delta_t;
			vel_y[i] += acc_y[i] * constants.delta_t;
			// Update the position
			pos_x[i] += vel_x[i] * constants.delta_t;
			pos_y[i] += vel_y[i] * constants.delta_t;
			// Reset the acceleration
			acc_x[i] = 0;
			acc_y[i] = 0;
		}

		// Destroy tree
		bh_tree.clear();
	}

	// Close output file
	output_file.close();
	// Free memory
	pos_x.clear();
	pos_y.clear();
	vel_x.clear();
	vel_y.clear();
	acc_x.clear();
	acc_y.clear();

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the time taken
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	// Print the time taken
	std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

	return 0;
}