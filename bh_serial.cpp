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

		// Loop over each particle to calculate th_ acceleration: TODO CUDA
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Get the acceleration for the particle
			add_node_acceleration(acc.x[i], acc.y[i], pos.x[i], pos.y[i], ROOT_INDEX, root.half_width, bh_tree, constants);
			acc.x[i] *= constants.gravity;
			acc.y[i] *= constants.gravity;
		}

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

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the time taken
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	// Print the time taken
	std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

	return 0;
}