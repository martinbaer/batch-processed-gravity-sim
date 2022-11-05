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

#include "bh_helpers.h"

#define USAGE "Usage: ./bh_serial [constants file]"
#define NUM_ARGS 2



/**
 * @brief 
 * Simulates a particle system under the given constants 
 * and saves the particle positions to a binary file.
 */
int main(int argc, char *argv[])
{
	// Check if the correct number of arguments were given
	if (argc != NUM_ARGS)
	{
		std::cout << USAGE << std::endl;
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
		std::cout << "Error opening file " << argv[2] << std::endl;
		return 1;
	}
	// Initialise physical vectors
	PhysicalVector2D pos;
	pos.x = std::vector<double>(constants.num_particles);
	pos.y = std::vector<double>(constants.num_particles);
	PhysicalVector2D vel;
	vel.x = std::vector<double>(constants.num_particles);
	vel.y = std::vector<double>(constants.num_particles);
	PhysicalVector2D acc;
	acc.x = std::vector<double>(constants.num_particles);
	acc.y = std::vector<double>(constants.num_particles);

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

	

	// Loop over the number of steps
	for (int step = 0; step < constants.num_steps; step++)
	{
		// Write the positions to the binary output file
		if (step % constants.write_interval == 0)
			write_positions(output_file, pos, constants);
		// Check the energy conservation
		if (constants.check_energy_conservation) 
			check_energy_conservation(pos, vel, constants);


		// Calculate the bounds of the simulation
		double min_x = pos.x[0];
		double max_x = pos.x[0];
		double min_y = pos.y[0];
		double max_y = pos.y[0];
		for (int i = 1; i < constants.num_particles; i++)
		{
			if (pos.x[i] < min_x)
				min_x = pos.x[i];
			if (pos.x[i] > max_x)
				max_x = pos.x[i];
			if (pos.y[i] < min_y)
				min_y = pos.y[i];
			if (pos.y[i] > max_y)
				max_y = pos.y[i];
		}
		std::cout << "HERE" << std::endl;
		// Create the root node based on the bounds
		QuadNode *root = new QuadNode;
		root->centre_x = (max_x + min_x) / 2;
		root->centre_y = (max_y + min_y) / 2;
		root->half_width = std::max(max_x - min_x, max_y - min_y) / 2;
		root->mass = 0;
		
		// Add each particle to the barnes-hut tree
		for (int i = 0; i < constants.num_particles; i++)
		{
			bh_tree_insert(root, pos.x[i], pos.y[i]);
		}

		// print here
		std::cout << "root->centre_x: " << root->centre_x << std::endl;
		// Loop over each particle to calculate the acceleration
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Get the acceleration for the particle
			calculate_acceleration(root, &acc.x[i], &acc.y[i], pos.x[i], pos.y[i], constants);
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
	pos.x.clear();
	pos.y.clear();
	vel.x.clear();
	vel.y.clear();
	acc.x.clear();
	acc.y.clear();
	return 0;
}