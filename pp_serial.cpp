/**
 * @file pp_serial.cpp
 * @author Martin Baer
 * @brief 
 * Particle-particle simulator serial implimentation
 * This code simulates a given particle system for a given number of steps, 
 * saving the particle positions of each step to a file.
 * It also prints the time it took to run the simulation.
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

#include "helpers.h"

#define USAGE "Usage: ./pp_serial [constants file] [output file]"
#define NUM_ARGS 3



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
	std::ofstream output_file(argv[2]);
	// Check if the file opened
	if (!output_file.is_open())
	{
		std::cout << "Error opening file " << argv[2] << std::endl;
		return 1;
	}
	// Initialise vectors
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

	// Set vector sizes
	size_t vector_size = pos.x.size();

	// Write the initial positions to the binary output file
	write_positions(output_file, pos, constants);
	// Loop over the number of steps
	for (int step = 0; step < constants.num_steps; step++)
	{
		// Loop over the particles
		for (int i = 0; i < constants.num_particles; i++)
		{
			// Loop over the other particles
			for (int j = 0; j < constants.num_particles; j++)
			{
				// Check if the particles are the same
				if (i == j)
				{
					continue;
				}
				// Calculate the distance between the particles
				double dx = pos.x[j] - pos.x[i];
				double dy = pos.y[j] - pos.y[i];
				double r = sqrt(dx * dx + dy * dy);
				// Calculate the acceleration
				double ax = constants.gravity * constants.particle_mass * dx / (r * r * r + constants.softening);
				double ay = constants.gravity * constants.particle_mass * dy / (r * r * r + constants.softening);
				// Add the acceleration to the total acceleration
				acc.x[i] += ax;
				acc.y[i] += ay;
			}
		}
		// Loop over the particles
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
		// Write the positions to the binary output file
		if (step % constants.write_interval == 0)
		{
			write_positions(output_file, pos, constants);
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