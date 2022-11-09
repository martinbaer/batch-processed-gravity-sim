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

#define USAGE "Usage: ./pp_serial [constants file]"
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
				double dx = pos_x[j] - pos_x[i];
				double dy = pos_y[j] - pos_y[i];
				double r = sqrt(dx * dx + dy * dy);
				// Calculate and cumulatively sum the acceleration
				acc_x[i] += constants.gravity * dx / (r * r * r + constants.softening);
				acc_y[i] += constants.gravity * dy / (r * r * r + constants.softening);
			}
		}
		// Loop over the particles
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

	// stop timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

	return 0;
}