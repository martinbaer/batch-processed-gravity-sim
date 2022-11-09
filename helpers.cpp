/**
 * @file parse_constants_file.cpp
 * @author Martin Baer
 * @brief Contains globally useful functions for the simulations
 * @version 0.1
 * @date 2022-11-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "helpers.h"

#include "bh_tree.h"

// TODO change everything to take input/* for input files and save everything to output (do this in the makefile?? and python script)

#define ENERGY_LOG_FILENAME "output/energy_log.txt"

#define RANDOM_SIZE 100

std::ofstream energy_log(ENERGY_LOG_FILENAME);
bool energy_log_open = false;



/**
 * @brief Checks for the conservation of energy in the system and logs the result
 * 
 * @param positions vector of particle positions
 * @param velocities 
 * @param constants 
 */
void log_energy_conservation(std::vector<double> pos_x, std::vector<double> pos_y, std::vector<double> vel_x, std::vector<double> vel_y, Constants constants)
{
	if (!energy_log_open)
	{
		energy_log << "total_energy, kinetic_energy, potential_energy" << std::endl;
		energy_log_open = true;
	}
	// Calculate the total kinetic energy
	// Formula: 1/2 * m * v^2
	double kinetic_energy = 0;
	for (int i = 0; i < constants.num_particles; i++)
	{
		kinetic_energy += 0.5 * (vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i]);
	}
	// Calculate the total potential energy
	// Formula: -G * m1 * m2 / r
	double potential_energy = 0;
	for (int i = 0; i < constants.num_particles; i++)
	{
		for (int j = i + 1; j < constants.num_particles; j++)
		{
			double dx = pos_x[i] - pos_x[j];
			double dy = pos_y[i] - pos_y[j];
			double r = sqrt(dx * dx + dy * dy);
			potential_energy -= constants.gravity / (r + constants.softening);
		}
	}
	// Calculate the total energy
	// Formula: E = K + U
	double total_energy = kinetic_energy + potential_energy;
	// Print the energy values to energy_log.txt
	energy_log << total_energy << "," << kinetic_energy << "," << potential_energy << std::endl;
}

/**
 * @brief 
 * Write the particle positions to the binary output file.
 * @param output_file given binary output file
 * @param positions vector of particle positions
 */
void write_positions(std::ofstream &output_file, std::vector<double> pos_x, std::vector<double> pos_y, Constants constants)
{
	output_file.write((char *)&(pos_x)[0], constants.num_particles * sizeof(double));
	output_file.write((char *)&(pos_y)[0], constants.num_particles * sizeof(double));
}

/**
 * @brief 
 * Parses the constants file and stores the constants in the constants struct.
 * Generates random initial positions for the particles if none are given.
 * 
 * Constants file format:
 * num_particles=[int]
 * num_steps=[int]
 * write_interval=[int]
 * delta_t=[double]
 * softening=[double]
 * gravity=[double]
 * particle_mass=[double]
 * init_pos=[string] (text file name containing initial positions seperated by commas and newlines)
 * @param filename name of file to parse
 * @param constants struct to store constants in
 */
void parse_constants(std::string filename, Constants &constants)
{
	// Open the file
	std::ifstream file(filename);
	// Check if the file opened
	if (!file.is_open())
	{
		std::cout << "Error opening file " << filename << std::endl;
		exit(1);
	}
	// Read the file
	std::string line;
	while (std::getline(file, line))
	{
		// Check if the line is a comment or blank
		if (line[0] == '#' || line[0] == '\n')
		{
			continue;
		}
		// Check if the line is a constant
		if (line.find("num_particles") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to an int
			constants.num_particles = std::stoi(value);
		}
		if (line.find("theta") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to an int
			constants.theta = std::stod(value);
		}
		else if (line.find("num_steps") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to an int
			constants.num_steps = std::stoi(value);
		}
		else if (line.find("delta_t") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to a double
			constants.delta_t = std::stod(value);
		}
		else if (line.find("softening") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to a double
			constants.softening = std::stod(value);
		}
		else if (line.find("gravity") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to a double
			constants.gravity = std::stod(value);
		}
		else if (line.find("write_interval") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to a double
			constants.write_interval = std::stoi(value);
		}
		else if (line.find("log_energy_conservation") != std::string::npos)
		{
			// Check if the value is true or false and set the constant
			constants.log_energy_conservation = line.find("true") != std::string::npos;
		}
		else if (line.find("log_tree_size") != std::string::npos)
		{
			// Check if the value is true or false and set the constant
			constants.log_tree_size = line.find("true") != std::string::npos;
		}
		else if (line.find("output_filename") != std::string::npos)
		{
			// Get the value of the constant
			constants.output_filename = line.substr(line.find('=') + 1);
		}
		else if (line.find("init_pos") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Check if the value is random
			if (value == "random")
			{
				// Generate random initial positions
				gen_random_points(constants);
			}
			else
			{
				// The value is the name of a file containing the initial positions
				// Open the file
				std::ifstream init_pos_file(value);
				// Check if the file opened
				if (!init_pos_file.is_open())
				{
					std::cout << "Error opening file " << value << std::endl;
					exit(1);
				}
				// Initialize the vectors
				constants.init_pos_x = std::vector<double>(constants.num_particles);
				constants.init_pos_y = std::vector<double>(constants.num_particles);
				// Read the file
				std::string init_pos_line;
				while (std::getline(init_pos_file, init_pos_line))
				{
					// Get the x and y values
					std::string x = init_pos_line.substr(0, init_pos_line.find(','));
					std::string y = init_pos_line.substr(init_pos_line.find(',') + 1);
					// Convert the values to doubles
					double x_d = std::stod(x);
					double y_d = std::stod(y);
					// Add the values to the vectors
					constants.init_pos_x.push_back(x_d);
					constants.init_pos_y.push_back(y_d);
				}
			}
		}
	}
}

void gen_random_points(Constants &constants)
{
	// Initialise vectors
	constants.init_pos_x = std::vector<double>(constants.num_particles);
	constants.init_pos_y = std::vector<double>(constants.num_particles);
	// Initialise random number generator
	std::uniform_real_distribution<double> unif(0, RANDOM_SIZE);
	std::default_random_engine re;
	// Generate random points
	for (int i = 0; i < constants.num_particles; i++)
	{
		constants.init_pos_x[i] = unif(re);
		constants.init_pos_y[i] = unif(re);
	}
}