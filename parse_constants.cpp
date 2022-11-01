/**
 * @file parse_constants_file.cpp
 * @author Martin Baer
 * @brief 
 * To parse the a file containing constants for the simulation.
 * 
 * Constants file format:
 * num_particles=[int]
 * box_size=[int]
 * num_steps=[int]
 * delta_t=[double]
 * init_pos=[string] (text file name containing initial positions seperated by commas and newlines)
 * @version 0.1
 * @date 2022-11-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <string>

#include "parse_constants.h"

/**
 * @brief 
 * Parses the constants file and stores the constants in the constants struct.
 * Generates random initial positions for the particles if none are given.
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
		else if (line.find("box_size") != std::string::npos)
		{
			// Get the value of the constant
			std::string value = line.substr(line.find('=') + 1);
			// Convert the value to an int
			constants.box_size = std::stoi(value);
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
				// Initialize the vector
				constants.init_pos = std::vector<Point2D>(constants.num_particles);
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
					// Create a point
					Point2D point;
					point.x = x_d;
					point.y = y_d;
					// Add the point to the vector
					constants.init_pos.push_back(point);
				}
			}
		}
	}
}

void gen_random_points(Constants &constants)
{
	// Initialise vector
	constants.init_pos = std::vector<Point2D>(constants.num_particles);
	// Generate random points
	for (int i = 0; i < constants.num_particles; i++)
	{
		constants.init_pos[i].x = rand() % constants.box_size;
		constants.init_pos[i].y = rand() % constants.box_size;
	}
}