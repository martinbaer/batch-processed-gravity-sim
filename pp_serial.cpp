/**
 * @file pp_serial.cpp
 * @author Martin Baer
 * @brief Particle-particle simulator serial implimentation
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

#include "parse_constants.h"

int main(int argc, char *argv[])
{
	Constants constants;
	// Check if the correct number of arguments were given
	if (argc != 2)
	{
		std::cout << "Usage: ./pp_serial <constants file>" << std::endl;
		return 1;
	}
	// Parse the constants file
	parse_constants(argv[1], constants);
	// print the constants
	std::cout << "num_particles: " << constants.num_particles << std::endl;
	std::cout << "box_size: " << constants.box_size << std::endl;
	std::cout << "num_steps: " << constants.num_steps << std::endl;
	std::cout << "delta_t: " << constants.delta_t << std::endl;
	// std::cout << "init_pos: " << constants.init_pos << std::endl;

	return 0;
}