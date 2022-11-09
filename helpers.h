#ifndef HELPERS_H_
#define HELPERS_H_

#include <string>
#include <vector>
#include <fstream>


typedef struct Constants {
	int num_particles;
	int box_size;
	int num_steps;
	int write_interval;
	double delta_t;
	double softening;
	double gravity; // Gravitational constant, multiplied by the mass of each particle squared
	bool log_energy_conservation;
	bool log_tree_size;
	std::vector<double> init_pos_x;
	std::vector<double> init_pos_y;
	std::string output_filename;
	double theta; // Used in Barnes-Hut algorithm
} Constants;

void parse_constants(std::string filename, Constants &constants);

void write_positions(std::ofstream &output_file, std::vector<double> pos_x, std::vector<double> pos_y, Constants constants);

void gen_random_points(Constants &constants);

void log_energy_conservation(std::vector<double> pos_x, std::vector<double> pos_y, std::vector<double> vel_x, std::vector<double> vel_y, Constants constants);

#endif