#ifndef HELPERS_H_
#define HELPERS_H_

#include <string>
#include <vector>
#include <fstream>

typedef struct ArrayVector2D {
	double *x;
	double *y;
} ArrayVector2D;

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
	ArrayVector2D init_pos;
	std::string output_filename;
	double theta; // Used in Barnes-Hut algorithm
} Constants;


void parse_constants(std::string filename, Constants &constants, bool get_init_pos = true);

void write_positions(std::ofstream &output_file, ArrayVector2D pos, Constants constants);

void gen_random_points(Constants &constants);

void log_energy_conservation(ArrayVector2D pos, ArrayVector2D vel, Constants constants);

#endif