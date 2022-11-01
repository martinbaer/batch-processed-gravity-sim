#ifndef PARSE_CONSTANTS_H_
#define PARSE_CONSTANTS_H_

#include <string>
#include <vector>
#include <fstream>

typedef struct Point2D {
	double x;
	double y;
} Point2D;

typedef struct Constants {
	int num_particles;
	int box_size;
	int num_steps;
	double delta_t;
	std::vector<Point2D> init_pos;
} Constants;

void parse_constants(std::string filename, Constants &constants);

void gen_random_points(Constants &constants);

#endif