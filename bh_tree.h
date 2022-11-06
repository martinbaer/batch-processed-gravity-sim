#include "helpers.h"

#ifndef BH_HELPERS_H_
#define BH_HELPERS_H_

typedef struct QuadNode {
	double centre_of_mass_x;
	double centre_of_mass_y;
	unsigned short mass;
	unsigned short top_left;
	unsigned short top_right;
	unsigned short bottom_left;
	unsigned short bottom_right;
} QuadNode;

// Additional information on each node used in the construction of the BH tree but not the the BH tree itself
typedef struct QuadNodeDesc {
	unsigned short index;
	double centre_x;
	double centre_y;
	double half_width;
} QuadNodeDesc;

void bh_tree_insert(double x, double y, std::vector<QuadNode> &tree, QuadNodeDesc node_desc);

void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, int node_index, double s, std::vector<QuadNode> bh_tree, Constants constants);

#endif