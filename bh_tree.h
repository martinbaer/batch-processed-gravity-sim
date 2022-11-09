#include "helpers.h"

#ifndef BH_HELPERS_H_
#define BH_HELPERS_H_

#define ROOT_INDEX 0

typedef struct Node {
	double centre_of_mass_x;
	double centre_of_mass_y;
	unsigned short mass;
	unsigned short top_left;
	unsigned short top_right;
	unsigned short bottom_left;
	unsigned short bottom_right;
} Node;

// Additional information on each node used in the construction of the BH tree but not the the BH tree itself
typedef struct NodeDescriber {
	unsigned short index;
	double centre_x;
	double centre_y;
	double half_width;
} NodeDescriber;

void bh_tree_insert(double x, double y, std::vector<Node> &tree, NodeDescriber node_desc);

void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, int node_index, double s, std::vector<Node> bh_tree, Constants constants);

void print_tree(int depth, std::vector<Node> tree, int node_index);

void log_tree_size(std::vector<Node> bh_tree, Constants constants);

#endif