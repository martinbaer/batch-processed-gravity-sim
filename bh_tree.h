#include "helpers.h"

#ifndef BH_HELPERS_H_
#define BH_HELPERS_H_

typedef struct QuadNode {
	double centre_x;
	double centre_y;
	double half_width;
	double centre_of_mass_x;
	double centre_of_mass_y;
	int mass;
	int top_left;
	int top_right;
	int bottom_left;
	int bottom_right;
	bool is_leaf;
} QuadNode;
// TODO: create an echo QuadNode strcuture that contains only info the bh builder needs

void bh_tree_insert(double x, double y, std::vector<QuadNode> &tree, int node_index);

void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, int node_index, std::vector<QuadNode> bh_tree, Constants constants);

#endif