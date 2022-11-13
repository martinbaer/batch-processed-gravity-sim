#include "helpers.h"

#ifndef BH_HELPERS_H_
#define BH_HELPERS_H_

#define ROOT_INDEX 0

typedef struct Node
{
	double centre_of_mass_x;
	double centre_of_mass_y;
	unsigned int mass;
	unsigned int top_left;
	unsigned int top_right;
	unsigned int bottom_left;
	unsigned int bottom_right;
} Node;

// Additional information on each node used in the construction of the BH tree but not the the BH tree itself
typedef struct NodeDescriber {
	unsigned int index;
	double centre_x;
	double centre_y;
	double half_width;
} NodeDescriber;

typedef struct BHTree {
	Node *nodes;
	unsigned int num_nodes;
	unsigned int max_nodes;
} BHTree;

void bh_tree_insert(double x, double y, BHTree &tree, NodeDescriber node_desc);

void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, unsigned int node_index, double s, BHTree bh_tree, Constants constants);

void print_tree(int depth, BHTree tree, unsigned int node_index);

void log_tree_size(BHTree bh_tree, Constants constants);

void reallocate_tree(BHTree &tree);

void zero_node(BHTree &bh_tree, unsigned int node_index);

void merge_trees(BHTree &dst_tree, Node* src_tree_nodes, unsigned short dst_node_index, unsigned short src_node_index);

#endif