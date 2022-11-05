#include "helpers.h"

#ifndef BH_HELPERS_H_
#define BH_HELPERS_H_

struct QuadNode {
	double centre_x;
	double centre_y;
	double half_width;
	double centre_of_mass_x;
	double centre_of_mass_y;
	int mass;
	struct QuadNode* top_left;
	struct QuadNode* top_right;
	struct QuadNode* bottom_left;
	struct QuadNode* bottom_right;
	bool is_leaf = true;
};
// TODO: create an echo QuadNode strcuture that contains only info the bh builder needs

typedef struct Point {
	double x;
	double y;
} Point;

void bh_tree_insert(struct QuadNode *node, double x, double y);

void calculate_acceleration(struct QuadNode *root, double *acc_x, double *acc_y, double x, double y, Constants constants);

#endif