/**
 * @file parse_constants_file.cpp
 * @author Martin Baer
 * @brief Contains globally useful functions for the simulations
 * @version 0.1
 * @date 2022-11-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "helpers.h"
#include "bh_tree.h"

#define ENERGY_LOG_FILENAME "energy_log.txt"

// Theta value for the Barnes-Hut algorithm
#define THETA 0.5

/**
 * @brief 
 * 
 * @param root 
 * @param acc_x 
 * @param acc_x 
 * @param x 
 * @param y 
 */
void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, int node_index, std::vector<QuadNode> bh_tree, Constants constants)
{
	QuadNode node = bh_tree[node_index];
	// Calculate the distance between the particles
	double dx = node.centre_x - x;
	double dy = node.centre_y - y;
	double r = sqrt(dx * dx + dy * dy);
	// If the node is a leaf, add the acceleration
	if (node.is_leaf)
	{
		// Calculate and add the acceleration (mass is 1)
		acc_x += constants.gravity * dx / (r * r * r + constants.softening);
		acc_y += constants.gravity * dy / (r * r * r + constants.softening);
	}
	// If the node is not a leaf, check if the node is far enough to take its centre of mass
	else
	{
		// Calculate the s/d ratio for the node
		// s is the width of the node, root->half_width
		// d is the distance between the particle and the centre of the node, r
		if (node.half_width / r < THETA)
		{
			// Calculate and cumulatively sum the acceleration
			acc_x += constants.gravity * node.mass * dx / (r * r * r + constants.softening);
			acc_y += constants.gravity * node.mass * dy / (r * r * r + constants.softening);
		}
		else
		{
			// Recursively calculate the acceleration
			if (node.bottom_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_left, bh_tree, constants);
			if (node.bottom_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_right, bh_tree, constants);
			if (node.top_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_left, bh_tree, constants);
			if (node.top_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_right, bh_tree, constants);
		}
	}
}

/**
 * @brief 
 * Return the index quadrant that the particle is in given the parent node.
 * Also creates the child node and adds it to the tree if it doesn't exist.
 * 
 * @param x The x position of the particle
 * @param y The y position of the particle
 * @param node_index The index of the parent node
 * @param tree The tree to add the child node to
 * @return int The index of the child node
 */
int get_child(double x, double y, std::vector<QuadNode> &tree, int node_index)
{

	QuadNode parent_node = tree[node_index];
	double parent_centre_x = parent_node.centre_x;
	double parent_centre_y = parent_node.centre_y;
	// Check if the particle is in the top left quadrant
	if (x <= parent_centre_x && y >= parent_centre_y)
	{
		int child_index = parent_node.top_left;
		// Check if the top left quadrant does not exists
		if (!parent_node.top_left)
		{
			// Create the top left quadrant
			// Set the top left quadrant of the parent node to the index of the top left quadrant
			child_index = tree.size();
			tree[node_index].top_left = child_index;
			tree.push_back(QuadNode());
			// Set values
			double half_width = parent_node.half_width / 2;
			tree[child_index].centre_x = parent_centre_x - half_width;
			tree[child_index].centre_y = parent_centre_y + half_width;
			tree[child_index].half_width = half_width;
			tree[child_index].is_leaf = true;
		}
		// Return the top left quadrant
		return child_index;
	}
	// Check if the particle is in the top right quadrant
	else if (x >= parent_centre_x && y >= parent_centre_y)
	{
		int child_index = parent_node.top_right;
		// Check if the top right quadrant does not exists
		if (!parent_node.top_right)
		{
			// Create the top right quadrant
			// Set the top right quadrant of the parent node to the index of the top right quadrant
			child_index = tree.size();
			tree[node_index].top_right = child_index;
			tree.push_back(QuadNode());
			// Set values
			double half_width = parent_node.half_width / 2;
			tree[child_index].centre_x = parent_centre_x + half_width;
			tree[child_index].centre_y = parent_centre_y + half_width;
			tree[child_index].half_width = half_width;
			tree[child_index].is_leaf = true;
		}
		// Return the top right quadrant
		return child_index;
	}
	// Check if the particle is in the bottom left quadrant
	else if (x <= parent_centre_x && y <= parent_centre_y)
	{
		int child_index = parent_node.bottom_left;
		// Check if the bottom left quadrant does not exists
		if (!parent_node.bottom_left)
		{
			// Create the bottom left quadrant
			// Set the bottom left quadrant of the parent node to the index of the bottom left quadrant
			child_index = tree.size();
			tree[node_index].bottom_left = child_index;
			tree.push_back(QuadNode());
			// Set values
			double half_width = parent_node.half_width / 2;
			tree[child_index].centre_x = parent_centre_x - half_width;
			tree[child_index].centre_y = parent_centre_y - half_width;
			tree[child_index].half_width = half_width;
			tree[child_index].is_leaf = true;
		}
		// Return the bottom left quadrant
		return child_index;
	}
	// Check if the particle is in the bottom right quadrant
	else
	{
		int child_index = parent_node.bottom_right;
		// Check if the bottom right quadrant does not exists
		if (!parent_node.bottom_right)
		{
			// Create the bottom right quadrant
			// Set the bottom right quadrant of the parent node to the index of the bottom right quadrant
			child_index = tree.size();
			tree[node_index].bottom_right = child_index;
			tree.push_back(QuadNode());
			// Set values
			double half_width = parent_node.half_width / 2;
			tree[child_index].centre_x = parent_centre_x + half_width;
			tree[child_index].centre_y = parent_centre_y - half_width;
			tree[child_index].half_width = half_width;
			tree[child_index].is_leaf = true;
		}
		// Return the bottom right quadrant
		return child_index;
	}
}

/**
 * @brief 
 * Recursively add a particle to the Barnes-Hutt tree, using the following steps
 * If node x does not contain a body, put the new body here.
 * If node x is an internal node, update the COM and mass of the node
 * If node x is an external node, create four children and recursively add the old and new bodies

 */
void bh_tree_insert(double x, double y, std::vector<QuadNode> &tree, int node_index)
{
// std::cout << "HERE" << std::endl;
	if (tree[node_index].is_leaf)
	{
		// Node is a leaf
		if (tree[node_index].mass == 0)
		{
			// Node is an empty leaf: add the particle here
			tree[node_index].mass = 1;
			tree[node_index].centre_of_mass_x = x;
			tree[node_index].centre_of_mass_y = y;
		}
		else
		{
			// Node is an occupied leaf: split the node and add the existing and new particles
			tree[node_index].is_leaf = false;
			// Add the old body to child of the tree
			int child_index = get_child(tree[node_index].centre_of_mass_x, tree[node_index].centre_of_mass_y, tree, node_index);
			bh_tree_insert(tree[node_index].centre_of_mass_x, tree[node_index].centre_of_mass_y, tree, child_index);
			// Re-add the new body to the tree
			bh_tree_insert(x, y, tree, node_index);
		}
	}
	else
	{
		// Node is not a leaf: update the COM and mass of the node and recursively add the new body to the appropriate child
		tree[node_index].centre_of_mass_x = (tree[node_index].centre_of_mass_x * tree[node_index].mass + x) / (tree[node_index].mass + 1);
		tree[node_index].centre_of_mass_y = (tree[node_index].centre_of_mass_y * tree[node_index].mass + y) / (tree[node_index].mass + 1);
		tree[node_index].mass++;
		int child_index = get_child(x, y, tree, node_index);
		bh_tree_insert(x, y, tree, child_index);
	}
}