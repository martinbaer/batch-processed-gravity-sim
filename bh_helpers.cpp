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
#include "bh_helpers.h"

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
void calculate_acceleration(struct QuadNode *root, double *acc_x, double *acc_y, double x, double y, Constants constants)
{
	// Check if the node exists
	if (root == NULL)
	{
		return;
	}
	// If the node is a leaf, add the acceleration
	if (root->is_leaf)
	{
		// Calculate the distance between the particles
		double dx = x - root->centre_of_mass_x;
		double dy = y - root->centre_of_mass_y;
		double r = sqrt(dx * dx + dy * dy);
		// Calculate the acceleration
		double ax = constants.gravity * dx / (r * r * r + constants.softening);
		double ay = constants.gravity * dy / (r * r * r + constants.softening);
		// Add the acceleration to the total acceleration
		*acc_x += ax;
		*acc_y += ay;
	}
	// If the node is not a leaf, check if the node is far enough to take its centre of mass
	else
	{
		// Calculate the s/d ratio for the node
		// Calculate the distance between the particles
		double dx = x - root->centre_x;
		double dy = y - root->centre_y;
		double r = sqrt(dx * dx + dy * dy);
		// s is the width of the node, root->half_width
		// d is the distance between the particle and the centre of the node, r
		if (root->half_width / r < THETA)
		{
			// Calculate and cumulatively sum the acceleration
			*acc_x += constants.gravity * dx / (r * r * r + constants.softening);
			*acc_y += constants.gravity * dy / (r * r * r + constants.softening);
		}
		else
		{
			// Recursively calculate the acceleration
			calculate_acceleration(root->top_left, acc_x, acc_y, x, y, constants);
			calculate_acceleration(root->top_right, acc_x, acc_y, x, y, constants);
			calculate_acceleration(root->bottom_left, acc_x, acc_y, x, y, constants);
			calculate_acceleration(root->bottom_right, acc_x, acc_y, x, y, constants);
		}
	}
}


/**
 * @brief 
 * Return the quadrant that the particle is in given the parent node.
 * Also creates the child node if it doesn't exist.
 * 
 * @param parant_node The parant node to add the particle to
 * @param x The x position of the particle
 * @param y The y position of the particle
 */
struct QuadNode *quadrant_of_particle(struct QuadNode *parent_node, double x, double y)
{
	double parent_centre_x = parent_node->centre_x;
	double parent_centre_y = parent_node->centre_y;
	// Check if the particle is in the top left quadrant
	if (x <= parent_centre_x && y >= parent_centre_y)
	{
		// Check if the top left quadrant does not exists
		if (parent_node->top_left == NULL)
		{
			// Create the top left quadrant
			QuadNode* new_node = new struct QuadNode();
			// Set the half width of the quadrant
			new_node->half_width = parent_node->half_width / 2;
			// Set the centre of the quadrant
			new_node->centre_x = parent_centre_x - new_node->half_width / 2;
			new_node->centre_y = parent_centre_y + new_node->half_width / 2;
			// Initialise the mass to 0
			new_node->mass = 0;
			// Initialise all the child nodes to NULL
			new_node->top_left = NULL;
			new_node->top_right = NULL;
			new_node->bottom_left = NULL;
			new_node->bottom_right = NULL;
			// Set the quadrant to be a leaf
			new_node->is_leaf = true;
			// Bind to the parent node
			parent_node->top_left = new_node;
		}
		// Return the top left quadrant
		return parent_node->top_left;
	}
	// Check if the particle is in the top right quadrant
	else if (x >= parent_centre_x && y >= parent_centre_y)
	{
		// Check if the top right quadrant does not exists
		if (parent_node->top_right == NULL)
		{
			// Create the top right quadrant
			QuadNode* new_node = new struct QuadNode();
			// Set the half width of the quadrant
			new_node->half_width = parent_node->half_width / 2;
			// Set the centre of the quadrant
			new_node->centre_x = parent_centre_x + new_node->half_width / 2;
			new_node->centre_y = parent_centre_y + new_node->half_width / 2;
			// Initialise the mass to 0
			new_node->mass = 0;
			// Initialise all the child nodes to NULL
			new_node->top_left = NULL;
			new_node->top_right = NULL;
			new_node->bottom_left = NULL;
			new_node->bottom_right = NULL;
			// Set the quadrant to be a leaf
			new_node->is_leaf = true;
			// Bind to the parent node
			parent_node->top_right = new_node;
		}
		// Return the top right quadrant
		return parent_node->top_right;
	}
	// Check if the particle is in the bottom left quadrant
	else if (x <= parent_centre_x && y <= parent_centre_y)
	{
		// Check if the bottom left quadrant does not exists
		if (parent_node->bottom_left == NULL)
		{
			// Create the bottom left quadrant
			QuadNode* new_node = new struct QuadNode();
			// Set the half width of the quadrant
			new_node->half_width = parent_node->half_width / 2;
			// Set the centre of the quadrant
			new_node->centre_x = parent_centre_x - new_node->half_width / 2;
			new_node->centre_y = parent_centre_y - new_node->half_width / 2;
			// Initialise the mass to 0
			new_node->mass = 0;
			// Initialise all the child nodes to NULL
			new_node->top_left = NULL;
			new_node->top_right = NULL;
			new_node->bottom_left = NULL;
			new_node->bottom_right = NULL;
			// Set the quadrant to be a leaf
			new_node->is_leaf = true;
			// Bind to the parent node
			parent_node->bottom_left = new_node;
		}
		// Return the bottom left quadrant
		return parent_node->bottom_left;
	}
	// Check if the particle is in the bottom right quadrant
	else if (x >= parent_centre_x && y <= parent_centre_y)
	{
		// Check if the bottom right quadrant does not exists
		if (parent_node->bottom_right == NULL)
		{
			// Create the bottom right quadrant
			QuadNode* new_node = new struct QuadNode();
			// Set the half width of the quadrant
			new_node->half_width = parent_node->half_width / 2;
			// Set the centre of the quadrant
			new_node->centre_x = parent_centre_x + new_node->half_width / 2;
			new_node->centre_y = parent_centre_y - new_node->half_width / 2;
			// Initialise the mass to 0
			new_node->mass = 0;
			// Initialise all the child nodes to NULL
			new_node->top_left = NULL;
			new_node->top_right = NULL;
			new_node->bottom_left = NULL;
			new_node->bottom_right = NULL;
			// Set the quadrant to be a leaf
			new_node->is_leaf = true;
			// Bind to the parent node
			parent_node->bottom_right = new_node;
		}
		// Return the bottom right quadrant
		return parent_node->bottom_right;
	}
}

/**
 * @brief 
 * Recursively add a particle to the Barnes-Hutt tree, using the following steps
 * If node x does not contain a body, put the new body here.
 * If node x is an internal node, update the COM and mass of the node
 * If node x is an external node, create four children and recursively add the old and new bodies
 * @param node The node to add the particle to
 * @param x The x position of the particle
 * @param y The y position of the particle
 */
void bh_tree_insert(QuadNode* node, double x, double y)
{
	// Check if node is empty
	if (!node->mass) 
	{
		// Set the node's centre of mass to the particle's position
		node->centre_of_mass_x = x;
		node->centre_of_mass_y = y;
		node->mass = 1;
	}
	else // Node has mass
	{
		// Check if node is a leaf
		if (node->is_leaf)
		{
			// Re-insert node's centre of mass as a particle into the tree
			node->is_leaf = false;
			QuadNode* quadrant = quadrant_of_particle(node, node->centre_of_mass_x, node->centre_of_mass_y);
			bh_tree_insert(quadrant, node->centre_of_mass_x, node->centre_of_mass_y);
		}
		// Add the new particle to the tree
		QuadNode* quadrant2 = quadrant_of_particle(node, x, y);
		bh_tree_insert(quadrant2, x, y);
		// Update the mass and centre of mass of the node
		node->centre_of_mass_x = (node->centre_of_mass_x * node->mass + x) / (node->mass + 1);
		node->centre_of_mass_y = (node->centre_of_mass_y * node->mass + y) / (node->mass + 1);
		node->mass++;
	}
}