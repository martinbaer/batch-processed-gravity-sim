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

#define TREE_SIZE_LOG_FILENAME "output/tree_size_log.txt"

#define REALLOCATE_FACTOR 1.2


std::ofstream tree_size_log(TREE_SIZE_LOG_FILENAME);
bool tree_size_log_open = false;

// increase the size of the nodes array in the BHTree by a factor of REALLOCATE_FACTOR
// This could be more efficient though it will only occur a few times
inline void reallocate_tree(BHTree &tree)
{
	Node* nodes = new Node[tree.max_nodes];
	for (int i = 0; i < tree.num_nodes; i++)
	{
		nodes[i] = tree.nodes[i];
	}
	delete[] tree.nodes;
	tree.max_nodes *= REALLOCATE_FACTOR;
	tree.nodes = new Node[tree.max_nodes];
	for (int i = 0; i < tree.num_nodes; i++)
	{
		tree.nodes[i] = nodes[i];
	}
}


/**
 * @brief 
 * 
 * @param root 
 * @param acc_x 
 * @param acc_x 
 * @param x 
 * @param y
 */
void add_node_acceleration(double &acc_x, double &acc_y, double x, double y, unsigned int node_index, double s, BHTree bh_tree, Constants constants)
{
	Node node = bh_tree.nodes[node_index];
	// Calculate the distance between the particle and the node
	double dx = node.centre_of_mass_x - x;
	double dy = node.centre_of_mass_y - y;
	double d = sqrt(dx * dx + dy * dy);
	// If the node is a leaf, add the acceleration
	if (node.mass == 1)
	{
		// Calculate and add the acceleration (mass is 1)
		acc_x += dx / (d * d * d + constants.softening);
		acc_y += dy / (d * d * d + constants.softening);
	}
	// If the node is not a leaf, check if the node is far enough to take its centre of mass
	else
	{
		// Check the s/d ratio for the node
		if (s / d < constants.theta)
		{
			// Calculate and add the acceleration (mass is >1)
			acc_x += node.mass * dx / (d * d * d + constants.softening);
			acc_y += node.mass * dy / (d * d * d + constants.softening);
		}
		else
		{
			// Recursively calculate the acceleration
			double new_s = s / 2;
			if (node.bottom_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_left, new_s, bh_tree, constants);
			if (node.bottom_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.bottom_right, new_s, bh_tree, constants);
			if (node.top_left)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_left, new_s, bh_tree, constants);
			if (node.top_right)
				add_node_acceleration(acc_x, acc_y, x, y, node.top_right, new_s, bh_tree, constants);
		}
	}
}




// find the maximum height of the tree using recursion
int tree_height(BHTree bh_tree, unsigned int node_index)
{
	Node node = bh_tree.nodes[node_index];
	if (!(node.top_left || node.top_right || node.bottom_left || node.bottom_right))
	{
		return 1;
	}
	else
	{
		int max_height = 0;
		if (node.top_left)
		{
			int height = tree_height(bh_tree, node.top_left);
			if (height > max_height)
			{
				max_height = height;
			}
		}
		if (node.top_right)
		{
			int height = tree_height(bh_tree, node.top_right);
			if (height > max_height)
			{
				max_height = height;
			}
		}
		if (node.bottom_left)
		{
			int height = tree_height(bh_tree, node.bottom_left);
			if (height > max_height)
			{
				max_height = height;
			}
		}
		if (node.bottom_right)
		{
			int height = tree_height(bh_tree, node.bottom_right);
			if (height > max_height)
			{
				max_height = height;
			}
		}
		return max_height + 1;
	}
}

void log_tree_size(BHTree bh_tree, Constants constants)
{
	if (!tree_size_log_open)
	{
		tree_size_log << "number_of_nodes,height" << std::endl;
		tree_size_log_open = true;
	}
	tree_size_log << bh_tree.num_nodes << "," << tree_height(bh_tree, ROOT_INDEX) << std::endl;
}

void print_tree(int depth, BHTree tree, unsigned int node_index)
{
	if (node_index == 0 && depth != 0)
	{
		return;
	}
	for (int i = 0; i < depth; i++)
	{
		std::cout << "  ";
	}
	Node node = tree.nodes[node_index];
	std::cout << "Node " << node_index << " with mass " << node.mass << ", centre of mass (" << node.centre_of_mass_x << ", " << node.centre_of_mass_y << ")" << std::endl;
	if (node.top_left)
		print_tree(depth + 1, tree, node.top_left);
	if (node.top_right)
		print_tree(depth + 1, tree, node.top_right);
	if (node.bottom_left)
		print_tree(depth + 1, tree, node.bottom_left);
	if (node.bottom_right)
		print_tree(depth + 1, tree, node.bottom_right);
}


void zero_node(BHTree &bh_tree, unsigned int node_index)
{
	bh_tree.nodes[node_index].mass = 0;
	bh_tree.nodes[node_index].centre_of_mass_x = 0;
	bh_tree.nodes[node_index].centre_of_mass_y = 0;
	bh_tree.nodes[node_index].top_left = 0;
	bh_tree.nodes[node_index].top_right = 0;
	bh_tree.nodes[node_index].bottom_left = 0;
	bh_tree.nodes[node_index].bottom_right = 0;
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
NodeDescriber get_child(double x, double y, BHTree &tree, NodeDescriber parent_node_desc)
{
	// Variables useful later
	int parent_index = parent_node_desc.index;
	Node parent_node = tree.nodes[parent_index];
	// Instantiate child node info struct
	NodeDescriber child_node_desc;
	child_node_desc.half_width = parent_node_desc.half_width / 2;
	// Check if the particle is in the top left quadrant
	if (x <= parent_node_desc.centre_x && y >= parent_node_desc.centre_y)
	{
		child_node_desc.centre_x = parent_node_desc.centre_x - child_node_desc.half_width;
		child_node_desc.centre_y = parent_node_desc.centre_y + child_node_desc.half_width;
		// Check if the top left quadrant does not exists
		if (!parent_node.top_left)
		{
			// Create the top left quadrant
			// Get the next available index
			child_node_desc.index = tree.num_nodes++;
			// Reallocate the tree to fit the new node
			if (tree.num_nodes > tree.max_nodes)
				reallocate_tree(tree);
			// Tell the parent node about the child node
			tree.nodes[parent_index].top_left = child_node_desc.index;
			// Initialise the child node
			zero_node(tree, child_node_desc.index);
		}
		else
		{
			// The top left quadrant already exists
			child_node_desc.index = parent_node.top_left;
		}
	}
	// Check if the particle is in the top right quadrant
	else if (x >= parent_node_desc.centre_x && y >= parent_node_desc.centre_y)
	{
		child_node_desc.centre_x = parent_node_desc.centre_x + child_node_desc.half_width;
		child_node_desc.centre_y = parent_node_desc.centre_y + child_node_desc.half_width;
		// Check if the top right quadrant does not exists
		if (!parent_node.top_right)
		{
			// Create the top right quadrant
			// Get the next available index
			child_node_desc.index = tree.num_nodes++;
			// Reallocate the tree to fit the new node
			if (tree.num_nodes > tree.max_nodes)
				reallocate_tree(tree);
			// Tell the parent node about the child node
			tree.nodes[parent_index].top_right = child_node_desc.index;
			// Initialise the child node
			zero_node(tree, child_node_desc.index);
		}
		else
		{
			// The top right quadrant already exists
			child_node_desc.index = parent_node.top_right;
		}
	}
	// Check if the particle is in the bottom left quadrant
	else if (x <= parent_node_desc.centre_x && y <= parent_node_desc.centre_y)
	{
		child_node_desc.centre_x = parent_node_desc.centre_x - child_node_desc.half_width;
		child_node_desc.centre_y = parent_node_desc.centre_y - child_node_desc.half_width;
		// Check if the bottom left quadrant does not exists
		if (!parent_node.bottom_left)
		{
			// Create the bottom left quadrant
			// Get the next available index
			child_node_desc.index = tree.num_nodes++;
			// Reallocate the tree to fit the new node
			if (tree.num_nodes > tree.max_nodes)
				reallocate_tree(tree);
			// Tell the parent node about the child node
			tree.nodes[parent_index].bottom_left = child_node_desc.index;
			// Initialise the child node
			zero_node(tree, child_node_desc.index);
		}
		else
		{
			// The bottom left quadrant already exists
			child_node_desc.index = parent_node.bottom_left;
		}
	}
	// Check if the particle is in the bottom right quadrant
	else
	{
		child_node_desc.centre_x = parent_node_desc.centre_x + child_node_desc.half_width;
		child_node_desc.centre_y = parent_node_desc.centre_y - child_node_desc.half_width;
		// Check if the bottom right quadrant does not exists
		if (!parent_node.bottom_right)
		{
			// Create the bottom right quadrant
			// Get the next available index
			child_node_desc.index = tree.num_nodes++;
			// Reallocate the tree to fit the new node
			if (tree.num_nodes > tree.max_nodes)
				reallocate_tree(tree);
			// Tell the parent node about the child node
			tree.nodes[parent_index].bottom_right = child_node_desc.index;
			// Initialise the child node
			zero_node(tree, child_node_desc.index);
		}
		else
		{
			// The bottom right quadrant already exists
			child_node_desc.index = parent_node.bottom_right;
		}
	}
	// Return the child node
	return child_node_desc;
}

/**
 * @brief 
 * Recursively add a particle to the Barnes-Hutt tree, using the following steps
 * If node x does not contain a body, put the new body here.
 * If node x is an internal node, update the COM and mass of the node
 * If node x is an external node, create four children and recursively add the old and new bodies

 */
void bh_tree_insert(double x, double y, BHTree &tree, NodeDescriber node_desc)
{
	int node_index = node_desc.index;
	Node node = tree.nodes[node_index];
	// Check if the node is a leaf
	if (!(node.top_left || node.top_right || node.bottom_left || node.bottom_right))
	{
		// Node is a leaf
		if (node.mass == 0)
		{
			// Node is an empty leaf: add the particle here
			tree.nodes[node_index].mass = 1;
			tree.nodes[node_index].centre_of_mass_x = x;
			tree.nodes[node_index].centre_of_mass_y = y;
		}
		else
		{
			// Node is an occupied leaf: split the node and add the existing and new particles
			// Add the existing particle to child of the tree
			NodeDescriber child_desc = get_child(node.centre_of_mass_x, node.centre_of_mass_y, tree, node_desc);
			bh_tree_insert(node.centre_of_mass_x, node.centre_of_mass_y, tree, child_desc);
			// Re-try adding new particle to the tree
			bh_tree_insert(x, y, tree, node_desc);
		}
	}
	else
	{
		// Node is not a leaf: update the mass and centre of mass of the node and add the node's child
		tree.nodes[node_index].centre_of_mass_x = (node.centre_of_mass_x * node.mass + x) / (node.mass + 1);
		tree.nodes[node_index].centre_of_mass_y = (node.centre_of_mass_y * node.mass + y) / (node.mass + 1);
		tree.nodes[node_index].mass++;
		// Add the particle to the appropriate child
		NodeDescriber child_desc = get_child(x, y, tree, node_desc);
		bh_tree_insert(x, y, tree, child_desc);
	}
}

// // Merge two trees based on the centre of mass of the nodes
// void merge_trees(std::vector<Node> &dst_tree, std::vector<Node> src_tree, unsigned short dst_node_index, unsigned short src_node_index)
// {
// 	Node dst_node = dst_tree[dst_node_index];
// 	Node src_node = src_tree[src_node_index];
// 	// Average the centre of mass of the two nodes
// 	dst_tree[dst_node_index].centre_of_mass_x = (dst_node.centre_of_mass_x * dst_node.mass + src_node.centre_of_mass_x * src_node.mass) / (dst_node.mass + src_node.mass);
// 	dst_tree[dst_node_index].centre_of_mass_y = (dst_node.centre_of_mass_y * dst_node.mass + src_node.centre_of_mass_y * src_node.mass) / (dst_node.mass + src_node.mass);
// 	// Add the mass of the two nodes
// 	dst_tree[dst_node_index].mass += src_node.mass;
// 	// Merge the src node's children with the dst node's children
// 	if (src_node.top_left)
// 	{
// 		// Create the top left child
// 		if (!dst_node.top_left)
// 		{
// 			// Create the top left child
// 			dst_tree[dst_node_index].top_left = dst_tree.size();
// 			dst_tree.push_back(Node());
// 		}
// 		// Merge the top left child
// 		merge_trees(dst_tree, src_tree, dst_tree[dst_node_index].top_left, src_node.top_left);
// 	}
// 	if (src_node.top_right)
// 	{
// 		// Create the top right child
// 		if (!dst_node.top_right)
// 		{
// 			// Create the top right child
// 			dst_tree[dst_node_index].top_right = dst_tree.size();
// 			dst_tree.push_back(Node());
// 		}
// 		// Merge the top right child
// 		merge_trees(dst_tree, src_tree, dst_tree[dst_node_index].top_right, src_node.top_right);
// 	}
// 	if (src_node.bottom_left)
// 	{
// 		// Create the bottom left child
// 		if (!dst_node.bottom_left)
// 		{
// 			// Create the bottom left child
// 			dst_tree[dst_node_index].bottom_left = dst_tree.size();
// 			dst_tree.push_back(Node());
// 		}
// 		// Merge the bottom left child
// 		merge_trees(dst_tree, src_tree, dst_tree[dst_node_index].bottom_left, src_node.bottom_left);
// 	}
// 	if (src_node.bottom_right)
// 	{
// 		// Create the bottom right child
// 		if (!dst_node.bottom_right)
// 		{
// 			// Create the bottom right child
// 			dst_tree[dst_node_index].bottom_right = dst_tree.size();
// 			dst_tree.push_back(Node());
// 		}
// 		// Merge the bottom right child
// 		merge_trees(dst_tree, src_tree, dst_tree[dst_node_index].bottom_right, src_node.bottom_right);
// 	}
// }