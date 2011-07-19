#include "Node.h"

#include <iostream>

Node::Node(Node *parent, int id, bool isLeaf)
{
	left = right = NULL;
	this->parent = parent;
	if (parent)
	{
		if (!parent->left)
			parent->left = this;
		else if (!parent->right)
			parent->right = this;
		else
		{
			if (!parent->parent) // this is the "root" in an unrooted tree
			{
				parent->parent = this;
			} else
			{
				cerr << "Error: parent node #" << parent->id << " already has two childs: #" << parent->left->id << "," << parent->right->id << endl;
			}
		}
	}

	this->id = id;
	this->isLeaf = isLeaf;
}

string Node::toString()
{
	if (isLeaf)
		return label;

	string result = "(";

	result+= left->toString();
	result+= ",";
	result+= right->toString();
	result+= ")";
	result+= label;

	return result;
}

Node::~Node()
{
}
