#include "Node.h"

#include <iostream>

Node::Node(Node *parent, int id, bool isLeaf)
{
	if (parent)
	{
		if (parent->neighbours.size() >= 3)
		{
			cerr << "Error: parent of node #" << id << " (node #" << parent->id << ") already has " << parent->neighbours.size() << " neighbours." << endl;
		}
		parent->neighbours.push_back(this);
		neighbours.push_back(parent);
	}

	this->id = id;
	this->isLeaf = isLeaf;
}

string Node::toString(Node *parent)
{
	if (isLeaf)
		return label;

	vector <Node*> list;
	for (unsigned int i=0; i<neighbours.size(); i++)
		if (neighbours[i] != parent)
			list.push_back(neighbours[i]);

	string result = "(";
	for (unsigned int i=0; i<list.size()-1; i++)
	{
		result+= list[i]->toString(this);
		result+= ",";
	}
	result+= list[list.size()-1]->toString(this);
	result+= ")";
	result+= label;

	return result;
}

Node::~Node()
{
}
