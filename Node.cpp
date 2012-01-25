#include "Node.h"
#include "Branch.h"
#include "globals.h"
#include <iostream>
#include <sstream>

Node::Node(Node *parent, int id, int seq)
{
	_id = id;
	if (seq >= 0)
		_isLeaf = true;
	else
		_isLeaf = false;
	_seq = seq;


	if (parent)
	{
		if (parent->_branches.size() >= 3)
		{
			cerr << "Error: parent of node #" << id << " (node #" << parent->_id << ") already has " << parent->_branches.size() << " neighbours." << endl;
		}
		Branch *branch = new Branch(id, parent, this);
		parent->_branches.push_back(branch);
		_branches.push_back(branch);
	}

	beta = 0.8;
	for (int i = 0; i < CharStates; i++)
		probs.push_back(1.0 / CharStates);
}



Node::~Node()
{
}


vector<Node*> Node::getTraversal(Node *parent)
{
	cerr << getIdent() << endl;
	vector<Node*> list;

//	if (!_isLeaf)
	{
		for (unsigned int i = 0; i < _branches.size(); i++)
		{
			Node *neighbour = _branches[i]->getNeighbour(this);
			if (neighbour != parent)
			{
				vector<Node*> childList = neighbour->getTraversal(this);
				list.insert(list.begin(), childList.begin(), childList.end());
			}
		}
	}

	return list;
}


string Node::toString(Node *parent)
{
	stringstream ss;
	Branch *parentBranch = NULL;

	if (_isLeaf)
	{
		parentBranch = _branches[0];
	} else
	{
		vector<Node*> list;
		for (unsigned int i = 0; i < _branches.size(); i++)
		{
			Node *neighbour = _branches[i]->getNeighbour(this);
			if (neighbour != parent)
				list.push_back(neighbour);
			else
				parentBranch = _branches[i];
		}

		ss << "(";
		for (unsigned int i = 0; i < list.size() - 1; i++)
		{
			ss << list[i]->toString(this);
			ss << ",";
		}
		ss << list[list.size() - 1]->toString(this);
		ss << ")";
	}

	ss << _label;
	if (parentBranch && parentBranch->getDistance() >= .0)
		ss << ":" << parentBranch->getDistance();

	return ss.str();
}



Branch* Node::getBranch(int id)
{
	if ((int) _branches.size() < id+1)
		return NULL;
	else
		return _branches[id];
}



Node* Node::getParent()
{
	if (_branches.empty())
		return NULL;
	else
		return _branches[0]->getNeighbour(this);
}



string Node::getIdent()
{
	stringstream ss;

	ss << _id;
	if (_isLeaf)
		ss << "[L]";
	else
		ss << "[I]";

	if (_label.length())
		ss << "(" << _label << ")";

	return ss.str();
}

