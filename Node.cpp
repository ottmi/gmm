#include "Node.h"
#include "Branch.h"
#include "globals.h"
#include <cmath>
#include <iostream>
#include <sstream>

Node::Node(Node *parent, int id)
{
	_id = id;
	_isLeaf = false; // this will be set to true in setSequence()

	if (parent)
	{
		if (parent->_branches.size() >= 3)
		{
			cerr << "Error: parent of node #" << _id << " (node #" << parent->_id << ") already has " << parent->_branches.size() << " neighbours." << endl;
		}
		Branch *branch = new Branch(_id - 1, parent, this);
		parent->_branches.push_back(branch);
		_branches.push_back(branch);
	}
}

Node::~Node()
{
}

void Node::setSequence(vector<unsigned int> &seq)
{
	_seq = seq;
	_isLeaf = true;
}

vector<Node*> Node::getTraversal(Node *parent)
{
	cerr << getIdent() << endl;
	vector<Node*> list;

	for (unsigned int i = 0; i < _branches.size(); i++)
	{
		Node *neighbour = _branches[i]->getNeighbour(this);
		if (neighbour != parent)
		{
			vector<Node*> childList = neighbour->getTraversal(this);
			list.insert(list.begin(), childList.begin(), childList.end());
		}
	}

	return list;
}

string Node::toString(Node *parent)
{
	stringstream ss;
	Branch *parentBranch = NULL;

	vector<Node*> list;
	for (unsigned int i = 0; i < _branches.size(); i++)
	{
		Node *neighbour = _branches[i]->getNeighbour(this);
		if (neighbour != parent)
			list.push_back(neighbour);
		else
			parentBranch = _branches[i];
	}

	if (list.size() >= 1)
	{
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
	if (parentBranch) ss << ":" << parentBranch->getLength();

	return ss.str();
}

Branch* Node::getBranch(int id)
{
	if ((int) _branches.size() > id)
		return _branches[id];
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getBranch(" << id << "): there are only " << _branches.size() << " branches";
		throw(ss.str());
	}
}

void Node::removeBranch(Branch* b)
{
	unsigned int i=0;
	while (_branches[i] != b)
			i++;

	if (_branches[i] == b)
			_branches.erase(_branches.begin() + i);
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::removeBranch(" << b->getId() << "): not found";
		throw(ss.str());
	}
}

void Node::addBranch(Branch* b)
{
	if ((_isLeaf && _branches.size() >= 1) || (!_isLeaf && _branches.size() >= 3))
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::addBranch(" << b->getId() << "): there are already " << _branches.size() << " branches";
		throw(ss.str());
	} else
	_branches.push_back(b);
}

void Node::reroot(Branch *branch)
{
	if (branch != NULL) // if this was the root node, branch would be NULL
	{
		if (branch->getNode(1) != this) // make _nodes[0] point towards the root
			branch->swapNodes();

		for (unsigned int i=0; i<_branches.size(); i++)
		{
			if (_branches[i] == branch)
			{
				_branches.erase(_branches.begin() + i); // remove branch from the list
				break;
			}
		}
		_branches.insert(_branches.begin(), branch); // re-insert at the front, so _branches[0] points towards the root
	}

	for (unsigned int i=0; i<_branches.size(); i++)
		if (_branches[i] != branch)
			_branches[i]->getNeighbour(this)->reroot(_branches[i]); // re-root all children recursively
}

Node* Node::getParent()
{
	if (!_branches.empty())
		return _branches[0]->getNeighbour(this);
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getParent(): there are no neighbours";
		throw(ss.str());
	}
}

Node* Node::getChild(int num)
{
	if ((int) _branches.size() > num)
		return _branches[num]->getNeighbour(this);
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getChild(" << num << "): there are only " << _branches.size() << " branches";
		throw(ss.str());
	}
}

string Node::getIdent()
{
	stringstream ss;

	if (_isLeaf)
		ss << "{";
	else
		ss << "[";

	ss << _id;

	if (_label.length()) ss << "|" << _label;

	if (_isLeaf)
		ss << "}";
	else
		ss << "]";

	return ss.str();
}

unsigned int Node::getBase(unsigned int site)
{
	if (site < _seq.size())
		return _seq[site];
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getBase(" << site << "): invalid site";
		throw(ss.str());
	}
}

vector<unsigned int>& Node::getSequence()
{
	if (_isLeaf)
		return _seq;
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getSequence(): is not a leaf";
		throw(ss.str());
	}

}
