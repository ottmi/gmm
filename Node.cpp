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
		Branch *branch = new Branch(_id, parent, this);
		parent->_branches.push_back(branch);
		_branches.push_back(branch);
	}

	_beta = 0.8;
	for (int i = 0; i < charStates; i++)
		probs.push_back(1.0 / charStates);
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
	if (parentBranch && parentBranch->getDistance() >= .0) ss << ":" << parentBranch->getDistance();

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

Node* Node::getParent()
{
	if (!_branches.empty())
		return _branches[0]->getNeighbour(this);
	else
	{
		stringstream ss;
		ss << "Node(" << getIdent() << ")::getParent(): there are neighbours";
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

double Node::computeValuesIntToInt(unsigned int numOfSites)
{
	cout << "computeValuesIntToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	Branch *parentBranch = _branches[0];

	vector<double> pG1 = parentBranch->pSiX2(numOfSites);
	vector<double> pG2 = parentBranch->pRiX1(numOfSites);
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			pG1[site * 4 + parentBase] /= parentBranch->getMarginalProbCol(parentBase);
			for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
			{
				siteProb += parentBranch->getProb(nodeBase, parentBase) * pG1[site * 4 + parentBase] * pG2[site * 4 + nodeBase];
			}
		}

		logLikelihood += log((1 - _beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}

double Node::computeValuesIntToLeaf(unsigned int numOfSites)
{
	cout << "computeValuesIntToLeaf() " << getIdent() << endl;
	double logLikelihood = 0;
	Branch *parentBranch = _branches[0];

	vector<double> pG1 = parentBranch->pSiX2(numOfSites);
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = parentBranch->getMarginalProbCol(parentBase);
			siteProb += parentBranch->getProb(parentBase, _seq[site]) * pG1[site * 4 + parentBase] / marginalProb;
		}

		logLikelihood += log((1 - _beta) * siteProb);
	}

	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}

double Node::computeValuesRootToInt(unsigned int numOfSites)
{
	cout << "computeValuesRootToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	Node *root = getParent();
	Branch *rootBranch = _branches[0];

	vector<double> pG2 = rootBranch->pRiX1(numOfSites);
	vector<unsigned int> rootSeq = root->getSequence();
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int j = 0; j < 4; j++)
		{
			siteProb += rootBranch->getProb(rootSeq[site], j) * pG2[site * 4 + j];
		}

		logLikelihood += log((1 - _beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}
