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
	if ((int) _branches.size() < id + 1)
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

Node* Node::getChild(int num)
{
	if (_branches.empty())
		return NULL;
	else
		return _branches[num]->getNeighbour(this);
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

double Node::pRiX1(int base, int site)
{
	Node *child1 = _branches[1]->getNeighbour(this);
	Node *child2 = _branches[2]->getNeighbour(this);

	double prob1;
	if (child1->isLeaf())
	{
		unsigned int base1 = child1->getBase(site);
		prob1 = _branches[1]->pX1X2(base, base1);
	} else
	{
		prob1 = 0;
		for (unsigned int base1 = 0; base1 < 4; base1++)
		{
			double pRiX1 = child1->pRiX1(base1, site); // TODO: maybe we want to store this
			prob1 += _branches[1]->pX1X2(base, base1) * pRiX1;
		}
	}

	double prob2;
	if (child2->isLeaf())
	{
		unsigned int base2 = child2->getBase(site);
		prob2 = _branches[2]->pX1X2(base, base2);
	} else
	{
		prob2 = 0;
		for (unsigned int base2 = 0; base2 < 4; base2++)
		{
			double pRiX1 = child2->pRiX1(base2, site); // TODO: maybe we want to store this
			prob2 += _branches[2]->pX1X2(base, base2) * pRiX1;
		}
	}

	return prob1 * prob2;
}

vector<double> Node::pRiX1(int site)
{
	Node *child1 = _branches[1]->getNeighbour(this);
	Node *child2 = _branches[2]->getNeighbour(this);

	double prob1;
	double prob2;
	vector<double> result(4);
	for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
	{
		if (child1->isLeaf())
			prob1 = _branches[1]->pX1X2(nodeBase, child1->getBase(site));
		else
		{
			prob1 = 0;
			vector<double> childProb = child1->pRiX1(site);
			for (unsigned int childBase = 0; childBase < 4; childBase++)
				prob1 += _branches[1]->pX1X2(nodeBase, childBase) * childProb[childBase];
		}

		if (child2->isLeaf())
			prob2 = _branches[2]->pX1X2(nodeBase, child2->getBase(site));
		else
		{
			prob2 = 0;
			vector<double> childProb = child2->pRiX1(site);
			for (unsigned int childBase = 0; childBase < 4; childBase++)
				prob2 += _branches[2]->pX1X2(nodeBase, childBase) * childProb[childBase];
		}

		result[nodeBase] = prob1 * prob2;
	}

	return result;
}

double Node::pSiX2(Node *blockedNode, unsigned int base, unsigned int site)
{
	Node *parent = getParent();
	Branch *childBranch;
	Node *child = getChild(1);
	if (child != blockedNode)
		childBranch = _branches[1];
	else
	{
		child = getChild(2);
		childBranch = _branches[2];
	}

	double childProb;
	if (child->isLeaf())
	{
		unsigned int base1 = child->getBase(site);
		childProb = childBranch->pX1X2(base, base1);
	} else
	{
		childProb = 0;
		vector<double> pRiX1 = child->pRiX1(site);
		for (unsigned int childBase = 0; childBase < 4; childBase++)
		{
			childProb += childBranch->pX1X2(base, childBase) * pRiX1[childBase];
		}
	}

	double parentProb;
	if (parent->isLeaf())
	{
		unsigned int parentBase = parent->getBase(site);
		parentProb = _branches[0]->getProb(parentBase, base);
	} else
	{
		parentProb = 0;
		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double pSiX2 = parent->pSiX2(this, parentBase, site); // TODO: maybe we want to store this
			parentProb += _branches[0]->pX1X2(parentBase, base) * pSiX2;
		}
	}

	return childProb * parentProb;
}

vector<double> Node::pSiX2(Node *blockedNode, unsigned int site)
{
	Node *parent = getParent();
	Branch *childBranch;
	Node *child = getChild(1);
	if (child != blockedNode)
		childBranch = _branches[1];
	else
	{
		child = getChild(2);
		childBranch = _branches[2];
	}

	double childProb;
	double parentProb;
	vector<double> result(4);
	for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
	{
		if (child->isLeaf())
		{
			childProb = childBranch->pX1X2(nodeBase, child->getBase(site));
		} else
		{
			childProb = 0;
			vector<double> pRiX1 = child->pRiX1(site);
			for (unsigned int childBase = 0; childBase < 4; childBase++)
			{
				childProb += childBranch->pX1X2(nodeBase, childBase) * pRiX1[childBase];
			}
		}

		if (parent->isLeaf())
		{
			parentProb = _branches[0]->getProb(parent->getBase(site), nodeBase);
		} else
		{
			parentProb = 0;
			vector <double> prob = parent->pSiX2(this, site);
			for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
			{
				parentProb += _branches[0]->pX1X2(parentBase, nodeBase) * prob[parentBase];
			}
		}

		result[nodeBase] = childProb * parentProb;
	}

	return result;
}

double Node::computeValuesIntToInt(unsigned int numOfSites)
{
	cout << "computeValuesIntToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	Node *parent = getParent();
	Branch *parentBranch = _branches[0];

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		vector<double> pG1j = parent->pSiX2(this, site);
		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			pG1j[parentBase] /= parentBranch->getMarginalProbCol(parentBase);
			vector<double> pG2j = pRiX1(site);
			for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
			{
				siteProb += parentBranch->getProb(nodeBase, parentBase) * pG1j[parentBase] * pG2j[nodeBase];
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
	Node *parent = getParent();
	Branch *parentBranch = _branches[0];

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;
		unsigned int base = getBase(site);

		/* For each site, consider all 4 values i.e. {A,C,G,T}
		 * Q(Xki,Xji)P(G1j|Xji) sum over all Xji
		 * G1j represents all end nodes connected to node j excluding
		 * the current leaf node
		 */

		//[0][base] is used as the direction of traversal is child -> parent
		vector<double> prob = parent->pSiX2(this, site);
		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = parentBranch->getMarginalProbCol(parentBase);
			siteProb += parentBranch->getProb(parentBase, base) * prob[parentBase] / marginalProb;
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

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;
		unsigned int base = root->getBase(site);

		vector<double> prob = pRiX1(site);
		for (unsigned int j = 0; j < 4; j++)
		{
			siteProb += rootBranch->getProb(base, j) * prob[j];
		}

		logLikelihood += log((1 - _beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}
