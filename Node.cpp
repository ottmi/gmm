#include "Node.h"
#include "Branch.h"
#include "globals.h"
#include <cmath>
#include <iostream>
#include <sstream>


Node::Node(Node *parent, int id) // internal node
{
	_id = id;
	_isLeaf = false;
	initialize(parent);
}


Node::Node(Node *parent, int id, const vector<unsigned int> &seq) // leaf node
{
	_id = id;
	_isLeaf = true;
	_seq = seq;
	initialize(parent);
}



Node::~Node()
{
}


void Node::initialize(Node *parent)
{
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

	ss << _id;
	if (_isLeaf)
		ss << "[L]";
	else
		ss << "[I]";

	if (_label.length())
		ss << "(" << _label << ") ";

	ss << _branches.size();

	return ss.str();
}


unsigned int Node::getBase(unsigned int position)
{
	if (position < _seq.size())
		return _seq[position];
	else
	{
		stringstream ss;
		ss << "Node::getBase(" << position << "): invalid position";
		throw(ss.str());
	}
}


double Node::pRiX1(int base, int position)
{
	Node *child1 = _branches[1]->getNeighbour(this);
	Node *child2 = _branches[2]->getNeighbour(this);

	double prob1;
	if (child1->isLeaf())
	{
		unsigned int base1 = child1->getBase(position);
		prob1 = _branches[1]->pX1X2(base, base1);
	} else
	{
		prob1 = 0;
		for (unsigned int base1 = 0; base1 < 4; base1++)
		{
			double pRiX1 = child1->pRiX1(base1, position); // TODO: maybe we want to store this
			prob1+= _branches[1]->pX1X2(base, base1) * pRiX1;
		}
	}

	double prob2;
	if (child2->isLeaf())
	{
		unsigned int base2 = child2->getBase(position);
		prob2 = _branches[2]->pX1X2(base, base2);
	} else
	{
		prob2 = 0;
		for (unsigned int base2 = 0; base2 < 4; base2++)
		{
			double pRiX1 = child2->pRiX1(base2, position); // TODO: maybe we want to store this
			prob2+= _branches[2]->pX1X2(base, base2) * pRiX1;
		}
	}

	return prob1 * prob2;
}


double Node::pSiX2(Node *blockedNode, unsigned int base, unsigned int position)
{
	Node *parent = getParent();
	Branch *childBranch;
	Node *child = getChild(1);
	if (child != blockedNode)
		childBranch = _branches[1];
	else {
		child = getChild(2);
		childBranch = _branches[2];
	}

	double childProb;
	if (child->isLeaf())
	{
		unsigned int base1 = child->getBase(position);
		childProb = childBranch->pX1X2(base, base1);
	} else
	{
		childProb = 0;
		for (unsigned int base1 = 0; base1 < 4; base1++)
		{
			double pRiX1 = child->pRiX1(base1, position); // TODO: maybe we want to store this
			childProb+= childBranch->pX1X2(base, base1) * pRiX1;
		}
	}

	double parentProb;
	if (parent->isLeaf())
	{
		unsigned int base1 = parent->getBase(position);
		parentProb = _branches[0]->getProb(base1, base);
	} else
	{
		parentProb = 0;
		for (unsigned int base1 = 0; base1 < 4; base1++)
		{
			double pSiX2 = parent->pSiX2(this, base1, position); // TODO: maybe we want to store this
			parentProb+= _branches[0]->pX1X2(base1, base) * pSiX2;
		}
	}

	return childProb * parentProb;
}


double Node::pG1jG2j(unsigned int nodeBase, int parentNodeBase, int position)
{

	double pG2j = pRiX1(nodeBase, position);
	Node *parent = getParent();
	double pG1j = parent->pSiX2(this, parentNodeBase, position);
	double sum = _branches[0]->getMarginalProbCol(parentNodeBase);
	pG1j = pG1j / sum;

	//cout << "pG1jG2j(" << nodeBase << "," << parentNodeBase << "," << position << ") pG1j=" << pG1j << " pG2j=" << pG2j << " pG1jG2j=" << pG1j * pG2j << endl;
	return pG1j * pG2j;
}


double Node::computeValuesIntToInt(unsigned int numOfSites)
{
	cout << getIdent() << endl;
	double logLikelihood = 0;
	Branch *parentBranch = _branches[0];

	/* Compute log likelihood for variable sites
	 * i represents the no of sites
	 * k represents {A,C,G,T} for parent
	 * j represents {A,C,G,T} for internal node
	 */
	for (unsigned int i = 0; i < numOfSites; i++)
	{
		double siteProb = 0.0;

		for (unsigned int j = 0; j < 4; j++) // base at parent node
		{
			for (unsigned int k=0; k<4; k++) // base at child node
			{
				//direction of traversal is from parent -> internal node
				//(site, parent node, child node)
				double val = pG1jG2j(k, j, i);
				siteProb+= parentBranch->getProb(j, k) * val;
			} //end of FOR loop for parent node
		} //end of FOR loop for internal node

//		probList.add(new Double(siteProb));
		logLikelihood += log((1-_beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}


double Node::computeValuesIntToLeaf(unsigned int numOfSites)
{
	cout << getIdent() << endl;
	double logLikelihood = 0;
	Node *parent = getParent();
	Branch *parentBranch = _branches[0];

	for(unsigned int i = 0; i < numOfSites; i++)
	{
		double siteProb = 0.0;
		unsigned int base = getBase(i);

		/* For each site, consider all 4 values i.e. {A,C,G,T}
		 * Q(Xki,Xji)P(G1j|Xji) sum over all Xji
		 * G1j represents all end nodes connected to node j excluding
		 * the current leaf node
		 */

		//[0][base] is used as the direction of traversal is child -> parent
		for (unsigned int j = 0; j < 4; j++)
		{
			double prob = parent->pSiX2(this, j, i);
			double marginalProb = parentBranch->getMarginalProbCol(j);
			siteProb+= parentBranch->getProb(j, base) * prob / marginalProb;
		}

		logLikelihood+= log((1-_beta) * siteProb);
	}

	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}


double Node::computeValuesRootToInt(unsigned int numOfSites)
{
	cout << getIdent() << endl;
	double logLikelihood = 0;
	Node *root = getParent();
	Branch *rootBranch = _branches[0];

	for(unsigned int i = 0; i < numOfSites; i++)
	{
		double siteProb = 0.0;
		unsigned int base = root->getBase(i);

		for (unsigned int j = 0; j < 4; j++)
		{
			double prob = pRiX1(j, i);
			siteProb+= rootBranch->getProb(base, j) * prob;
		}

		logLikelihood+= log((1-_beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}
