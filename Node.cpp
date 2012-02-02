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

vector<double> Node::pRiX1(unsigned int numOfSites)
{
	Node *child1 = _branches[1]->getNeighbour(this);
	Node *child2 = _branches[2]->getNeighbour(this);

	cout << "pRiX1 node=" << getIdent() << " child1=" << child1->getIdent() << " child2=" << child2->getIdent() << endl;

	double prob1;
	double prob2;
	vector<double> result(4 * numOfSites);

	vector<unsigned int> childSeq1;
	vector<double> childProb1;
	if (child1->isLeaf())
		childSeq1 = child1->getSequence();
	else
		childProb1 = child1->pRiX1(numOfSites);

	vector<unsigned int> childSeq2;
	vector<double> childProb2;
	if (child2->isLeaf())
		childSeq2 = child2->getSequence();
	else
		childProb2 = child2->pRiX1(numOfSites);

	for (unsigned int site = 0; site < numOfSites; site++)
	{

		for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
		{
			if (child1->isLeaf())
				prob1 = _branches[1]->pX1X2(nodeBase, childSeq1[site]);
			else
			{
				prob1 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob1 += _branches[1]->pX1X2(nodeBase, childBase) * childProb1[site * 4 + childBase];
			}

			if (child2->isLeaf())
				prob2 = _branches[2]->pX1X2(nodeBase, childSeq2[site]);
			else
			{
				prob2 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob2 += _branches[2]->pX1X2(nodeBase, childBase) * childProb2[site * 4 + childBase];
			}

			result[site * 4 + nodeBase] = prob1 * prob2;
		}
	}

	return result;
}

vector<double> Node::pSiX2(Node *blockedNode, unsigned int numOfSites)
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

	cout << "pSiX2 node=" << getIdent() << " parent=" << parent->getIdent() << " child=" << child->getIdent() << endl;

	vector<unsigned int> childSeq;
	vector<double> childProbRiX1;
	if (child->isLeaf())
		childSeq = child->getSequence();
	else
		childProbRiX1 = child->pRiX1(numOfSites);

	vector<unsigned int> parentSeq;
	vector<double> parentProbSiX2;
	if (parent->isLeaf())
		parentSeq = parent->getSequence();
	else
		parentProbSiX2 = parent->pSiX2(this, numOfSites);

	double childProb;
	double parentProb;
	vector<double> result(4 * numOfSites);
	for (unsigned int site = 0; site < numOfSites; site++)

		for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
		{
			if (child->isLeaf())
			{
				childProb = childBranch->pX1X2(nodeBase, childSeq[site]);
			} else
			{
				childProb = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
				{
					childProb += childBranch->pX1X2(nodeBase, childBase) * childProbRiX1[site * 4 + childBase];
				}
			}

			if (parent->isLeaf())
			{
				parentProb = _branches[0]->getProb(parentSeq[site], nodeBase);
			} else
			{
				parentProb = 0;
				for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
				{
					parentProb += _branches[0]->pX1X2(parentBase, nodeBase) * parentProbSiX2[site * 4 + parentBase];
				}
			}

			result[site * 4 + nodeBase] = childProb * parentProb;
		}

	return result;
}

double Node::computeValuesIntToInt(unsigned int numOfSites)
{
	cout << "computeValuesIntToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	Node *parent = getParent();
	Branch *parentBranch = _branches[0];

	vector<double> pG1j = parent->pSiX2(this, numOfSites);
	vector<double> pG2j = pRiX1(numOfSites);
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			pG1j[site * 4 + parentBase] /= parentBranch->getMarginalProbCol(parentBase);
			for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
			{
				siteProb += parentBranch->getProb(nodeBase, parentBase) * pG1j[site * 4 + parentBase] * pG2j[site * 4 + nodeBase];
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

	vector<double> prob = parent->pSiX2(this, numOfSites);
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = parentBranch->getMarginalProbCol(parentBase);
			siteProb += parentBranch->getProb(parentBase, _seq[site]) * prob[site * 4 + parentBase] / marginalProb;
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

	vector<double> prob = pRiX1(numOfSites);
	vector<unsigned int> rootSeq = root->getSequence();
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int j = 0; j < 4; j++)
		{
			siteProb += rootBranch->getProb(rootSeq[site], j) * prob[site * 4 + j];
		}

		logLikelihood += log((1 - _beta) * siteProb);
	}
	cout << "logLH=" << logLikelihood << endl;

	return logLikelihood;
}
