#include "Branch.h"
#include "globals.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

Branch::Branch(int id, Node *n1, Node *n2)
{
	_id = id;
	_nodes.push_back(n1);
	_nodes.push_back(n2);
	_distance = -1.0;

	_beta = 0.8;
	_q = new Matrix(charStates);
	_q->setDiag(1.0 / 8.0);
	_q->setOffDiag(1.0 / 24.0);

	_qVersion = 0;
	_pRiX1Version = 0;
	_pSiX2Version = 0;
}

Branch::~Branch()
{
	// TODO Auto-generated destructor stub
}

Node* Branch::getNeighbour(Node *node)
{
	if (_nodes[0] == node)
		return _nodes[1];
	else if (_nodes[1] == node)
		return _nodes[0];
	else
		return NULL;
}

void Branch::setDistance(double &distance)
{
	this->_distance = distance;
	distance = -1.0;
}

void Branch::print()
{
	cout << getIdent() << endl;
	_q->print();
}

string Branch::getIdent()
{
	stringstream ss;
	ss << _nodes[0]->getIdent() << "<---(" << _distance << ")--->" << _nodes[1]->getIdent();
	return ss.str();
}

double Branch::computeLH(unsigned int numOfSites)
{
	double lh;

	if (_nodes[0]->isLeaf())
		lh = computeValuesRootToInt(numOfSites);
	else if (_nodes[1]->isLeaf())
		lh = computeValuesIntToLeaf(numOfSites);
	else
		lh = computeValuesIntToInt(numOfSites);

	cout << "logLH=" << fixed << setprecision(10) << lh << endl << endl;
	return lh;
}

void Branch::updateQ()
{
	_q->update(*_updatedQ);
	_qVersion++;
	free(_updatedQ);
}

double Branch::computeValuesIntToInt(unsigned int numOfSites)
{
	cout << "computeValuesIntToInt() " << getIdent() << endl;
	double logLikelihood = 0;

	if (verbose >= 3) _q->print();

	vector<double> pG1 = pSiX2(numOfSites);
	vector<double> pG2 = pRiX1(numOfSites);
	if (_siteProb.empty()) _siteProb = vector<double>(numOfSites);
	Branch *grandParentBranch = _nodes[0]->getBranch(0);

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			for (unsigned int childBase = 0; childBase < 4; childBase++)
			{
				siteProb += getProb(parentBase, childBase) * pG1[site * 4 + parentBase] * pG2[site * 4 + childBase] / marginalProb;
			}
		}
		_siteProb[site] = siteProb;
		logLikelihood += log((1 - _beta) * siteProb);
	}
	updateQIntToInt(numOfSites);

	return logLikelihood;
}

void Branch::updateQIntToInt(unsigned int numOfSites)
{
	vector<double> &pG1 = _pSiX2;
	vector<double> &pG2 = _pRiX1;
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	vector<vector<double> > sum(4, vector<double>(4, 0.0));

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			for (unsigned int childBase = 0; childBase < 4; childBase++)
			{
				double siteProb = getProb(parentBase, childBase) * pG1[site * 4 + parentBase] * pG2[site * 4 + childBase] / marginalProb;
				double denominator = _siteProb[site];
				sum[parentBase][childBase] += siteProb / denominator;
			}
		}
	}

	_updatedQ = new Matrix(4);
	for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		for (unsigned int childBase = 0; childBase < 4; childBase++)
			_updatedQ->setEntry(parentBase, childBase, sum[parentBase][childBase] / numOfSites);
}

double Branch::computeValuesIntToLeaf(unsigned int numOfSites)
{
	cout << "computeValuesIntToLeaf() " << getIdent() << endl;
	double logLikelihood = 0;

	if (verbose >= 3) _q->print();

	vector<double> pG1 = pSiX2(numOfSites);
	vector<unsigned int> leafSeq = _nodes[1]->getSequence();
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	if (_siteProb.empty()) _siteProb = vector<double>(numOfSites);

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			siteProb += getProb(parentBase, leafSeq[site]) * pG1[site * 4 + parentBase] / marginalProb;
		}
		_siteProb[site] = siteProb;
		logLikelihood += log((1 - _beta) * siteProb);
	}
	updateQIntToLeaf(numOfSites);

	return logLikelihood;
}

void Branch::updateQIntToLeaf(unsigned int numOfSites)
{
	vector<double> &pG1 = _pSiX2;
	vector<unsigned int> leafSeq = _nodes[1]->getSequence();
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	vector<vector<double> > sum(4, vector<double>(4, 0));

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			unsigned int childBase = leafSeq[site];
			double siteProb = getProb(parentBase, childBase) * pG1[site * 4 + parentBase] / marginalProb;
			double denominator = _siteProb[site];
			sum[parentBase][childBase] += siteProb / denominator;
		}
	}

	_updatedQ = new Matrix(4);
	for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		for (unsigned int childBase = 0; childBase < 4; childBase++)
			_updatedQ->setEntry(parentBase, childBase, sum[parentBase][childBase] / numOfSites);
}

double Branch::computeValuesRootToInt(unsigned int numOfSites)
{
	cout << "computeValuesRootToInt() " << getIdent() << endl;
	double logLikelihood = 0;

	if (verbose >= 3) _q->print();

	vector<double> pG2 = pRiX1(numOfSites);
	vector<unsigned int> rootSeq = _nodes[0]->getSequence();
	if (_siteProb.empty()) _siteProb = vector<double>(numOfSites);

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int j = 0; j < 4; j++)
		{
			siteProb += getProb(rootSeq[site], j) * pG2[site * 4 + j];
		}
		_siteProb[site] = siteProb;
		logLikelihood += log((1 - _beta) * siteProb);
	}

	updateQRootToInt(numOfSites);
	return logLikelihood;
}

void Branch::updateQRootToInt(unsigned int numOfSites)
{
	vector<unsigned int> rootSeq = _nodes[0]->getSequence();
	vector<double> &pG2 = _pRiX1;
	vector<vector<double> > sum(4, vector<double>(4, 0));

	for (unsigned int site = 0; site < numOfSites; site++)
	{
		unsigned int rootBase = rootSeq[site];
		for (unsigned int childBase = 0; childBase < 4; childBase++)
		{
			double siteProb = getProb(rootBase, childBase) * pG2[site * 4 + childBase];
			double denominator = _siteProb[site];
			sum[rootBase][childBase] += siteProb / denominator;
		}
	}

	_updatedQ = new Matrix(4);
	for (unsigned int rootBase = 0; rootBase < 4; rootBase++)
		for (unsigned int childBase = 0; childBase < 4; childBase++)
			_updatedQ->setEntry(rootBase, childBase, sum[rootBase][childBase] / numOfSites);
}

/* This method computes the conditional probability P(x2|x1) */
double Branch::pX1X2(unsigned int parent, unsigned int child)
{
	double marginalProb = _q->getRowSum(parent);
	double condProb = _q->getEntry(parent, child) / marginalProb;

	return condProb;
}

// probability away from root
vector<double>& Branch::pRiX1(unsigned int numOfSites)
{
	Node *parent = _nodes[1]; // this should be the node away from root
	Branch *childBranch1 = parent->getBranch(1);
	Branch *childBranch2 = parent->getBranch(2);
	Node *child1 = childBranch1->getNeighbour(parent);
	Node *child2 = childBranch2->getNeighbour(parent);

	if (verbose >= 5)
		cout << "Branch::pRiX1 parent=" << parent->getIdent() << " child1=" << child1->getIdent() << " child2=" << child2->getIdent() << " qVer=" << _qVersion
				<< " ownVer=" << _pRiX1Version << endl;

	if (_pRiX1.empty())
		_pRiX1 = vector<double>(4 * numOfSites);
	else if (_qVersion < _pRiX1Version) return _pRiX1;

	vector<unsigned int> childSeq1;
	vector<double> childProb1;
	if (child1->isLeaf())
		childSeq1 = child1->getSequence();
	else
		childProb1 = childBranch1->pRiX1(numOfSites);

	vector<unsigned int> childSeq2;
	vector<double> childProb2;
	if (child2->isLeaf())
		childSeq2 = child2->getSequence();
	else
		childProb2 = childBranch2->pRiX1(numOfSites);

	double prob1;
	double prob2;
	for (unsigned int site = 0; site < numOfSites; site++)
	{

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			if (child1->isLeaf())
				prob1 = childBranch1->pX1X2(parentBase, childSeq1[site]);
			else
			{
				prob1 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob1 += childBranch1->pX1X2(parentBase, childBase) * childProb1[site * 4 + childBase];
			}

			if (child2->isLeaf())
				prob2 = childBranch2->pX1X2(parentBase, childSeq2[site]);
			else
			{
				prob2 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob2 += childBranch2->pX1X2(parentBase, childBase) * childProb2[site * 4 + childBase];
			}

			_pRiX1[site * 4 + parentBase] = prob1 * prob2;
		}
	}

	_pRiX1Version++;
	return _pRiX1;
}

// probability towards root
vector<double>& Branch::pSiX2(unsigned int numOfSites)
{
	Node *parent = _nodes[0]; // this should be the node towards the root
	Node *grandParent = parent->getParent();
	Branch *grandParentBranch = parent->getBranch(0);
	Branch *siblingBranch;
	Node *sibling = parent->getChild(1);
	if (sibling != _nodes[1])
		siblingBranch = parent->getBranch(1);
	else
	{
		sibling = parent->getChild(2);
		siblingBranch = parent->getBranch(2);
	}

	if (verbose >= 5)
		cout << "Branch::pSiX2 parent=" << parent->getIdent() << " grandParent=" << grandParent->getIdent() << " sibling=" << sibling->getIdent() << " qVer="
				<< _qVersion << " ownVer=" << _pSiX2Version << endl;

	if (_pSiX2.empty())
		_pSiX2 = vector<double>(4 * numOfSites);
	else if (_qVersion < _pSiX2Version) return _pSiX2;

	vector<unsigned int> siblingSeq;
	vector<double> siblingProbRiX1;
	if (sibling->isLeaf())
		siblingSeq = sibling->getSequence();
	else
		siblingProbRiX1 = siblingBranch->pRiX1(numOfSites);

	vector<unsigned int> grandParentSeq;
	vector<double> grandParentProbSiX2;
	if (grandParent->isLeaf())
		grandParentSeq = grandParent->getSequence();
	else
		grandParentProbSiX2 = grandParentBranch->pSiX2(numOfSites);

	double siblingProb;
	double grandParentProb;
	for (unsigned int site = 0; site < numOfSites; site++)

		for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
		{
			if (sibling->isLeaf())
			{
				siblingProb = siblingBranch->pX1X2(parentBase, siblingSeq[site]);
			} else
			{
				siblingProb = 0;
				for (unsigned int siblingBase = 0; siblingBase < 4; siblingBase++)
				{
					siblingProb += siblingBranch->pX1X2(parentBase, siblingBase) * siblingProbRiX1[site * 4 + siblingBase];
				}
			}

			if (grandParent->isLeaf())
			{
				grandParentProb = grandParentBranch->getProb(grandParentSeq[site], parentBase);
			} else
			{
				grandParentProb = 0;
				for (unsigned int grandParentBase = 0; grandParentBase < 4; grandParentBase++)
				{
					grandParentProb += grandParentBranch->pX1X2(grandParentBase, parentBase) * grandParentProbSiX2[site * 4 + grandParentBase];
				}
			}

			_pSiX2[site * 4 + parentBase] = siblingProb * grandParentProb;
		}

	_pSiX2Version++;
	return _pSiX2;
}

double Branch::getProb(unsigned int from, unsigned int to)
{
	return _q->getEntry(from, to);
}

double Branch::getMarginalProbRow(unsigned int row)
{
	return _q->getRowSum(row);
}

double Branch::getMarginalProbCol(unsigned int col)
{
	return _q->getColSum(col);
}
