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
	_q->print();

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
			for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
			{
				siteProb += getProb(parentBase, nodeBase) * pG1[site * 4 + parentBase] * pG2[site * 4 + nodeBase] / marginalProb;
			}
		}
		_siteProb[site] = siteProb;
//		cout << "i=" << site << " siteProb=" << siteProb << endl;

		logLikelihood += log((1 - _beta) * siteProb);
	}
	updateQIntToInt(numOfSites);

	return logLikelihood;
}

void Branch::updateQIntToInt(unsigned int numOfSites)
{
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	_updatedQ = new Matrix(4);

	double sum;
	for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
	{
		double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
		for (unsigned int childBase = 0; childBase < 4; childBase++)
		{
			sum = 0.0;
			for (unsigned int site = 0; site < numOfSites; site++)
			{
					double siteProb = getProb(parentBase, childBase) * _pSiX2[site * 4 + parentBase] * _pRiX1[site * 4 + childBase] / marginalProb;
					double denominator = _siteProb[site];
					sum += siteProb / denominator;
					//cout << "i=" << parentBase << " j=" << childBase << " k=" << site << " siteProb=" << siteProb << " denom=" << denominator << " sum=" << sum << endl;
			}

			_updatedQ->setEntry(parentBase, childBase, sum/numOfSites);
		} //end of FOR loop for child node
	} //end of FOR loop for parent node
	//_updatedQ->print();
	//_q->update(newQ);
}

double Branch::computeValuesIntToLeaf(unsigned int numOfSites)
{
	cout << "computeValuesIntToLeaf() " << getIdent() << endl;
	double logLikelihood = 0;
	_q->print();

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
	vector<unsigned int> leafSeq = _nodes[1]->getSequence();
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	_updatedQ = new Matrix(4);

	double sum;
	for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
	{
		double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
		for (unsigned int childBase = 0; childBase < 4; childBase++)
		{
			sum = 0.0;
			for (unsigned int site = 0; site < numOfSites; site++)
			{
				if (leafSeq[site] == childBase)
				{
					double siteProb = getProb(parentBase, childBase) * _pSiX2[site * 4 + parentBase] / marginalProb;
					double denominator = _siteProb[site];
					sum += siteProb / denominator;
					//cout << "i=" << parentBase << " j=" << childBase << " k=" << site << " siteProb=" << siteProb << " denom=" << denominator << " sum=" << sum << endl;
				}
			}

			_updatedQ->setEntry(parentBase, childBase, sum/numOfSites);
		} //end of FOR loop for child node
	} //end of FOR loop for parent node
	//_updatedQ->print();
	//_q->update(newQ);
}


double Branch::computeValuesRootToInt(unsigned int numOfSites)
{
	cout << "computeValuesRootToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	_q->print();

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
//		cout << "i=" << site << " siteProb=" << siteProb << endl;
		logLikelihood += log((1 - _beta) * siteProb);
	}

	updateQRootToInt(numOfSites);
	return logLikelihood;
}

void Branch::updateQRootToInt(unsigned int numOfSites)
{
	vector<unsigned int> rootSeq = _nodes[0]->getSequence();
	_updatedQ = new Matrix(4);

	double sum;
	for (unsigned int rootBase = 0; rootBase < 4; rootBase++)
	{
		for (unsigned int childBase = 0; childBase < 4; childBase++)
		{
			sum = 0.0;
			for (unsigned int site = 0; site < numOfSites; site++)
			{
				if (rootSeq[site] == rootBase)
				{
					double siteProb = getProb(rootBase, childBase) * _pRiX1[site * 4 + childBase];
					double denominator = _siteProb[site];
					sum += siteProb / denominator;
//					cout << "i=" << rootBase << " j=" << childBase << " k=" << site << " siteProb=" << siteProb << " denom=" << denominator << " sum=" << sum << endl;
				}
			}

			_updatedQ->setEntry(rootBase, childBase, sum/numOfSites);
		} //end of FOR loop for child node
	} //end of FOR loop for parent node
	//_updatedQ->print();
	//_q->update(newQ);
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
	Node *node = _nodes[1]; // this should be the node away from root
	Branch *childBranch1 = node->getBranch(1);
	Branch *childBranch2 = node->getBranch(2);
	Node *child1 = childBranch1->getNeighbour(node);
	Node *child2 = childBranch2->getNeighbour(node);

	cout << "Branch::pRiX1 node=" << node->getIdent() << " child1=" << child1->getIdent() << " child2=" << child2->getIdent() << " qVer=" << _qVersion << " ownVer=" << _pRiX1Version << endl;

	if (_pRiX1.empty())
		_pRiX1 = vector<double>(4 * numOfSites);
	else if (_qVersion < _pRiX1Version)
		return _pRiX1;

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

		for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
		{
			if (child1->isLeaf())
				prob1 = childBranch1->pX1X2(nodeBase, childSeq1[site]);
			else
			{
				prob1 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob1 += childBranch1->pX1X2(nodeBase, childBase) * childProb1[site * 4 + childBase];
			}

			if (child2->isLeaf())
				prob2 = childBranch2->pX1X2(nodeBase, childSeq2[site]);
			else
			{
				prob2 = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
					prob2 += childBranch2->pX1X2(nodeBase, childBase) * childProb2[site * 4 + childBase];
			}

			_pRiX1[site * 4 + nodeBase] = prob1 * prob2;
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

	cout << "Branch::pSiX2 parent=" << parent->getIdent() << " grandParent=" << grandParent->getIdent() << " sibling=" << sibling->getIdent() << " qVer=" << _qVersion << " ownVer=" << _pSiX2Version << endl;

	if (_pSiX2.empty())
		_pSiX2 = vector<double>(4 * numOfSites);
	else if (_qVersion < _pSiX2Version)
		return _pSiX2;


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

		for (unsigned int nodeBase = 0; nodeBase < 4; nodeBase++)
		{
			if (sibling->isLeaf())
			{
				siblingProb = siblingBranch->pX1X2(nodeBase, siblingSeq[site]);
			} else
			{
				siblingProb = 0;
				for (unsigned int childBase = 0; childBase < 4; childBase++)
				{
					siblingProb += siblingBranch->pX1X2(nodeBase, childBase) * siblingProbRiX1[site * 4 + childBase];
				}
			}

			if (grandParent->isLeaf())
			{
				grandParentProb = grandParentBranch->getProb(grandParentSeq[site], nodeBase);
			} else
			{
				grandParentProb = 0;
				for (unsigned int parentBase = 0; parentBase < 4; parentBase++)
				{
					grandParentProb += grandParentBranch->pX1X2(parentBase, nodeBase) * grandParentProbSiX2[site * 4 + parentBase];
				}
			}

			_pSiX2[site * 4 + nodeBase] = siblingProb * grandParentProb;
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
