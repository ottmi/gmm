#include "Branch.h"
#include "globals.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

Branch::Branch(int id, Node *n1, Node *n2)
{
	_id = id;
	_nodes.push_back(n1);
	_nodes.push_back(n2);

	_beta = 0.8;
	for (unsigned int i = 0; i < charStates; i++)
		_invar.push_back(0.25);

	_q = new Matrix(charStates);
	_q->setDiag(1.0 / (2 * charStates));
	_q->setOffDiag(1.0 / (6 * charStates));

	_qVersion = 0;
	_pRiX1 = NULL;
	_pRiX1Version = 0;
	_pSiX2 = NULL;
	_pSiX2Version = 0;
	_siteProb = NULL;
}

Branch::~Branch()
{
	free(_pRiX1);
	free(_pSiX2);
	free(_siteProb);
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

Node* Branch::getNode(unsigned int num)
{
	if (num > _nodes.size() - 1)
	{
		stringstream ss;
		ss << "Branch(" << _id << ")::getNode(" << num << "): there are only " << _nodes.size() << " adjacent nodes";
	}
	return _nodes[num];
}

void Branch::print()
{
	cout << getIdent() << endl;
	_q->print();
}

string Branch::getIdent()
{
	stringstream ss;
	ss << _nodes[0]->getIdent() << "<---(" << _id << ")---" << _nodes[1]->getIdent();
	return ss.str();
}

double Branch::getLength()
{
	return -0.25 * log(_q->determinant());
}

void Branch::unlinkNode(Node *node)
{
	if (_nodes[0] == node)
		_nodes.erase(_nodes.begin());
	else if (_nodes[1] == node)
		_nodes.erase(_nodes.begin() + 1);
	else
	{
		stringstream ss;
		ss << "Branch(" << _id << ")::unlinkNode(" << node->getIdent() << "): not found";
		throw(ss.str());
	}
	node->removeBranch(this);
	if (verbose >= 5)
		cout << "Unlinked branch " <<  getId() << " from node " << node->getIdent() << endl;

}

void Branch::linkNode(Node *node)
{
	if (_nodes.size() == 2)
	{
		stringstream ss;
		ss << "Branch(" << _id << ")::linkNode(" << node->getIdent() << "): there are already " << _nodes.size() << " adjacent nodes";
		throw(ss.str());
	} else
	{
		_nodes.push_back(node);
		node->addBranch(this);
		if (verbose >= 5)
			cout << "Linked branch " <<  getId() << " to node " << node->getIdent() << endl;
	}
}

void Branch::swapNodes()
{
	Node* temp = _nodes[0];
	_nodes[0] = _nodes[1];
	_nodes[1] = temp;
}

double Branch::computeLH(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	double lh;

	if (_nodes[0]->isLeaf())
		lh = computeValuesRootToInt(patternCount, invarSites, invarStart);
	else if (_nodes[1]->isLeaf())
		lh = computeValuesIntToLeaf(patternCount, invarSites, invarStart);
	else
		lh = computeValuesIntToInt(patternCount, invarSites, invarStart);

	cout << "logLH=" << fixed << setprecision(10) << lh << endl << endl;
	return lh;
}

void Branch::computeUpdatedQ(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	if (_nodes[0]->isLeaf())
		updateQRootToInt(numOfSites, patternCount, invarSites, invarStart);
	else if (_nodes[1]->isLeaf())
		updateQIntToLeaf(numOfSites, patternCount, invarSites, invarStart);
	else
		updateQIntToInt(numOfSites, patternCount, invarSites, invarStart);

	double sum = 0;
	for (unsigned int row = 0; row < charStates; row++)
		for (unsigned int col = 0; col < charStates; col++)
			sum += _updatedQ->getEntry(row, col) * numOfSites;

	for (unsigned int row = 0; row < charStates; row++)
		for (unsigned int col = 0; col < charStates; col++)
			_updatedQ->setEntry(row, col, _updatedQ->getEntry(row, col) * numOfSites / sum);

}

bool Branch::updateQ(double qDelta)
{
	double sum = 0;
	for (unsigned int row = 0; row < charStates; row++)
		for (unsigned int col = 0; col < charStates; col++)
			sum += pow(_q->getEntry(row, col) - _updatedQ->getEntry(row, col), 2);

	free(_q);
	_q = _updatedQ;
	_updatedQ = NULL;
	_qVersion++;

	return (bool) (sum > qDelta);
}

bool Branch::updateParameters(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart,
		double betaDelta)
{
	unsigned int numOfUniqueSites = patternCount.size();
	unsigned int siteCount = numOfSites;
	for (unsigned int site = invarStart; site < numOfUniqueSites; site++)
		siteCount -= patternCount[site];

	double alphaSum = 0;
	double betaSum = 0;
	double invarSum = 0;
	vector<double> invar(charStates, 0.0);
	for (unsigned int site = invarStart; site < numOfUniqueSites; site++)
	{
		unsigned int invarChar = invarSites[site - invarStart];
		double denominator = (1 - _beta) * _siteProb[site] + _beta * _invar[invarChar];
		alphaSum += patternCount[site] * _siteProb[site] / denominator;
		betaSum += patternCount[site] * _invar[invarChar] / denominator;
		invar[invarChar] = _invar[invarChar] * patternCount[site] * _beta / denominator;
		invarSum += invar[invarChar];
	}
	double alpha = (1 - _beta) * (siteCount / (1 - _beta) + alphaSum);
	double beta = _beta * betaSum;
	beta = beta / (alpha + beta);

	bool updateRequired = (fabs(_beta - beta) > betaDelta);

	_beta = beta;
	for (unsigned int i = 0; i < charStates; i++)
		_invar[i] = invar[i] / invarSum;

	return updateRequired;
}

void Branch::NNI(int branch)
{
	// branch indicates which branch of the left node will be swapped, i.e. 0 for the 1st branch, 1 for the 2nd

	if (_nodes[0]->isLeaf())
		throw("NNI: Node " + _nodes[0]->getIdent() + " is a leaf node");
	else if (_nodes[1]->isLeaf()) throw("NNI: Node " + _nodes[1]->getIdent() + " is a leaf node");

	Branch *branch0; // this is the outgoing branch from the left node we are going to swap
	for (unsigned int i = 0; i < 3; i++)
	{
		branch0 = _nodes[0]->getBranch(i);
		if (branch0 != this && !branch--) break;
	}

	Branch *branch1; // this is the outgoing branch from the right node we are going to swap
	for (unsigned int i = 0; i < 3; i++)
	{
		branch1 = _nodes[1]->getBranch(i);
		if (branch1 != this) break;
	}

	//unlink the branch from the left node
	branch0->unlinkNode(_nodes[0]);

	//unlink the branch from the right node
	branch1->unlinkNode(_nodes[1]);

	//link the right node's branch to the left node
	branch1->linkNode(_nodes[0]);

	//link the left node's branch to the right node
	branch0->linkNode(_nodes[1]);
}

void Branch::SPR(Node *parentNode, Branch *insertBranch)
{
	if (verbose >= 2)
		cout << "Branch(" << getId() << ")::SPR(" << parentNode->getIdent() << ", " << insertBranch->getId() << ")" << endl;

	Node *childNode = getNeighbour(parentNode);

	if (parentNode->isLeaf()) throw("Branch::SPR() parentNode is a leaf node !");

	/* Find the nearest internal node towards the root and the branch leading to it
	 * This will usually be our grandparent, but it has to be an internal node
	 * So if it's the actual root (which is a leaf) we choose the next internal node instead */
	Node *grandParentNode = parentNode->getParent();
	if (grandParentNode->isLeaf())
		grandParentNode = parentNode->getNeighbour(childNode, grandParentNode);
	Branch *grandParentBranch = parentNode->getNeighbourBranch(grandParentNode);

	cout << "Grandparent: " << grandParentNode->getIdent() << " Branch: " << grandParentBranch->getIdent() << endl;

	// Identify the branch that leads to our sibling
	Branch *siblingBranch = parentNode->getNeighbourBranch(childNode, grandParentNode);

	// Unlink our sibling including its branch from our parent
	siblingBranch->unlinkNode(parentNode);

	// Unlink our parent including its branch from its parent.
	grandParentBranch->unlinkNode(grandParentNode);

	// Link our sibling there instead
	siblingBranch->linkNode(grandParentNode);

	// That's the parent node of the branch where we're going to insert our subtree
	Node *insertParentNode = insertBranch->getNode(0);

	// Unlink the insertion branch from its parent
	insertBranch->unlinkNode(insertParentNode);

	// Link the branch leading into our parent there instead
	grandParentBranch->linkNode(insertParentNode);

	// Link the insertion branch as sibling to our parent node
	insertBranch->linkNode(parentNode);
}


/*
 * Private functions start here
 */

double Branch::computeValuesIntToInt(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	cout << "computeValuesIntToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	unsigned int numOfUniqueSites = patternCount.size();

	if (verbose >= 3)
	{
		_q->print();
		cout << "beta=" << _beta << " invar[0]=" << _invar[0] << " invar[1]=" << _invar[1] << " invar[2]=" << _invar[2] << " invar[3]=" << _invar[3] << endl;
	}

	double* pG1 = pSiX2(numOfUniqueSites);
	double* pG2 = pRiX1(numOfUniqueSites);
	if (_siteProb == NULL) _siteProb = new double[numOfUniqueSites];
	Branch *grandParentBranch = _nodes[0]->getBranch(0);

#ifdef _OPENMP
#pragma omp parallel for reduction(+:logLikelihood)
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			for (unsigned int childBase = 0; childBase < charStates; childBase++)
			{
				siteProb += getProb(parentBase, childBase) * pG1[site * charStates + parentBase] * pG2[site * charStates + childBase] / marginalProb;
			}
		}
		_siteProb[site] = siteProb;

		if (site < invarStart)
			logLikelihood += patternCount[site] * log((1 - _beta) * siteProb);
		else
		{
			unsigned int invarChar = invarSites[site - invarStart];
			logLikelihood += patternCount[site] * log(_beta * _invar[invarChar] + (1 - _beta) * siteProb);
		}
	}

	return logLikelihood;
}

void Branch::updateQIntToInt(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	unsigned int numOfUniqueSites = patternCount.size();
	double *pG1 = _pSiX2;
	double *pG2 = _pRiX1;
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	double sum[omp_get_max_threads()][charStates][charStates];
	memset(sum, 0.0, sizeof(sum));

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			for (unsigned int childBase = 0; childBase < charStates; childBase++)
			{
				double siteProb = getProb(parentBase, childBase) * pG1[site * charStates + parentBase] * pG2[site * charStates + childBase] / marginalProb;
				double denominator = _siteProb[site];

				if (site < invarStart)
					sum[omp_get_thread_num()][parentBase][childBase] += patternCount[site] * siteProb / denominator;
				else
				{
					unsigned int invarChar = invarSites[site - invarStart];
					sum[omp_get_thread_num()][parentBase][childBase] += patternCount[site] * (1 - _beta) * siteProb
							/ (_beta * _invar[invarChar] + (1 - _beta) * denominator);
				}
			}
		}
	}

#ifdef _OPENMP
	for (int i = 1; i < omp_get_max_threads(); i++)
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
			for (unsigned int childBase = 0; childBase < charStates; childBase++)
				sum[0][parentBase][childBase] += sum[i][parentBase][childBase];
#endif

	_updatedQ = new Matrix(charStates);
	for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		for (unsigned int childBase = 0; childBase < charStates; childBase++)
			_updatedQ->setEntry(parentBase, childBase, sum[0][parentBase][childBase] / numOfSites);
}

double Branch::computeValuesIntToLeaf(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	cout << "computeValuesIntToLeaf() " << getIdent() << endl;
	double logLikelihood = 0;
	unsigned int numOfUniqueSites = patternCount.size();

	if (verbose >= 3)
	{
		_q->print();
		cout << "beta=" << _beta << " invar[0]=" << _invar[0] << " invar[1]=" << _invar[1] << " invar[2]=" << _invar[2] << " invar[3]=" << _invar[3] << endl;
	}

	double* pG1 = pSiX2(numOfUniqueSites);
	vector<unsigned int> leafSeq = _nodes[1]->getSequence();
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	if (_siteProb == NULL) _siteProb = new double[numOfUniqueSites];

#ifdef _OPENMP
#pragma omp parallel for reduction(+:logLikelihood)
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			siteProb += getProb(parentBase, leafSeq[site]) * pG1[site * charStates + parentBase] / marginalProb;
		}
		_siteProb[site] = siteProb;

		if (site < invarStart)
			logLikelihood += patternCount[site] * log((1 - _beta) * siteProb);
		else
		{
			unsigned int invarChar = invarSites[site - invarStart];
			logLikelihood += patternCount[site] * log(_beta * _invar[invarChar] + (1 - _beta) * siteProb);
		}
	}

	return logLikelihood;
}

void Branch::updateQIntToLeaf(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	unsigned int numOfUniqueSites = patternCount.size();
	double* pG1 = _pSiX2;
	vector<unsigned int> leafSeq = _nodes[1]->getSequence();
	Branch *grandParentBranch = _nodes[0]->getBranch(0);
	double sum[omp_get_max_threads()][charStates][charStates];
	memset(sum, 0.0, sizeof(sum));

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double marginalProb = grandParentBranch->getMarginalProbCol(parentBase);
			unsigned int childBase = leafSeq[site];
			double siteProb = getProb(parentBase, childBase) * pG1[site * charStates + parentBase] / marginalProb;
			double denominator = _siteProb[site];

			if (site < invarStart)
				sum[omp_get_thread_num()][parentBase][childBase] += patternCount[site] * siteProb / denominator;
			else
			{
				unsigned int invarChar = invarSites[site - invarStart];
				sum[omp_get_thread_num()][parentBase][childBase] += patternCount[site] * (1 - _beta) * siteProb
						/ (_beta * _invar[invarChar] + (1 - _beta) * denominator);
			}
		}
	}

#ifdef _OPENMP
	for (int i = 1; i < omp_get_max_threads(); i++)
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
			for (unsigned int childBase = 0; childBase < charStates; childBase++)
				sum[0][parentBase][childBase] += sum[i][parentBase][childBase];
#endif

	_updatedQ = new Matrix(charStates);
	for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		for (unsigned int childBase = 0; childBase < charStates; childBase++)
			_updatedQ->setEntry(parentBase, childBase, sum[0][parentBase][childBase] / numOfSites);
}

double Branch::computeValuesRootToInt(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	cout << "computeValuesRootToInt() " << getIdent() << endl;
	double logLikelihood = 0;
	unsigned int numOfUniqueSites = patternCount.size();

	if (verbose >= 3)
	{
		_q->print();
		cout << "beta=" << _beta << " invar[0]=" << _invar[0] << " invar[1]=" << _invar[1] << " invar[2]=" << _invar[2] << " invar[3]=" << _invar[3] << endl;
	}

	double* pG2 = pRiX1(numOfUniqueSites);
	vector<unsigned int> rootSeq = _nodes[0]->getSequence();
	if (_siteProb == NULL) _siteProb = new double[numOfUniqueSites];

#ifdef _OPENMP
#pragma omp parallel for reduction(+:logLikelihood)
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		double siteProb = 0.0;

		for (unsigned int j = 0; j < charStates; j++)
		{
			siteProb += getProb(rootSeq[site], j) * pG2[site * charStates + j];
		}
		_siteProb[site] = siteProb;

		if (site < invarStart)
			logLikelihood += patternCount[site] * log((1 - _beta) * siteProb);
		else
		{
			unsigned int invarChar = invarSites[site - invarStart];
			logLikelihood += patternCount[site] * log(_beta * _invar[invarChar] + (1 - _beta) * siteProb);
		}
	}

	return logLikelihood;
}

void Branch::updateQRootToInt(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart)
{
	unsigned int numOfUniqueSites = patternCount.size();
	vector<unsigned int> rootSeq = _nodes[0]->getSequence();
	double* pG2 = _pRiX1;
	double sum[omp_get_max_threads()][charStates][charStates];
	memset(sum, 0.0, sizeof(sum));

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < numOfUniqueSites; site++)
	{
		unsigned int rootBase = rootSeq[site];
		for (unsigned int childBase = 0; childBase < charStates; childBase++)
		{
			double siteProb = getProb(rootBase, childBase) * pG2[site * charStates + childBase];
			double denominator = _siteProb[site];

			if (site < invarStart)
				sum[omp_get_thread_num()][rootBase][childBase] += patternCount[site] * siteProb / denominator;
			else
			{
				unsigned int invarChar = invarSites[site - invarStart];
				sum[omp_get_thread_num()][rootBase][childBase] += patternCount[site] * (1 - _beta) * siteProb / (_beta * _invar[invarChar] + (1 - _beta) * denominator);
			}
		}
	}

	#ifdef _OPENMP
	for (int i = 1; i < omp_get_max_threads(); i++)
		for (unsigned int rootBase = 0; rootBase < charStates; rootBase++)
			for (unsigned int childBase = 0; childBase < charStates; childBase++)
				sum[0][rootBase][childBase] += sum[i][rootBase][childBase];
#endif

	_updatedQ = new Matrix(charStates);
	for (unsigned int rootBase = 0; rootBase < charStates; rootBase++)
		for (unsigned int childBase = 0; childBase < charStates; childBase++)
			_updatedQ->setEntry(rootBase, childBase, sum[0][rootBase][childBase] / numOfSites);
}

/* This method computes the conditional probability P(x2|x1) */
double Branch::pX1X2(unsigned int parent, unsigned int child)
{
	double marginalProb = _q->getRowSum(parent);
	double condProb = _q->getEntry(parent, child) / marginalProb;

	return condProb;
}

// probability away from root
double* Branch::pRiX1(unsigned int numOfSites)
{
	Node *parent = _nodes[1]; // this should be the node away from root
	Branch *childBranch1 = parent->getBranch(1);
	Branch *childBranch2 = parent->getBranch(2);
	Node *child1 = childBranch1->getNeighbour(parent);
	Node *child2 = childBranch2->getNeighbour(parent);

	if (verbose >= 5)
		cout << "Branch::pRiX1 parent=" << parent->getIdent() << " child1=" << child1->getIdent() << " child2=" << child2->getIdent() << " qVer=" << _qVersion
				<< " ownVer=" << _pRiX1Version << endl;

	if (_pRiX1 == NULL)
		_pRiX1 = new double[charStates * numOfSites];
	else if (_qVersion < _pRiX1Version) return _pRiX1;

	vector<unsigned int> childSeq1;
	double* childProb1;
	if (child1->isLeaf())
		childSeq1 = child1->getSequence();
	else
		childProb1 = childBranch1->pRiX1(numOfSites);

	vector<unsigned int> childSeq2;
	double* childProb2;
	if (child2->isLeaf())
		childSeq2 = child2->getSequence();
	else
		childProb2 = childBranch2->pRiX1(numOfSites);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < numOfSites; site++)
	{
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double prob1;
			double prob2;
			if (child1->isLeaf())
				prob1 = childBranch1->pX1X2(parentBase, childSeq1[site]);
			else
			{
				prob1 = 0;
				for (unsigned int childBase = 0; childBase < charStates; childBase++)
					prob1 += childBranch1->pX1X2(parentBase, childBase) * childProb1[site * charStates + childBase];
			}

			if (child2->isLeaf())
				prob2 = childBranch2->pX1X2(parentBase, childSeq2[site]);
			else
			{
				prob2 = 0;
				for (unsigned int childBase = 0; childBase < charStates; childBase++)
					prob2 += childBranch2->pX1X2(parentBase, childBase) * childProb2[site * charStates + childBase];
			}

			_pRiX1[site * charStates + parentBase] = prob1 * prob2;
		}
	}

	_pRiX1Version++;
	return _pRiX1;
}

// probability towards root
double* Branch::pSiX2(unsigned int numOfSites)
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

	if (_pSiX2 == NULL)
		_pSiX2 = new double[charStates * numOfSites];
	else if (_qVersion < _pSiX2Version) return _pSiX2;

	vector<unsigned int> siblingSeq;
	double* siblingProbRiX1;
	if (sibling->isLeaf())
		siblingSeq = sibling->getSequence();
	else
		siblingProbRiX1 = siblingBranch->pRiX1(numOfSites);

	vector<unsigned int> grandParentSeq;
	double* grandParentProbSiX2;
	if (grandParent->isLeaf())
		grandParentSeq = grandParent->getSequence();
	else
		grandParentProbSiX2 = grandParentBranch->pSiX2(numOfSites);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < numOfSites; site++)
		for (unsigned int parentBase = 0; parentBase < charStates; parentBase++)
		{
			double siblingProb;
			double grandParentProb;

			if (sibling->isLeaf())
			{
				siblingProb = siblingBranch->pX1X2(parentBase, siblingSeq[site]);
			} else
			{
				siblingProb = 0;
				for (unsigned int siblingBase = 0; siblingBase < charStates; siblingBase++)
				{
					siblingProb += siblingBranch->pX1X2(parentBase, siblingBase) * siblingProbRiX1[site * charStates + siblingBase];
				}
			}

			if (grandParent->isLeaf())
			{
				grandParentProb = grandParentBranch->getProb(grandParentSeq[site], parentBase);
			} else
			{
				grandParentProb = 0;
				for (unsigned int grandParentBase = 0; grandParentBase < charStates; grandParentBase++)
				{
					grandParentProb += grandParentBranch->pX1X2(grandParentBase, parentBase) * grandParentProbSiX2[site * charStates + grandParentBase];
				}
			}

			_pSiX2[charStates * site + parentBase] = siblingProb * grandParentProb;
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
