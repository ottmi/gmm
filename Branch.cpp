#include "Branch.h"
#include "globals.h"
#include <iostream>


Branch::Branch(int id, Node *n1, Node *n2)
{
	_id = id;
	_nodes.push_back(n1);
	_nodes.push_back(n2);
	_distance = -1.0;

	_q = new Matrix(charStates);
	_q->setDiag(1.0/8.0);
	_q->setOffDiag(1.0/24.0);
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
	cout << "{" << _nodes[0]->getIdent() << "}<---(" << _distance << ")--->{" << _nodes[1]->getIdent() << "}" << endl;
	_q->print();
}


/* This method computes the conditional probability P(x2|x1) */
double Branch::pX1X2(unsigned int parent, unsigned int child)
{
	double marginalProb = _q->getRowSum(parent);
	double condProb = _q->getEntry(parent, child) / marginalProb;

	return condProb;
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
