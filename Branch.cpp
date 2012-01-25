#include "Branch.h"
#include "globals.h"
#include <iostream>


Branch::Branch(int id, Node *n1, Node *n2)
{
	_id = id;
	_nodes.push_back(n1);
	_nodes.push_back(n2);
	_distance = -1.0;

	_q = new Matrix(CharStates);
	_q->setDiag(1.0/8.0);
	_q->setOffDiag(1.0/24.0);

	if (Verbose >= 2)
	{
		cout << "Branch #" << _id << "|";
		if (n1)
			cout << n1->getIdent();
		cout << "|";
		if (n2)
			cout << n2->getIdent();
		cout << endl;
		_q->print();
		cout << endl;
	}
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
