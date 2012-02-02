#ifndef BRANCH_H_
#define BRANCH_H_

#include <vector>
#include "Node.h"
#include "Matrix.h"
using namespace std;

class Branch
{
public:
	Branch(int id, Node *n1, Node *n2);
	virtual
	~Branch();
	Node* getNeighbour(Node *node);
	void setDistance(double &distance);
	double getDistance() { return _distance; };

	double pX1X2(unsigned int parent, unsigned int child);
	vector<double> pRiX1(unsigned int numOfSites);
	vector<double> pSiX2(unsigned int numOfSites);

	double getProb(unsigned int from, unsigned int to);
	double getMarginalProbRow(unsigned int row);
	double getMarginalProbCol(unsigned int col);

	void print();

	void initializeQ();

private:
	vector <Node*> _nodes;
	double _distance;
	Matrix *_q;
	int _id;
};

#endif /* BRANCH_H_ */
