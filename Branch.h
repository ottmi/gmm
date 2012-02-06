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

	double computeLH(unsigned int numOfSites, vector<unsigned int> &invarSites);
	void updateQ();

	void print();
	string getIdent();

	void initializeQ();

private:
	vector <Node*> _nodes;
	double _distance;
	Matrix* _q;
	Matrix* _updatedQ;
	unsigned int _qVersion;
	double _beta;
	vector<double> _invar;
	vector<double> _pRiX1;
	unsigned int _pRiX1Version;
	vector<double> _pSiX2;
	unsigned int _pSiX2Version;
	vector<double> _siteProb;
	int _id;

	double computeValuesIntToInt(unsigned int numOfSites, vector<unsigned int> &invarSites);
	double computeValuesIntToLeaf(unsigned int numOfSites, vector<unsigned int> &invarSites);
	double computeValuesRootToInt(unsigned int numOfSites, vector<unsigned int> &invarSites);
	void updateQIntToInt(unsigned int numOfSites, vector<unsigned int> &invarSites);
	void updateQIntToLeaf(unsigned int numOfSites, vector<unsigned int> &invarSites);
	void updateQRootToInt(unsigned int numOfSites, vector<unsigned int> &invarSites);

	double pX1X2(unsigned int parent, unsigned int child);
	vector<double>& pRiX1(unsigned int numOfSites);
	vector<double>& pSiX2(unsigned int numOfSites);

	double getProb(unsigned int from, unsigned int to);
	double getMarginalProbRow(unsigned int row);
	double getMarginalProbCol(unsigned int col);
};

#endif /* BRANCH_H_ */
