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

	double computeLH(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	void computeUpdatedQ(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	bool updateQ(double qDelta);
	bool updateParameters(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart, double betaDelta);

	void print();
	string getIdent();
	double getLength();

	void initializeQ();

private:
	vector <Node*> _nodes;
	Matrix* _q;
	Matrix* _updatedQ;
	unsigned int _qVersion;
	double _beta;
	vector<double> _invar;
	double* _pRiX1;
	unsigned int _pRiX1Version;
	double* _pSiX2;
	unsigned int _pSiX2Version;
	double* _siteProb;
	int _id;

	double computeValuesIntToInt(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	double computeValuesIntToLeaf(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	double computeValuesRootToInt(vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	void updateQIntToInt(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	void updateQIntToLeaf(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);
	void updateQRootToInt(unsigned int numOfSites, vector<unsigned int> &patternCount, vector<unsigned int> &invarSites, unsigned int invarStart);

	double pX1X2(unsigned int parent, unsigned int child);
	double* pRiX1(unsigned int numOfSites);
	double* pSiX2(unsigned int numOfSites);

	double getProb(unsigned int from, unsigned int to);
	double getMarginalProbRow(unsigned int row);
	double getMarginalProbCol(unsigned int col);
};

#endif /* BRANCH_H_ */
