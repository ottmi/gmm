#ifndef TREE_H_
#define TREE_H_
#include <string>
#include <vector>
#include "Node.h"
#include "Alignment.h"

using namespace std;

class Tree
{
public:
	Tree(Alignment &alignment);
	virtual
	~Tree();

	void readNewick(string &treeString);
	void readNewickFromFile(string &fileName);
	void computeLH();
	void printBranches();
	void print();
private:
	Node *_root;
	vector<Node*> _leaves;
	vector<Node*> _internalNodes;
	vector<Branch*> _branches;
	Alignment _alignment;
	unsigned int _nodeCount;
	unsigned int _numOfSites;

};

#endif /* TREE_H_ */
