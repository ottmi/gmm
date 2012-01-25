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
	Tree();
	virtual
	~Tree();

	void readNewick(string &treeString, Alignment &alignment);
	void readNewickFromFile(string &fileName, Alignment &alignment);
	void computeLH();
	void print();
private:
	Node *_root;
	vector<Node*> _leaves;
	vector<Node*> _internalNodes;
	vector<Branch*> _branches;
	unsigned int _nodeCount;

};

#endif /* TREE_H_ */
