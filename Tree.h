/*
 * Tree.h
 *
 *  Created on: 15/07/2011
 *      Author: ott029
 */

#ifndef TREE_H_
#define TREE_H_
#include <string>
#include <vector>
#include "Node.h"

using namespace std;

class Tree
{
public:
	Tree();
	virtual
	~Tree();

	void readNewick(string fileName);
	void print();
private:
	Node *root;
	vector<Node*> leaves;
	vector<Node*> internalNodes;

};

#endif /* TREE_H_ */
