#ifndef TREE_H_
#define TREE_H_
#include <string>
#include <vector>
#include "globals.h"
#include "Node.h"
#include "Alignment.h"

using namespace std;

class Tree
{
	friend class Optimizer;

public:
	Tree();
	Tree(Tree const &tree);
	virtual
	~Tree();

	bool operator==(Tree const &tree) const;
	bool operator>(Tree const &tree) const;
	bool operator<(Tree const &tree) const;
	bool operator>=(Tree const &tree) const;
	bool operator<=(Tree const &tree) const;
	Tree& operator= (Tree const &tree);

	void readNewick(Alignment *alignment, string treeString, Options &options);
	double getLogLH();
	void computeLH();
	void updateModel(double qDelta, double betaDelta, bool thorough = false);

	void printBranches();
	void printNodes();
	void print();

	vector<Node*> const& getLeaves() const { return _leaves; };
	vector<Node*> const& getInternalNodes() const { return _internalNodes; };
	vector<Branch*> const& getBranches() const { return _branches; };

private:
	void copy(Tree const &tree);

	Node *_root;
	bool _unrooted;
	vector<Node*> _leaves;
	vector<Node*> _internalNodes;
	vector<Branch*> _branches;
	Alignment* _alignment;
	double _logLH;
};

#endif /* TREE_H_ */
