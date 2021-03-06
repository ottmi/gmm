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
	bool operator!=(Tree const &tree) const;
	bool operator>(Tree const &tree) const;
	bool operator<(Tree const &tree) const;
	Tree& operator= (Tree const &tree);

    void removeNode(Node *node);
    void removeBranch(Branch *branch);
    Branch* findBranch(int id);
  
	void readNewick(Alignment *alignment, string treeString, Options &options);
	double getLogLH() const;
	void computeLH();
	void updateModel(double qDelta, double betaDelta, bool thorough = false);
	void clearModel() ;

	void printBranches();
	void printNodes();
    string toString(const bool topologyOnly=true) const;
	void print();

	vector<Node*> const& getLeaves() const { return _leaves; };
	vector<Node*> const& getInternalNodes() const { return _internalNodes; };
	vector<Branch*> const& getBranches() const { return _branches; };
	
	void setComment(const string &comment) { this->_comment = comment; };
	string getComment() const { return _comment; };

private:
  void clean();
	void copy(Tree const &tree);

	Node *_root;
	bool _unrooted;
	vector<Node*> _leaves;
	vector<Node*> _internalNodes;
	vector<Branch*> _branches;
	Alignment* _alignment;
	double _logLH;
	string _comment;
};

#endif /* TREE_H_ */
