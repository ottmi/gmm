#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>
using namespace std;

class Branch;

class Node
{
public:
	Node(int id);
	virtual
	~Node();

	void setSequence(unsigned int* seq);
	void setLabel(string &label) { this->_label = label; };
	string& getLabel() { return _label; };
	int getId() { return _id; };
	string getIdent();
	Branch* getBranch(int id);
	vector<Branch*>& getBranches() { return _branches; };

	void removeBranch(Branch* b);
	void addBranch(Branch* b);
	void reroot(Branch *branch);

	string toString(Node *parent = NULL);
	Node* getParent();
	Node *getChild(int num);

	Branch* getNeighbourBranch(Node *neighbour);
	Branch* getNeighbourBranch(Node *exclude1, Node *exclude2);
	Node* getNeighbour(Node *neighbour);
	Node* getNeighbour(Node *exclude1, Node *exclude2);

	vector<Node*> getTraversal(Node *parent = NULL);
	void getDescendantBranches(Node *parent, vector<int> &branches);

	bool isLeaf() { return _isLeaf; };

	unsigned int getBase(unsigned int position);
	unsigned int* getSequence();

private:
	string _label;
	vector <Branch*> _branches;
	bool _isLeaf;
	int _id;
	unsigned int* _seq;
};

#endif /* NODE_H_ */