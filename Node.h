#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>
using namespace std;

class Branch;

class Node
{
public:
	Node(Node *parent, int id);
	virtual
	~Node();

	void setSequence(vector<unsigned int> &seq);
	void setLabel(string &label) { this->_label = label; label.clear(); };
	string& getLabel() { return _label; };
	int getId() { return _id; };
	string getIdent();
	Branch* getBranch(int id);
	vector<Branch*>& getBranches() { return _branches; };
	string toString(Node *parent = NULL);
	Node* getParent();
	Node *getChild(int num);
	vector<Node*> getTraversal(Node *parent = NULL);

	bool isLeaf() { return _isLeaf; };

	unsigned int getBase(unsigned int position);
	vector<unsigned int>& getSequence();

private:
	string _label;
	vector <Branch*> _branches;
	vector <double> probs;
	bool _isLeaf;
	int _id;
	vector<unsigned int> _seq;
};

#endif /* NODE_H_ */
