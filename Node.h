#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>
using namespace std;

class Branch;

class Node
{
public:
	Node(Node *parent, int id, int seq);
	virtual
	~Node();

	void setLabel(string &label) { this->_label = label; label.clear(); };
	string& getLabel() { return _label; };
	int getId() { return _id; };
	string getIdent();
	Branch* getBranch(int id);
	string toString(Node *parent = NULL);
	Node* getParent();
	vector<Node*> getTraversal(Node *parent = NULL);

private:
	string _label;
	vector <Branch*> _branches;
	vector <double> probs;
	double beta;
	bool _isLeaf;
	int _id;
	int _seq;
};

#endif /* NODE_H_ */
