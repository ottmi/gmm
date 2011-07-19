#ifndef NODE_H_
#define NODE_H_

#include <string>
using namespace std;

class Node
{
public:
	Node(Node *parent, int id, bool isLeaf);
	virtual
	~Node();

	void setLabel(string label) { this->label = label; };
	string& getLabel() { return label; };
	string toString();

	Node *left;
	Node *right;
	Node *parent;
	bool isLeaf;
	int id;

private:
	string label;
};

#endif /* NODE_H_ */
