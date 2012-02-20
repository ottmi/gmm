#include <iostream>
#include "Optimizer.h"
#include "globals.h"

Optimizer::Optimizer()
{
	// TODO Auto-generated constructor stub

}

Optimizer::~Optimizer()
{
	// TODO Auto-generated destructor stub
}

void Optimizer::rearrange(Tree &tree)
{
	vector<Branch*> branches = tree.getBranches();
	vector<Node*> internalNodes = tree.getInternalNodes();

	SPR(branches[1], internalNodes[0], branches[5], internalNodes[2]);
}

void Optimizer::NNI(Branch* branch, int swap)
{
	// swap indicates which branch of the left node will be swapped, i.e. 0 for the 1st branch, 1 for the 2nd

	Node *leftNode = branch->getNode(0);
	Node *rightNode = branch->getNode(1);
	if (leftNode->isLeaf())
		throw("NNI: Node " + leftNode->getIdent() + " is a leaf node");
	else if (rightNode->isLeaf()) throw("NNI: Node " + rightNode->getIdent() + " is a leaf node");

	Branch *leftBranch; // this is the outgoing branch from the left node we are going to swap
	for (unsigned int i = 0; i < 3; i++)
	{
		leftBranch = leftNode->getBranch(i);
		if (leftBranch != branch && !swap--) break;
	}

	Branch *rightBranch; // this is the outgoing branch from the right node we are going to swap
	for (unsigned int i = 0; i < 3; i++)
	{
		rightBranch = rightNode->getBranch(i);
		if (rightBranch != branch) break;
	}

	//unlink the branch from the left node
	leftBranch->unlinkNode(leftNode);

	//unlink the branch from the right node
	rightBranch->unlinkNode(rightNode);

	//link the right node's branch to the left node
	rightBranch->linkNode(leftNode);

	//link the left node's branch to the right node
	leftBranch->linkNode(rightNode);
}


void Optimizer::SPR(Branch *fromBranch, Node *fromParent, Branch *toBranch, Node *toParent)
{
	if (verbose >= 2)
		cout << "SPR: fromBranch=" << fromBranch->getId() << " fromParent=" << fromParent->getIdent() << " toBranch=" << toBranch->getId() << " toParent=" << toParent->getIdent() << endl;

	Node *fromChild = fromBranch->getNeighbour(fromParent);

	if (fromParent->isLeaf()) throw("SPR: fromParent is a leaf node !");

	/* Find the nearest internal node towards the root and the branch leading to it
	 * This will usually be our grandparent, but it has to be an internal node
	 * So if it's the actual root (which is a leaf) we choose the next internal node instead */
	Node *fromGrandParent = fromParent->getParent();
	if (fromGrandParent->isLeaf())
		fromGrandParent = fromParent->getNeighbour(fromChild, fromGrandParent);
	Branch *fromParentBranch = fromParent->getNeighbourBranch(fromGrandParent);

	// Identify the branch that leads to our sibling
	Branch *fromSiblingBranch = fromParent->getNeighbourBranch(fromChild, fromGrandParent);

	// Unlink our sibling including its branch from our parent
	fromSiblingBranch->unlinkNode(fromParent);

	// Unlink our parent including its branch from its parent.
	fromParentBranch->unlinkNode(fromGrandParent);

	// Link our sibling there instead
	fromSiblingBranch->linkNode(fromGrandParent);

	// Unlink the insertion branch from its parent
	toBranch->unlinkNode(toParent);

	// Link the branch leading into our parent there instead
	fromParentBranch->linkNode(toParent);

	// Link the insertion branch as sibling to our parent node
	toBranch->linkNode(fromParent);
}


