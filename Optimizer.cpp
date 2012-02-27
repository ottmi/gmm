#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "Optimizer.h"
#include "helper.h"
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
	Tree savedTree = tree;
	for (unsigned int i = 0; i < tree._branches.size(); i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			if (!tree._branches[i]->getNode(j)->isLeaf())
			{
				Tree reducedTree = tree;
				Branch *fromCandidate = reducedTree._branches[tree._branches[i]->getId()];

				vector<int> toCandidates;
				subtreePrune(fromCandidate, fromCandidate->getNode(j), toCandidates);
				//reducedTree.printBranches();

				for (unsigned int k = 0; k < toCandidates.size(); k++)
				{
					for (unsigned int l = 0; l < 2; l++)
					{
						Tree t = reducedTree;
						Branch *fromBranch = t._branches[fromCandidate->getId()];
						Branch *toBranch = t._branches[toCandidates[k]];

						if (!toBranch->getNode(l)->isLeaf())
						{
							subtreeRegraft(fromBranch, fromBranch->getNode(j), toBranch, toBranch->getNode(l), t._root);
							t.print();
							toBranch->computeLH(t._alignment->getPatternCount(), t._alignment->getInvarSites(), t._alignment->getInvarStart());
						}
					}
				}
			}
		}
	}
}

void Optimizer::subtreePrune(Branch *fromBranch, Node *fromParent, vector<int>& insertCandidates)
{
	if (fromBranch->getNode(0) != fromParent && fromBranch->getNode(1) != fromParent)
		throw("SPR: fromParent " + fromParent->getIdent() + " doesn't link to fromBranch " + str(fromBranch->getId()));
	if (fromParent->isLeaf()) throw("SPR: fromParent " + fromParent->getIdent() + " is a leaf node !");

	Node *fromChild = fromBranch->getNeighbour(fromParent);

	// Find the nearest internal node linking into fromParent, starting with the grandparent
	Node *fromGrandParent = fromParent->getParent();
	// If fromChild is closer to the root than fromParent, fromParent->getParent() will link to fromChild
	// We don't want this, so we pick a child node instead
	// Also, it has to be an internal node
	if (fromGrandParent->isLeaf() || fromGrandParent == fromChild) fromGrandParent = fromParent->getChild(1);
	// Make sure it's an internal node, otherwise choose the other link.
	// If both children are internal nodes, it doesn't matter which one we choose
	if (fromGrandParent->isLeaf() || fromGrandParent == fromChild) fromGrandParent = fromParent->getChild(2);
	// If both children are leaf nodes, we can't perform an SPR move on that branch
	if (fromGrandParent->isLeaf()) return;

	Branch *fromParentBranch = fromParent->getNeighbourBranch(fromGrandParent);

	if (verbose >= 2) cout << "subtreePrune: fromBranch=" << fromBranch->getId() << " fromParent=" << fromParent->getIdent() << " fromParentBranch=" << fromParentBranch->getId() << " fromGrandParent=" << fromGrandParent->getIdent() << endl;

	// Identify the branch that leads to our sibling
	Branch *fromSiblingBranch = fromParent->getNeighbourBranch(fromChild, fromGrandParent);
	Node *fromSibling = fromSiblingBranch->getNeighbour(fromParent);

	// Unlink our sibling including its branch from our parent
	fromSiblingBranch->unlinkNode(fromParent);

	// Unlink our parent including its branch from its parent.
	fromParentBranch->unlinkNode(fromGrandParent);

	// Link our sibling there instead
	fromSiblingBranch->linkNode(fromGrandParent);

	// Let's enumerate all possible insertion branches
	fromGrandParent->getDescendantBranches(fromSibling, insertCandidates);
	fromSibling->getDescendantBranches(fromGrandParent, insertCandidates);
}

void Optimizer::subtreeRegraft(Branch *fromBranch, Node *fromParent, Branch *toBranch, Node *toParent, Node *root)
{
	if (verbose >= 2)
		cout << "subtreeRegraft: fromBranch=" << fromBranch->getId() << " fromParent=" << fromParent->getIdent() << " toBranch=" << toBranch->getId()
				<< " toParent=" << toParent->getIdent() << endl;
	// fromParent is an internal node with only 2 branches. Find the branch that leads away from the parent
	Branch *fromParentBranch = fromParent->getBranch(0);
	if (fromParentBranch == fromBranch) fromParentBranch = fromParent->getBranch(1);

	// Unlink the insertion branch from its parent
	toBranch->unlinkNode(toParent);

	// Link the branch leading away from our parent there instead
	fromParentBranch->linkNode(toParent);

	// Link the insertion branch as sibling to our parent node
	toBranch->linkNode(fromParent);

	// Make sure that all Node::_branches[0] and Branch::_nodes[0] point towards the root
	root->reroot(NULL);

	// Reset all _pRiX1 and _pSiX2 vectors towards the root
	fromParentBranch->resetVectors();
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
