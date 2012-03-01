#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <set>
#include <ctime>
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

#define MAX_BEST_TREES 5
void Optimizer::rearrange(Tree &tree, Options &options)
{
	double currentCutOff = 0.01;
	bool improved = true;

	tree.updateModel(currentCutOff, currentCutOff);
	tree.computeLH();
	bestTree = tree;
	set<Tree> bestTrees;
	bestTrees.insert(tree);

	unsigned int round = 0;
	while (improved)
	{
		cout << endl << "Starting round #" << round << endl;
		improved = false;

		if (round % 2 == 0)
			optimizeNNI(bestTree, currentCutOff, bestTrees);
		else
			optimizeSPR(bestTree, currentCutOff, bestTrees);

		for (set<Tree>::reverse_iterator it = bestTrees.rbegin(); it != bestTrees.rend(); it++)
		{
			Tree t = *it;
			t.updateModel(options.cutOff, options.cutOff);
			t.computeLH();
			if (t > bestTree)
			{
				if (t.getLogLH() - bestTree.getLogLH() > options.cutOff) improved = true;
				bestTree = t;
			}
		}
		if (currentCutOff > options.cutOff) currentCutOff /= 2;
		if (currentCutOff < options.cutOff) currentCutOff = options.cutOff;

		round++;
		cout << "\r  Best tree: " << fixed << setprecision(6) << bestTree.getLogLH() << "                                                       " << endl;
	}
	bestTree.updateModel(options.cutOff, options.cutOff);
	bestTree.computeLH();

	cout << endl << "Best tree found: " << fixed << setprecision(10) << bestTree.getLogLH() << endl;
	tree = bestTree;
}

void Optimizer::optimizeSPR(Tree &tree, double cutOff, set<Tree> &bestTrees)
{
	cout << "Performing SPR moves" << endl;
	time_t t = time(NULL);
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
							assessTree(t, cutOff, bestTrees);
						}
					}
				}
			}
			time_t elapsed = time(NULL) -t ;
			cout << "\r" << (i * 100) / tree._branches.size() << "%  " << "\tTime elapsed: " << printTime(elapsed) << flush;
			if (i - elapsed > 0)
			{
				time_t eta = (elapsed * tree._branches.size()) / i - elapsed;
				cout << "\tETA: " << printTime(eta) << "  " << "\t[" << fixed << setprecision(4) << bestTrees.begin()->_logLH << ","	<< bestTrees.rbegin()->_logLH << "]" << flush;
			}
		}
	}
}

void Optimizer::optimizeNNI(Tree &tree, double cutOff, set<Tree> &bestTrees)
{
	cout << "Performing NNI moves" << endl;
	time_t t = time(NULL);
	for (unsigned int i = 0; i < tree._branches.size(); i++)
	{
		cout << "\r" << (i * 100) / tree._branches.size() << "%   " << flush;
		if (!tree._branches[i]->getNode(0)->isLeaf() && !tree._branches[i]->getNode(1)->isLeaf())
		{
			for (unsigned int j = 0; j < 2; j++)
			{
				Tree t = tree;
				Branch *b = t._branches[tree._branches[i]->getId()];
				NNI(b, j);
				t._root->reroot(NULL);
				assessTree(t, cutOff, bestTrees);
			}
		}
		time_t elapsed = time(NULL) -t ;
		cout << "\r" << (i * 100) / tree._branches.size() << "%  " << "\tTime elapsed: " << printTime(elapsed) << flush;
		if (i - elapsed > 0)
		{
			time_t eta = (elapsed * tree._branches.size()) / i - elapsed;
			cout << "\tETA: " << printTime(eta) << "  " << "\t[" << fixed << setprecision(4) << bestTrees.begin()->_logLH << ","	<< bestTrees.rbegin()->_logLH << "]" << flush;
		}
	}
}

void Optimizer::assessTree(Tree &tree, double cutOff, set<Tree> &bestTrees)
{
	for (unsigned int m = 0; m < tree._branches.size(); m++)
	{
		tree._branches[m]->resetVectors();
		tree._branches[m]->resetQ();
	}

	tree.updateModel(cutOff, cutOff);
	tree.computeLH();
	if (tree > *bestTrees.begin())
	{
		bestTrees.insert(tree);
		if (bestTrees.size() > MAX_BEST_TREES) bestTrees.erase(bestTrees.begin());
		if (verbose >= 2)
		{
			cout << "\rlogLH: " << tree.getLogLH() << "                                                       " << endl;
			tree.print();
			cout << endl;
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
	// Also, it has to be an internal node and we need to make sure not to get back to the original parent
	if (fromGrandParent->isLeaf() || fromGrandParent == fromChild || fromGrandParent == fromParent) fromGrandParent = fromParent->getChild(1);
	// Make sure it's an internal node, otherwise choose the other link.
	// If both children are internal nodes, it doesn't matter which one we choose
	if (fromGrandParent->isLeaf() || fromGrandParent == fromChild || fromGrandParent == fromParent) fromGrandParent = fromParent->getChild(2);
	// If both children are leaf nodes, we can't perform an SPR move on that branch
	if (fromGrandParent->isLeaf() || fromGrandParent == fromChild || fromGrandParent == fromParent) return;

	Branch *fromParentBranch = fromParent->getNeighbourBranch(fromGrandParent);

	if (verbose >= 3)
		cout << "subtreePrune: fromBranch=" << fromBranch->getId() << " fromParent=" << fromParent->getIdent() << " fromParentBranch=" << fromParentBranch->getId()
				<< " fromGrandParent=" << fromGrandParent->getIdent() << endl;

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
	if (verbose >= 3)
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
