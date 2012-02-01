#include "Tree.h"
#include "Branch.h"
#include "helper.h"
#include "globals.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <sstream>

Tree::Tree(Alignment &alignment)
{
	_nodeCount = 0;
	_alignment = alignment;
	_numOfSites = _alignment.getCols();
}

Tree::~Tree()
{
	// TODO Auto-generated destructor stub
}

void Tree::readNewick(string &tree)
{
	string treeString;
	if (tree.find(';') == string::npos)
	{
		ifstream fileReader;
		fileReader.open(tree.c_str());
		if (!fileReader.is_open()) throw("Cannot open file " + tree);

		safeGetline(fileReader, treeString);
		fileReader.close();
	} else
		treeString = tree;

	cout << "Tree: " << treeString << endl;

	unsigned int i = 0;
	Node *prevInternalNode, *currentNode;
	prevInternalNode = currentNode = NULL;
	string label;
	double distance = -1.0;
	bool nextCouldBeLeaf = false;

	while (i < treeString.length() && treeString[i] != ';')
	{
		if (treeString[i] == '(') // internal node starts
		{
			if (verbose >= 5) cout << "Internal Node #" << _internalNodes.size() << " starts" << endl;
			currentNode = new Node(currentNode, _nodeCount++);
			_internalNodes.push_back(currentNode);
			Branch *branch = currentNode->getBranch(0);
			if (branch) _branches.push_back(branch);
			nextCouldBeLeaf = true;
			i++;
		} else if (treeString[i] == ')' || treeString[i] == ',') // node ends, could be internal or leaf
		{
			if (nextCouldBeLeaf)
			{
				if (verbose >= 5) cout << "  Leaf #" << _leaves.size() << " (" << label << ") " << distance << endl;

				Node *leaf = new Node(currentNode, _nodeCount++);
				leaf->setLabel(label);
				Branch *branch = leaf->getBranch(0);
				branch->setDistance(distance);
				_branches.push_back(branch);
				_leaves.push_back(leaf);
			} else
			{
				prevInternalNode->setLabel(label);
				Branch *branch = prevInternalNode->getBranch(0);
				branch->setDistance(distance);
			}

			if (treeString[i] == ')') // internal node
			{
				if (verbose >= 5) cout << "Internal Node #" << currentNode->getId() << " ends " << endl;
				prevInternalNode = currentNode;
				currentNode = currentNode->getParent();
				nextCouldBeLeaf = false;
			} else // node will follow, could be a leaf
			{
				nextCouldBeLeaf = true;
			}
			i++;
		} else if (treeString[i] == ':') // distance for last node
		{
			int j = treeString.find_first_of(",():;", i + 1);
			stringstream ss(treeString.substr(i + 1, j - i - 1));
			ss >> distance;
			i = j;
		}

		if (isalpha(treeString[i])) // this is a label
		{
			int j = treeString.find_first_of(",():;", i + 1);
			label = treeString.substr(i, j - i);
			i = j;
		}
	}

	_root = prevInternalNode;
	_root->setLabel(label);
	if (_root->getBranches().size() == 1)
	{
		cout << "The root at " << _root->getLabel() << " is a leaf" << endl;
		_leaves.push_back(_root);
		if (_internalNodes.front() == _root) _internalNodes.erase(_internalNodes.begin());
	}

	for (unsigned int i = 0; i < _leaves.size(); i++)
	{
		label = _leaves[i]->getLabel();
		int a = _alignment.find(label);
		if (a < 0)
			throw("The alignment contains no sequence \"" + label + "\"");
		else
		{
			vector<unsigned int> seq = _alignment.getNumericalSeq(a);
			_leaves[i]->setSequence(seq);
		}
	}

	cout << "The tree has " << _internalNodes.size() << " internal nodes, " << _leaves.size() << " leaves and " << _branches.size() << " branches." << endl;
	if (_alignment.getRows() != (int) _leaves.size())
	{
		stringstream ss;
		ss << "The number of leaves (" << _leaves.size() << ") and the number of sequences in the alignment (" << _alignment.getRows() << ") do not match";
		throw(ss.str());
	}
}

void Tree::computeLH()
{
	/*
	 vector<Node*> traversal = _leaves[0]->getTraversal();
	 for (unsigned int i=0; i<traversal.size(); i++)
	 cout << traversal[i]->getIdent() << endl;
	 */
//	_internalNodes[1]->computeValuesIntToInt(_numOfSites);
	cout << endl << "Root" << endl;
	_root->computeValuesRootToInt(_numOfSites);

	cout << endl << "Internal Nodes" << endl;
	for (unsigned int i = 0; i < _internalNodes.size(); i++)
		if (_internalNodes[i] != _root) _internalNodes[i]->computeValuesIntToInt(_numOfSites);

	cout << endl << "Leaves" << endl;
	for (unsigned int i = 0; i < _leaves.size(); i++)
		if (_leaves[i]->getParent() != _root) _leaves[i]->computeValuesIntToLeaf(_numOfSites);
}

void Tree::printBranches()
{
	for (unsigned int i = 0; i < _branches.size(); i++)
	{
		_branches[i]->print();
	}
}

void Tree::print()
{
	for (unsigned int i = 0; i < _internalNodes.size(); i++)
	{
		cout << _internalNodes[i]->toString() << ";" << endl;
	}
}
