#include "Tree.h"
#include "Branch.h"
#include "helper.h"
#include "globals.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <sstream>
#include <map>

Tree::Tree(Alignment* alignment)
{
	_alignment = alignment;
}

Tree::Tree(Tree const &tree)
{
	copy(tree);
}

Tree::~Tree()
{
	for (unsigned int i = 0; i < _leaves.size(); i++)
		delete _leaves[i];

	for (unsigned int i = 0; i < _internalNodes.size(); i++)
		delete _internalNodes[i];

	for (unsigned int i = 0; i < _branches.size(); i++)
		delete _branches[i];
}

Tree& Tree::operator=(Tree const &tree)
{
	if (this != &tree) copy(tree);

	return *this;
}

void Tree::copy(Tree const &tree)
{
	if (verbose >= 5) cout << "Tree::copy start" << endl;
	vector<Node*> nodes(tree._internalNodes.size() + tree._leaves.size(), NULL);

	_alignment = tree._alignment;
	_internalNodes.resize(tree._internalNodes.size(), NULL);
	for (unsigned int i = 0; i < tree._internalNodes.size(); i++)
	{
		Node *node = new Node(tree._internalNodes[i]->getId());
		_internalNodes[i] = node;
		_internalNodes[i]->setLabel(tree._internalNodes[i]->getLabel());
		nodes[node->getId()] = node;
	}

	_leaves.resize(tree._leaves.size(), NULL);
	for (unsigned int i = 0; i < tree._leaves.size(); i++)
	{
		Node *node = new Node(tree._leaves[i]->getId());
		node->setSequence(tree._leaves[i]->getSequence());
		node->setLabel(tree._leaves[i]->getLabel());
		_leaves[i] = node;
		nodes[node->getId()] = node;
	}

	_root = nodes[tree._root->getId()];

	_branches.resize(tree._branches.size(), NULL);
	for (unsigned int i = 0; i < tree._branches.size(); i++)
	{
		int id = tree._branches[i]->getId();
		unsigned int numOfNodes = tree._branches[i]->getNumOfNodes();
		Node *node0 = NULL;
		Node *node1 = NULL;
		if (numOfNodes >= 1) node0 = nodes[tree._branches[i]->getNode(0)->getId()];
		if (numOfNodes >= 2) node1 = nodes[tree._branches[i]->getNode(1)->getId()];
		Branch *branch = new Branch(tree._branches[i], node0, node1, _alignment->getNumOfUniqueSites());
		_branches[id] = branch;
	}
	if (verbose >= 5) cout << "Tree::copy finish" << endl;
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
	unsigned int nodeCount = 0;
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
			if (currentNode)
			{
				Node *newNode = new Node(nodeCount++);
				Branch *branch = new Branch(nodeCount - 2, currentNode, newNode);
				currentNode = newNode;
				_branches.push_back(branch);
			} else // there was no current node, so this is the first node, hence no branch leading into it
			{
				currentNode = new Node(nodeCount++);
			}
			_internalNodes.push_back(currentNode);
			nextCouldBeLeaf = true;
			i++;
		} else if (treeString[i] == ')' || treeString[i] == ',') // node ends, could be internal or leaf
		{
			if (nextCouldBeLeaf)
			{
				if (verbose >= 5) cout << "  Leaf #" << _leaves.size() << " (" << label << ") " << distance << endl;

				Node *leaf;
				if (currentNode)
				{
					leaf = new Node(nodeCount++);
					Branch *branch = new Branch(nodeCount - 2, currentNode, leaf);
					_branches.push_back(branch);
				} else // same as above, no node and branch yet
				{
					leaf = new Node(nodeCount++);
				}
				leaf->setLabel(label);
				label.clear();
				_leaves.push_back(leaf);
			} else
			{
				prevInternalNode->setLabel(label);
				label.clear();
//				Branch *branch = prevInternalNode->getBranch(0);
//				branch->setDistance(distance);
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
			int j = treeString.find_first_of(",():;[", i + 1);
			stringstream ss(treeString.substr(i + 1, j - i - 1));
			ss >> distance;
			i = j;
		}

		if (isalpha(treeString[i])) // this is a label
		{
			int j = treeString.find_first_of(",():;[", i + 1);
			label = treeString.substr(i, j - i);
			i = j;
		}

		if (treeString[i] == '[') // this is a comment
		{
			size_t j = treeString.find_first_of(']', i + 1);
			if (j == string::npos)
				throw("Error paarsing Newick tree: a square bracket opened at position " + str(i) + " but no closing bracket could be found.");
			else
			{
				cout << "Comment: " << treeString.substr(i+1, j-i-1) << endl;
				i = j+1;
			}
		}
	}

	_root = prevInternalNode;
	_root->setLabel(label);
	label.clear();
	if (_root->getBranches().size() == 1)
	{
		cout << "The root at " << _root->getLabel() << " is a leaf" << endl;
		_leaves.push_back(_root);
		if (_internalNodes.front() == _root) _internalNodes.erase(_internalNodes.begin());
	}

	for (unsigned int i = 0; i < _leaves.size(); i++)
	{
		label = _leaves[i]->getLabel();
		int a = _alignment->find(label);
		if (a < 0)
			throw("The alignment contains no sequence \"" + label + "\"");
		else
		{
			unsigned int* seq = _alignment->getNumericalSeq(a);
			_leaves[i]->setSequence(seq);
		}
	}

	cout << "The tree has " << _internalNodes.size() << " internal nodes, " << _leaves.size() << " leaves and " << _branches.size() << " branches." << endl;
	if (_alignment->getNumOfSequences() != _leaves.size())
	{
		stringstream ss;
		ss << "The number of leaves (" << _leaves.size() << ") and the number of sequences in the alignment (" << _alignment->getNumOfSequences()
				<< ") do not match";
		throw(ss.str());
	}
}

void Tree::computeLH()
{
	for (unsigned int i = 0; i < _branches.size(); i++)
	{
		_branches[i]->computeLH(_alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart());
	}
}

bool Tree::updateModel(double qDelta, double betaDelta)
{
	for (unsigned int i = 0; i < _branches.size(); i++)
		_branches[i]->computeUpdatedQ(_alignment->getNumOfSites(), _alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart());

	unsigned int updates = 0;
	for (unsigned int i = 0; i < _branches.size(); i++)
	{
		if (_branches[i]->updateQ(qDelta)) updates++;

		if (_branches[i]->updateParameters(_alignment->getNumOfSites(), _alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart(),
				betaDelta)) updates++;
	}

	return (updates > 0);
}

void Tree::printNodes()
{
	cout << "Internal Nodes: " << endl;
	for (unsigned int i = 0; i < _internalNodes.size(); i++)
	{
		cout << _internalNodes[i]->getIdent() << " parent=" << _internalNodes[i]->getParent()->getIdent() << " child1="
				<< _internalNodes[i]->getChild(1)->getIdent() << " child2=" << _internalNodes[i]->getChild(2)->getIdent() << endl;
	}

	cout << "Leaves: " << endl;
	for (unsigned int i = 0; i < _leaves.size(); i++)
	{
		cout << _leaves[i]->getIdent() << " parent=" << _leaves[i]->getParent()->getIdent() << endl;
	}
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
	cout << _root->toString() << ";" << endl;
}
