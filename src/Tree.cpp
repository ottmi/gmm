#include "Tree.h"
#include "Branch.h"
#include "helper.h"
#include "globals.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>

Tree::Tree()
{
	_alignment = NULL;
	_root = NULL;
	_logLH = 0;
}

Tree::Tree(Tree const &tree)
{
	copy(tree);
}

Tree::~Tree()
{
    clean();
}

bool Tree::operator==(Tree const &tree) const {
	return toString(true) == tree.toString(true);
}

bool Tree::operator!=(Tree const &tree) const {
	return !(*this == tree);
}

bool Tree::operator>(Tree const &tree) const {
	if (_logLH == 0)
      throw(string("Tree::operator>() LHS logLH hasn't been computed yet"));
	if (tree._logLH == 0)
      throw(string("Tree::operator>() RHS logLH hasn't been computed yet"));
	return (_logLH > tree._logLH);
}

bool Tree::operator<(Tree const &tree) const {
	if (_logLH == 0)
      throw(string("Tree::operator<() LHS logLH hasn't been computed yet"));
	if (tree._logLH == 0)
      throw(string("Tree::operator<() RHS logLH hasn't been computed yet"));
	return (_logLH < tree._logLH);
}

Tree& Tree::operator=(Tree const &tree)
{
  if (this != &tree) {
      clean();
      copy(tree);
    }

	return *this;
}

void Tree::clean() {
  _root = NULL;
  _unrooted = false;

  for (unsigned int i = 0; i < _leaves.size(); i++)
    delete _leaves[i];
  _leaves.clear();

  for (unsigned int i = 0; i < _internalNodes.size(); i++)
    delete _internalNodes[i];
  _internalNodes.clear();

  for (unsigned int i = 0; i < _branches.size(); i++)
    delete _branches[i];
  _branches.clear();

  _alignment = NULL;
  _logLH = 0;
}

void Tree::copy(Tree const &tree)
{
	if (verbose >= 5) cout << "Tree::copy start" << endl;
    map<int,Node*> nodes;
  
	_alignment = tree._alignment;
	for (unsigned int i = 0; i < tree._internalNodes.size(); i++)
	{
		Node *node = new Node(tree._internalNodes[i]->getId());
		node->setLabel(tree._internalNodes[i]->getLabel());
        _internalNodes.push_back(node);
		nodes.insert(pair<int,Node*>(tree._internalNodes[i]->getId(), node));
	}

	for (unsigned int i = 0; i < tree._leaves.size(); i++)
	{
		Node *node = new Node(tree._leaves[i]->getId());
		node->setSequence(tree._leaves[i]->getSequence());
		node->setLabel(tree._leaves[i]->getLabel());
        node->setLeaf();
		_leaves.push_back(node);
        nodes.insert(pair<int,Node*>(tree._leaves[i]->getId(), node));
	}

	unsigned int links = 0;
	for (unsigned int i = 0; i < tree._branches.size(); i++)
	{
		unsigned int numOfNodes = tree._branches[i]->getNumOfNodes();
		Node *node0 = NULL;
		Node *node1 = NULL;
		if (numOfNodes >= 1)
		{
			node0 = nodes[tree._branches[i]->getNode(0)->getId()];
			links++;
		}
		if (numOfNodes >= 2)
		{
			node1 = nodes[tree._branches[i]->getNode(1)->getId()];
			links++;
		}
		Branch *branch = new Branch(tree._branches[i], node0, node1, _alignment->getNumOfUniqueSites());
		_branches.push_back(branch);
	}
	_root = nodes[tree._root->getId()];
	if (links == _branches.size() * 2) _root->reroot(NULL);

    _unrooted = tree._unrooted;
	_logLH = tree._logLH;

	if (verbose >= 5) cout << "Tree::copy finish" << endl;
}

void Tree::removeNode(Node *node) {
  if (node->isLeaf()) {
    vector<Node*>::iterator it = _leaves.begin();
    while (it != _leaves.end() && *it != node) {
      it++;
    }
    if (*it == node) {
      _leaves.erase(it);
      delete node;
    } else {
      cout << "Tried to delete leaf node " << node->getIdent() << " but could not find it in the list of leaves." << endl;
    }
  } else {
    vector<Node*>::iterator it = _internalNodes.begin();
    while (it != _internalNodes.end() && *it != node) {
      it++;
    }
    if (*it == node) {
      _internalNodes.erase(it);
      delete node;
    } else {
      cout << "Tried to delete internal node " << node->getIdent() << " but could not find it in the list of internal nodes." << endl;
    }
  }
}

void Tree::removeBranch(Branch *branch) {
  vector<Branch*>::iterator it = _branches.begin();
  while (it != _branches.end() && *it != branch) {
    it++;
  }
  if (*it == branch) {
    _branches.erase(it);
    delete branch;
  } else {
    cout << "Tried to delete branch " << branch->getIdent() << " but could not find it in the list of branches." << endl;
  }
}

Branch* Tree::findBranch(int id) {
  vector<Branch*>::iterator it = _branches.begin();
  while (it != _branches.end() && (*it)->getId() != id) {
    it++;
  }
  if ((*it)->getId() == id) {
    return *it;
  } else {
    cout << "Could not find branch " << id << " in the list of branches." << endl;
    return NULL;
  }
}

void Tree::readNewick(Alignment *alignment, string treeString, Options &options)
{
	cout << "Tree: " << treeString << endl;

	_alignment = alignment;
	unsigned int i = 0;
	unsigned int nodeCount = 0;
	Node *prevInternalNode, *currentNode;
	prevInternalNode = currentNode = NULL;
	string label;
	string comment;
	double distance = -1.0;
	bool nextCouldBeLeaf = false;

	while (i < treeString.length() && treeString[i] != ';')
	{
        if (treeString[i] == '(') { // internal node starts
          Node *newNode = new Node(nodeCount++);
          if (verbose >= 5)
            cout << "Internal node " << newNode->getIdent() << " " << _internalNodes.size() + 1 << "/" << nodeCount << ") starts" << endl;
          if (currentNode) { // current node exists, so let's link it to the new node
                Branch *branch = new Branch(nodeCount - 2, currentNode, newNode);
                _branches.push_back(branch);
            }
            currentNode = newNode;
			_internalNodes.push_back(currentNode);
			nextCouldBeLeaf = true;
			i++;
        } else if (treeString[i] == ')' || treeString[i] == ',') { // node ends, could be internal or leaf
			if (nextCouldBeLeaf) {
                Node *leaf = new Node(nodeCount++);
                leaf->setLabel(label);
                label.clear();
				leaf->setComment(comment);
				comment.clear();
                leaf->setLeaf();
                if (verbose >= 5)
                    cout << "Leaf node " << leaf->getIdent() << " " << _leaves.size() + 1 << "/" << nodeCount << ") starts" << endl;
                if (currentNode) { // current node exists, so let's link it to the new node
                    Branch *branch = new Branch(nodeCount - 2, currentNode, leaf);
                    _branches.push_back(branch);
                }
				_leaves.push_back(leaf);
			} else {
				prevInternalNode->setLabel(label);
				label.clear();
				prevInternalNode->setComment(comment);
				comment.clear();
//				Branch *branch = prevInternalNode->getBranch(0);
//				branch->setDistance(distance);
			}

            if (treeString[i] == ')') { // internal node ends
                if (verbose >= 5)
                    cout << "Internal node " << currentNode->getIdent() << " " << _internalNodes.size() << "/" << nodeCount << ") ends" << endl;
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
				throw("Error parsing Newick tree: a square bracket opened at position " + str(i) + " but no closing bracket could be found.");
			else
			{
				comment = treeString.substr(i + 1, j - i - 1);
				i = j + 1;
				//if the tree is yet empty, the comment applies to the tree
				if ((_branches.size() == 0) && (_internalNodes.size() == 0) && (_leaves.size() == 0)) {
					setComment(comment);
					comment.clear();
				}
			}
		}
	}

	prevInternalNode->setLabel(label);
	label.clear();
	
	if (comment.size() > 0) {
		setComment(comment);
		comment.clear();
	}

	if (prevInternalNode->getBranches().size() == 1)
	{
		_root = prevInternalNode;
		cout << "This Newick representation appears to be a rooted tree, rooted at the leaf node " << _root->getLabel() << ". Perfect." << endl;
		_leaves.push_back(_root);
		if (_internalNodes.front() == _root) _internalNodes.erase(_internalNodes.begin());
		_unrooted = false;
    } else {
      if (prevInternalNode->getBranches().size() == 2) {
		cout << "This Newick representation appears to be a rooted tree, rooted at the internal node " << prevInternalNode->getIdent() << " with two descendants " << prevInternalNode->getChild(0)->getIdent() << " and " << prevInternalNode->getChild(1)->getIdent() << ". Removing internal node " << prevInternalNode->getIdent() << endl;
        
        Node* child0 = prevInternalNode->getChild(0);
        Node* child1 = prevInternalNode->getChild(1);
        Branch* child0Branch = prevInternalNode->getBranch(0);
        Branch* child1Branch = prevInternalNode->getBranch(1);
        child0Branch->unlinkNode(prevInternalNode); // Remove internal node from both outgoing branches
        child1Branch->unlinkNode(prevInternalNode);
        child1Branch->unlinkNode(child1); // Remove child node from outgoing branch
        child0Branch->linkNode(child1); // Link child node to other branch
        removeNode(prevInternalNode); // Remove node from the list
        removeBranch(child1Branch); // Remove branch from list
        if (child0->isLeaf()) {
          _root = child0;
        } else if (child1->isLeaf()) {
          _root = child1;
        } else {
          prevInternalNode = child0;
        }
        if (_root != NULL) {
          _root->reroot(NULL);
        }
        _unrooted = false;
      } else if (prevInternalNode->getBranches().size() == 3) {
		cout << "This Newick representation appears to be an unrooted tree, rooted at the internal node " << prevInternalNode->getIdent() << "." << endl;
      }
      if (_root == NULL && !options.rootNode.length()) {
        vector<Branch*> branches = prevInternalNode->getBranches();
		for (unsigned int i = 0; i < branches.size(); i++) {
          _root = branches[i]->getNeighbour(prevInternalNode);
          if (_root->getBranches().size() == 1)
            break;
        }
        if (_root->getBranches().size() == 1) {
          cout << "The BH+I model requires the root to be placed at a leaf node, picking " << _root->getIdent() << " as root." << endl;
          _root->reroot(NULL);
        } else {
          cout << "Arbitrarily picking leaf node " << _leaves[0]->getIdent() << " as root." << endl;
          _root = _leaves[0];
          _root->reroot(NULL);
//				throw(string("The BH+I model requires the root to placed at a leaf node, but this node has no leaves as children."));
        }
      }
      _unrooted = true;
	}

	if (options.rootNode.length())
	{
		_root = NULL;
		for (unsigned int i = 0; i < _leaves.size(); i++)
			if (_leaves[i]->getLabel() == options.rootNode)
			{
				_root = _leaves[i];
				_root->reroot(NULL);
				break;
			}
		if (_root != NULL)
			cout << "Placing root at the user-specified leaf node " << _root->getIdent() << endl;
		else
			throw("The user-specified root \"" + options.rootNode + "\" does not exist in the tree or is not a leaf node");
	}

	for (unsigned int i = 0; i < _leaves.size(); i++)
	{
		label = _leaves[i]->getLabel();
		int a = _alignment->find(label);
		if (a < 0)
			throw("The alignment contains no sequence \"" + label + "\"");
		else
		{
			unsigned int* seq = _alignment->getCompressedSeq(a);
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

double Tree::getLogLH() const
{
	return _logLH;
}

void Tree::computeLH()
{
	_logLH = _branches[0]->computeLH(_alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart());
}

void Tree::updateModel(double qDelta, double betaDelta, bool thorough)
{
	unsigned int updates;
	double prevLH;

	do
	{
		prevLH  = getLogLH();
		for (unsigned int i = 0; i < _branches.size(); i++)
			_branches[i]->computeUpdatedQ(_alignment->getNumOfSites(), _alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart());

		updates = 0;
		for (unsigned int i = 0; i < _branches.size(); i++)
		{
			if (_branches[i]->updateQ(qDelta)) updates++;

			if (_branches[i]->updateParameters(_alignment->getNumOfSites(), _alignment->getPatternCount(), _alignment->getInvarSites(), _alignment->getInvarStart(),
					betaDelta)) updates++;
		}
		computeLH();
	} while (updates > 0 || (thorough && fabs(prevLH - getLogLH()) > 0.1));
}

void Tree::clearModel() {
	for (unsigned int m = 0; m < _branches.size(); m++)
	{
		_branches[m]->resetVectors();
		_branches[m]->resetQ();
	}
}

void Tree::printNodes()
{
	cout << "Internal Nodes: " << endl;
	for (unsigned int i = 0; i < _internalNodes.size(); i++)
	{
      cout << _internalNodes[i]->getIdent() << " parent=" << _internalNodes[i]->getParent()->getIdent();
      if (_internalNodes[i]->getChildCount() >= 1) {
        cout << " child1=" << _internalNodes[i]->getChild(1)->getIdent();
        if (_internalNodes[i]->getChildCount() == 2) {
          cout << " child2=" << _internalNodes[i]->getChild(2)->getIdent();
        }
      }
      cout << endl;
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

string Tree::toString(const bool topologyOnly) const
{
	if (_root == NULL) {
		throw("Tree::toString() _root is NULL");
	}
	
    stringstream ss;
	if (!topologyOnly) {
		ss << "[" << fixed << setprecision(6) << getLogLH() << "]";
	}
	if (_unrooted) {
		Node* neighbour = _root->getBranch(0)->getNeighbour(_root);
		if (neighbour) {
			ss << _root->getBranch(0)->getNeighbour(_root)->toString(NULL, topologyOnly);
		} else {
			throw("Tree::toString() Unable to locate neighbour");
		}
	} else {
      ss << _root->toString(NULL, topologyOnly);
	}
	
	ss << ";";
	return ss.str();
}

void Tree::print()
{
  cout << toString(false) << endl;
}
