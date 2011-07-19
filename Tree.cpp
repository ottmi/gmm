#include "Tree.h"
#include "Node.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>

Tree::Tree()
{
	// TODO Auto-generated constructor stub

}

Tree::~Tree()
{
	// TODO Auto-generated destructor stub
}

istream& safeGetline(istream& is, string& t)
{
	/* Courtesy of http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf */
    t.clear();
    istream::sentry se(is);
    streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        default:
            t += (char)c;
        }
    }
}

void Tree::readNewick(string fileName)
{
/*
	cout << "readNewick(" << fileName << ")" << endl;

	ifstream _fileReader;
	_fileReader.open(fileName.c_str());
	if (! _fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName );

	string str;
	safeGetline(_fileReader, str);
	_fileReader.close();
*/
	string treeStr = fileName;
	cout << "Tree: " << treeStr << endl;

	unsigned int i = 0;
	unsigned int nodeCount = 0;
	Node *lastInternalNode, *currentNode;
	lastInternalNode = currentNode = NULL;
	string name;
	while (i < treeStr.length() && treeStr[i] != ';')
	{
		if (treeStr[i] == '(') // internal node starts
		{
			cout << "Internal Node #" << nodeCount << " starts" << endl;
			currentNode = new Node(currentNode, nodeCount, false);
			lastInternalNode = currentNode;
			nodeCount++;
			i++;
		} else
		if (treeStr[i] == ')') // internal node ends
		{
			if (!currentNode->right) // the right child must have been a leaf
			{
				cout << "  Leaf #" << nodeCount << " (" << name << ")" << endl;
				Node *leaf = new Node(currentNode,  nodeCount, true);
				nodeCount++;
				leaf->setLabel(name);
				leaves.push_back(leaf);
			} else // that label belongs to the previous internal node
			{
				lastInternalNode->setLabel(name);
				name = "";
			}

			cout << "Internal Node #" << currentNode->id << " ends" << endl;

			if (!currentNode->left)
				cerr << "Error: Internal node has no left child !" << endl;
			internalNodes.push_back(currentNode);
			lastInternalNode = currentNode;
			currentNode = currentNode->parent;
			i++;
		} else
		if (treeStr[i] == ',') // the 2nd or 3rd child starts here
		{
			if (!currentNode->left) // the 1st child must have been a leaf
			{
				cout << "  Leaf #" << nodeCount << " (" << name << ")" << endl;
				Node *leaf = new Node(currentNode,  nodeCount, true);
				nodeCount++;
				leaf->setLabel(name);
				leaves.push_back(leaf);
			} else if (!currentNode->right) // the 2nd child must have been a leaf and this is the root
			{
					cout << "  Leaf #" << nodeCount << " (" << name << ")" << endl;
					Node *leaf = new Node(currentNode,  nodeCount, true);
					nodeCount++;
					leaf->setLabel(name);
					leaves.push_back(leaf);
			}
			i++;
		}

		if (isalpha(treeStr[i])) // this is a label
		{
			int j=treeStr.find_first_of(",():;", i+1);
			name = treeStr.substr(i, j-i);
			i = j;
		}
	}

	lastInternalNode->setLabel(name);
	root = lastInternalNode;

	cout << "The tree has " << internalNodes.size() << " internal nodes and " << leaves.size() << " leaves." << endl;
}

void Tree::print()
{
	for (unsigned int i=0; i<internalNodes.size(); i++)
	{
		Node *root = internalNodes[i];
		string children = root->toString();
		string parent = root->parent->toString();

		cout << children.substr(0, children.find_last_of(')')) << "," << parent << ")" << root->getLabel() << ";" << endl;
	}
}
