#include "Tree.h"
#include "Node.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <sstream>

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
	Node *prevInternalNode, *currentNode;
	prevInternalNode = currentNode = NULL;
	string label;
	double distance = -1.0;
	bool nextCouldBeLeaf = false;
	while (i < treeStr.length() && treeStr[i] != ';')
	{
		if (treeStr[i] == '(') // internal node starts
		{
			cout << "Internal Node #" << nodeCount << " starts" << endl;
			currentNode = new Node(currentNode, nodeCount, false);
			internalNodes.push_back(currentNode);
			nodeCount++;
			nextCouldBeLeaf = true;
			i++;
		} else if (treeStr[i] == ')' || treeStr[i] == ',') // node ends, could be internal or leaf
		{
			if (nextCouldBeLeaf)
			{
				cout << "  Leaf #" << nodeCount << " (" << label << ") " << distance << endl;
				Node *leaf = new Node(currentNode,  nodeCount, true);
				leaf->setLabel(label);
				label = "";
				leaves.push_back(leaf);
				nodeCount++;
			} else
			{
				prevInternalNode->setLabel(label);
				cout << label << " " << distance << endl;
				label = "";
			}

			if  (treeStr[i] == ')') // internal node
			{
				cout << "Internal Node #" << currentNode->id << " ends " << endl;
				prevInternalNode = currentNode;
				currentNode = currentNode->neighbours[0];
				nextCouldBeLeaf=false;
			} else // node will follow, could be a leaf
			{
				nextCouldBeLeaf=true;
			}
			i++;
		} else if (treeStr[i] == ':') // distance for last node
		{
			int j=treeStr.find_first_of(",():;", i+1);
			stringstream ss(treeStr.substr(i+1, j-i-1));
			ss >> distance;
			i = j;
		}

		if (isalpha(treeStr[i])) // this is a label
		{
			int j=treeStr.find_first_of(",():;", i+1);
			label = treeStr.substr(i, j-i);
			i = j;
		}
	}

	root = prevInternalNode;
	root->setLabel(label);

	cout << "The tree has " << internalNodes.size() << " internal nodes and " << leaves.size() << " leaves." << endl;
}

void Tree::print()
{
	for (unsigned int i=0; i<internalNodes.size(); i++)
	{
		Node *root = internalNodes[i];
		cout << root->toString(NULL) << ";" << endl;
/*
		string children = root->toString();
		string parent = root->parent->toString();
		cout << children.substr(0, children.find_last_of(')')) << "," << parent << ")" << root->getLabel() << ";" << endl;
*/

	}
}
