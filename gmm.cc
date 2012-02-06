#include "globals.h"
#include "Tree.h"
#include "Alignment.h"
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;

int verbose = 0;
int charStates = 4;


void printSyntax()
{
	cout << "Syntax:" << endl;
	cout << "  gmm -s <FILE> -t<FILE> [-v[NUM]]" << endl;
	cout << "  gmm -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -s\tSequence alignment" << endl;
	cout << "  -t\tInput tree" << endl;
	cout << "  -v\tBe increasingly verbose" << endl;
	cout << "  -h\tThis help page" << endl;
	cout << endl;
}


int parseArguments(int argc, char** argv, Options *options)
{
	char c;

	options->help = false;

	while ( (c = getopt(argc, argv, "s:t:v::h")) != -1)
	{
		switch (c)
		{
			case 's':
				options->alignment = optarg;
				break;
			case 't':
				options->inputTree = optarg;
				break;
			case 'v':
				if (optarg)
					verbose = atoi(optarg);
				else
					verbose = 1;
				break;
			case 'h':
				options->help = true;
				break;
			default:
				if (c != '?')
					cerr << "Unknown parameter: " << c << endl;
				return 1;
		}
	}

	if (argc == 1 || options->help)
	{
		printSyntax();
		return 255;
	}

	if (options->alignment.empty())
	{
		cerr << "Please specify a sequence alignment with -s" << endl;
		return 254;
	}

	if (options->inputTree.empty())
	{
		cerr << "Please specify an input tree with -t" << endl;
		return 253;
	}

	return 0;
}


int main(int argc, char **argv)
{
	Options options;

	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;

	Alignment alignment;
	try {
		alignment.read(options.alignment);

		Tree tree(alignment);
		tree.readNewick(options.inputTree);
		if (verbose >= 2)
			tree.printNodes();

		tree.computeLH();
		tree.updateQ();
		tree.computeLH();
		tree.updateQ();
		tree.computeLH();
	}
	catch (string& s)
	{
		cerr << s << endl;
		return(255);
	}

	return 0;
}
