#include "globals.h"
#include "Tree.h"
#include "Alignment.h"
#include "Optimizer.h"
#include "helper.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <csignal>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int verbose = 0;
unsigned int charStates = 4;
Tree bestTree;

void cancellationHandler(int parameter) {
	cout << endl << endl;
	cout << "Program has been canceled. Current best tree:" << endl;
	cout << "logLH: " << fixed << setprecision(10) << bestTree.getLogLH() << endl;
	bestTree.print();
	exit(1);
}

void printSyntax() {
	cout << "Syntax:" << endl;
	cout << "  gmm -s<FILE> [-d|-c] -t<FILE> [-r<STRING>] [-e] [-x<NUM>] [-v[NUM]]" << endl;
	cout << "  gmm -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -s\tSequence alignment" << endl;
	cout << "  -d\tTreat alignment sequences as dicodons" << endl;
	cout << "  -c\tTreat alignment sequences as codons" << endl;
	cout << "  -t\tInput tree" << endl;
	cout << "  -r\tRoot the tree at this node (has to be a leaf)" << endl;
	cout << "  -e\tOnly evaluate given tree topology" << endl;
	cout << "  -x\tOptimization cutoff [default: 0.0001]" << endl;
	cout << "  -v\tBe increasingly verbose" << endl;
	cout << "  -h\tThis help page" << endl;
	cout << endl;
}

int parseArguments(int argc, char** argv, Options *options) {
	char c;

	options->help = false;
	options->alignmentGrouping = 1;
	options->evaluateOnly = false;
	options->cutOff = 0.0001;

	while ((c = getopt(argc, argv, "s:dct:r:ex:v::h")) != -1) {
		switch (c) {
		case 's':
			options->alignment = optarg;
			break;
		case 'd':
			options->alignmentGrouping = 2;
			break;
		case 'c':
			options->alignmentGrouping = 3;
			break;
		case 't':
			options->inputTree = optarg;
			break;
		case 'r':
			options->rootNode = optarg;
			break;
		case 'e':
			options->evaluateOnly = true;
			break;
		case 'x': {
			stringstream ss(optarg);
			ss >> options->cutOff;
			break;
		}
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

	if (argc == 1 || options->help) {
		printSyntax();
		return 255;
	}

	if (options->alignment.empty()) {
		cerr << "Please specify a sequence alignment with -s" << endl;
		return 254;
	}

	if (options->inputTree.empty()) {
		cerr << "Please specify an input tree with -t" << endl;
		return 253;
	}

	return 0;
}

int main(int argc, char **argv) {
	cout << PROGNAME << " " << VERSION << "|";
#ifdef _OPENMP
	cout << "OpenMP|";
#endif
#ifdef _DEBUG
	cout << "DEBUG|";
#endif

	cout << PROGDATE << endl << endl;

	Options options;
	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;

	signal(SIGINT, cancellationHandler);
	signal(SIGTERM, cancellationHandler);

#ifdef _OPENMP
	cout << "Running in parallel on " << omp_get_max_threads() << " Threads" << endl << endl;
#endif


	ifstream treeFileReader;
	vector<string> treeStrings;
	if (options.inputTree.find(';') == string::npos) {
		treeFileReader.open(options.inputTree.c_str());
		if (!treeFileReader.is_open())
			throw("Cannot open file " + options.inputTree);
		while (!treeFileReader.eof()) {
			string s;
			safeGetline(treeFileReader, s);
			if (s.find(';') != string::npos)
				treeStrings.push_back(s);
		}
		treeFileReader.close();
	} else
		treeStrings.push_back(options.inputTree);

	Alignment alignment;
	try {
		alignment.read(options.alignment, options.alignmentGrouping);
	} catch (string& s) {
		cerr << s << endl;
		return (255);
	}

	for (unsigned int i = 0; i < treeStrings.size(); i++) {
		cout << endl << "Processing tree " << i + 1 << "/" << treeStrings.size() << endl;
		try {
			Tree tree;
			tree.readNewick(&alignment, treeStrings[i], options);
			if (verbose >= 2) {
				tree.printNodes();
				tree.printBranches();
			}

			if (options.evaluateOnly) {
				tree.updateModel(options.cutOff, options.cutOff, true);
				cout << "logLH: " << fixed << setprecision(10) << tree.getLogLH() << endl;
			} else {
				Optimizer optimizer;
				optimizer.rearrange(tree, options);
			}

			tree.print();
		} catch (string& s) {
			cerr << s << endl;
			if (treeStrings.size() > 1) {
				cout << "Skipping tree" << endl;
			}
		}
	}

	return 0;
}
