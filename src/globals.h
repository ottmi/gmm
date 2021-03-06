#ifndef GLOBALS_H_
#define GLOBALS_H_
#include <string>
using namespace std;

#define PROGNAME "gmm"
#define VERSION "0.1.47"
#define PROGDATE "2020-12-22"

#define _DNA_DATA				0
#define	_AA_DATA				1
#define _DICODON_DATA		2

#define _DNA_MAP      "ACGTURYKMSWBDHVN?-"
#define _AA_MAP       "ACDEFGHIKLMNPQRSTVWYBJXZ?-"

typedef struct opt_struct
{
	string alignment;
	unsigned int alignmentGrouping;
	string inputTree;
	string rootNode;
	bool evaluateOnly;
	double cutOff;
	bool restart;
	unsigned int maxBestTrees;
	bool help;
} Options;

extern int verbose;
extern unsigned int charStates;
//extern class Tree bestTree;

#endif /* GLOBALS_H_ */
