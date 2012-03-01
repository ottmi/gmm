#ifndef GLOBALS_H_
#define GLOBALS_H_
#include <string>
using namespace std;

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
	bool evaluateOnly;
	double cutOff;
	bool help;
} Options;

extern int verbose;
extern unsigned int charStates;

#endif /* GLOBALS_H_ */
