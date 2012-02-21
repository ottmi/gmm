#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <string>
#include <vector>
#include <fstream>

using namespace std;

class Alignment
{
public:
	Alignment();
	virtual
	~Alignment();
	void read(string fileName, unsigned int grouping);
	int find(string name);
	unsigned int getNumOfSites() { return _numOfSites; };
	unsigned int getNumOfUniqueSites() { return _numOfUniqueSites; };
	unsigned int getNumOfSequences() { return _numOfSequences; };
	int getDataType() { return _dataType; };

	unsigned int* getNumericalSeq(unsigned int row);
	vector <unsigned int>& getPatternCount() { return _patternCount; };
	unsigned int getInvarStart() {return _invarStart; };
	vector <unsigned int>& getInvarSites() {return _invarSites; };

private:
	int _dataType;
	vector <string> _names;
	vector <string> _sequences;
	vector <unsigned int*> _numericalSequences;
	vector <unsigned int*> _compressedSequences;
	vector <unsigned int> _patternCount;
	vector <unsigned int> _invarSites;
	unsigned int _invarStart;
	unsigned int _numOfSites;
	unsigned int _numOfUniqueSites;
	unsigned int _numOfSequences;

	void readPhylip(string fileName);
	void readFasta(string fileName);
	void identifyDataTpe(unsigned int grouping);
	void translateToNumerical(unsigned int grouping);
	void compress();
};

#endif /* ALIGNMENT_H_ */
