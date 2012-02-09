#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <string>
#include <vector>
#include <fstream>

using namespace std;

typedef vector <unsigned int> unsignedIntVec_t;

class Alignment
{
public:
	Alignment();
	virtual
	~Alignment();
	void read(string fileName, unsigned int grouping);
	int find(string name);
	unsigned int getNumOfSites() { return _numericalSequences[0].size(); };
	unsigned int getNumOfUniqueSites() { return _compressedSequences[0].size(); };
	unsigned int getNumOfSequences() { return _numericalSequences.size(); };
	int getDataType() { return _dataType; };

	unsignedIntVec_t& getNumericalSeq(unsigned int row);
	unsignedIntVec_t& getPatternCount() { return _patternCount; };
	unsigned int getInvarStart() {return _invarStart; };
	unsignedIntVec_t& getInvarSites() {return _invarSites; };

private:
	int _dataType;
	vector <string> _names;
	vector <string> _sequences;
	vector <unsignedIntVec_t> _numericalSequences;
	vector <unsignedIntVec_t> _compressedSequences;
	unsignedIntVec_t _patternCount;
	unsigned int _invarStart;
	unsignedIntVec_t _invarSites;

	void readPhylip(string fileName);
	void readFasta(string fileName);
	void identifyDataTpe(unsigned int grouping);
	void translateToNumerical(unsigned int grouping);
	void compress();
};

#endif /* ALIGNMENT_H_ */
