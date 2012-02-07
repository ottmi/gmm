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
	void read(string fileName);
	int find(string name);
	unsigned int getNumOfSites() { return _sequences[0].size(); };
	unsigned int getNumOfUniqueSites() { return _compressedSequences[0].size(); };
	unsigned int getNumOfSequences() { return _sequences.size(); };
	int getDataType() { return _dataType; };

	vector<unsigned int> getNumericalSeq(unsigned int row);
	vector<unsigned int>& getPatternCount() { return _patternCount; };
	unsigned int getInvarStart() {return _invarStart; };
	vector<unsigned int>& getInvarSites() {return _invarSites; };

private:
	int _dataType;
	vector <string> _names;
	vector <string> _sequences;
	vector <string> _compressedSequences;
	vector<unsigned int> _patternCount;
	unsigned int _invarStart;
	vector<unsigned int> _invarSites;

	void readPhylip(string fileName);
	void readFasta(string fileName);
	void identifyDataTpe();
	void compress();
	void identifyInvarSites();
};

#endif /* ALIGNMENT_H_ */
