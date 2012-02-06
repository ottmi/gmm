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
	int getCols() { return _cols; };
	int getRows() { return _rows; };
	int getDataType() { return _dataType; };

	vector<unsigned int> getNumericalSeq(unsigned int row);
	vector<unsigned int>& getInvarSites() {return _invarSites; };

private:
	int _cols;
	int _rows;
	int _dataType;
	vector <string> _sequences;
	vector <string> _names;
	vector<unsigned int> _invarSites;

	void readPhylip(string fileName);
	void readFasta(string fileName);
	void identifyDataTpe();
	void identifyInvarSites();
};

#endif /* ALIGNMENT_H_ */
