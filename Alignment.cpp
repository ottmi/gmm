#include "Alignment.h"
#include "globals.h"
#include "helper.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>

Alignment::Alignment()
{
}

Alignment::~Alignment()
{
	// TODO Auto-generated destructor stub
}

void Alignment::read(string fileName, unsigned int grouping)
{
	string ext = fileName.substr(fileName.find_last_of('.') + 1);

	if (!ext.compare("phy") || !ext.compare("phylip"))
		readPhylip(fileName);
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
		readFasta(fileName);
	else
	{
		cerr << "Unknown input alignment format" << endl;
		exit(255);
	}
	cout << "The alignment contains " << _sequences.size() << " sequences with " << _sequences[0].size() << " characters each." << endl;

	identifyDataTpe(grouping);
	string dataTypeDesc[] = { "DNA", "AA" };
	cout << "The data appears to be " << dataTypeDesc[_dataType] << ", using " << charStates << "-state system." << endl;

	translateToNumerical(grouping);
	compress();
	cout << "There are " << getNumOfUniqueSites() << " unique sites, " << _invarSites.size() << " of which are invariant." << endl << endl;
}

int Alignment::find(string name)
{
	for (unsigned int i = 0; i < _names.size(); i++)
	{
		if (_names[i] == name) return i;
	}

	return -1;
}

unsignedIntVec_t& Alignment::getNumericalSeq(unsigned int row)
{
	if (row < getNumOfSequences())
		return _compressedSequences[row];
	else
	{
		stringstream ss;
		ss << "Alignment::getNumericalSeq(" << row << "): There are only " << getNumOfSequences() << " sequences in the alignment";
		throw(ss.str());
	}
}

void Alignment::readPhylip(string fileName)
{
	ifstream fileReader;
	fileReader.open(fileName.c_str());
	if (!fileReader.is_open()) throw("\n\nError, cannot open file " + fileName);

	string str;
	safeGetline(fileReader, str);

	string whiteSpace = " \t";
	str = str.substr(str.find_first_not_of(whiteSpace));
	string rowsStr = str.substr(0, str.find_first_of(whiteSpace));

	str = str.substr(str.find_first_of(whiteSpace));
	str = str.substr(str.find_first_not_of(whiteSpace));
	string colsStr = str.substr(0, str.find_first_of(whiteSpace));

	unsigned int cols = atoi(colsStr.c_str());
	unsigned int rows = atoi(rowsStr.c_str());

	while (!fileReader.eof())
	{
		safeGetline(fileReader, str);
		if (str.length())
		{
			str = str.substr(str.find_first_not_of(whiteSpace));
			string name = str.substr(0, str.find_first_of(whiteSpace));
			str = str.substr(str.find_first_of(whiteSpace));
			string seq = str.substr(str.find_first_not_of(whiteSpace));

			if (seq.length() < cols) cerr << "Sequence #" << _sequences.size() + 1 << " (" << name << ") consists only of" << seq.length() << " characters." << endl;
			seq = adjustString(seq);
			if (!name.empty() && !seq.empty())
			{
				_names.push_back(name);
				_sequences.push_back(seq);
			}
		}
	}
	if ( _sequences.size() < rows) cerr << "The alignment consists only of " << _sequences.size() << " rows." << endl;
}

void Alignment::readFasta(string fileName)
{
	ifstream fileReader;
	fileReader.open(fileName.c_str());
	if (!fileReader.is_open()) throw("\n\nError, cannot open file " + fileName);

	string s;
	safeGetline(fileReader, s);
	while ((!fileReader.eof()) && s[0] != '>')
		safeGetline(fileReader, s);

	while (!fileReader.eof())
	{
		string header;
		string seq;
		header = s;
		s.clear();
		while (!fileReader.eof() && s[0] != '>')
		{
			safeGetline(fileReader, s);
			if (s[0] != '>') seq += s;
		}
		seq = adjustString(seq, true);
		header = adjustString(header.substr(1), false);
		if (header.length() > 1 && seq.length())
		{
			_names.push_back(header);
			_sequences.push_back(seq);
		}
	}
}

void Alignment::identifyDataTpe(unsigned int grouping)
{
	map<char, unsigned long> baseOccurences;
	for (unsigned int i = 0; i < _sequences.size(); i++)
	{
		string s = _sequences[i];
		for (unsigned int j = 0; j < s.length(); j++)
			baseOccurences[s[j]]++;
	}

	string maps[] = { _DNA_MAP, _AA_MAP };
	unsigned long counts[2];
	for (unsigned int i = 0; i < 2; i++)
	{
		counts[i] = 0;
		string map = maps[i];
		for (unsigned j = 0; j < map.length(); j++)
			counts[i] += baseOccurences[map[j]];
	}

	if (verbose) cout << counts[0] << " DNA characters and " << counts[1] << " AA characters." << endl;
	if (counts[0] >= counts[1])
	{
		_dataType = _DNA_DATA;
		if (grouping == 1)
			charStates = 4;
		else if (grouping == 2)
			charStates = 16;
		else if (grouping == 3)
			charStates = 64;
	} else
	{
		_dataType = _AA_DATA;
		charStates = 20;
	}
}

void Alignment::translateToNumerical(unsigned int grouping)
{
	for (unsigned int i = 0; i<_sequences.size(); i++)
	{
		unsignedIntVec_t numSeq;
		for (unsigned int j = 0; j<_sequences[i].size(); j+= grouping)
		{
			if (_dataType == _DNA_DATA)
				numSeq.push_back(mapDNAToNum(_sequences[i].substr(j, grouping)));
			else
				numSeq.push_back(mapAAToNum(_sequences[i][j]));
		}
		_numericalSequences.push_back(numSeq);
	}
}

void Alignment::compress()
{
	map<unsignedIntVec_t, unsigned int> patterns;
	for (unsigned int col = 0; col < getNumOfSites(); col++)
	{
		unsignedIntVec_t site;
		for (unsigned int row = 0; row < getNumOfSequences(); row++)
			site.push_back(_numericalSequences[row][col]);
		patterns[site]++;
	}

	_compressedSequences = vector<unsignedIntVec_t>(getNumOfSequences());
	vector<unsignedIntVec_t> invar(getNumOfSequences());
	unsignedIntVec_t invarCount;
	for (map<unsignedIntVec_t, unsigned int>::iterator it = patterns.begin(); it != patterns.end(); it++)
	{
		unsigned int row = 1;
		while (row < getNumOfSequences() && it->first[0] == it->first[row])
			row++;

		if (row == getNumOfSequences())	// this site contains only identical characters
		{
			for (unsigned int i = 0; i < getNumOfSequences(); i++)
				invar[i].push_back(it->first[i]);
			invarCount.push_back(it->second);
		} else
		{
			for (unsigned int i = 0; i < getNumOfSequences(); i++)
				_compressedSequences[i].push_back(it->first[i]);
			_patternCount.push_back(it->second);
		}
	}

	_invarStart = _patternCount.size();
	for (unsigned col = 0; col < invarCount.size(); col++) // append the invariable sites to the variable sites
	{
		for (unsigned int row = 0; row < getNumOfSequences(); row++)
			_compressedSequences[row].push_back(invar[row][col]);
		_patternCount.push_back(invarCount[col]);
		_invarSites.push_back(invar[0][col]);
	}


	if (verbose >= 10)
	{
		cout << "\t";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << col / 10;
		cout << endl << "\t";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << col % 10;
		cout << endl;

		for (unsigned int row = 0; row < getNumOfSequences(); row++)
		{
			cout << _names[row] << "\t";
			for (unsigned int col= 0; col < getNumOfUniqueSites(); col++)
				 cout << mapNumToDNA(_compressedSequences[row][col], 1);
			cout << endl;

		}

		cout << "\t";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << _patternCount[col];
		cout << endl;
	}

}
