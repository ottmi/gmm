#include "Alignment.h"
#include "globals.h"
#include "helper.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>

Alignment::Alignment()
{
}

Alignment::~Alignment()
{
	for (unsigned int i = 0; i < _compressedSequences.size(); i++)
		delete[] _compressedSequences[i];
}

void Alignment::read(string fileName, unsigned int grouping)
{
	string ext = fileName.substr(fileName.find_last_of('.') + 1);

	if (!ext.compare("phy") || !ext.compare("phylip"))
		readPhylip(fileName);
	else if (!ext.compare("fsa") || !ext.compare("fst") || !ext.compare("fasta"))
		readFasta(fileName);
	else
		throw(string("Unknown input alignment format, only PHYLIP (.phy|.phylip) and FASTA (.fsa|.fst|.fasta) are supported."));
	cout << "The alignment contains " << _sequences.size() << " sequences with " << _sequences[0].size() << " characters each." << endl;
	if (_sequences[0].size() % grouping != 0)
		throw("The alignment is supposed to be grouped into sites of " + str(grouping) + " columns each, but " + str(_sequences[0].size()) + " is not divisible by "
				+ str(grouping));
	_numOfSites = _sequences[0].size() / grouping;

	identifyDataTpe(grouping);
	string dataTypeDesc[] = { "DNA", "AA" };
	cout << "The data appears to be " << dataTypeDesc[_dataType] << ", using " << charStates << "-state system." << endl;

	compress(grouping);
	cout << "There are " << getNumOfUniqueSites() << " unique sites, " << _invarSites.size() << " of which are invariant." << endl;
	if (_ambiguousSites) cout << _ambiguousSites << " sites contain invalid or ambiguous characters and have therefore been ignored." << endl;
}

int Alignment::find(string name)
{
	for (unsigned int i = 0; i < _names.size(); i++)
	{
		if (_names[i] == name) return i;
	}

	return -1;
}

unsigned int* Alignment::getCompressedSeq(unsigned int row)
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

	string s;
	safeGetline(fileReader, s);

	string whiteSpace = " \t";
	s = s.substr(s.find_first_not_of(whiteSpace));
	string rowsStr = s.substr(0, s.find_first_of(whiteSpace));

	s = s.substr(s.find_first_of(whiteSpace));
	s = s.substr(s.find_first_not_of(whiteSpace));
	string colsStr = s.substr(0, s.find_first_of(whiteSpace));

	unsigned int cols = atoi(colsStr.c_str());
	unsigned int rows = atoi(rowsStr.c_str());

	while (!fileReader.eof())
	{
		safeGetline(fileReader, s);
		if (s.length())
		{
			s = s.substr(s.find_first_not_of(whiteSpace));
			string name = s.substr(0, s.find_first_of(whiteSpace));
			s = s.substr(s.find_first_of(whiteSpace));
			string seq = s.substr(s.find_first_not_of(whiteSpace));
			if (seq.length() < cols) throw("Sequence #" + str(_sequences.size() + 1) + " (" + name + ") consists only of " + str(seq.length()) + " characters.");
			seq = adjustString(seq);
			if (!name.empty() && !seq.empty())
			{
				_names.push_back(name);
				_sequences.push_back(seq);
			}
		}
	}
	if (_sequences.size() < rows) throw("The alignment consists only of " + str(_sequences.size()) + " rows.");
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
		else if (grouping == 3) charStates = 64;
	} else
	{
		_dataType = _AA_DATA;
		charStates = 20;
	}
}

void Alignment::compress(unsigned int grouping)
{
	_ambiguousSites = 0;
	map<vector<unsigned int>, unsigned int> patterns;
	for (unsigned int col = 0; col < _numOfSites * grouping; col += grouping)
	{
		vector<unsigned int> site;
		bool ambiguous = false;
		for (unsigned int row = 0; row < getNumOfSequences(); row++)
		{
			unsigned int val;
			if (_dataType == _DNA_DATA)
				val = mapDNAToNum(_sequences[row].substr(col, grouping));
			else
				val = mapAAToNum(_sequences[row][col]);
			if (val != 0xffffffff)
				site.push_back(val);
			else
			{
				ambiguous = true;
				break;
			}
		}
		if (!ambiguous)
			patterns[site]++;
		else
		{
			if (verbose >= 3)
			{
				cout << "Site " << col + 1;
				if (grouping > 1) cout << " (columns " << col * grouping + 1 << "-" << col * grouping + grouping << ")";
				cout << " contains invalid or ambiguous characters and will therefore be ignored." << endl;
			}
			_ambiguousSites++;
		}
	}

	_compressedSequences = vector<unsigned int*>(getNumOfSequences());
	for (unsigned int i = 0; i < getNumOfSequences(); i++)
		_compressedSequences[i] = new unsigned int[patterns.size()];
	_patternCount.resize(patterns.size(), 0);

	unsigned int k = 0;
	unsigned int l = patterns.size();
	for (map<vector<unsigned int>, unsigned int>::iterator it = patterns.begin(); it != patterns.end(); it++)
	{
		unsigned int row = 1;
		while (row < getNumOfSequences() && it->first[0] == it->first[row])
			row++;

		if (row == getNumOfSequences()) // this site contains only identical characters
		{
			_patternCount[--l] = it->second;
			_invarSites.push_back(it->first[0]);
			for (unsigned int i = 0; i < getNumOfSequences(); i++)
				_compressedSequences[i][l] = it->first[i];
		} else
		{
			for (unsigned int i = 0; i < getNumOfSequences(); i++)
				_compressedSequences[i][k] = it->first[i];
			_patternCount[k++] = it->second;
		}
	}
	_invarStart = l;

	if (verbose >= 10)
	{
		cout << setw(11) << " ";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << setw(grouping) << col / 10 << " ";
		cout << endl;
		cout << setw(11) << " ";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << setw(grouping) << col % 10 << " ";
		cout << endl;

		for (unsigned int row = 0; row < getNumOfSequences(); row++)
		{
			cout.flags(ios::left);
			cout << setw(10) << _names[row] << " ";
			for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
				cout << mapNumToDNA(_compressedSequences[row][col], grouping) << " ";
			cout << endl;

		}
		cout.flags(ios::right);

		cout << setw(11) << " ";
		for (unsigned int col = 0; col < getNumOfUniqueSites(); col++)
			cout << setw(grouping) << _patternCount[col] << " ";
		cout << endl;
	}

}
