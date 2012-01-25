#include "Alignment.h"
#include "globals.h"
#include "helper.h"
#include <stdlib.h>
#include <iostream>
#include <map>



Alignment::Alignment() {
	_cols = _rows = 0;
}



Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
}



void Alignment::read(string fileName) {
	string ext = fileName.substr(fileName.find_last_of('.') + 1);

	if (!ext.compare("phy") || !ext.compare("phylip"))
		readPhylip(fileName);
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
		readFasta(fileName);
	else {
		cerr << "Unknown input alignment format" << endl;
		exit(255);
	}

	cout << "The alignment contains " << _rows << " sequences with " << _cols
			<< " sites each." << endl;

	string dataTypeDesc[] = { "DNA", "AA" };
	map<char, unsigned long> baseOccurences;
	for (unsigned int i = 0; i < _sequences.size(); i++) {
		string s = _sequences[i];
		for (unsigned int j = 0; j < s.length(); j++)
			baseOccurences[s[j]]++;
	}

	string maps[] = { _DNA_MAP, _AA_MAP };
	unsigned long counts[2];
	for (unsigned int i = 0; i < 2; i++) {
		counts[i] = 0;
		string map = maps[i];
		for (unsigned j = 0; j < map.length(); j++)
			counts[i] += baseOccurences[map[j]];
	}

	if (Verbose)
		cout << counts[0] << " DNA characters and " << counts[1]
				<< " AA characters." << endl;
	if (counts[0] >= counts[1])
	{
		_dataType = _DNA_DATA;
		CharStates = 4;
	}
	else
	{
		_dataType = _AA_DATA;
		CharStates = 20;
	}

	cout << "The data appears to be " << dataTypeDesc[_dataType] << "." << endl;
}



void Alignment::readPhylip(string fileName) {
	ifstream fileReader;
	fileReader.open(fileName.c_str());
	if (!fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName);

	string str;
	safeGetline(fileReader, str);

	string whiteSpace = " \t";
	str = str.substr(str.find_first_not_of(whiteSpace));
	string rows = str.substr(0, str.find_first_of(whiteSpace));

	str = str.substr(str.find_first_of(whiteSpace));
	str = str.substr(str.find_first_not_of(whiteSpace));
	string cols = str.substr(0, str.find_first_of(whiteSpace));

	_cols = atoi(cols.c_str());
	_rows = atoi(rows.c_str());

	while (!fileReader.eof()) {
		safeGetline(fileReader, str);
		if (str.length()) {
			str = str.substr(str.find_first_not_of(whiteSpace));
			string name = str.substr(0, str.find_first_of(whiteSpace));
			str = str.substr(str.find_first_of(whiteSpace));
			string seq = str.substr(str.find_first_not_of(whiteSpace));

			if ((int) seq.length() < _cols)
				cerr << "Sequence #" << _sequences.size() + 1 << " (" << name
						<< ") has only " << seq.length() << " characters."
						<< endl;
			seq = adjustString(seq);
			if (!name.empty() && !seq.empty()) {
				_names.push_back(name);
				_sequences.push_back(seq);
			}
		}
	}
	if ((int) _sequences.size() < _rows)
		cerr << "The alignment has only " << _sequences.size() << " rows."
				<< endl;
}



void Alignment::readFasta(string fileName) {
	ifstream fileReader;
	fileReader.open(fileName.c_str());
	if (!fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName);

	string s;
	safeGetline(fileReader, s);
	while ((!fileReader.eof()) && s[0] != '>')
		safeGetline(fileReader, s);

	while (!fileReader.eof()) {
		string header;
		string seq;
		header = s;
		s.clear();
		while (!fileReader.eof() && s[0] != '>') {
			safeGetline(fileReader, s);
			if (s[0] != '>')
				seq += s;
		}
		seq = adjustString(seq, true);
		header = adjustString(header.substr(1), false);
		if (header.length() > 1 && seq.length()) {
			_names.push_back(header);
			_sequences.push_back(seq);
		}
	}
	_rows=_sequences.size();
	_cols=_sequences[0].length();
}



int Alignment::find(string name) {
	for (unsigned int i = 0; i < _names.size(); i++) {
		if (_names[i] == name)
			return i;
	}

	return -1;
}
