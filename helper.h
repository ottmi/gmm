#ifndef HELPER_H_
#define HELPER_H_

#include <string>
#include <iostream>
using namespace std;

string adjustString(string s, bool upercase=false);
istream& safeGetline(istream& is, string& t);

unsigned int mapDNAToNum(char c);

#endif /* HELPER_H_ */
