#include "helper.h"
#include "globals.h"
#include <fstream>
#include <ctime>
#include <iomanip>

string adjustString(string s, bool upercase)
{
	string r = "";

	for (unsigned int i = 0; i < s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != '\r' && c != ' ')
		{
			if (upercase)
				r += toupper(c);
			else
				r += c;
		}
	}

	return (r);
}

istream& safeGetline(istream& is, string& t)
{
	/* Courtesy of http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf */
	t.clear();
	istream::sentry se(is);
	streambuf* sb = is.rdbuf();

	for (;;)
	{
		int c = sb->sbumpc();
		switch (c)
		{
			case '\r':
				c = sb->sgetc();
				if (c == '\n') sb->sbumpc();
				return is;
			case '\n':
			case EOF:
				return is;
			default:
				t += (char) c;
		}
	}
}

unsigned int mapDNAToNum(string s)
{
	unsigned int d = 0;
	for (unsigned int i = 0; i < s.size(); i++)
	{
		d = d << 2;
		switch (s[i])
		{
			case 'A':
			case 'a':
				d += 0x00;
				break;

			case 'C':
			case 'c':
				d += 0x01;
				break;

			case 'G':
			case 'g':
				d += 0x02;
				break;

			case 'T':
			case 't':
				d += 0x03;
				break;

			default:
				throw("Unsupported character: " + s[i]);
				break;
		}
	}

	return d;
}

string mapNumToDNA(unsigned int val, unsigned int len)
{
	string map = _DNA_MAP;
	string s;
	for (unsigned int i = 0; i < len; i++)
	{
		s = map[val & 0x3] + s;
		val = val >> 2;
	}

	return s;
}

unsigned int mapAAToNum(char c)
{
	unsigned int d = 0;
	switch (c)
	{
		case 'A':
		case 'a':
			d = 0x00;
			break;

		case 'C':
		case 'c':
			d = 0x01;
			break;

		case 'D':
		case 'd':
			d = 0x02;
			break;

		case 'E':
		case 'e':
			d = 0x03;
			break;

		case 'F':
		case 'f':
			d = 0x04;
			break;

		case 'G':
		case 'g':
			d = 0x05;
			break;

		case 'H':
		case 'h':
			d = 0x06;
			break;

		case 'I':
		case 'i':
			d = 0x07;
			break;

		case 'K':
		case 'k':
			d = 0x08;
			break;

		case 'L':
		case 'l':
			d = 0x09;
			break;

		case 'M':
		case 'm':
			d = 0x0A;
			break;

		case 'N':
		case 'n':
			d = 0x0B;
			break;

		case 'P':
		case 'p':
			d = 0x0C;
			break;

		case 'Q':
		case 'q':
			d = 0x0D;
			break;

		case 'R':
		case 'r':
			d = 0x0E;
			break;

		case 'S':
		case 's':
			d = 0x0F;
			break;

		case 'T':
		case 't':
			d = 0x10;
			break;

		case 'V':
		case 'v':
			d = 0x11;
			break;

		case 'W':
		case 'w':
			d = 0x12;
			break;

		case 'Y':
		case 'y':
			d = 0x13;
			break;

		default:
			string s = "Unsupported character: " + c;
			throw(s);
			break;
	}
	return d;
}

string printTime(time_t t)
{
	stringstream s;
	s << t / 3600 << ":" << setfill('0') << setw(2);
	t = t % 3600;
	s << t / 60 << ":" << setfill('0') << setw(2) << t % 60;

	return s.str();
}
