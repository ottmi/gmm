#include "helper.h"
#include <fstream>



string adjustString(string s, bool upercase)
{
	string r = "";

	for (unsigned int i=0; i<s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != '\r' && c != ' ')
	  {
			if (upercase)
				r+= toupper(c);
			else
				r+= c;
	  }
	}

	return(r);
}



istream& safeGetline(istream& is, string& t)
{
	/* Courtesy of http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf */
    t.clear();
    istream::sentry se(is);
    streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        default:
            t += (char)c;
        }
    }
}

