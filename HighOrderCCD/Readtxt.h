#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>

using namespace std;

class Readtxt
{
	static void read_txt(const string& pathname, vector<vector<double>>& res)
	{
		res.clear();
		ifstream infile;
		infile.open(pathname.data());
		assert(infile.is_open());
		vector<double> suanz;
		string s;
		while (getline(infile, s))
		{
			istringstream is(s);
			double d;
			while (!is.eof())
			{
				is >> d;
				suanz.push_back(d);
			}
			res.push_back(suanz);
			suanz.clear();
			s.clear();
		}
		infile.close();
	}
};