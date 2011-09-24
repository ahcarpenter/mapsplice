#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
using namespace std;

int
main(int argc, const char** argv)
{
	
	ifstream ifs(argv[1]);

	ofstream ofs(argv[2]);

	if (ifs.is_open())
	{
		//string skipline;
		//getline(mapfile,skipline);

		string line;
		
		while (getline(ifs,line))
		{
			if (line.empty())
				continue;

			char readsname[1000];

			sscanf(line.c_str(), "%s\t", readsname);

			string readsnamestr = readsname;

			if (readsnamestr.find("~") != string::npos)
				line = line.substr(line.find("~") + 1);

			ofs << line<<endl/*<<'\t' <<strand<<'\t' << chromnamechr<<'\t' <<start<<'\t' << skip1<<'\t' <<spliceway<<'\t' <<skip2<<'\t' <<skip3<<'\t' <<skip4<<'\t' <<readseq<<'\t' <<qualseq<<"\tNM:i:" <<cur_mismatch*/;
			
		}
	}
}