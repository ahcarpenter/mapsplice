//#define VS
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string.h>

#ifdef VS
#include <hash_map> //vc only
#else
#include <ext/hash_map> //g++ only
#endif

using namespace std;
#ifdef VS
using namespace stdext;
#endif

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;

namespace __gnu_cxx
{
	template<typename Traits, typename Allocator>
	struct hash<std::basic_string<char, Traits, Allocator> >
	{
		size_t operator()(const std::basic_string<char, Traits, Allocator>& __s) const
		{
			return __stl_hash_string(__s.c_str());
		}
	};
}
#endif

bool paircomp ( const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs)
{
	return lhs.first < rhs.first;
}

int
readisland(const char* islands, hash_map<string, vector<pair<size_t, size_t> > >& islandmap)
{
	ifstream ifs(islands);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs, skipline);
		while (!ifs.eof() )
		{
			string line;
			getline(ifs, line);

			if (line == "")
				continue;

			char chrom[1000], group[1000], islandname[1000];
			size_t start, end;

			sscanf(line.c_str(), "%s\t%s\t%s\t%llu\t%llu", chrom, group, islandname, &start, &end);

			string chromstr = chrom;

			if (islandmap[chromstr].empty() || islandmap[chromstr].back().second + 1 < start)
				islandmap[chromstr].push_back(make_pair(start, end));
			else if(islandmap[chromstr].back().first > start)
			{
				cout << islandmap[chromstr].back().first << '\t' << islandmap[chromstr].back().second<<endl;
				cout << line << endl;
				//getchar();
				cout << "islands not ordered!"<<endl;
				exit(0);
			}
			else if(islandmap[chromstr].back().second < end)
				islandmap[chromstr].back().second = end;
			else
			{
			}
		}

		ifs.close();
	}

	hash_map<string, vector<pair<size_t, size_t> > >::iterator m_iter;

	for (m_iter = islandmap.begin(); m_iter != islandmap.end(); ++m_iter)
	{
		sort(m_iter->second.begin(), m_iter->second.end(), paircomp);
	}

	return 0;
}


int
filternotinisland(hash_map<string, vector<pair<size_t, size_t> > > &islandmap, const char* BWTfile, const char* notinisland, const char* inisland)
{
	ifstream ifs(BWTfile);

	ofstream ofs(notinisland);

	ofstream inofs(inisland);

	if (ifs.is_open())
	{
		int ncount = 0;
		int icount = 0;

		while (!ifs.eof() )
		{
			string line;
			getline(ifs, line);

			if (line.empty())
				continue;

			char chrom[1000], readname[1000], strand, readseq[1000];
			size_t start, end;

			sscanf(line.c_str(), "%s\t%c\t%s\t%llu\t%s", readname, &strand, chrom, &start, readseq);

			if (islandmap.find(chrom) == islandmap.end())
			{
				//cout << "can't find chromosome "<<endl;
				//cout << chrom << endl;
				//cout << line << endl;
				inofs<<line<<endl;
				continue;
			}

			vector<pair<size_t, size_t> >::iterator vpit;

			size_t read_len = strlen(readseq);

			end = start + read_len;

			string chromstr = chrom;

			//1~3~HWI-EAS217:4:1:11:986#0/1   +       chr1    157043645       GATGAAGCTTACATGAACAAGGTAN       `aab_`_aaaaa_a_^aaBBBBBBB       2       24:G>N

			pair<size_t, size_t> current_island = make_pair(start, end);

			vector<pair<size_t, size_t> >::iterator lowerbound = lower_bound(islandmap[chromstr].begin(), islandmap[chromstr].end(), current_island, paircomp);

			if (lowerbound != islandmap[chromstr].begin())
				--lowerbound;

			for (vpit = lowerbound; vpit != islandmap[chromstr].end(); ++vpit)
			{
				if (start + read_len >= vpit->first && start <= vpit->second + read_len)
					break;
				else if (start < vpit->first)
				{
					vpit = islandmap[chromstr].end();
					break;
				}
			}

			if (vpit == islandmap[chromstr].end())
			{
				ofs<<line<<endl;
				++ncount;
				continue;
			}

			inofs<<line<<endl;
			++icount;
		}

		ifs.close();
		ofs.close();
		inofs.close();
		cout << "not in island"<<ncount<<endl;
		cout << "in island"<<icount<<endl;
	}
	return 0;
}

int
main(int argc, const char** argv)
{
	if (argc < 5)
	{
		cout << " need  three files\n";
		return 0;
	}

	const char* islands = argv[1];
	const char* BWAfile = argv[2];
	const char* notinisland = argv[3];
	const char* inisland = argv[4];

	hash_map<string, vector<pair<size_t, size_t> > > islandmap;
	
	readisland(islands, islandmap);

	filternotinisland(islandmap, BWAfile, notinisland, inisland);

}