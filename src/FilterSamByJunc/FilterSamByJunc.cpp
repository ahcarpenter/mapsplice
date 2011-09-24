//#define VS

#include <iostream>
#include <vector>

#include <string>

#ifdef VS
#include <hash_map> //vc only
#include <hash_set> //vc only
#else
#include <ext/hash_map> //g++ only
#include <ext/hash_set> //g++ only
#endif

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>

using namespace std;
#ifdef VS
using namespace stdext;
#endif

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_set;
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

const size_t THIRTY_TWO = 32;
const size_t ALL_BITS_ON = static_cast<size_t>(-1);
const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

void read_juncbed(hash_map<string, hash_set<size_t> >& junc_bed_hash, string junction_file)
{
	ifstream junc_ifs(junction_file.c_str());

	if (junc_ifs.is_open())
	{
		string skipline;
		getline(junc_ifs,skipline);

		while (!junc_ifs.eof() )
		{
			string line;
			getline(junc_ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000];
			size_t juncst, juncend, prefixst, suffixend, kinds, hits, prefixlen, suffixlen;

			//chr1    14757   15039   JUNC_2  521     +       14829   14970   255,0,0 2 72,70   0,212			
			sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%llu\t%s\t%llu\t%llu\t%s\t%llu\t%llu,%llu\t%s", chromname,  &prefixst, &suffixend, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, &prefixlen, &suffixlen, blocksoffset);

			string chromstr = chromname;

			size_t comb_offset = ((prefixst + prefixlen) << THIRTY_TWO) + (suffixend - suffixlen + 1);

			junc_bed_hash[chromstr].insert(comb_offset);
		}
		junc_ifs.close();
	}
}

void
FilterSamByJunc(hash_map<string, hash_set<size_t> >& junc_bed_hash, vector<string>& m_mapreads_files, ofstream& remained_sam_ofs, ofstream& filtered_sam_ofs)
{
	for (size_t mi = 0; mi < m_mapreads_files.size(); ++mi)
	{
		ifstream ifs(m_mapreads_files[mi].c_str());		

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "")
					continue;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				size_t /*prim, */ start, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/, mate_offset;
					
				int mate_diff_dist;

				unsigned short mis_match;

				char mate_match;

//				TRAN00000074719:37	1	chr1	971958	0	21M84N54M	GCCAGTGGGGGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAG
//TRAN00000074719:100	1	chr1	972105	0	75M	CCACCTGTCGAGAGGGCTTCCTGCTACAACTCCGCGTTGGGCTGCTGCTCTGATGGGAAGACGCCCTCGCTGGAC
//TRAN00000074719:46	1	chr1	971967	0	12M84N63M	GGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAGAGGGCTTCC

				                                       //%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%llu\t%s\t%llu\t%llu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu", 
					readname, &strand_t, chromname, &start, &incorrect, spliceway, &mate_match, &mate_offset, &mate_diff_dist, chromseq, qualseq, &mis_match);

				//cout << line<<endl;
				vector<pair<size_t, int> > spliceway_vec;

				string splicewaystr = spliceway;
				size_t index = 0;

				bool exceed = false;

				bool junc_bed_matched = true, unspliced = true;

				string flag_str = " ";

				while (true)
				{
					if (index >= splicewaystr.length())
						break;

					int maplen;

					char flag;

					sscanf(splicewaystr.c_str() + index, "%d%c", &maplen, &flag);

					if (flag_str[0] == ' ')
					{
						if (flag == 'I')
							spliceway_vec.push_back(make_pair(start, -maplen));
						else if (flag == 'M')
							spliceway_vec.push_back(make_pair(start, maplen));
						else if (flag == 'N')
						{
							cout<<"start with N?"<<endl;
							spliceway_vec.push_back(make_pair(start + maplen, 0));
						}
					}
					else if (flag_str[0] == 'M')
					{
						if (flag == 'I')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
						else if (flag == 'M')
						{
							cout << "continue Ms?"<<endl;
							spliceway_vec.back().second += maplen;
						}
						else if (flag == 'N')
						{
							unspliced = false;
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
						}
					}
					else if (flag_str[0] == 'N')
					{
						if (flag == 'I')
							spliceway_vec.back().second = -maplen;
						else if (flag == 'M')
							spliceway_vec.back().second = maplen;
						else if (flag == 'N')
						{
							cout << "continue Ns?"<<endl;
							spliceway_vec.back().first += maplen;
						}
					}
					else if (flag_str[0] == 'I')
					{
						if (flag == 'I')
						{
							cout << "continue Is?"<<endl;
							spliceway_vec.back().second += -maplen;
						}
						else if (flag == 'M')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));
						else if (flag == 'N')
						{
							unspliced = false;
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
						}
					}

					flag_str[0] = flag;

					index = splicewaystr.find(flag_str, index) + 1;
				}

				if (spliceway_vec.size() > 1)
				{
					vector<pair<size_t, int> >::iterator vp_iter;

					for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end() - 1; ++vp_iter)
					{
						if (vp_iter->second < 0 || (vp_iter + 1)->second < 0)
						{
							continue;
						}

						size_t prefixend = vp_iter->first + vp_iter->second - 1;

						size_t suffixst = (vp_iter + 1)->first;// + spliceoutlen;

						size_t comb_offset = (prefixend << THIRTY_TWO) + suffixst;

						if (junc_bed_hash[chromname].find(comb_offset) == junc_bed_hash[chromname].end())
						{
							junc_bed_matched = false;
							break;
						}
					}
				}

				if (unspliced || junc_bed_matched)
					remained_sam_ofs << line << endl;
				else
					filtered_sam_ofs << line << endl;
			}
			ifs.close();
		}
		else cout << "Unable to open file";
	}
}

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		cout << "mapreads_file junction_file read_width min_intron max_intron "<<endl;
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}	

	string junc_bed_file = argv[1];

	string remained_sam = argv[2];

	string filtered_sam = argv[3];

	vector<string> comb_mapreads_files;

	cout << argv[4] << endl;
	for (int i = 4; i < argc; ++i)
		comb_mapreads_files.push_back(argv[i]);
	
	cout << junc_bed_file << endl;

	hash_map<string, hash_set<size_t> > junc_bed_hash;

	read_juncbed(junc_bed_hash, junc_bed_file);

	ofstream remained_sam_ofs(remained_sam.c_str());

	ofstream filtered_sam_ofs(filtered_sam.c_str());

	FilterSamByJunc(junc_bed_hash, comb_mapreads_files, remained_sam_ofs, filtered_sam_ofs);

	//Covert2JuncComb(juncfile.c_str(), comb_mapreads_files, read_width, chrom_dir, max_rank, min_intron, max_intron, min_anchor);

	//Covert2Junc(juncfile.c_str(), comb_mapreads_files, read_width, chrom_dir, max_rank, min_intron, max_intron, min_anchor);
}