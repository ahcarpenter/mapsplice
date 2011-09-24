/*    
 *    fusionsam2junc.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//#define VS

#include <iostream>
#include <vector>

#include <string>

#ifdef VS
#include <hash_map> //vc only
#else
#include <ext/hash_map> //g++ only
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
using __gnu_cxx::hash_map;

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

struct JuncInfo{
	JuncInfo(/*int pm, *//*const string fs, */int loc, int suffix_len, size_t rw, size_t tagidx) : /*prim(pm), *//*flankstr(fs),*/ p(rw - 1, 0) 
	{
		++p[loc-1];
		max_prefix_len = loc;
		max_suffix_len = suffix_len;

		m[tagidx] = 1;
	}
	bool inc_hits(int idx, int suffix_len, size_t tagidx)
	{
		if (tagidx != -1 && m.find(tagidx) != m.end())
			return false;

		++p[idx-1];

		m[tagidx] = 1;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;
		return true;
	}
	//int prim;
	//string flankstr;
	size_t max_prefix_len;
	size_t max_suffix_len;
	vector<unsigned short> p;

	map<size_t, int> m;
};

typedef hash_map<size_t, JuncInfo> JUNC_SUFFIX;
typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;

typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
typedef JUNC_HASH::iterator JUNC_HASH_ITER;
typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef map<string, JUNC_HASH > CONJ_HASH_MAP;

struct JuncInfoFusion{
	JuncInfoFusion( size_t loc, size_t suffix_len, size_t rw, const string& ss, char sd1, char sd2, unsigned short mis, const string& fl) : p(rw - 1, 0), max_prefix_len(loc), 
		max_suffix_len(suffix_len),	synstr(ss), strand1(sd1), strand2(sd2), min_mismatch(mis), max_mismatch(mis), sum_mismatch(mis), flankstr(fl)
	{
		prim = 0;

		if (flankstr == "ATAC")
			prim = 1;
		else if (flankstr == "CTAC")
			prim = 6;
		else if (flankstr == "CTGC")
			prim = 3;
		else if (flankstr == "GCAG")
			prim = 4;
		else if (flankstr == "GTAG")
			prim = 5;
		else if (flankstr == "GTAT")
			prim = 2;

		min_mismatch = mis;

		max_mismatch = mis;

		sum_mismatch = mis;

		++p[loc - 1];
	}

	bool inc_hits(int idx, size_t suffix_len, unsigned short mis)
	{
		++p[idx - 1];

		sum_mismatch += mis;

		if (mis > max_mismatch)
			max_mismatch = mis;

		if (mis < min_mismatch)
			min_mismatch = mis;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;

		return true;
	}
	size_t max_prefix_len, max_suffix_len;
	int prim;
	string flankstr;
	string synstr;
	char strand1, strand2;
	vector<unsigned short> p;
	size_t min_mismatch;
	size_t max_mismatch;
	size_t sum_mismatch;	
};

static const size_t ALL_BITS_ON = static_cast<size_t>(-1);
static const size_t THIRTY_TWO = 32;
//static const size_t ALL_BITS_ON = static_cast<size_t>(-1);
static const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
static const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

typedef hash_map<size_t, JuncInfoFusion> JUNCFusion_SUFFIX;
typedef JUNCFusion_SUFFIX::iterator JUNCFusion_SUFFIX_ITER;
typedef JUNCFusion_SUFFIX::const_iterator JUNCFusion_SUFFIX_CITER;

typedef hash_map<size_t, JUNCFusion_SUFFIX > JUNCFusion_HASH;
typedef JUNCFusion_HASH::iterator JUNCFusion_HASH_ITER;
typedef JUNCFusion_HASH::const_iterator JUNCFusion_HASH_CITER;

typedef map<string, JUNCFusion_HASH > CONJFusion_HASH_MAP;

typedef hash_map<size_t, JuncInfoFusion> FUSION_JUNC_HASH_COMB;

typedef FUSION_JUNC_HASH_COMB::iterator FUSION_JUNC_HASH_COMB_ITER ;

typedef hash_map<string, FUSION_JUNC_HASH_COMB> CHROM_FUSION_JUNC_HASH_COMB;

typedef CHROM_FUSION_JUNC_HASH_COMB::iterator CHROM_FUSION_JUNC_HASH_COMB_ITER;

struct JuncForSortFusion{
	int juncst;
	int juncend;
	int hits;
	int kinds;
	string blocks;
	string blocksoffset;
	string rank;
	string flankcase;
	string flankchr;
	string strand;
	string synstr;

	unsigned short min_mismatch;
	unsigned short max_mismatch;
	double ave_mismatch;

	JuncForSortFusion(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& fc, const string& fcr, const string& sd, const string& st, unsigned short min_mis, unsigned short max_mis, double ave_mis) : 
	    juncst(jst), juncend(jend), hits(hts), kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), flankcase(fc), flankchr(fcr), strand(sd), synstr(st),
		min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis){}
};

bool
compjuncfusion(const JuncForSortFusion& lhs, const JuncForSortFusion& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

struct CoverageBlock{
	int blockst;
	int blockend;
	int hits;
	CoverageBlock(const int& blst, const int& blend, const int& ht) : blockst(blst), blockend(blend), hits(ht) {}
};

void
SeprateMapreads(vector<string>& comb_mapreads_files, vector<string>& m_mapreads_files, string juncfile)
{
	vector<string>::iterator cm_iter;
	string prevchrom = "";
	ofstream ofs;
	map<string, int> mapreads_filemap;

	string seprate_path = juncfile;
	//seprate_path.append(comb_mapreads_files.front());

	seprate_path.append(".sepdir");

	string mkdir_cmdstr = "mkdir ";
	mkdir_cmdstr.append(seprate_path);
	system(mkdir_cmdstr.c_str());

	seprate_path.append("/");

	for (cm_iter = comb_mapreads_files.begin(); cm_iter != comb_mapreads_files.end(); ++cm_iter)
	{
		ifstream ifs((*cm_iter).c_str());

		if (ifs.is_open())
		{
			string line;
			while (getline(ifs,line))
			{
				if (line == "")
					continue;
				char /*chromname[1000], */readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
				char chromname1[1000], chromname2[1000];
				char strand = '+';
				int prefixst, strand_t, incorrect, prefixend;
//HWI-EAS217:4:100:1001:1129#0/1_6458158	chr11	+	10486490	10486546	chr3	+	97819005	97819047	CGGGGTAAGATTTGCCGAGTTCCTTTTACTTTTTTTAACCTTTCCTTGTGGGCATGCCTGTGTTGGGTTAA
//CAGTGAGGGTAATAATGGCTTGTTGGTTG	CTGC	NM:i:4
				sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d", 
					readname, chromname1, &strand, &prefixst, &prefixend, chromname2, &strand, &prefixst, &prefixend);

				string chromname = chromname1;

				chromname.append(chromname2);
				if (prevchrom != string(chromname))
				{
					char mapreads_file[1000];

					sprintf(mapreads_file, "%smapreads_%s.txt", seprate_path.c_str(), chromname.c_str());

					if (ofs.is_open())
						ofs.close();

					if (mapreads_filemap.find(mapreads_file) != mapreads_filemap.end())
						ofs.open(mapreads_file, fstream::app);
					else
						ofs.open(mapreads_file);

					prevchrom = chromname;

					mapreads_filemap[mapreads_file] = 1;

				}

				ofs<< line<<endl;
			}

			ifs.close();
		}
	}

	map<string, int>::iterator mf_iter;
	for (mf_iter = mapreads_filemap.begin(); mf_iter != mapreads_filemap.end(); ++mf_iter)
	{
		m_mapreads_files.push_back(mf_iter->first);
	}
}

struct JuncForSort{
	int juncst;
	int juncend;
	int hits;
	int kinds;
	string blocks;
	string blocksoffset;
	string rank;
	string lpq;

	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& l) : juncst(jst), juncend(jend), hits(hts), 
		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l) {}

};

bool
compjunc(const JuncForSort& lhs, const JuncForSort& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

void
readchrom(const char* filename, string& longseq)
{
	long size;  

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);

		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";
}

void
SortJunc(const char* juncfile, string chrom_dir, size_t min_intron, size_t max_intron)
{
	map<string, vector<JuncForSort> > m_p;

    ifstream ifs(juncfile);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], lpq[100];
			int juncst, juncend, prefixend, suffixst, kinds, hits;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, lpq);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset,	rank, lpq));
		}

		ifs.close();
	}
	map<string, vector<JuncForSort> >::iterator mpit;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		sort(mpit->second.begin(), mpit->second.end(), compjunc);
	}

	ofstream ofs(juncfile);

	string juncidstr = "JUNC_";
	size_t juncid = 1;

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		string chromfile = chrom_dir + mpit->first;

		chromfile.append(".fa");

		string chromseq;

		readchrom(chromfile.c_str(), chromseq);

		if (chromseq == "")
		{
			cout <<"empty chrom: "<<chromfile<<endl;
			exit(1);
		}

		vector<JuncForSort>::iterator vit;
		for (vit = mpit->second.begin(); vit != mpit->second.end(); ++vit)
		{
			string flankstr = chromseq.substr(vit->juncst, 2) + chromseq.substr(vit->juncend - 3, 2);

			for (size_t i = 0; i < flankstr.length(); ++i)
			{
				if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
					flankstr[i] = flankstr[i] + 'A' - 'a';
			}

			int flankcase = 0;

			if (flankstr == "ATAC")
				flankcase = 1;
			else if (flankstr == "CTAC")
				flankcase = 6;
			else if (flankstr == "CTGC")
				 flankcase = 3;
			else if (flankstr == "GCAG")
				 flankcase = 4;
			else if (flankstr == "GTAG")
				 flankcase = 5;
			else if (flankstr == "GTAT")
				 flankcase = 2;

			double il = 1.0 - (((double) (vit->juncend - vit->juncst - min_intron + 1)) / (double (max_intron - min_intron + 2))); 
			ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<'+'<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< flankcase <<'\t'<<flankstr <<'\t'<<il<<'\t'<<vit->lpq<<endl;
			++juncid;
		}
	} 
}

void
WriteJunc(const CONJ_HASH_MAP& conj_hash_map, ofstream& ofs, size_t m_read_width, double m_max_rank, string chrom_dir, size_t max_intron)
{
	
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CONJ_HASH_MAP::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		JUNC_HASH_CITER iter_conj = chm_iter->second.begin();

		string chromfile = chrom_dir + chm_iter->first;

		chromfile.append(".fa");

		string chromseq;

		readchrom(chromfile.c_str(), chromseq);

		if (chromseq == "")
		{
			cout <<"empty chrom: "<<chromfile<<endl;
			exit(1);
		}

		size_t chrom_size = chromseq.size() - 1;

		int j=0;
		while (iter_conj != chm_iter->second.end())
		{
			JUNC_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			while (iter_conj_pre!=iter_conj->second.end())
			{
				int kinds = 0;
				int hits = 0;
				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						++kinds;
						hits += iter_conj_pre->second.p[i];
					}
				}

				double rank = 0;

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						double pi = (double)iter_conj_pre->second.p[i] / (double)hits;
						rank += pi * log(pi);
					}
				}

				if (rank!=0)
					rank = -rank;

				if (rank < m_max_rank)
				{
					iter_conj_pre++;
					continue;
				}

				ofs <<chm_iter->first <<'\t'<< iter_conj_pre->first << '\t' << iter_conj->first<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;
				
				//int juncidlen = (int) log10((double)juncid);

				//string suffix0(8 - juncidlen

				//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';

				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
				//		ofs << k << ',';
				//		break;
				//	}
				//}
				ofs << iter_conj_pre->second.max_prefix_len << ',';

				//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				//{
				//	if (iter_conj_pre->second.p[i] > 0)
				//	{
				//		ofs << m_read_width - i - 1 << ',';
				//		break;
				//	}
				//}
				ofs << iter_conj_pre->second.max_suffix_len << ',';
				ofs << '\t';

				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
				//		ofs << 0 << ',';
				//		break;
				//	}
				//}
				ofs << 0 << ',';

				//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				//{
				//	if (iter_conj_pre->second.p[i] > 0)
				//	{
				//		ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
				//		break;
				//	}
				//}
				ofs << int (iter_conj->first + iter_conj_pre->second.max_suffix_len - iter_conj_pre->first + 1 )<<',';
				ofs << '\t';

				//double ppower = pow(double(iter_conj_pre->second.max_prefix_len), 4);

				////double pNpower = pow(1.0 - ppower, (double)chrom_size);

				//double qpower = pow(double(iter_conj_pre->second.max_suffix_len), 4);

				////double pDpower = pow(1.0 - qpower, (double)max_intron);

				//double lpq = ((double)chrom_size / ppower) + ((double)max_intron / qpower);

				double ppower = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower = pow(1.0 - ppower, (double)chrom_size);

				double qpower = pow(0.25, double(iter_conj_pre->second.max_suffix_len));

				double pDpower = pow(1.0 - qpower, (double)max_intron);

				double lpq = 1.0 - (pNpower * pDpower);

				double ppower2 = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower2 = pow(1.0 - ppower2, (double)max_intron );

				double qpower2 = pow(0.25, double(iter_conj_pre->second.max_suffix_len));

				double pDpower2 = pow(1.0 - qpower2, (double)chrom_size);

				double lpq2 = 1.0 - (pNpower2 * pDpower2);

				double lpqave = 1.0 - (lpq + lpq2) / 2;

				//cout <<"max_prefix_len\t" << iter_conj_pre->second.max_prefix_len << endl;

				//cout <<"ppower\t" << ppower << endl;

				//cout <<"chrom_size\t" << chrom_size << endl;

				//cout <<"pNpower\t" << pNpower << endl;

				//cout <<"max_suffix_len\t" << iter_conj_pre->second.max_suffix_len << endl;

				//cout <<"qpower\t" << qpower << endl;

				//cout <<"max_intron\t" << max_intron << endl;

				//cout <<"pDpower\t" << pDpower << endl;

				//cout <<"lpq\t" << lpq << endl;

				//getchar();

				ofs << rank << '\t'<< lpqave/*iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr */<<endl;

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}


void
Covert2Junc(const char* junc_filename, vector<string>& m_mapreads_files, size_t m_read_width, string chrom_dir, double m_max_rank, size_t min_intron, size_t max_intron)
{
	//char junc_filename[1000];
	//char wig_filename[1000];

	//sprintf(junc_filename, "%sjunctions.txt", m_path.c_str());
	ofstream ofs(junc_filename);

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	//sprintf(wig_filename, "%scoverage.wig.txt", m_path.c_str());
	//ofstream ofs_wig(wig_filename);

	//string headline2 = "track type=bedGraph name=\"Mapsplice - read coverage\"";

	//ofs_wig << headline2 << endl;

	for (size_t mi = 0; mi < m_mapreads_files.size(); ++mi)
	{
		ifstream ifs(m_mapreads_files[mi].c_str());

		CONJ_HASH_MAP conj_hash_map;

		//VEC_COVER_BLOCK v_coverage_block;

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "")
					continue;

				//cout <<line<<endl;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				int /*prim, */ prefixst, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/;
				//sscanf(line.c_str(), "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d", chromname, readname, &prim, flankseq, strand, readseq, prefixseq, suffixseq, &prefixlen, &prefixend, &suffixst);

//				TRAN00000074719:37	1	chr1	971958	0	21M84N54M	GCCAGTGGGGGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAG
//TRAN00000074719:100	1	chr1	972105	0	75M	CCACCTGTCGAGAGGGCTTCCTGCTACAACTCCGCGTTGGGCTGCTGCTCTGATGGGAAGACGCCCTCGCTGGAC
//TRAN00000074719:46	1	chr1	971967	0	12M84N63M	GGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAGAGGGCTTCC

				                                       //%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s"/*\t*\t0\t0\t%s\t%s"*/, 
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway/*, &prefixlen, &spliceoutlen, &suffixlen,*//* chromseq, qualseq*/);

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_last_of("_");
				size_t tagidx = -1;

				if (last_idx != string::npos)
				{
					string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
					tagidx = atoi(tagidxstr.c_str()) - 1;
				}

				//cout << line<<endl;
				vector<pair<size_t, size_t> > spliceway_vec;

				string splicewaystr = spliceway;
				size_t index = 0;

				bool exceed = false;
				while (true)
				{
					int maplen, intron;

					if (index == 0)
					{
						sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
						spliceway_vec.push_back(make_pair(prefixst, maplen));

						if (maplen <=0)
							exceed =true;

						//if (splicewaystr.find("M", index) == string::npos)
						//	break;

						index = splicewaystr.find("M", index) + 1;

						//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

						//getchar();

						
						//if (index == splicewaystr.length())
						//	break;
						//if (index == string::npos)
						//	break;
					}
					else
					{
						if (index == splicewaystr.length())
							break;
						sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);

						if (intron != 0 && maplen != 0)
						{
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));
						}

						if (maplen <=0)
							exceed =true;

						if (intron < min_intron || intron > max_intron)
							exceed = true;
			
						if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
						{
							//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;
							//if (splicewaystr.find("M", index) == splicewaystr.length() - 1)
							//	getchar();
							break;
						}

						index = splicewaystr.find("M", index) + 1;

						//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

						//getchar();
						
					}
				}

				if (spliceway_vec.size() > 1 && !exceed)
				{
					vector<pair<size_t, size_t> >::iterator vp_iter;
					for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end() - 1; ++vp_iter)
					{
						//suffixend = prefixst + spliceoutlen + (int)m_read_width;

						prefixst = vp_iter->first;
						size_t prefixend = vp_iter->first + vp_iter->second - 1;

						size_t suffixst = (vp_iter + 1)->first;// + spliceoutlen;

						size_t prefixlen = vp_iter->second;

						size_t suffixlen = (vp_iter + 1)->second;

						//string flankstr = flankseq;

						//coverage block
						//v_coverage_block.push_back(CoverageBlock(prefixst + 1, prefixst + prefixlen, 1));
						//v_coverage_block.push_back(CoverageBlock(suffixst, suffixlen, 1));

						string chrom_ID = chromname;

						//if (chrom_ID == "chr1" &&  prefixend == 4862 && suffixst == 5002)
						//{
						//	cout <<line<<endl;
						//	for (size_t i = 0; i < spliceway_vec.size(); ++i)
						//		cout << spliceway_vec[i].first<<'\t'<<spliceway_vec[i].second<<endl;
						//	getchar();
						//}
						JUNC_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

						if (iter_conj == conj_hash_map[chrom_ID].end())
						{
							JUNC_SUFFIX match_hash_temp;
							conj_hash_map[chrom_ID].insert(JUNC_HASH::value_type(suffixst, match_hash_temp));
							iter_conj = conj_hash_map[chrom_ID].find(suffixst);
						}

						JUNC_SUFFIX_ITER iter_conj_pre; 

						if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
						{
							//if (iter_conj_pre->second.flankstr != flankstr)
							//{
							//	cout <<"flank does not match"<<endl;
							//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
							//}
							iter_conj_pre->second.inc_hits(prefixlen, suffixlen, tagidx);
						}	  
						else 
						{
							//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
							(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(/*prim, *//*flankseq, */prefixlen, suffixlen, m_read_width, tagidx)));
						}
					}
				}
			}
			ifs.close();
		}
		else cout << "Unable to open file";

		cout <<"write junc"<<endl;
		WriteJunc(conj_hash_map, ofs, m_read_width, m_max_rank, chrom_dir, max_intron);

		//cout <<"sort v_coverage_block"<<endl;
		//sort(v_coverage_block.begin(), v_coverage_block.end(), compblk);

		//cout <<"WriteCoverage"<<endl;
		//WriteCoverage(v_coverage_block, conj_hash_map.begin()->first, ofs_wig, mi);
	}

	ofs.close();

	cout <<"sort junc"<<endl;

	SortJunc(junc_filename, chrom_dir, min_intron, max_intron);
}

void
WriteJuncFusion(const CONJFusion_HASH_MAP& conj_hash_map, ofstream& ofs,  size_t m_read_width, double m_max_rank)
{
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CONJFusion_HASH_MAP::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		JUNCFusion_HASH_CITER iter_conj = chm_iter->second.begin();

		int j=0;
		while (iter_conj != chm_iter->second.end())
		{
			JUNCFusion_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			while (iter_conj_pre!=iter_conj->second.end())
			{
				int kinds = 0;
				int hits = 0;
				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						++kinds;
						hits += iter_conj_pre->second.p[i];
					}
				}

				double rank = 0;

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						double pi = (double)iter_conj_pre->second.p[i] / (double)hits;
						rank += pi * log(pi);
					}
				}

				if (rank!=0)
					rank = -rank;

				if (rank < m_max_rank)
				{
					iter_conj_pre++;
					continue;
				}

				ofs <<chm_iter->first <<'\t'<< iter_conj_pre->first << '\t' << iter_conj->first<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;
				
				//int juncidlen = (int) log10((double)juncid);

				//string suffix0(8 - juncidlen

				//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';

				for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				{
					if (iter_conj_pre->second.p[k - 1] > 0)
					{
						ofs << k << ',';
						break;
					}
				}

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						ofs << m_read_width - i - 1 << ',';
						break;
					}
				}
				ofs << '\t';

				for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				{
					if (iter_conj_pre->second.p[k - 1] > 0)
					{
						ofs << 0 << ',';
						break;
					}
				}

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
						break;
					}
				}
				ofs << '\t';

				ofs << rank << '\t'<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'<< iter_conj_pre->second.strand1<<iter_conj_pre->second.strand2 <<'\t'<< iter_conj_pre->second.synstr 
					<<'\t'<<iter_conj_pre->second.min_mismatch << '\t'<<iter_conj_pre->second.max_mismatch << '\t'<<(double)iter_conj_pre->second.sum_mismatch / (double)hits <<endl;

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}

void
WriteJuncFusion(const CHROM_FUSION_JUNC_HASH_COMB& conj_hash_map, ofstream& ofs, size_t read_width, double m_max_rank)
{
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CHROM_FUSION_JUNC_HASH_COMB::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		FUSION_JUNC_HASH_COMB::const_iterator iter_conj;// = chm_iter->second.begin();

		int j=0;
		for (iter_conj = chm_iter->second.begin();iter_conj != chm_iter->second.end(); ++iter_conj)
		{
			size_t comb_offset = iter_conj->first;

			size_t prefix_end = comb_offset >> THIRTY_TWO;

			size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

			//JUNCFusion_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			//while (iter_conj_pre!=iter_conj->second.end())
			//{
			int kinds = 0;
			int hits = 0;
			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					++kinds;
					hits += iter_conj->second.p[i];
				}
			}

			double rank = 0;

			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					double pi = (double)iter_conj->second.p[i] / (double)hits;
					rank += pi * log(pi);
				}
			}

			if (rank!=0)
				rank = -rank;

			if (rank < m_max_rank)
			{
				//iter_conj_pre++;
				continue;
			}

			ofs <<chm_iter->first <<'\t'<< prefix_end << '\t' << suffix_st<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;

			//int juncidlen = (int) log10((double)juncid);

			//string suffix0(8 - juncidlen

			//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
			ofs << hits << "\t"<<iter_conj->second.strand1 <<iter_conj->second.strand2 <</*"\t"<<prefix_end << '\t' << suffix_st<<*/ "\t255,0,0\t"<<2<<'\t' << 
				iter_conj->second.max_prefix_len << ',' << iter_conj->second.max_suffix_len <<"\t0," << static_cast<int> (suffix_st - prefix_end)<<'\t';

			//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
			//{
			//	if (iter_conj_pre->second.p[k - 1] > 0)
			//	{
			//		ofs << k << ',';
			//		break;
			//	}
			//}

			//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
			//{
			//	if (iter_conj_pre->second.p[i] > 0)
			//	{
			//		ofs << m_read_width - i - 1 << ',';
			//		break;
			//	}
			//}
			//ofs << '\t';

			//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
			//{
			//	if (iter_conj_pre->second.p[k - 1] > 0)
			//	{
			//		ofs << 0 << ',';
			//		break;
			//	}
			//}

			//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
			//{
			//	if (iter_conj_pre->second.p[i] > 0)
			//	{
			//		ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
			//		break;
			//	}
			//}
			//ofs << '\t';

			//string flankstr = FlankString(iter_conj->second.flank_bits, iter_conj->second.bads);

			//int flankcase = 0;

			//if (flankstr == "ATAC")
			//	flankcase = 1;
			//else if (flankstr == "CTAC")
			//	flankcase = 2;
			//else if (flankstr == "CTGC")
			//	 flankcase = 3;
			//else if (flankstr == "GCAG")
			//	 flankcase = 4;
			//else if (flankstr == "GTAG")
			//	 flankcase = 5;
			//else if (flankstr == "GTAT")
			//	 flankcase = 6;

			ofs << rank << '\t'<< iter_conj->second.prim<<'\t'<< iter_conj->second.flankstr <<'\t'/*<< iter_conj_pre->second.strand <<'\t'*/<< iter_conj->second.synstr 
				<<'\t'<<iter_conj->second.min_mismatch<<'\t'<<iter_conj->second.max_mismatch<<'\t'
				<<(double)iter_conj->second.sum_mismatch / (double)hits<<endl;

			++juncid;
			j++;
				//iter_conj_pre++;
			//}

			//iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}

void 
SortJuncFusion(const char* juncfile)
{
	map<string, vector<JuncForSortFusion> > m_p;

    ifstream ifs(juncfile);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], flankcase[100], flankchr[100], strand12[10], synstr[1000];
			int juncst, juncend, prefixend, suffixst, kinds, hits;

			unsigned short min_mis, max_mis;

			double ave_mis;

//			chr2_chr10      170331943       27842899        JUNC_1  32      ++      255,0,0 2       58      74
//2.92163 0       TGCA    CCCCTGCTCTTGAGGTTCACAGATGATACTTTTGATCCAGAACTTGCAGCAACAATAGCAGTTGGGGAGGAGATGACTGC
//ATTGTGAAAAACCTCCAGACTGCCTGTCTGCATTGGCTGT        0       3       1.8125

			//sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
			//	juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, flankcase, flankchr, strand12, synstr, &min_mis, &max_mis, &ave_mis);

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, rgb, &kinds, blocks, blocksoffset, rank, flankcase, flankchr, synstr, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSortFusion(prefixend, suffixst, hits, kinds, blocks, blocksoffset, rank, flankcase, flankchr, strand, synstr, min_mis, max_mis, ave_mis));
		}

		ifs.close();
	}
	map<string, vector<JuncForSortFusion> >::iterator mpit;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		sort(mpit->second.begin(), mpit->second.end(), compjuncfusion);
	}

	ofstream ofs(juncfile);

	string juncidstr = "JUNC_";
	size_t juncid = 1;

	string headline = "track name=junctions description=\"Mapsplice fusion junctions\"";

	ofs << headline << endl;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		vector<JuncForSortFusion>::iterator vit;
		for (vit = mpit->second.begin(); vit != mpit->second.end(); ++vit)
		{
			ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<vit->strand <</*'\t'<<vit->juncst<< '\t'<<vit->juncend<<*/"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank <<'\t'<< vit->flankcase<<'\t'<< vit->flankchr<</*'\t'<< vit->strand<<*/'\t'<< vit->synstr
				<<'\t'<<vit->min_mismatch << '\t'<< vit->max_mismatch<< '\t'<< vit->ave_mismatch<<endl;
			++juncid;
		}
	} 
}


void
CovertFusion2Junc(const char* junc_filename, vector<string>& m_mapreads_files, size_t m_read_width, double m_max_rank)
{
	ofstream ofs(junc_filename);

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	for (size_t mi = 0; mi < m_mapreads_files.size(); ++mi)
	{
		ifstream ifs(m_mapreads_files[mi].c_str());

		CONJFusion_HASH_MAP conj_hash_map;

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "")
					continue;

				//cout <<line<<endl;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000], synchr[1000], flankchr[100];
				char chromname1[1000], chromname2[1000];
				char strand = '+', strand1, strand2;
				int /*prim, */ prefixst, prefixend, suffixst, suffixend, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/;
				//sscanf(line.c_str(), "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d", chromname, readname, &prim, flankseq, strand, readseq, prefixseq, suffixseq, &prefixlen, &prefixend, &suffixst);
				unsigned short mis_match;

//				TRAN00000074719:37	1	chr1	971958	0	21M84N54M	GCCAGTGGGGGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAG
//TRAN00000074719:100	1	chr1	972105	0	75M	CCACCTGTCGAGAGGGCTTCCTGCTACAACTCCGCGTTGGGCTGCTGCTCTGATGGGAAGACGCCCTCGCTGGAC
//TRAN00000074719:46	1	chr1	971967	0	12M84N63M	GGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAGAGGGCTTCC

//HWI-EAS217:4:52:346:1875#0/1_3287705	chr1	-	18451	18525	chr1	+	88603914	88603938	GCGATACACATGTACAAGTTTTGTCACTTAAAAAGAATTTGAAAAAACTTCATAGATGCAAAAAAAAAAAAAAACAAAAACCACCTA
//TTTTTAAAGATTA	TGTT	NM:i:2

//1012966~TRAN00000110580:1400:334        chr10   +       66128298        66128322        chr1    +       241548206       24154
//8255    tgaggcttatcctaggtaaactaacgcaagaacagaaaaccaaatactgcatattcttacgtataagtgggagctaaatgataagaacttatgagcacaaagtaggaaaccacagac
//agt     TTAC    NM:i:2
				                                       //%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\tNM:i:%hu"/*\t*\t0\t0\t%s\t%s"*/, 
					readname, chromname1, &strand1, &prefixst, &prefixend, chromname2, &strand2, &suffixst, &suffixend, synchr, flankchr, &mis_match/*, &prefixlen, &spliceoutlen, &suffixlen,*//* chromseq, qualseq*/);

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_last_of("_");
				size_t tagidx = -1;

				if (last_idx != string::npos)
				{
					string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
					tagidx = atoi(tagidxstr.c_str()) - 1;
				}

				size_t max_loc = prefixend - prefixst + 1;

				string chr1 = chromname1, chr2 = chromname2;

				string chrom_ID = chr1;
				chrom_ID.append("_");
				chrom_ID.append(chr2);

				if (strand1 == '-')
				{
					size_t temp = prefixst;
					prefixst = prefixend;
					prefixend = temp;					
				}

				if (strand2 == '-')
				{		
					size_t temp = suffixst;
					suffixst = suffixend;
					suffixend = temp;
				}

				//cout << "find(suffixst) "<<endl;
				JUNCFusion_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

				if (iter_conj == conj_hash_map[chrom_ID].end())
				{
					//cout << "insert(JUNC_HASH: "<<endl;
					JUNCFusion_SUFFIX match_hash_temp;
					conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
					iter_conj = conj_hash_map[chrom_ID].find(suffixst);
				}

				JUNCFusion_SUFFIX_ITER iter_conj_pre;


				if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
				{
					//if (iter_conj_pre->second.flankstr != flankstr)
					//{
					//	cout <<"flank does not match"<<endl;
					//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
					//}
					//cout << "inc_hits(max_loc) "<<endl;
					iter_conj_pre->second.inc_hits(max_loc, 0, mis_match);	
				}	  
				else 
				{
					//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
					//cout << "insert(JUNC_SUFFIX "<<endl;
					//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

					//string synstr = prefixsyn;
					//synstr.append(suffixsyn);

					string strand12 = "++";
					strand12[0] = strand1;
					strand12[1] = strand2;
					//(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankchr, max_loc, m_read_width, synchr, strand12, mis_match)));
				}

				////cout << line<<endl;
				//vector<pair<size_t, size_t> > spliceway_vec;

				//string splicewaystr = spliceway;
				//size_t index = 0;

				//bool exceed = false;
				//while (true)
				//{
				//	int maplen, intron;

				//	if (index == 0)
				//	{
				//		sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
				//		spliceway_vec.push_back(make_pair(prefixst, maplen));

				//		if (maplen <=0)
				//			exceed =true;

				//		//if (splicewaystr.find("M", index) == string::npos)
				//		//	break;

				//		index = splicewaystr.find("M", index) + 1;

				//		//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

				//		//getchar();

				//		
				//		//if (index == splicewaystr.length())
				//		//	break;
				//		//if (index == string::npos)
				//		//	break;
				//	}
				//	else
				//	{
				//		if (index == splicewaystr.length())
				//			break;
				//		sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);

				//		if (intron != 0 && maplen != 0)
				//		{
				//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));
				//		}

				//		if (maplen <=0)
				//			exceed =true;

				//		if (intron < min_intron || intron > max_intron)
				//			exceed = true;
			
				//		if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
				//		{
				//			//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;
				//			//if (splicewaystr.find("M", index) == splicewaystr.length() - 1)
				//			//	getchar();
				//			break;
				//		}

				//		index = splicewaystr.find("M", index) + 1;

				//		//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

				//		//getchar();
				//		
				//	}
				//}

				/*if (spliceway_vec.size() > 1 && !exceed)
				{
					vector<pair<size_t, size_t> >::iterator vp_iter;
					for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end() - 1; ++vp_iter)
					{
						//suffixend = prefixst + spliceoutlen + (int)m_read_width;

						prefixst = vp_iter->first;
						size_t prefixend = vp_iter->first + vp_iter->second - 1;

						size_t suffixst = (vp_iter + 1)->first;// + spliceoutlen;

						size_t prefixlen = vp_iter->second;

						size_t suffixlen = (vp_iter + 1)->second;

						//string flankstr = flankseq;

						//coverage block
						//v_coverage_block.push_back(CoverageBlock(prefixst + 1, prefixst + prefixlen, 1));
						//v_coverage_block.push_back(CoverageBlock(suffixst, suffixlen, 1));

						string chrom_ID = chromname;

						//if (chrom_ID == "chr1" &&  prefixend == 4862 && suffixst == 5002)
						//{
						//	cout <<line<<endl;
						//	for (size_t i = 0; i < spliceway_vec.size(); ++i)
						//		cout << spliceway_vec[i].first<<'\t'<<spliceway_vec[i].second<<endl;
						//	getchar();
						//}
						JUNC_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

						if (iter_conj == conj_hash_map[chrom_ID].end())
						{
							JUNC_SUFFIX match_hash_temp;
							conj_hash_map[chrom_ID].insert(JUNC_HASH::value_type(suffixst, match_hash_temp));
							iter_conj = conj_hash_map[chrom_ID].find(suffixst);
						}

						JUNC_SUFFIX_ITER iter_conj_pre; 

						if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
						{
							//if (iter_conj_pre->second.flankstr != flankstr)
							//{
							//	cout <<"flank does not match"<<endl;
							//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
							//}
							iter_conj_pre->second.inc_hits(prefixlen, suffixlen, tagidx);
						}	  
						else 
						{
							//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
							(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(/*prim, *//*flankseq, *///prefixlen, suffixlen, m_read_width, tagidx)));
						/*}
					}
				}*/
			}
			ifs.close();
		}
		else cout << "Unable to open file";

		cout <<"write junc"<<endl;

		WriteJuncFusion(conj_hash_map, ofs, m_read_width, m_max_rank);
		//WriteJunc(conj_hash_map, ofs, m_read_width, m_max_rank, chrom_dir, max_intron);

		//cout <<"sort v_coverage_block"<<endl;
		//sort(v_coverage_block.begin(), v_coverage_block.end(), compblk);

		//cout <<"WriteCoverage"<<endl;
		//WriteCoverage(v_coverage_block, conj_hash_map.begin()->first, ofs_wig, mi);
	}

	ofs.close();

	//getchar();

	cout <<"sort junc"<<endl;

	SortJuncFusion(junc_filename);
}

void
CovertFusion2JuncComb(const char* junc_filename, vector<string>& m_mapreads_files, size_t m_read_width, double m_max_rank)
{
	ofstream ofs(junc_filename);

	string headline = "track name=junctions description=\"Mapsplice fusion junctions\"";

	ofs << headline << endl;

	CHROM_FUSION_JUNC_HASH_COMB conj_hash_map;

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

				//cout <<line<<endl;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000], synchr[1000], flankchr[100];
				char chromname1[1000], chromname2[1000];
				char strand = '+', strand1, strand2;
				size_t /*prim, */ prefixst, prefixend, suffixst, suffixend, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/;
				//sscanf(line.c_str(), "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d", chromname, readname, &prim, flankseq, strand, readseq, prefixseq, suffixseq, &prefixlen, &prefixend, &suffixst);
				unsigned short mis_match;

//				TRAN00000074719:37	1	chr1	971958	0	21M84N54M	GCCAGTGGGGGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAG
//TRAN00000074719:100	1	chr1	972105	0	75M	CCACCTGTCGAGAGGGCTTCCTGCTACAACTCCGCGTTGGGCTGCTGCTCTGATGGGAAGACGCCCTCGCTGGAC
//TRAN00000074719:46	1	chr1	971967	0	12M84N63M	GGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAGAGGGCTTCC

//HWI-EAS217:4:52:346:1875#0/1_3287705	chr1	-	18451	18525	chr1	+	88603914	88603938	GCGATACACATGTACAAGTTTTGTCACTTAAAAAGAATTTGAAAAAACTTCATAGATGCAAAAAAAAAAAAAAACAAAAACCACCTA
//TTTTTAAAGATTA	TGTT	NM:i:2

//1012966~TRAN00000110580:1400:334        chr10   +       66128298        66128322        chr1    +       241548206       24154
//8255    tgaggcttatcctaggtaaactaacgcaagaacagaaaaccaaatactgcatattcttacgtataagtgggagctaaatgataagaacttatgagcacaaagtaggaaaccacagac
//agt     TTAC    NM:i:2
				                                       //%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu"/*\t*\t0\t0\t%s\t%s"*/, 
					readname, chromname1, &strand1, &prefixst, &prefixend, chromname2, &strand2, &suffixst, &suffixend, synchr, flankchr, &mis_match/*, &prefixlen, &spliceoutlen, &suffixlen,*//* chromseq, qualseq*/);

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_first_of("~");
				size_t tagidx = -1;

				if (last_idx != string::npos)
				{
					string tagidxstr = tagnamestr.substr(0, last_idx);
					tagidx = atoi(tagidxstr.c_str()) - 1;
				}

				size_t max_loc = prefixend - prefixst + 1;

				string chr1 = chromname1, chr2 = chromname2;

				string chrom_ID = chr1;
				chrom_ID.append("_");
				chrom_ID.append(chr2);

				size_t prefixlen = prefixend - prefixst + 1, suffixlen = suffixend - suffixst + 1;

				if (strand1 == '-')
				{
					size_t temp = prefixst;
					prefixst = prefixend;
					prefixend = temp;					
				}

				if (strand2 == '-')
				{		
					size_t temp = suffixst;
					suffixst = suffixend;
					suffixend = temp;
				}

				//cout << "find(suffixst) "<<endl;
				size_t combined_offset = (prefixend << THIRTY_TWO) + suffixst;

				CHROM_FUSION_JUNC_HASH_COMB_ITER chrom_junc_hash_iter = conj_hash_map.find(chrom_ID);

				if (chrom_junc_hash_iter == conj_hash_map.end())
				{
					FUSION_JUNC_HASH_COMB junc_hash_comb;

					chrom_junc_hash_iter = (conj_hash_map.insert(CHROM_FUSION_JUNC_HASH_COMB::value_type(chrom_ID, junc_hash_comb))).first;							
				}

				FUSION_JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

				FUSION_JUNC_HASH_COMB_ITER junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

				if (junc_hash_comb_iter != junc_hash_comb.end())
				{
					junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, mis_match);
				}
				else
				{
					//string synstr = prefixsyn;
					//synstr.append(suffixsyn);

					//string strand12 = "++";
					//strand12[0] = strand1;
					//strand12[1] = strand2;

					//bool isbads = m_matched_bads ? true : false;

					//unsigned short pm = static_cast<unsigned short>(prim);

					chrom_junc_hash_iter->second.insert(FUSION_JUNC_HASH_COMB::value_type(combined_offset, JuncInfoFusion(prefixlen, suffixlen, m_read_width, synchr, strand1, strand2, mis_match, flankchr)));
				}


				//FUSION_JUNC_HASH_COMB_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

				//if (iter_conj == conj_hash_map[chrom_ID].end())
				//{
				//	//cout << "insert(JUNC_HASH: "<<endl;
				//	JUNCFusion_SUFFIX match_hash_temp;
				//	conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
				//	iter_conj = conj_hash_map[chrom_ID].find(suffixst);
				//}

				//JUNCFusion_SUFFIX_ITER iter_conj_pre;


				//if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
				//{
				//	//if (iter_conj_pre->second.flankstr != flankstr)
				//	//{
				//	//	cout <<"flank does not match"<<endl;
				//	//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
				//	//}
				//	//cout << "inc_hits(max_loc) "<<endl;
				//	iter_conj_pre->second.inc_hits(max_loc, mis_match);	
				//}	  
				//else 
				//{
				//	//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
				//	//cout << "insert(JUNC_SUFFIX "<<endl;
				//	//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

				//	//string synstr = prefixsyn;
				//	//synstr.append(suffixsyn);

				//	string strand12 = "++";
				//	strand12[0] = strand1;
				//	strand12[1] = strand2;
				//	(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankchr, max_loc, m_read_width, synchr, strand12, mis_match)));
				//}

				////cout << line<<endl;
				//vector<pair<size_t, size_t> > spliceway_vec;

				//string splicewaystr = spliceway;
				//size_t index = 0;

				//bool exceed = false;
				//while (true)
				//{
				//	int maplen, intron;

				//	if (index == 0)
				//	{
				//		sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
				//		spliceway_vec.push_back(make_pair(prefixst, maplen));

				//		if (maplen <=0)
				//			exceed =true;

				//		//if (splicewaystr.find("M", index) == string::npos)
				//		//	break;

				//		index = splicewaystr.find("M", index) + 1;

				//		//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

				//		//getchar();

				//		
				//		//if (index == splicewaystr.length())
				//		//	break;
				//		//if (index == string::npos)
				//		//	break;
				//	}
				//	else
				//	{
				//		if (index == splicewaystr.length())
				//			break;
				//		sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);

				//		if (intron != 0 && maplen != 0)
				//		{
				//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));
				//		}

				//		if (maplen <=0)
				//			exceed =true;

				//		if (intron < min_intron || intron > max_intron)
				//			exceed = true;
			
				//		if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
				//		{
				//			//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;
				//			//if (splicewaystr.find("M", index) == splicewaystr.length() - 1)
				//			//	getchar();
				//			break;
				//		}

				//		index = splicewaystr.find("M", index) + 1;

				//		//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

				//		//getchar();
				//		
				//	}
				//}

				/*if (spliceway_vec.size() > 1 && !exceed)
				{
					vector<pair<size_t, size_t> >::iterator vp_iter;
					for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end() - 1; ++vp_iter)
					{
						//suffixend = prefixst + spliceoutlen + (int)m_read_width;

						prefixst = vp_iter->first;
						size_t prefixend = vp_iter->first + vp_iter->second - 1;

						size_t suffixst = (vp_iter + 1)->first;// + spliceoutlen;

						size_t prefixlen = vp_iter->second;

						size_t suffixlen = (vp_iter + 1)->second;

						//string flankstr = flankseq;

						//coverage block
						//v_coverage_block.push_back(CoverageBlock(prefixst + 1, prefixst + prefixlen, 1));
						//v_coverage_block.push_back(CoverageBlock(suffixst, suffixlen, 1));

						string chrom_ID = chromname;

						//if (chrom_ID == "chr1" &&  prefixend == 4862 && suffixst == 5002)
						//{
						//	cout <<line<<endl;
						//	for (size_t i = 0; i < spliceway_vec.size(); ++i)
						//		cout << spliceway_vec[i].first<<'\t'<<spliceway_vec[i].second<<endl;
						//	getchar();
						//}
						JUNC_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

						if (iter_conj == conj_hash_map[chrom_ID].end())
						{
							JUNC_SUFFIX match_hash_temp;
							conj_hash_map[chrom_ID].insert(JUNC_HASH::value_type(suffixst, match_hash_temp));
							iter_conj = conj_hash_map[chrom_ID].find(suffixst);
						}

						JUNC_SUFFIX_ITER iter_conj_pre; 

						if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
						{
							//if (iter_conj_pre->second.flankstr != flankstr)
							//{
							//	cout <<"flank does not match"<<endl;
							//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
							//}
							iter_conj_pre->second.inc_hits(prefixlen, suffixlen, tagidx);
						}	  
						else 
						{
							//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
							(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(/*prim, *//*flankseq, *///prefixlen, suffixlen, m_read_width, tagidx)));
						/*}
					}
				}*/
			}
			ifs.close();
		}
		else cout << "Unable to open file";
		
		//WriteJunc(conj_hash_map, ofs, m_read_width, m_max_rank, chrom_dir, max_intron);

		//cout <<"sort v_coverage_block"<<endl;
		//sort(v_coverage_block.begin(), v_coverage_block.end(), compblk);

		//cout <<"WriteCoverage"<<endl;
		//WriteCoverage(v_coverage_block, conj_hash_map.begin()->first, ofs_wig, mi);
	}

	cout <<"write junc"<<endl;

	WriteJuncFusion(conj_hash_map, ofs, m_read_width, m_max_rank);

	ofs.close();

	//getchar();

	cout <<"sort fusion junc"<<endl;

	SortJuncFusion(junc_filename);
}

void
RemoveDupMapreads(vector<string>& comb_mapreads_files)
{
	map<string, int> mapped_reads;
	vector<string>::iterator cm_iter;
	for (cm_iter = comb_mapreads_files.begin(); cm_iter != comb_mapreads_files.end(); ++cm_iter)
	{
		ifstream ifs((*cm_iter).c_str());

		if (ifs.is_open())
		{
			string line;
			while (getline(ifs,line))
			{
				if (line == "")
					continue;

				if (mapped_reads.find(line) == mapped_reads.end())
					mapped_reads[line] = 1;
				else
					mapped_reads[line]++;
			}

			ifs.close();

			ofstream ofs((*cm_iter).c_str());

			map<string, int>::iterator msi_iter;

			for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
			{
				ofs<<msi_iter->first<<endl;
			}

			ofs.close();
		}
		else
		{
			cout << "can't open file "<< (*cm_iter) <<endl;
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		cout << "mapreads_file junction_file read_width min_intron max_intron "<<endl;
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	string juncfile = argv[1];

	int read_width = atoi(argv[2]);

	//string chrom_dir = argv[3];

	//size_t min_intron = atoi(argv[4]);
	//
	//size_t max_intron = atoi(argv[5]);

	//chrom_dir.append("/");

	double max_rank = 0;

	vector<string> comb_mapreads_files;

	cout << argv[3] << endl;
	for (int i = 3; i < argc; ++i)
		comb_mapreads_files.push_back(argv[i]);

	//vector<string> m_mapreads_files;

	//cout << "separate fusion sam"<<endl;
	
	//SeprateMapreads(comb_mapreads_files, m_mapreads_files, juncfile.c_str());

	//RemoveDupMapreads(m_mapreads_files);

	CovertFusion2JuncComb(juncfile.c_str(), comb_mapreads_files, read_width, max_rank);
}