/*    
 *    newsam2junc_.cpp		
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
#endif

struct JuncInfo{
	JuncInfo(/*int pm, *//*const string fs, */int loc, int suffix_len, size_t rw, size_t tagidx, unsigned short mis) : /*prim(pm), *//*flankstr(fs),*/ p(rw - 1, 0) 
	{
		++p[loc-1];
		max_prefix_len = loc;
		max_suffix_len = suffix_len;

		max_mismatch = mis;

		min_mismatch = mis;

		sum_mismatch = mis;

		m[tagidx] = 1;
	}
	bool inc_hits(int idx, int suffix_len, size_t tagidx, unsigned short mis)
	{
		if (tagidx != -1 && m.find(tagidx) != m.end())
			return false;

		++p[idx-1];

		m[tagidx] = 1;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;

		if (mis > max_mismatch)
			max_mismatch = mis;
		
		if (mis < min_mismatch)
			min_mismatch = mis;

		sum_mismatch += mis;

		return true;
	}
	//int prim;
	//string flankstr;
	size_t max_prefix_len;
	size_t max_suffix_len;
	vector<unsigned short> p;

	unsigned short max_mismatch;

	unsigned short min_mismatch;

	unsigned short sum_mismatch;

	map<size_t, int> m;
};

typedef hash_map<size_t, JuncInfo> JUNC_SUFFIX;
typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;

typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
typedef JUNC_HASH::iterator JUNC_HASH_ITER;
typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef map<string, JUNC_HASH > CONJ_HASH_MAP;

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

	string seprate_path = juncfile;
	//seprate_path.append(comb_mapreads_files.front());

	seprate_path.append(".sepdir");

	string mkdir_cmdstr = "mkdir ";
	mkdir_cmdstr.append(seprate_path);
	system(mkdir_cmdstr.c_str());

	seprate_path.append("/");

	map<string, ofstream* > chrom_ofs_map; 

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
				char chromname[1000], readname[1000];

				int strand_t;

				sscanf(line.c_str(), "%s\t%d\t%s", 
					readname, &strand_t, chromname);

				char mapreads_file[1000];

				sprintf(mapreads_file, "%smapreads_%s.txt", seprate_path.c_str(), chromname);

				if (chrom_ofs_map.find(mapreads_file) == chrom_ofs_map.end())
				{
					ofstream* ofs = new ofstream;

					ofs->open(mapreads_file);

					chrom_ofs_map[mapreads_file] = ofs;
				}

				(*chrom_ofs_map[mapreads_file])<< line<<endl;
			}

			ifs.close();
		}
	}

	map<string, ofstream*>::iterator mf_iter;
	for (mf_iter = chrom_ofs_map.begin(); mf_iter != chrom_ofs_map.end(); ++mf_iter)
	{
		m_mapreads_files.push_back(mf_iter->first);

		mf_iter->second->close();

		delete mf_iter->second;
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
	unsigned short min_mismatch;
	unsigned short max_mismatch;
	double ave_mismatch;

	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& l, unsigned short min_mis, unsigned short max_mis, double ave_mis) : juncst(jst), juncend(jend), hits(hts), 
		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}

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

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, lpq, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset,	rank, lpq, min_mis, max_mis, ave_mis));
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
				flankcase = 2;
			else if (flankstr == "CTGC")
				 flankcase = 3;
			else if (flankstr == "GCAG")
				 flankcase = 4;
			else if (flankstr == "GTAG")
				 flankcase = 5;
			else if (flankstr == "GTAT")
				 flankcase = 6;

			double il = 1.0 - (((double) (vit->juncend - vit->juncst - min_intron + 1)) / (double (max_intron - min_intron + 2))); 
			ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<'+'<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< flankcase <<'\t'<<flankstr <<'\t'<<il<<'\t'<<vit->lpq
				<<'\t' << vit->min_mismatch<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<<endl;
			++juncid;
		}
	} 

	cout << juncid <<" junctions"<<endl;
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
#ifdef DEBUG
				cout <<chm_iter->first <<'\t'<< iter_conj_pre->first << '\t' << iter_conj->first<< '\t'<<juncidstr<<juncid<<'\t'; 
#endif
				//int juncidlen = (int) log10((double)juncid);

				//string suffix0(8 - juncidlen

				//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';
#ifdef DEBUG
				cout<< hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';
#endif
				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
				//		ofs << k << ',';
				//		break;
				//	}
				//}
				ofs << iter_conj_pre->second.max_prefix_len << ',';
#ifdef DEBUG
				cout << iter_conj_pre->second.max_prefix_len << ',';
#endif
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
#ifdef DEBUG
				cout << iter_conj_pre->second.max_suffix_len << ',';
				cout <<  '\t';
#endif
				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
				//		ofs << 0 << ',';
				//		break;
				//	}
				//}
				ofs << 0 << ',';
#ifdef DEBUG
				cout<< 0 << ',';
#endif

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
#ifdef DEBUG
				cout << int (iter_conj->first + iter_conj_pre->second.max_suffix_len - iter_conj_pre->first + 1 )<<',';
				cout << '\t';
#endif


				size_t intron_len = iter_conj->first - iter_conj_pre->first - 1;

				//double ppower = pow(double(iter_conj_pre->second.max_prefix_len), 4);

				////double pNpower = pow(1.0 - ppower, (double)chrom_size);

				//double qpower = pow(double(iter_conj_pre->second.max_suffix_len), 4);

				////double pDpower = pow(1.0 - qpower, (double)max_intron);

				//double lpq = ((double)chrom_size / ppower) + ((double)max_intron / qpower);

				double ppower = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower = pow(1.0 - ppower, (double)chrom_size);

				double qpower = pow(0.25, double(iter_conj_pre->second.max_suffix_len));

				double pDpower = pow(1.0 - qpower, (double)intron_len);

				double lpq = 1.0 - (pNpower * pDpower);

				double ppower2 = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower2 = pow(1.0 - ppower2, (double)intron_len );

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

				ofs << rank << '\t'<< lpqave/*iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr */
					<<'\t'<<iter_conj_pre->second.min_mismatch<<'\t' <<iter_conj_pre->second.max_mismatch<<'\t'
					<< (double)iter_conj_pre->second.sum_mismatch / (double)hits<<endl;

#ifdef DEBUG
				cout << rank << '\t'<< lpqave/*iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr */
					<<'\t'<<iter_conj_pre->second.min_mismatch<<'\t' <<iter_conj_pre->second.max_mismatch<<'\t'
					<< (double)iter_conj_pre->second.sum_mismatch / (double)hits<<endl;

				cout <<"p q chrom_size intron_len"<<endl;
				cout << iter_conj_pre->second.max_prefix_len <<'\t' << iter_conj_pre->second.max_suffix_len << '\t' << chrom_size << '\t'<<intron_len<<endl;
				cout << "lpqave: "<<lpqave<<endl;
#endif

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}

		//cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}


void
Covert2Junc(const char* junc_filename, vector<string>& m_mapreads_files, size_t m_read_width, string chrom_dir, double m_max_rank, size_t min_intron, size_t max_intron, size_t min_anchor)
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

				unsigned short mis_match;
				//sscanf(line.c_str(), "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d", chromname, readname, &prim, flankseq, strand, readseq, prefixseq, suffixseq, &prefixlen, &prefixend, &suffixst);

//				TRAN00000074719:37	1	chr1	971958	0	21M84N54M	GCCAGTGGGGGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAG
//TRAN00000074719:100	1	chr1	972105	0	75M	CCACCTGTCGAGAGGGCTTCCTGCTACAACTCCGCGTTGGGCTGCTGCTCTGATGGGAAGACGCCCTCGCTGGAC
//TRAN00000074719:46	1	chr1	971967	0	12M84N63M	GGTGGCTCTGGGGGGCTCGAGCCCTTGGAGGGCAGCAGCGTGGCCACCCCTGGGCCACCTGTCGAGAGGGCTTCC

				                                       //%dM\t%dN\t%dM
				sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%hu", 
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway/*, &prefixlen, &spliceoutlen, &suffixlen,*/ ,chromseq, qualseq, &mis_match);

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_first_of("~");
				size_t tagidx = -1;

				if (last_idx != string::npos)
				{
					string tagidxstr = tagnamestr.substr(0, last_idx);
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

						if (maplen < 0/* || maplen < min_anchor*/)
						{
							exceed =true;
							//cout << splicewaystr << endl;
						}

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

						if (intron >= 0 && maplen >= 0)
						{
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));
						}

						if (maplen < 0/* || maplen < min_anchor*/)
						{
							exceed =true;
							//cout << splicewaystr << endl;
						}

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

						if (prefixlen <= 0 || suffixlen <= 0)
							continue;

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
							iter_conj_pre->second.inc_hits(prefixlen, suffixlen, tagidx, mis_match);
						}	  
						else 
						{
							//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
							(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(/*prim, *//*flankseq, */prefixlen, suffixlen, m_read_width, tagidx, mis_match)));
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

	//getchar();

	SortJunc(junc_filename, chrom_dir, min_intron, max_intron);
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
	if (argc < 6)
	{
		cout << "mapreads_file junction_file read_width min_intron max_intron "<<endl;
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	string juncfile = argv[1];

	int read_width = atoi(argv[2]);

	string chrom_dir = argv[3];

	size_t min_intron = atoi(argv[4]);
	
	size_t max_intron = atoi(argv[5]);

	size_t min_anchor = atoi(argv[6]);

	chrom_dir.append("/");

	double max_rank = 0;

	vector<string> comb_mapreads_files;

	cout << argv[7] << endl;
	for (int i = 7; i < argc; ++i)
		comb_mapreads_files.push_back(argv[i]);

	vector<string> m_mapreads_files;
	
	SeprateMapreads(comb_mapreads_files, m_mapreads_files, juncfile.c_str());

	//RemoveDupMapreads(m_mapreads_files);

	cout << juncfile << endl;

	Covert2Junc(juncfile.c_str(), m_mapreads_files, read_width, chrom_dir, max_rank, min_intron, max_intron, min_anchor);
}