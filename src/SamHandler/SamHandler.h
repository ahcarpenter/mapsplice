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
	JuncInfo(int loc, size_t suffix_len, size_t rw, size_t tagidx) : p(rw - 1, 0) 
	{
		++p[loc-1];
		max_prefix_len = loc;
		max_suffix_len = suffix_len;

		m[tagidx] = 1;
	}
	bool inc_hits(size_t idx, size_t suffix_len, size_t tagidx)
	{
		if (tagidx != size_t(-1) && m.find(tagidx) != m.end())
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

bool 
comp_str_less_than_linux(const string& lhs, const string& rhs)
{
	if (lhs == rhs)
		return lhs < rhs;

	if (lhs.find(rhs) != string::npos)
		return true;

	if (rhs.find(lhs) != string::npos)
		return false;

	return lhs < rhs;
}

string
base_name(string filename) {
	const string s(filename.substr(0, filename.find_last_of(".")));
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
		return s.substr(final_slash + 1);
	else return s;
}


typedef hash_map<size_t, JuncInfo> JUNC_SUFFIX;
typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;

typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
typedef JUNC_HASH::iterator JUNC_HASH_ITER;
typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef map<string, JUNC_HASH > CONJ_HASH_MAP;

static const string Is(500, 'I');

struct SamRec {
	string tag_name;
	unsigned short strand_t;
	string chrom_name;
	size_t start, end;
	unsigned short confid;
	string splice_way;
	string mapped_seq;
	unsigned short mis_match;

	vector<pair<size_t, int> > spliceway_vec;
	bool wrong_format;

	size_t tagidx;

	double filter_score;

	size_t ave_intron_len;

	double ave_junc_mis;

	vector<string> junc_id;

	int best;

	string alters;

	string qual_str;

	bool isunique, isexonic, isspliced, issmallins, issmalldel;
	
	size_t canon_count, noncanon_count;

	double canon_rate;

	int matched_id;

	char mate_match;

	size_t mate_offest;

	int mate_diff;

	string cur_line;

	bool is_insert;

	size_t mappedlen;

	//bool is_spliced;

	SamRec() {}

	SamRec(const string& tname, size_t strand, const string& cname, size_t st, size_t conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
		char matematch, size_t mateoffest, int matediff, string line = "") : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
		mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offest(mateoffest), mate_diff(matediff), wrong_format(false), tagidx(tidx), filter_score(0), best(0), alters(alt), qual_str(qualstr), 
		isunique(false), isexonic(false), isspliced(false), ave_intron_len(0), ave_junc_mis(-1), canon_count(0), noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), cur_line(line),
		is_insert(false)/*, is_spliced(false)*/, mappedlen(0)
	{
		size_t index = 0;

		if ((spliceway.find("I") != string::npos))
			issmallins = true;
		else if ((spliceway.find("D") != string::npos))
			issmalldel = true;
		else if (spliceway.find("N") != string::npos)
			isspliced = true;
		else
			isexonic = true;		

		string flag_str = " ";

		//isexonic = true;

		//isspliced = false;

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
					spliceway_vec.push_back(make_pair(start, -maplen));
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));

					mappedlen += maplen;
				}
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

					mappedlen += maplen;
				}
				else if (flag == 'N')
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
					spliceway_vec.back().second = -maplen;
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;

					mappedlen += maplen;
				}
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
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

					mappedlen += maplen;
				}
				else if (flag == 'N')
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			//if (

			//if (index == 0)
			//{
			//	sscanf(splice_way.c_str() + index, "%dM", &maplen);

			//	spliceway_vec.push_back(make_pair(start, maplen));

			//	isexonic = true;

			//	isspliced = false;

			//	if (maplen <=0)
			//	{
			//		wrong_format =true;
			//		break;
			//	}

			//	index = splice_way.find("M", index) + 1;				
			//}
			//else
			//{
			//	if (index == splice_way.length())
			//		break;

			//	sscanf(splice_way.c_str() + index, "%dN%dM", &intron, &maplen);

			//	if (maplen <=0 || intron <= 0)
			//	{
			//		wrong_format = true;
			//		break;
			//	}

			//	spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));

			//	isexonic = false;

			//	isspliced = true;

			//	if (splice_way.find("M", index) == string::npos || splice_way.find("M", index) == splice_way.length() - 1)
			//		break;

			//	index = splice_way.find("M", index) + 1;

			//}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;
	}

	string tostring()
	{
		if (!cur_line.empty())
			return cur_line;

		char sam_rec_char[5000];

		//cout << tag_name.c_str()<<'\t' << strand_t<<'\t' <<  chrom_name.c_str()<<'\t' <<  start<<'\t' <<  confid<<'\t' <<  splice_way.c_str()<<"\t*\t0\t0\t" <<  mapped_seq.c_str()<<'\t' <<  Is.substr(0, mapped_seq.length()).c_str()<<'\t' <<  mis_match << endl;

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s"/*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
			chrom_name.c_str(), start, confid, splice_way.c_str(), mate_match, mate_offest, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, alters.c_str()/*, best, filter_score*/);

		return sam_rec_char;
	}
};

extern bool comp_offset(const SamRec* lhs, const SamRec* rhs);

extern bool comp_tagname_offset(const SamRec& lhs, const SamRec& rhs);

extern bool comp_filterscore(const SamRec* lhs, const SamRec* rhs);

struct JuncForSort{
	size_t juncst;
	size_t juncend;
	unsigned short hits;
	unsigned short kinds;
	string blocks;
	string blocksoffset;
	double rank;
	double lpq;
	string flankstring;
	unsigned short flankcase;
	double intronlen;
	string juncname;

	unsigned short min_mismatch, max_mismatch;

	double ave_mismatch;

	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const double& rk, const double& l, const string& jn) : juncst(jst), juncend(jend), hits(hts), 
		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), juncname(jn) {}

	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const double& rk, const double& l, const string& fs, const int fc, const double i, const string& jn, unsigned short min_mis, unsigned short max_mis, double ave_mis) : juncst(jst), juncend(jend), hits(hts), 
		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), flankstring(fs), flankcase(fc), intronlen(i), juncname(jn), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}
};

class SamHandler {
private:

	size_t FilterByMisMatch(map<string, map<size_t, vector<SamRec*> > >& sam_rec_map_filter_mismatch, size_t maxhits, const char* filtered_file);

	size_t FilterByMisMatchByRead(map<size_t, vector<SamRec> >& sam_rec_map_filter_mismatch, ofstream* ofs_notminmis, ofstream* ofs_exceed);

	size_t FilterByQualScoreByRead(vector<SamRec>& sam_rec_map_filter_mismatch, ofstream* ofs_notminqualscore, ofstream* ofs_exceed);

	size_t CalcFilterScore(map<string, map<size_t, vector<SamRec*> >  >& sam_rec_map_filter_mismatch, map<string, vector<JuncForSort> >& junc_map, 
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	size_t CalcFilterScoreByRead(vector<SamRec>& sam_rec_map_filter_mismatch, map<string, vector<JuncForSort> >& junc_map, 
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);
	//double entrpy_weight = 0.22419185, double pqlen_weight = 0.80924761

	map<string, vector<SamRec> > sam_rec_map;

	map<string, vector<SamRec> > sam_rec_map_next;
	
	string prev_tag;

	map<size_t, size_t> multiple_mapped_counts;

public:

	ifstream cur_ifs;

	string file_name;

	string base_file_name;

	size_t ReadSortedJunction(const char* juncfile, map<string, vector<JuncForSort> >& junc_map);

	size_t MarkCanonNoncanon(map<string, vector<JuncForSort> >& junc_map);

	void MarkCanonNoncanonByReads(vector<SamRec>& read_sam, map<string, vector<JuncForSort> >& junc_map);

	size_t exonic_reads, exonic_unique, exonic_multiple, spliced_reads, spliced_unique, spliced_multiple, unique_reads, multiple_reads, total_mapped_reads, exonic_shared, exonic_notshared, exonic_same_start_or_end, exonic_only,
		spliced_shared, spliced_notshared, spliced_same_start_or_end, spliced_only, small_del, small_insert, small_del_multiple, small_insert_multiple, small_del_unique, small_insert_unique;

	size_t unique_shared, unique_notshared, unique_same_start_or_end, unique_only, multiple_shared, multiple_notshared, multiple_same_start_or_end, multiple_only;

	SamHandler(): exonic_reads(0), spliced_reads(0), unique_reads(0), multiple_reads(0), exonic_shared(0), exonic_notshared(0), spliced_shared(0), spliced_notshared(0), small_del(0), small_insert(0),
		unique_shared(0), unique_notshared(0), multiple_shared(0), multiple_notshared(0), total_mapped_reads(0), exonic_unique(0), exonic_multiple(0), spliced_unique(0), spliced_multiple(0),
		small_del_multiple(0), small_insert_multiple(0), small_del_unique(0), small_insert_unique(0), exonic_same_start_or_end(0), spliced_same_start_or_end(0), unique_same_start_or_end(0), multiple_same_start_or_end(0),
		exonic_only(0), spliced_only(0), unique_only(0), multiple_only(0) {}

	bool OpenSam(const char* sam_file)
	{
		if (cur_ifs.is_open())
		{
			cur_ifs.close();
		}

		cur_ifs.open(sam_file, ios_base::in);

		file_name = sam_file;

		base_file_name = base_name(file_name);

		return cur_ifs.is_open();
	}

	size_t ReadNextTag();

	size_t ReadAppendSamFile(const char* samfile);

	size_t ReadAppendSamFileSimp(const char* samfile);

	size_t CollectStat();

	size_t CollectStat2();

	const void WriteCompStat(const SamHandler& synthetic_sam, const char* output_file);

	const void WriteStat(const char* output_file, int unmapped);

	const void WriteStat(const char* output_file);

	size_t ReadNewSamFile(const char* samfile)
	{
		Clear();

		return ReadAppendSamFile(samfile);
	}

	size_t WriteSamFile(const char* samfile);

	size_t Sam2Junc(const char* juncfile, size_t read_width, string chrom_dir, size_t min_intron, size_t max_intron);

	bool Empty()
	{
		if (sam_rec_map.size())
			return false;
		else
			return true;
	}

	void WriteSamByChromoOffset(const char* samfile);

	void SortSamByChromoOffset();

	void SortSamByTagnameOffset();

	void WriteSamByTagnameOffset(const char* samfile);

	size_t WriteUniqueMappedTagSam(const char* unique_samfile);

	size_t WriteMultipleMappedTagSam(const char* unique_samfile);

	size_t FilterMultipleMappedTag(const char* filtered_samfile, const char* junc_file, size_t maxhits, double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	size_t FilterMultipleMappedTagByRead(const char* sam_file, const char* unique_samfile, const char* filtered_samfile, const char* junc_file, size_t maxhits, const char* stat_file, 
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	size_t FilterMultipleMappedTagByReadFilterByRepeat(const char* sam_file, const char* unique_samfile, const char* filtered_samfile, const char* junc_file, size_t maxhits, const char* stat_file, 
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	size_t FilterCanonNoncaonTag(const char* sam_file, const char* filtered_canon, const char* filtered_noncanon, const char* fitlered_noncanon_canon);

	size_t FilterCanonNoncaonByReads(const char* sam_file, const char* filtered_canon, const char* filtered_noncanon, 
		const char* fitlered_noncanon_canon, const char* filtered_ins_file, const char* stat_file, map<string, vector<JuncForSort> >& junc_map);

	size_t CompareToSyntheticMappedReads(SamHandler& synthetic_sam, const char* output_file, const char* same_offest_output, size_t anchor_width, size_t range = 0);

	size_t Compare2Sam(SamHandler& compare_sam, const char* output_file, const char* same_offest_output, size_t anchor_width, size_t range = 0);

	size_t Compare2Sam(SamHandler& compare_sam, string output_dir);

	pair<size_t, size_t> SepSplicedUnspliced(const char* spliced, const char* unspliced);

	size_t FilterByAnchor(const char* sam_file, size_t min_anchor);

	size_t FilterByRepeated(const char* repeated_reads, size_t fq, const char* sam_file);

	size_t FilterByUnmapped(const char* unmapped_reads, size_t fq, const char* sam_file);

	size_t FilterByIsland(const char* islandfile, const char* sam_file);

	size_t FilterByIntronLenHmer(const char* sam_file, size_t intronlen, size_t seg_len);

	size_t FilterBySmallExon(const char* sam_file, size_t seg_len);

	size_t FilterByUnsplicedMapped(SamHandler& unsplice_sam, char* sam_file);

	//size_t SeparateSplicedUnspliced(const char* output_file, const char* same_offest_output, size_t range = 0);

	bool Clear()
	{
		sam_rec_map.clear();

		return true;
	}

	size_t CompareTwoSamByOffset(const char* cur_file, const char* syn_file, const char* output, size_t anchor_width, size_t range = 0);

};

extern string base_name(string filename);

extern void readchrom(const char* filename, string& longseq);

extern size_t WriteJunc(const CONJ_HASH_MAP& conj_hash_map, ofstream& ofs, size_t m_read_width, double m_max_rank, string chrom_dir);

extern bool compjunc(const JuncForSort& lhs, const JuncForSort& rhs);

extern size_t SortJunc(const char* juncfile, string chrom_dir, size_t min_intron, size_t max_intron);

bool comp_offset(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->chrom_name == rhs->chrom_name)
		return lhs->start < rhs->start;
	else
		return lhs->chrom_name < rhs->chrom_name;
}

bool comp_chromname_offset_spliceway(const SamRec& lhs, const SamRec& rhs)
{
	if (lhs.chrom_name == rhs.chrom_name)
	{
		if (lhs.start == rhs.start)
		{
			if (lhs.splice_way.length() == rhs.splice_way.length())
				return lhs.splice_way < rhs.splice_way;
			else
				return lhs.splice_way.length() < rhs.splice_way.length();
		}
		else
			return lhs.start < rhs.start;
	}
	else
		return lhs.chrom_name < rhs.chrom_name;
}

bool comp_tagname_offset(const SamRec& lhs, const SamRec& rhs)
{
	if (lhs.chrom_name == rhs.chrom_name)
	{
		if (lhs.start == rhs.start)
		{
			return lhs.splice_way < rhs.splice_way;
		}
		else
			return lhs.start < rhs.start;
	}
	else
		return lhs.chrom_name < rhs.chrom_name;
}


void
readchrom(const char* filename, string& longseq)
{
	size_t size;  

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

size_t
WriteJunc(const CONJ_HASH_MAP& conj_hash_map, ofstream& ofs, size_t m_read_width, double m_max_rank, string chrom_dir)
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
				
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';

				ofs << iter_conj_pre->second.max_prefix_len << ',';

				ofs << iter_conj_pre->second.max_suffix_len << ',';
				ofs << '\t';

				ofs << 0 << ',';

				ofs << int (iter_conj->first + iter_conj_pre->second.max_suffix_len - iter_conj_pre->first + 1 )<<',';
				ofs << '\t';

				size_t intron_len = iter_conj->first - iter_conj_pre->first - 1;

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

				ofs << rank << '\t'<< lpqave<<endl;

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}
	}

	return juncid;
}

bool
compjunc(const JuncForSort& lhs, const JuncForSort& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

size_t
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

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000];
			int juncst, juncend, prefixend, suffixst, kinds, hits;
			double rank, lpq;

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%lf\t%lf\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, &rank, &lpq, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset,	rank, lpq, juncname));
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
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< flankcase <<'\t'<<flankstr <<'\t'<<il<<'\t'<<vit->lpq<<'\t'<<vit->min_mismatch
				<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<<endl;
			++juncid;
		}
	} 

	cout << juncid <<" junctions"<<endl;

	return juncid;
}

bool
comp_filterscore(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->filter_score == rhs->filter_score)
	{
		if (lhs->ave_intron_len == rhs->ave_intron_len)
		{
			return lhs->ave_junc_mis < rhs->ave_junc_mis;
		}
		else
			return lhs->ave_intron_len < rhs->ave_intron_len;
	}
	else
		return lhs->filter_score > rhs->filter_score;
}

bool
comp_filterscore_ref(const SamRec& lhs, const SamRec& rhs)
{
	if (lhs.filter_score == rhs.filter_score)
	{
		if (lhs.ave_intron_len == rhs.ave_intron_len)
		{
			return lhs.ave_junc_mis < rhs.ave_junc_mis;
		}
		else
			return lhs.ave_intron_len < rhs.ave_intron_len;
	}
	else
		return lhs.filter_score > rhs.filter_score;
}

bool
comp_qualscore_ref(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.confid < rhs.confid;
	//{
	//	if (lhs.ave_intron_len == rhs.ave_intron_len)
	//	{
	//		return lhs.ave_junc_mis < rhs.ave_junc_mis;
	//	}
	//	else
	//		return lhs.ave_intron_len < rhs.ave_intron_len;
	//}
	//else
	//	return lhs.filter_score > rhs.filter_score;
}

size_t
SamHandler::ReadNextTag()
{
	string line;

	int count = 0;

	sam_rec_map.clear();

	if (!sam_rec_map_next.empty())
	{
		sam_rec_map = sam_rec_map_next;

		sam_rec_map_next.clear();

		++count;
	}

	while(getline(cur_ifs, line))
	{
		if (line.empty())
			continue;

		char chromname[1000], readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
		char strand = '+';
		size_t prefixst,strand_t, incorrect, mate_offest;

		string alterstr;

		int mate_diff;

		unsigned short mis_match;

		char mate_match;

		//4803~61KTMAAXX:1:1:2785:2670#0/1	16	chr13	114542764	34	60M6676N15M	*	0	0	CGTGGCAGGCTGGGGTCCCCTTCCTATCGCAGGGGTTGGGCGTGCACTGGTCAGTCAGGTTTTGCACGCAGGTGG	BBBBBBBBBBBBBBBBa\\a[XaaY^^a^a^^^^a^aJaaa\`_b\`\Y]_^b^Tc]Y[]`Yc`cbdaccdefcf	NM:i:3
		sscanf(line.c_str(), "%s\t%llu\t%s\t%llu\t%llu\t%s\t%c\t%llu\t%d\t%s\t%s", 
			readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, &mate_match, &mate_offest, &mate_diff, chromseq, qualseq);

		size_t NMi_loc = line.find("NM:i:");
		if (NMi_loc != string::npos)
			sscanf(line.c_str() + NMi_loc, "NM:i:%hu", &mis_match);

		string tagnamestr = readname;
		size_t last_idx = tagnamestr.find_first_of("~");
		size_t tagidx = -1;

		if (last_idx != string::npos)
		{
			string tagidxstr = tagnamestr.substr(0, last_idx);

			tagidx = atoi(tagidxstr.c_str()) - 1;

			tagnamestr = tagnamestr.substr(tagnamestr.find_first_of("~") + 1);
		}

		if (tagnamestr.find_last_of("/") != string::npos)
			tagnamestr = tagnamestr.substr(0, tagnamestr.find_last_of("/"));

		string tagnamestr_temp;

		//cout << "rem non-character"<<endl;

		//for (size_t i = 0; i < tagnamestr.length(); ++i)
		//	if ((tagnamestr[i] >= '0' && tagnamestr[i] <= '9') || (tagnamestr[i] >= 'a' && tagnamestr[i] <= 'z') || (tagnamestr[i] >= 'A' && tagnamestr[i] <= 'Z'))
		//		tagnamestr_temp.push_back(tagnamestr[i]);

		//tagnamestr = tagnamestr_temp;

		//cout << "rem non-character"<<endl;

		if (strcmp(prev_tag.c_str(), tagnamestr.c_str()) > 0)
		//if (comp_str_less_than_linux(tagnamestr, prev_tag))
			cout <<"tag name not sorted\n" <<prev_tag << endl << tagnamestr << endl;

		if (prev_tag.empty())
		{
			sam_rec_map[tagnamestr].push_back(SamRec(readname, strand, chromname, prefixst, incorrect, spliceway, chromseq, mis_match, tagidx, alterstr, qualseq, mate_match, mate_offest, mate_diff, line));

			prev_tag = tagnamestr;

			++count;
		}
		else if (prev_tag == tagnamestr)
		{
			sam_rec_map[tagnamestr].push_back(SamRec(readname, strand, chromname, prefixst, incorrect, spliceway, chromseq, mis_match, tagidx, alterstr, qualseq, mate_match, mate_offest, mate_diff, line));

			++count;
		}
		else
		{
			sam_rec_map_next[tagnamestr].push_back(SamRec(readname, strand, chromname, prefixst, incorrect, spliceway, chromseq, mis_match, tagidx, alterstr, qualseq, mate_match, mate_offest, mate_diff, line));

			prev_tag = tagnamestr;

			CollectStat2();

			return count;
		}
	}

	CollectStat2();

	return count;
}

size_t
SamHandler::ReadAppendSamFile(const char *samfile)
{
	ifstream input_fs(samfile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n"); cerr << samfile<<endl; return 0;
	}

	string line;
	int count=0;
	while(getline(input_fs, line))
	{
		if (line == "")
			continue;

		char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
		unsigned short strand, something, mismatch;
		size_t offset;

		string alterstr = "";

		string qualstr = "I";

		char mate_match;

		size_t mate_offest, mate_diff;

		if (line.find("\t*\t") != string::npos || line.find("\t=\t") != string::npos)
		{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A
			size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

			if (read_count == 10)
				alterstr = alters;

			qualstr = qual_str;
		}
		else
		{
			sscanf(line.c_str(), "%s\t%hu\t%s\t%lu\t%hu\t%s\t%s\t%s\tNM:i:%hu", tagname, &strand, chrom, &offset, &something, mapped, seq, qual_str, &mismatch);

			mismatch = 0;
		}

		//TRAN00000031091:400:43	1	chr1	12489651	0	75M	AAGGTGTATGTGCAGGTGACCAAGAAAGCCGTGAGCACGAGCAGTGGAGTGTCCATCCCCGGCCCCTCCCACCAG


		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("~");
		size_t tagidx = -1;

		if (last_idx != string::npos)
		{
			string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
			tagidx = atoi(tagidxstr.c_str()) - 1;

			//name = name.substr(name.find_last_of("~") + 1);
			tagnamestr = tagnamestr.substr(tagnamestr.find_last_of("~") + 1);
		}

		if (tagnamestr.find_last_of("/"))
		{
			tagnamestr = tagnamestr.substr(0, tagnamestr.find_last_of("/"));
		}

		//cout << line << endl;

		sam_rec_map[tagnamestr].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff, line));

		//cout << sam_rec_map[tagnamestr].back().tostring() << endl;

		++count;
	}

	input_fs.close();

	return sam_rec_map.size();
}

size_t
SamHandler::ReadAppendSamFileSimp(const char *samfile)
{
	ifstream input_fs(samfile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n"); cerr << samfile<<endl; return 0;
	}

	string line;
	int count=0, unmapped = 0;
	while(getline(input_fs, line))
	{
		if (line == "" || line[0] == '@')
			continue;

		char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
		unsigned short strand, something, mismatch(0);
		size_t offset;

		string alterstr = "";

		string qualstr = "I";

		char mate_match;

		size_t mate_offest, mate_diff;

		//if (line.find("\t*\t") != string::npos || line.find("\t=\t") != string::npos)
		//{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A
		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str);

			//if (read_count == 10)
			//	alterstr = alters;

			//qualstr = qual_str;
		//}
		//else
		//{
		//	sscanf(line.c_str(), "%s\t%hu\t%s\t%lu\t%hu\t%s\t%s\t%s\tNM:i:%hu", tagname, &strand, chrom, &offset, &something, mapped, seq, qual_str, &mismatch);

		//	mismatch = 0;
		//}

		//TRAN00000031091:400:43	1	chr1	12489651	0	75M	AAGGTGTATGTGCAGGTGACCAAGAAAGCCGTGAGCACGAGCAGTGGAGTGTCCATCCCCGGCCCCTCCCACCAG

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("~");
		size_t tagidx = -1;

		if (last_idx != string::npos)
		{
			string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
			tagidx = atoi(tagidxstr.c_str()) - 1;

			//name = name.substr(name.find_last_of("~") + 1);
			tagnamestr = tagnamestr.substr(tagnamestr.find_last_of("~") + 1);
		}

		if (tagnamestr.find_last_of("/"))
		{
			tagnamestr = tagnamestr.substr(0, tagnamestr.find_last_of("/"));
		}

		//cout << line << endl;

		string spliceway = mapped;

		if (spliceway.find("*") != string::npos || strand == 4)
		{
			++unmapped;
			continue;
		}

		sam_rec_map[tagnamestr].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, "", mismatch, tagidx, "", "", mate_match, mate_offest, mate_diff, line));

		//cout << sam_rec_map[tagnamestr].back().tostring() << endl;

		++count;

		//if (count % 10000 == 0)
		//	cout << count << endl;

		//cout << line << endl;
	}

	input_fs.close();

	//cout << unmapped << endl;

	return unmapped;
}

size_t
SamHandler::CollectStat()
{
	map<string, vector<SamRec> >::iterator tag_iter;

	total_mapped_reads = sam_rec_map.size();

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		bool unique = false;
		if (tag_iter->second.size() == 1)
		{
			++unique_reads;
			tag_iter->second.front().isunique = true;

			unique = true;
		}
		else
			++multiple_reads;

		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->isexonic)
			{
				++exonic_reads;

				for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
				{
					samrec_iter->isexonic = true;
				}

				if (unique)
					++exonic_unique;
				else
					++exonic_multiple;
				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->isspliced)
			{
				++spliced_reads;

				for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
				{
					samrec_iter->isspliced = true;
				}

				if (unique)
					++spliced_unique;
				else
					++spliced_multiple;

				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->issmalldel)
			{
				++small_del;

				for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
				{
					samrec_iter->issmalldel = true;
				}

				if (unique)
					++small_del_unique;
				else
					++small_del_multiple;

				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->issmallins)
			{
				++small_insert;

				for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
				{
					samrec_iter->issmallins = true;
				}

				if (unique)
					++small_insert_unique;
				else
					++small_insert_multiple;

				break;
			}
		}

		//small_del_multiple, small_insert_multiple, small_del_unique, small_insert_unique
	}

	return sam_rec_map.size();
}

size_t
SamHandler::CollectStat2()
{
	map<string, vector<SamRec> >::iterator tag_iter;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		bool unique = false;
		if (tag_iter->second.size() == 1)
		{
			++unique_reads;

			unique = true;
		}
		else
		{
			++multiple_reads;

			map<size_t, size_t>::iterator multiple_mapped_counts_iter;

			multiple_mapped_counts_iter = multiple_mapped_counts.find(tag_iter->second.size());

			if (multiple_mapped_counts_iter == multiple_mapped_counts.end())
				multiple_mapped_counts[tag_iter->second.size()] = 1;
			else
				++multiple_mapped_counts[tag_iter->second.size()];
		}

		vector<SamRec>::iterator samrec_iter;

		bool isexonic = false;

		bool isspliced = false;

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->isexonic)
			{
				++exonic_reads;

				if (unique)
					++exonic_unique;
				else
					++exonic_multiple;
				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->isspliced)
			{
				++spliced_reads;

				if (unique)
					++spliced_unique;
				else
					++spliced_multiple;

				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->issmalldel)
			{
				++small_del;

				if (unique)
					++small_del_unique;
				else
					++small_del_multiple;

				break;
			}
		}

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->issmallins)
			{
				++small_insert;

				if (unique)
					++small_insert_unique;
				else
					++small_insert_multiple;

				break;
			}
		}

		//small_del_multiple, small_insert_multiple, small_del_unique, small_insert_unique
	}

	return sam_rec_map.size();
}

const void
SamHandler::WriteCompStat(const SamHandler& synthetic_sam, const char* output_file)
{
	ofstream ofs(output_file);

	ofs<<'\t'<<base_file_name<<"_only"<<'\t'<<base_file_name<<"_notmatched"<<'\t'<<base_file_name<<"_same_start_or_end"<<'\t'<<base_file_name<<"_matched"<<'\t'
		<<synthetic_sam.base_file_name<<"_only"<<'\t'<<synthetic_sam.base_file_name<<"_notmatched"<<'\t'<<synthetic_sam.base_file_name<<"_same_start_or_end"<<'\t'<<synthetic_sam.base_file_name<<"_matched"<<endl;

	ofs<<"exonic\t"<<exonic_only<<'\t'<<exonic_notshared<<'\t'<<exonic_same_start_or_end<<'\t'<<exonic_shared<<'\t'<< synthetic_sam.exonic_only<<'\t'<< synthetic_sam.exonic_notshared<<'\t'<< synthetic_sam.exonic_same_start_or_end<<'\t'<<synthetic_sam.exonic_shared<<endl;

	ofs<<"spliced\t"<<spliced_only<<'\t'<<spliced_notshared<<'\t'<<spliced_same_start_or_end<<'\t'<<spliced_shared<<'\t'<< synthetic_sam.spliced_only<<'\t'<< synthetic_sam.spliced_notshared<<'\t'<< synthetic_sam.spliced_same_start_or_end<<'\t'<<synthetic_sam.spliced_shared<<endl;

	ofs<<"unique\t"<<unique_only<<'\t'<<unique_notshared<<'\t'<<unique_same_start_or_end<<'\t'<<unique_shared<<'\t'<< synthetic_sam.unique_only<<'\t'<< synthetic_sam.unique_notshared<<'\t'<< synthetic_sam.unique_same_start_or_end<<'\t'<<synthetic_sam.unique_shared<<endl;

	ofs<<"multiple\t"<<multiple_only<<'\t'<<multiple_notshared<<'\t'<<multiple_same_start_or_end<<'\t'<<multiple_shared<<'\t'<< synthetic_sam.multiple_only<<'\t'<< synthetic_sam.multiple_notshared<<'\t'<< synthetic_sam.multiple_same_start_or_end<<'\t'<<synthetic_sam.multiple_shared<<endl;

	map<size_t, size_t>::iterator multiple_mapped_counts_iter;

	for (multiple_mapped_counts_iter = multiple_mapped_counts.begin(); multiple_mapped_counts_iter != multiple_mapped_counts.end(); ++multiple_mapped_counts_iter)
		ofs << multiple_mapped_counts_iter->first  <<'\t'<<multiple_mapped_counts_iter->second  <<endl;

	//multiple_mapped_counts_iter = multiple_mapped_counts.find(tag_iter->second.size());
}

const void
SamHandler::WriteStat(const char* output_file, int unmapped)
{
	ofstream ofs(output_file);

	ofs<<"exonic\t"<<exonic_reads<<endl;

	ofs <<"exonic unique\t"<<exonic_unique<<endl;

	ofs <<"exonic multiple\t"<<exonic_multiple<<endl; 

	ofs<<"spliced\t"<<spliced_reads << endl;

	ofs<<"spliced unique\t"<<spliced_unique << endl;

	ofs<<"spliced multiple\t"<<spliced_multiple << endl;

	ofs<<"small insert\t"<<small_insert << endl;

	ofs<<"small insert unique\t"<<small_insert_unique << endl;

	ofs<<"small insert multiple\t"<<small_insert_multiple << endl;

	ofs<<"small delete\t"<<small_del << endl;

	ofs<<"small delete unique\t"<<small_del_unique << endl;

	ofs<<"small delete multiple\t"<<small_del_multiple << endl;

	ofs<<"unique\t"<<unique_reads<<endl;

	ofs<<"multiple\t"<<multiple_reads<<endl;

	ofs << "unmapped\t"<<unmapped << endl;

	ofs<<"total\t"<<total_mapped_reads<<endl;
}

const void
SamHandler::WriteStat(const char* output_file)
{
	ofstream ofs(output_file);

	ofs<<"exonic\t"<<exonic_reads<<endl;

	ofs <<"exonic unique\t"<<exonic_unique<<endl;

	ofs <<"exonic multiple\t"<<exonic_multiple<<endl; 

	ofs<<"spliced\t"<<spliced_reads << endl;

	ofs<<"spliced unique\t"<<spliced_unique << endl;

	ofs<<"spliced multiple\t"<<spliced_multiple << endl;

	ofs<<"small insert\t"<<small_insert << endl;

	ofs<<"small insert unique\t"<<small_insert_unique << endl;

	ofs<<"small insert multiple\t"<<small_insert_multiple << endl;

	ofs<<"small delete\t"<<small_del << endl;

	ofs<<"small delete unique\t"<<small_del_unique << endl;

	ofs<<"small delete multiple\t"<<small_del_multiple << endl;

	ofs<<"unique\t"<<unique_reads<<endl;

	ofs<<"multiple\t"<<multiple_reads<<endl;

	//ofs << "unmapped\t"<<unmapped << endl;

	ofs<<"total\t"<<total_mapped_reads<<endl;
}

size_t
SamHandler::WriteSamFile(const char* samfile)
{
	ofstream ofs(samfile);
	map<string, vector<SamRec> >::iterator tag_iter;

	size_t count = 0;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			ofs << samrec_iter->tostring() << endl;
			++count;
		}
	}

	ofs.close();

	return count;
}

size_t
SamHandler::Sam2Junc(const char* juncfile, size_t read_width, string chrom_dir, size_t min_intron, size_t max_intron)
{
	ofstream ofs(juncfile);

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	CONJ_HASH_MAP conj_hash_map;

	map<string, vector<SamRec> >::iterator tag_iter;

	cout << "insert junction"<<endl;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		vector<SamRec>::iterator sam_iter;

		for (sam_iter = tag_iter->second.begin(); sam_iter != tag_iter->second.end(); ++sam_iter)
		{
			if (sam_iter->spliceway_vec.size() > 1 && sam_iter->wrong_format == false)
			{
				vector<pair<size_t, int> >::iterator vp_iter;

				for (vp_iter =  sam_iter->spliceway_vec.begin(); vp_iter != sam_iter->spliceway_vec.end() - 1; ++vp_iter)
				{
					size_t prefixst = vp_iter->first;

					size_t prefixend = vp_iter->first + vp_iter->second - 1;

					size_t suffixst = (vp_iter + 1)->first;

					size_t prefixlen = vp_iter->second;

					size_t suffixlen = (vp_iter + 1)->second;

					string chrom_ID = sam_iter->chrom_name;

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
						iter_conj_pre->second.inc_hits(prefixlen, suffixlen, sam_iter->tagidx);
					}	  
					else 
					{
						(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(prefixlen, suffixlen, read_width, sam_iter->tagidx)));
					}
				}
			}
		}
	}

	size_t junc_count;

	cout << "write junction"<<endl;

	WriteJunc(conj_hash_map, ofs, read_width, 0, chrom_dir);

	ofs.close();

	cout << "sort junction"<<endl;

	junc_count = SortJunc(juncfile, chrom_dir, min_intron, max_intron);

	return junc_count;
}

size_t
SamHandler::WriteUniqueMappedTagSam(const char* unique_samfile)
{
	ofstream ofs(unique_samfile);
	map<string, vector<SamRec> >::iterator tag_iter;

	size_t count = 0;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		if (tag_iter->second.size() == 1)
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
			{
				ofs << samrec_iter->tostring() << endl;
				++count;
			}
		}
	}

	ofs.close();

	return count;
}

size_t
SamHandler::WriteMultipleMappedTagSam(const char* multiple_samfile)
{
	ofstream ofs(multiple_samfile);
	map<string, vector<SamRec> >::iterator tag_iter;

	size_t count = 0;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		if (tag_iter->second.size() > 1)
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
			{
				ofs << samrec_iter->tostring() << endl;
				++count;
			}
		}
	}

	ofs.close();

	return count;
}

size_t
SamHandler::FilterByMisMatchByRead(map<size_t, vector<SamRec> >& sam_rec_map_filter_mismatch, ofstream* ofs_notminmis, ofstream* ofs_exceed)
{
	map<size_t, vector<SamRec> >::iterator mismatch_iter;

	vector<SamRec>::iterator sam_iter;

	size_t count = 0;

	if (sam_rec_map_filter_mismatch.size() > 1)
	{
		map<size_t, vector<SamRec> >::iterator sam_rec_iter = sam_rec_map_filter_mismatch.begin();

		++sam_rec_iter;

		map<size_t, vector<SamRec> >::iterator cur_sam_rec_iter;

		for (cur_sam_rec_iter = sam_rec_iter; cur_sam_rec_iter != sam_rec_map_filter_mismatch.end(); ++cur_sam_rec_iter)
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = cur_sam_rec_iter->second.begin(); samrec_iter != cur_sam_rec_iter->second.end(); ++samrec_iter)
			{
				(*ofs_notminmis)<<samrec_iter->tostring()<<endl;
			}
		}

		sam_rec_map_filter_mismatch.erase(sam_rec_iter, sam_rec_map_filter_mismatch.end());
	}

	return count;
}

size_t
SamHandler::FilterByQualScoreByRead(vector<SamRec>& sam_rec_map_filter_mismatch, ofstream* ofs_notminqualscore, ofstream* ofs_exceed)
{
	sort(sam_rec_map_filter_mismatch.begin(), sam_rec_map_filter_mismatch.end(), comp_qualscore_ref);

	//vector<SamRec>::iterator sam_iter;

	size_t count = 0;

	//if (sam_rec_map_filter_mismatch.size() > 1)
	//{
	//	map<size_t, vector<SamRec> >::iterator sam_rec_iter = sam_rec_map_filter_mismatch.begin();

	//	++sam_rec_iter;

		//map<size_t, vector<SamRec> >::iterator cur_sam_rec_iter;

	//for (cur_sam_rec_iter = sam_rec_iter; cur_sam_rec_iter != sam_rec_map_filter_mismatch.end(); ++cur_sam_rec_iter)
	//{
	//	vector<SamRec>::iterator samrec_iter;

	for (size_t i = 1; i < sam_rec_map_filter_mismatch.size(); ++i)
	{
		if (sam_rec_map_filter_mismatch[i].confid > sam_rec_map_filter_mismatch[i - 1].confid)
		{
			for (size_t j = i; j < sam_rec_map_filter_mismatch.size(); ++j)
				(*ofs_notminqualscore)<<sam_rec_map_filter_mismatch[j].tostring()<<endl;

			sam_rec_map_filter_mismatch.resize(i);

			return i;
		}		
	}
	//}

	//	sam_rec_map_filter_mismatch.erase(sam_rec_iter, sam_rec_map_filter_mismatch.end());
	//}

	return sam_rec_map_filter_mismatch.size();
}

size_t
SamHandler::FilterByMisMatch(map<string, map<size_t, vector<SamRec*> > >& sam_rec_map_filter_mismatch, size_t maxhits, const char* filtered_file)
{
	string filtered_by_notminmis = filtered_file;filtered_by_notminmis.append(".notminmis");ofstream ofs_notminmis(filtered_by_notminmis.c_str());

	string filtered_by_exceed = filtered_file;filtered_by_exceed.append(".exceed");ofstream ofs_exceed (filtered_by_exceed.c_str());

	map<string, vector<SamRec> >::iterator tag_iter;

	size_t count = 0;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		if (tag_iter->second.size() > 1 && tag_iter->second.size() <= maxhits)
		{
			vector<SamRec>::iterator samrec_iter;

			bool inserted = false;
			for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
			{
				if (samrec_iter->spliceway_vec.size() > 0)
				{
					sam_rec_map_filter_mismatch[tag_iter->first][samrec_iter->mis_match].push_back(&(*samrec_iter));

					inserted = true;
				}
			}

			if (inserted && sam_rec_map_filter_mismatch[tag_iter->first].size() > 1)
			{
				map<size_t, vector<SamRec*> >::iterator sam_rec_iter = sam_rec_map_filter_mismatch[tag_iter->first].begin();

				++sam_rec_iter;

				map<size_t, vector<SamRec*> >::iterator cur_sam_rec_iter;

				for (cur_sam_rec_iter = sam_rec_iter; cur_sam_rec_iter != sam_rec_map_filter_mismatch[tag_iter->first].end(); ++cur_sam_rec_iter)
				{
					vector<SamRec*>::iterator samrec_iter;

					for (samrec_iter = cur_sam_rec_iter->second.begin(); samrec_iter != cur_sam_rec_iter->second.end(); ++samrec_iter)
					{
						ofs_notminmis<<(*samrec_iter)->tostring()<<endl;
					}
				}

				sam_rec_map_filter_mismatch[tag_iter->first].erase(sam_rec_iter, sam_rec_map_filter_mismatch[tag_iter->first].end());
			}

			if (inserted)
				count += sam_rec_map_filter_mismatch[tag_iter->first].begin()->second.size();
		}
		else if (tag_iter->second.size() > maxhits)
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
			{
				ofs_exceed<<samrec_iter->tostring()<<endl;
			}
		}
	}

	return count;
}

size_t
SamHandler::ReadSortedJunction(const char* juncfile, map<string, vector<JuncForSort> >& junc_map)
{
	ifstream ifs(juncfile);

	size_t count = 0;

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

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], flankchr[100];
			int juncst, juncend, prefixend, suffixst, kinds, hits, flankcase;

			float rank, lpq, il;

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%f\t%d\t%s\t%f\t%f\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, &rank, &flankcase, flankchr, &il, &lpq, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			junc_map[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset, rank, lpq, flankchr, flankcase, il, juncname, min_mis, max_mis, ave_mis));

			++count;
		}

		ifs.close();
	}

	return count;
}

bool comp_sam_rec_maplen(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.mappedlen > rhs.mappedlen;
}

void
SamHandler::MarkCanonNoncanonByReads(vector<SamRec>& read_sam, map<string, vector<JuncForSort> >& junc_map)
{
	sort(read_sam.begin(), read_sam.end(), comp_sam_rec_maplen);

	vector<SamRec >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin() + 1; sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		if (sam_rec_iter->mappedlen != (sam_rec_iter - 1)->mappedlen)
		{
			read_sam.resize(sam_rec_iter - read_sam.begin());

			break;
		}
	}

	for (sam_rec_iter = read_sam.begin(); sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		//cout << (*sam_rec_iter)->tostring()<<endl;
		if ((sam_rec_iter)->spliceway_vec.size() == 1)
			continue;

		vector<pair<size_t, int> >::iterator oft_mpl_iter;

		//cout << "calculating "<<endl;

		for (oft_mpl_iter = (sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
		{
			if (oft_mpl_iter->second < 0 || (oft_mpl_iter + 1)->second < 0)
			{
				sam_rec_iter->is_insert = true;

				continue;
			}

			JuncForSort jfs(oft_mpl_iter->first + oft_mpl_iter->second - 1, (oft_mpl_iter + 1)->first, 0, 0, "0", "0", 0, 0, "0", 0, 0, "0", 0, 0, 0);

			vector<JuncForSort>::iterator lower_bound_iter, upper_bound_iter;

			lower_bound_iter = lower_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

			upper_bound_iter = upper_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

			if (lower_bound_iter == junc_map[(sam_rec_iter)->chrom_name].end() || lower_bound_iter->juncst != jfs.juncst || lower_bound_iter->juncend != jfs.juncend)
			{
				cout << "not same junction"<<endl;

				cout << (sam_rec_iter)->tostring() << endl;

				cout << lower_bound_iter->juncname<<endl;

				continue;

				getchar();
			}

			if (upper_bound_iter - lower_bound_iter != 1)
			{
				cout << "more than 1 junc"<<endl;

				cout << (sam_rec_iter)->tostring() << endl;

				cout << lower_bound_iter->juncname<<endl;

				continue;

				getchar();
			}

			(sam_rec_iter)->junc_id.push_back(lower_bound_iter->juncname);

			if (lower_bound_iter->flankcase >= 5)
				(sam_rec_iter)->canon_count++;
			else
				(sam_rec_iter)->noncanon_count++;
		}

		(sam_rec_iter)->canon_rate = (double)(sam_rec_iter)->canon_count / (double)((sam_rec_iter)->canon_count + (sam_rec_iter)->noncanon_count);

		//++count;
	}
}

size_t
SamHandler::MarkCanonNoncanon(map<string, vector<JuncForSort> >& junc_map)
{
	map<string, vector<SamRec> >::iterator sm_iter; 

	size_t count = 0;

	for (sm_iter = sam_rec_map.begin(); sm_iter != sam_rec_map.end(); ++sm_iter)
	{
		vector<SamRec >::iterator sam_rec_iter;

		//cout << "process sam"<<endl;

		for (sam_rec_iter = sm_iter->second.begin(); sam_rec_iter != sm_iter->second.end(); ++sam_rec_iter)
		{
			//cout << (*sam_rec_iter)->tostring()<<endl;
			if ((sam_rec_iter)->spliceway_vec.size() == 1)
				continue;

			vector<pair<size_t, int> >::iterator oft_mpl_iter;

			//cout << "calculating "<<endl;

			for (oft_mpl_iter = (sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
			{
				JuncForSort jfs(oft_mpl_iter->first + oft_mpl_iter->second - 1, (oft_mpl_iter + 1)->first, 0, 0, "0", "0", 0, 0, "0", 0, 0, "0", 0, 0, 0);

				vector<JuncForSort>::iterator lower_bound_iter, upper_bound_iter;

				lower_bound_iter = lower_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

				upper_bound_iter = upper_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

				if (lower_bound_iter->juncst != jfs.juncst || lower_bound_iter->juncend != jfs.juncend)
				{
					cout << "not same junction"<<endl;

					cout << (sam_rec_iter)->tostring() << endl;

					cout << lower_bound_iter->juncname<<endl;

					continue;

					getchar();
				}

				if (upper_bound_iter - lower_bound_iter != 1)
				{
					cout << "more than 1 junc"<<endl;

					cout << (sam_rec_iter)->tostring() << endl;

					cout << lower_bound_iter->juncname<<endl;

					continue;

					getchar();
				}

				(sam_rec_iter)->junc_id.push_back(lower_bound_iter->juncname);

				if (lower_bound_iter->flankcase)
					(sam_rec_iter)->canon_count++;
				else
					(sam_rec_iter)->noncanon_count++;
			}

			(sam_rec_iter)->canon_rate = (double)(sam_rec_iter)->canon_count / (double)((sam_rec_iter)->canon_count + (sam_rec_iter)->noncanon_count);

			++count;
		}

		//cout << "sorting"<<endl;
		//sort(sm_iter->second.begin()->second.begin(), sm_iter->second.begin()->second.end(), comp_filterscore);

		//cout << "sorting finshed"<<endl;

		
	}

	return count;
}

size_t 
SamHandler::CalcFilterScoreByRead(vector<SamRec>& sam_rec_map_filter_mismatch, map<string, vector<JuncForSort> >& junc_map, 
		double entrpy_weight, double pqlen_weight, double ave_mis_weight)
{
	size_t count = 0;

	vector<SamRec>::iterator sam_rec_iter;

	//cout << "process sam"<<endl;

	for (sam_rec_iter = sam_rec_map_filter_mismatch.begin(); sam_rec_iter != sam_rec_map_filter_mismatch.end(); ++sam_rec_iter)
	{
		//cout << (*sam_rec_iter)->tostring()<<endl;
		if ((sam_rec_iter)->spliceway_vec.size() == 1)
			continue;

		vector<pair<size_t, int> >::iterator oft_mpl_iter;

		double sumscore = 0;

		size_t sum_intron_len = 0;

		double sum_ave_mis = 0;

		//cout << "calculating "<<endl;

		for (oft_mpl_iter = (sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
		{
			if (oft_mpl_iter->second <=  0 || (oft_mpl_iter + 1)->second <= 0)
				continue;

			(sam_rec_iter)->isspliced = true;

			JuncForSort jfs(oft_mpl_iter->first + oft_mpl_iter->second - 1, (oft_mpl_iter + 1)->first, 0, 0, "0", "0", 0, 0, "0", 0, 0, "0", 0, 0, 0);

			vector<JuncForSort>::iterator lower_bound_iter, upper_bound_iter;

			lower_bound_iter = lower_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

			upper_bound_iter = upper_bound(junc_map[(sam_rec_iter)->chrom_name].begin(), junc_map[(sam_rec_iter)->chrom_name].end(), jfs, compjunc);

			if (lower_bound_iter->juncst != jfs.juncst || lower_bound_iter->juncend != jfs.juncend)
				cout << "not same junction"<<endl;
			if (upper_bound_iter - lower_bound_iter != 1)
			{
				cout << "more than 1 junc"<<endl;

				cout << (sam_rec_iter)->tostring() << endl;

				cout << lower_bound_iter->juncname<<endl;

				continue;

				getchar();
			}

			sumscore += (lower_bound_iter->rank * entrpy_weight) + (lower_bound_iter->lpq * pqlen_weight) + (lower_bound_iter->ave_mismatch * ave_mis_weight);

			sum_intron_len += lower_bound_iter->juncend - lower_bound_iter->juncst;

			sum_ave_mis += lower_bound_iter->ave_mismatch;

			(sam_rec_iter)->junc_id.push_back(lower_bound_iter->juncname);
		}

		sumscore = sumscore / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_intron_len =  sum_intron_len / ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_ave_mis = sum_ave_mis / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		//cout << "assigning "<<endl;
		(sam_rec_iter)->filter_score = sumscore;

		(sam_rec_iter)->ave_intron_len = sum_intron_len;

		(sam_rec_iter)->ave_junc_mis = sum_ave_mis;

		++count;
	}

	//cout << "sorting"<<endl;
	sort(sam_rec_map_filter_mismatch.begin(), sam_rec_map_filter_mismatch.end(), comp_filterscore_ref);

	//cout << "sorting finshed"<<endl;

	return count;
}

size_t 
SamHandler::CalcFilterScore(map<string, map<size_t, vector<SamRec*> >  >& sam_rec_map_filter_mismatch, map<string, vector<JuncForSort> >& junc_map, 
		double entrpy_weight, double pqlen_weight, double ave_mis_weight)
{
	map<string, map<size_t, vector<SamRec* > > >::iterator sm_iter; 

	size_t count = 0;

	for (sm_iter = sam_rec_map_filter_mismatch.begin(); sm_iter != sam_rec_map_filter_mismatch.end(); ++sm_iter)
	{
		if (sm_iter->second.size() > 1)
		{
			cout << "not filtered "<<sm_iter->second.size() <<endl;

		}
		else if (sm_iter->second.size() == 0)
		{
			//cout <<"0 hits"<<endl;
			cout << sm_iter->first<<endl;
			continue;
		}
		

		vector<SamRec* >::iterator sam_rec_iter;

		//cout << "process sam"<<endl;

		for (sam_rec_iter = sm_iter->second.begin()->second.begin(); sam_rec_iter != sm_iter->second.begin()->second.end(); ++sam_rec_iter)
		{
			//cout << (*sam_rec_iter)->tostring()<<endl;
			if ((*sam_rec_iter)->spliceway_vec.size() == 1)
				continue;

			vector<pair<size_t, int> >::iterator oft_mpl_iter;

			double sumscore = 0;

			size_t sum_intron_len = 0;

			double sum_ave_mis = 0;

			//cout << "calculating "<<endl;

			for (oft_mpl_iter = (*sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (*sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
			{
				if (oft_mpl_iter->second <=  0 || (oft_mpl_iter + 1)->second <= 0)
					continue;

				JuncForSort jfs(oft_mpl_iter->first + oft_mpl_iter->second - 1, (oft_mpl_iter + 1)->first, 0, 0, "0", "0", 0, 0, "0", 0, 0, "0", 0, 0, 0);

				vector<JuncForSort>::iterator lower_bound_iter, upper_bound_iter;

				lower_bound_iter = lower_bound(junc_map[(*sam_rec_iter)->chrom_name].begin(), junc_map[(*sam_rec_iter)->chrom_name].end(), jfs, compjunc);

				upper_bound_iter = upper_bound(junc_map[(*sam_rec_iter)->chrom_name].begin(), junc_map[(*sam_rec_iter)->chrom_name].end(), jfs, compjunc);

				if (lower_bound_iter->juncst != jfs.juncst || lower_bound_iter->juncend != jfs.juncend)
					cout << "not same junction"<<endl;
				if (upper_bound_iter - lower_bound_iter != 1)
				{
					cout << "more than 1 junc"<<endl;

					cout << (*sam_rec_iter)->tostring() << endl;

					cout << lower_bound_iter->juncname<<endl;

					continue;

					getchar();
				}

				sumscore += (lower_bound_iter->rank * entrpy_weight) + (lower_bound_iter->lpq * pqlen_weight) + (lower_bound_iter->ave_mismatch * ave_mis_weight);

				sum_intron_len += lower_bound_iter->juncend - lower_bound_iter->juncst;

				sum_ave_mis += lower_bound_iter->ave_mismatch;

				(*sam_rec_iter)->junc_id.push_back(lower_bound_iter->juncname);
			}

			sumscore = sumscore / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_intron_len =  sum_intron_len / ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_ave_mis = sum_ave_mis / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//cout << "assigning "<<endl;
			(*sam_rec_iter)->filter_score = sumscore;

			(*sam_rec_iter)->ave_intron_len = sum_intron_len;

			(*sam_rec_iter)->ave_junc_mis = sum_ave_mis;

			++count;
		}

		//cout << "sorting"<<endl;
		sort(sm_iter->second.begin()->second.begin(), sm_iter->second.begin()->second.end(), comp_filterscore);

		//cout << "sorting finshed"<<endl;

		
	}

	return count;
}

size_t
SamHandler::FilterCanonNoncaonByReads(const char* sam_file, const char* filtered_canon_file, const char* filtered_noncanon_file, const char* fitlered_noncanon_canon_file,
									  const char* filtered_ins_file, const char* stat_file, map<string, vector<JuncForSort> >& junc_map)
{
	ofstream ofs_filtered_canon(filtered_canon_file);

	ofstream ofs_filtered_noncanon(filtered_noncanon_file);

	ofstream ofs_fitlered_noncanon_canon(fitlered_noncanon_canon_file);

	ofstream ofs_fitlered_ins(filtered_ins_file);

	string unspliced_reads = sam_file; unspliced_reads.append(".unspliced.sam");

	ofstream& unspliced_ofs = ofs_filtered_canon;//(unspliced_reads.c_str());

	ifstream input_fs(sam_file);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n"); cerr << sam_file<<endl; return 0;
	}

	string line;
	int count=0, tag_count = 0, unspliced = 0, filtered_canon = 0, filtered_noncanon = 0, filtered_canon_noncanon = 0;

	vector<SamRec> cur_read_sam;

	string prev_tagname = "", prev_unspliced = "";

	while(getline(input_fs, line))
	{
		if (line == "")
			continue;

		char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
		unsigned short strand, something, mismatch;
		size_t offset;

		char mate_match;

		size_t mate_offest, mate_diff;

		string alterstr = "";

		string qualstr = "I";

		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

		++count;

		string spliceway = mapped;

		//if (spliceway.find("N") == string::npos)
		//{
		//	unspliced_ofs << line << endl;

		//	if (prev_unspliced != tagname)
		//	{
		//		++unspliced;
		//		prev_unspliced = tagname;				
		//	}

		//	continue;
		//}

		if (read_count == 10)
			alterstr = alters;

		qualstr = qual_str;

		size_t tagidx = 0;

		if (prev_tagname.empty() || prev_tagname == tagname)
		{
			cur_read_sam.push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));

			if (prev_tagname.empty())
				prev_tagname = tagname;

			tag_count++;
		}
		else
		{
			//filter
			if (tag_count == 1)
			{
				ofs_filtered_canon << cur_read_sam.front().tostring()<<endl;

				++filtered_canon;

				//++count;
			}
			else
			{
				MarkCanonNoncanonByReads(cur_read_sam, junc_map);

				vector<SamRec>::iterator samrec_iter;

				bool canon = false, noncanon = false, is_insert = false, not_insert = false;

				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (samrec_iter->canon_count)
						canon =true;

					if (samrec_iter->noncanon_count)
						noncanon =true;

					if (samrec_iter->is_insert)
						is_insert = true;
					else
						not_insert = true;
				}

				int insert_model;

				if (is_insert && not_insert)
					insert_model = 0;
				else if (!is_insert && not_insert)
					insert_model = 1;
				else if (is_insert && !not_insert)
					insert_model = -1;
				else
				{
					insert_model = -2;

					cout << "not possible"<<endl;
				}

				if (canon && !noncanon)
				{
					if (insert_model == 0)
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							if (samrec_iter->is_insert == false)
								ofs_filtered_canon << samrec_iter->tostring()<<endl;
							else 
								ofs_fitlered_ins << samrec_iter->tostring()<<endl;
						}
					}
					else
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							ofs_filtered_canon << samrec_iter->tostring()<<endl;
						}
					}

					++filtered_canon;
				}
				else if (!canon && noncanon)
				{
					if (insert_model == 0)
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							if (samrec_iter->is_insert == false)
								ofs_filtered_canon << samrec_iter->tostring()<<endl;
							else 
								ofs_fitlered_ins << samrec_iter->tostring()<<endl;
						}
					}
					else
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							ofs_filtered_canon << samrec_iter->tostring()<<endl;
						}
					}

					++filtered_canon;
				}
				else if (canon && noncanon)
				{
					double max_canon_rate = 0;

					bool noncan_can = false, noncan = false, unspliced = false;

					if (insert_model == 0)
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							if (max_canon_rate < samrec_iter->canon_rate && samrec_iter->is_insert == false)
								max_canon_rate = samrec_iter->canon_rate;
						}

						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							if (samrec_iter->is_insert)
							{
								ofs_fitlered_ins << samrec_iter->tostring()<<endl;
							}
							else if (max_canon_rate == samrec_iter->canon_rate)
							{
								ofs_filtered_canon << samrec_iter->tostring()<<endl;
								//++count;
							}
							else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
							{
								ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl;

								noncan_can = true;
							}
							else if (samrec_iter->noncanon_count)
							{
								ofs_filtered_noncanon << samrec_iter->tostring()<<endl;

								noncan = true;
							}
							else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
							{
								//cout << "should not be here one reads"<<endl;
								unspliced_ofs << samrec_iter->tostring()<<endl;

								unspliced = true;
							}

						}
					}
					else
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							if (max_canon_rate < samrec_iter->canon_rate)
								max_canon_rate = samrec_iter->canon_rate;
						}

						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{

							if (max_canon_rate == samrec_iter->canon_rate)
							{
								ofs_filtered_canon << samrec_iter->tostring()<<endl;
								//++count;
							}
							else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
							{
								ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl;

								noncan_can = true;
							}
							else if (samrec_iter->noncanon_count)
							{
								ofs_filtered_noncanon << samrec_iter->tostring()<<endl;

								noncan = true;
							}
							else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
							{
								//cout << "should not be here one reads"<<endl;
								unspliced_ofs << samrec_iter->tostring()<<endl;

								unspliced = true;
							}

						}
					}

					++filtered_canon;

					if (noncan_can)
						++filtered_canon_noncanon;

					if (noncan)
						++filtered_noncanon;

					if (unspliced)
						++unspliced;
				}
				else
				{
					if (insert_model == 0)
					{

						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							//cout << "should not be here"<<endl;
							if (!samrec_iter->is_insert)
							{
								unspliced_ofs << samrec_iter->tostring()<<endl;
							}
							else
								ofs_fitlered_ins << samrec_iter->tostring()<<endl;
						}
					}
					else
					{
						for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
						{
							unspliced_ofs << samrec_iter->tostring()<<endl;
						}
					}

					++unspliced;
				}
			}

			cur_read_sam.clear();
			cur_read_sam.push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));
			prev_tagname = tagname;
			tag_count = 1;
		}
	}

	if (tag_count == 1)
	{
		ofs_filtered_canon << cur_read_sam.front().tostring()<<endl;

		++filtered_canon;

		//++count;
	}
	else
	{
		MarkCanonNoncanonByReads(cur_read_sam, junc_map);

		vector<SamRec>::iterator samrec_iter;

		bool canon = false, noncanon = false, is_insert = false, not_insert = false;

		for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
		{
			if (samrec_iter->canon_count)
				canon =true;

			if (samrec_iter->noncanon_count)
				noncanon =true;

			if (samrec_iter->is_insert)
				is_insert = true;
			else
				not_insert = true;
		}

		int insert_model;

		if (is_insert && not_insert)
			insert_model = 0;
		else if (!is_insert && not_insert)
			insert_model = 1;
		else if (is_insert && !not_insert)
			insert_model = -1;
		else
		{
			insert_model = -2;

			cout << "not possible"<<endl;
		}

		if (canon && !noncanon)
		{
			if (insert_model == 0)
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (samrec_iter->is_insert == false)
						ofs_filtered_canon << samrec_iter->tostring()<<endl;
					else 
						ofs_fitlered_ins << samrec_iter->tostring()<<endl;
				}
			}
			else
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					ofs_filtered_canon << samrec_iter->tostring()<<endl;
				}
			}

			++filtered_canon;
		}
		else if (!canon && noncanon)
		{
			if (insert_model == 0)
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (samrec_iter->is_insert == false)
						ofs_filtered_canon << samrec_iter->tostring()<<endl;
					else 
						ofs_fitlered_ins << samrec_iter->tostring()<<endl;
				}
			}
			else
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					ofs_filtered_canon << samrec_iter->tostring()<<endl;
				}
			}

			++filtered_canon;
		}
		else if (canon && noncanon)
		{
			double max_canon_rate = 0;

			bool noncan_can = false, noncan = false, unspliced = false;

			if (insert_model == 0)
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (max_canon_rate < samrec_iter->canon_rate && samrec_iter->is_insert == false)
						max_canon_rate = samrec_iter->canon_rate;
				}

				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (samrec_iter->is_insert)
					{
						ofs_fitlered_ins << samrec_iter->tostring()<<endl;
					}
					else if (max_canon_rate == samrec_iter->canon_rate)
					{
						ofs_filtered_canon << samrec_iter->tostring()<<endl;
						//++count;
					}
					else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
					{
						ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl;

						noncan_can = true;
					}
					else if (samrec_iter->noncanon_count)
					{
						ofs_filtered_noncanon << samrec_iter->tostring()<<endl;

						noncan = true;
					}
					else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
					{
						//cout << "should not be here one reads"<<endl;
						unspliced_ofs << samrec_iter->tostring()<<endl;

						unspliced = true;
					}

				}
			}
			else
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					if (max_canon_rate < samrec_iter->canon_rate)
						max_canon_rate = samrec_iter->canon_rate;
				}

				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{

					if (max_canon_rate == samrec_iter->canon_rate)
					{
						ofs_filtered_canon << samrec_iter->tostring()<<endl;
						//++count;
					}
					else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
					{
						ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl;

						noncan_can = true;
					}
					else if (samrec_iter->noncanon_count)
					{
						ofs_filtered_noncanon << samrec_iter->tostring()<<endl;

						noncan = true;
					}
					else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
					{
						//cout << "should not be here one reads"<<endl;
						unspliced_ofs << samrec_iter->tostring()<<endl;

						unspliced = true;
					}

				}
			}

			++filtered_canon;

			if (noncan_can)
				++filtered_canon_noncanon;

			if (noncan)
				++filtered_noncanon;

			if (unspliced)
				++unspliced;
		}
		else
		{
			if (insert_model == 0)
			{

				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					//cout << "should not be here"<<endl;
					if (!samrec_iter->is_insert)
					{
						unspliced_ofs << samrec_iter->tostring()<<endl;
					}
					else
						ofs_fitlered_ins << samrec_iter->tostring()<<endl;
				}
			}
			else
			{
				for (samrec_iter = cur_read_sam.begin(); samrec_iter != cur_read_sam.end(); ++samrec_iter)
				{
					unspliced_ofs << samrec_iter->tostring()<<endl;
				}
			}

			++unspliced;
		}
	}

	ofstream stat_ofs(stat_file, ios_base::app);
	
	stat_ofs << "FilterCanonNoncanon:"<<endl << "total_mapped: "<< count<<endl<< "total_mapped_reads: "<<unspliced + filtered_canon + filtered_noncanon + filtered_canon_noncanon<<endl
		<< "unspliced_reads: "<<unspliced<<endl <<"filtered_canon_reads: " <<  filtered_canon<<endl<<"filtered_noncanon_reads: "<< filtered_noncanon<<endl
		<<"filtered_canon_noncanon_reads: " << filtered_canon_noncanon<<endl;

	return count;
}

size_t
SamHandler::FilterCanonNoncaonTag(const char* sam_file, const char* filtered_canon, const char* filtered_noncanon, const char* fitlered_noncanon_canon)
{
	ofstream ofs_filtered_canon(filtered_canon);

	ofstream ofs_filtered_noncanon(filtered_noncanon);

	ofstream ofs_fitlered_noncanon_canon(fitlered_noncanon_canon);

	string unspliced_reads = sam_file; unspliced_reads.append(".unspliced.sam");

	ofstream unspliced_ofs(unspliced_reads.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		bool canon = false, noncanon = false;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->canon_count)
				canon =true;

			if (samrec_iter->noncanon_count)
				noncanon =true;
		}

		if (canon && !noncanon)
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				if (samrec_iter->spliceway_vec.size() == 1)
					unspliced_ofs << samrec_iter->tostring()<<endl;
				else
					ofs_filtered_canon << samrec_iter->tostring()<<endl;
				++count;
			}
		}
		else if (!canon && noncanon)
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				if (samrec_iter->spliceway_vec.size() == 1)
					unspliced_ofs << samrec_iter->tostring()<<endl;
				else
					ofs_filtered_canon << samrec_iter->tostring()<<endl;
				++count;
			}
		}
		else if (canon && noncanon)
		{
			double max_canon_rate = 0;

			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				if (max_canon_rate < samrec_iter->canon_rate)
					max_canon_rate = samrec_iter->canon_rate;
			}

			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				if (max_canon_rate == samrec_iter->canon_rate)
				{
					ofs_filtered_canon << samrec_iter->tostring()<<endl;
					++count;
				}
				else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
				{
					ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl;
				}
				else if (samrec_iter->noncanon_count)
				{
					ofs_filtered_noncanon << samrec_iter->tostring()<<endl;
				}
				else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
				{
					unspliced_ofs << samrec_iter->tostring()<<endl;
					//cout << "should not be here one reads"<<endl;
				}
			}
		}
		else
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				unspliced_ofs << samrec_iter->tostring()<<endl;
				//cout << "should not be here"<<endl;
			}
		}

		++lhs_tag_iter;

	}

	return count;
}

size_t
SamHandler::FilterMultipleMappedTagByRead(const char* sam_file, const char* unique_samfile, const char* filtered_samfile, const char* junc_file, 
										  size_t maxhits, const char* stat_file, double entrpy_weight, double pqlen_weight, double ave_mis_weight)
{
	map<size_t, vector<SamRec> > sam_rec_map_filter_mismatch;

	map<string, vector<JuncForSort> > junc_map;

	cout << "ReadSortedJunction"<<endl; 

	ReadSortedJunction(junc_file, junc_map);

	string filtered_by_samescore = filtered_samfile;filtered_by_samescore.append(".same_score");ofstream ofs_same_score(filtered_by_samescore.c_str());

	string filtered_by_notminmis = filtered_samfile;filtered_by_notminmis.append(".notminmis");ofstream ofs_notminmis (filtered_by_notminmis.c_str());

	string filtered_by_notminqual = filtered_samfile;filtered_by_notminqual.append(".notminqual");ofstream ofs_notminqual (filtered_by_notminqual.c_str());

	string filtered_by_maxhits = filtered_samfile;filtered_by_maxhits.append(".maxhits");ofstream ofs_maxhits (filtered_by_maxhits.c_str());

	ofstream ofs_unique (unique_samfile);

	ofstream ofs_filtered_multiple(filtered_samfile);

	ifstream input_fs(sam_file);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n"); cerr << sam_file<<endl; return 0;
	}

	string line;
	int count=0, tag_count = 0, multipled_mapped = 0, exceed_maxhits = 0, unique_mapped = 0;

	string prev_tagname = "";

	while(getline(input_fs, line))
	{
		if (line == "")
			continue;

		char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
		unsigned short strand, something, mismatch;
		size_t offset;

		string alterstr = "";

		string qualstr = "I";

		char mate_match;

		size_t mate_offest, mate_diff;

		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

		++count;

		if (read_count == 10)
			alterstr = alters;

		qualstr = qual_str;

		size_t tagidx = 0;

		if (prev_tagname.empty() || prev_tagname == tagname)
		{
			sam_rec_map_filter_mismatch[mismatch].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));

			if (prev_tagname.empty())
				prev_tagname = tagname;

			tag_count++;
		}
		else
		{
			//filter

			if (tag_count == 1)
			{
				map<size_t, vector<SamRec> >::iterator mismatch_iter;
				
				for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
				{
					vector<SamRec>::iterator sam_iter;

					for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
					{
						ofs_unique << sam_iter->tostring() << endl;
					}
				}

				++unique_mapped;
			}
			else if (tag_count > maxhits)
			{
				map<size_t, vector<SamRec> >::iterator mismatch_iter;
				
				for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
				{
					vector<SamRec>::iterator sam_iter;

					for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
					{
						ofs_maxhits << sam_iter->tostring() << endl;
					}
				}

				++exceed_maxhits;
			}
			else
			{
				//cout << "FilterByMisMatch"<<endl; 

				FilterByMisMatchByRead(sam_rec_map_filter_mismatch, &ofs_notminmis, &ofs_maxhits);

				//cout << "CalcFilterScore"<<endl; 

				vector<SamRec>& vec_sam_rec = sam_rec_map_filter_mismatch.begin()->second;

				//if (vec_sam_rec.size() > 1)
				//	FilterByQualScoreByRead(vec_sam_rec, &ofs_notminqual, &ofs_maxhits);

				CalcFilterScoreByRead(vec_sam_rec, junc_map, entrpy_weight, pqlen_weight, ave_mis_weight);

				//cout << "start select "<<endl;

				if (vec_sam_rec.size() > 1)
				{
					vec_sam_rec[0].best = 1;
					ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

					//++count;

					for (size_t idx = 1; idx < vec_sam_rec.size(); ++idx)
					{
						if (vec_sam_rec[idx - 1].filter_score == vec_sam_rec[idx].filter_score &&
							vec_sam_rec[idx - 1].ave_intron_len == vec_sam_rec[idx].ave_intron_len &&
							vec_sam_rec[idx - 1].ave_junc_mis == vec_sam_rec[idx].ave_junc_mis)
						{
							//vec_sam_rec[idx - 1]->best = 1;
							ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

							//++count;
						}
						else
						{
							//tag_iter->second.begin()->second.front()->best = 1;
							//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

							//++count;

							for (; idx < vec_sam_rec.size(); ++idx)
							{
								if (vec_sam_rec[idx].isexonic)
									ofs_filtered_multiple << vec_sam_rec[idx].tostring() << endl;
								else
								{
									ofs_same_score<< vec_sam_rec[idx].tostring()<<'\t'<<vec_sam_rec[idx].filter_score<<'\t' 
										<< vec_sam_rec[idx].ave_intron_len<<'\t' <<vec_sam_rec[idx].ave_junc_mis<<'\t';

									vector<string>::iterator junc_id_iter;

									for (junc_id_iter = vec_sam_rec[idx].junc_id.begin(); junc_id_iter != vec_sam_rec[idx].junc_id.end(); ++junc_id_iter)
										ofs_same_score<<*junc_id_iter<<'\t';

									ofs_same_score <<endl;
								}
							}

							break;
						}
					}
				}
				else
				{
					vec_sam_rec[0].best = 1;
					ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

					//++count;
				}

				++multipled_mapped;
			}

			sam_rec_map_filter_mismatch.clear();
			sam_rec_map_filter_mismatch[mismatch].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));
			prev_tagname = tagname;
			tag_count = 1;
			
		}

		//++count;
	}

	//filter

	if (tag_count == 1)
	{
		map<size_t, vector<SamRec> >::iterator mismatch_iter;

		for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
		{
			vector<SamRec>::iterator sam_iter;

			for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
			{
				ofs_unique << sam_iter->tostring() << endl;
			}
		}

		++unique_mapped;
	}
	else if (tag_count > maxhits)
	{
		map<size_t, vector<SamRec> >::iterator mismatch_iter;

		for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
		{
			vector<SamRec>::iterator sam_iter;

			for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
			{
				ofs_maxhits << sam_iter->tostring() << endl;
			}
		}

		++exceed_maxhits;
	}
	else
	{
		cout << "FilterByMisMatch"<<endl; 

		FilterByMisMatchByRead(sam_rec_map_filter_mismatch, &ofs_notminmis, &ofs_maxhits);

		cout << "CalcFilterScore"<<endl; 

		vector<SamRec>& vec_sam_rec = sam_rec_map_filter_mismatch.begin()->second;

		CalcFilterScoreByRead(vec_sam_rec, junc_map, entrpy_weight, pqlen_weight, ave_mis_weight);

		cout << "start select "<<endl;

		if (vec_sam_rec.size() > 1)
		{
			vec_sam_rec[0].best = 1;
			ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

			//++count;

			for (size_t idx = 1; idx < vec_sam_rec.size(); ++idx)
			{
				if (vec_sam_rec[idx - 1].filter_score == vec_sam_rec[idx].filter_score &&
					vec_sam_rec[idx - 1].ave_intron_len == vec_sam_rec[idx].ave_intron_len &&
					vec_sam_rec[idx - 1].ave_junc_mis == vec_sam_rec[idx].ave_junc_mis)
				{
					//vec_sam_rec[idx - 1]->best = 1;
					ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

					//++count;
				}
				else
				{
					//tag_iter->second.begin()->second.front()->best = 1;
					//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

					//++count;

					for (; idx < vec_sam_rec.size(); ++idx)
					{
						if (vec_sam_rec[idx].isexonic)
							ofs_filtered_multiple << vec_sam_rec[idx].tostring();
						else
						{

							ofs_same_score<< vec_sam_rec[idx].tostring()<<'\t'<<vec_sam_rec[idx].filter_score<<'\t' 
								<< vec_sam_rec[idx].ave_intron_len<<'\t' <<vec_sam_rec[idx].ave_junc_mis<<'\t';

							vector<string>::iterator junc_id_iter;

							for (junc_id_iter = vec_sam_rec[idx].junc_id.begin(); junc_id_iter != vec_sam_rec[idx].junc_id.end(); ++junc_id_iter)
								ofs_same_score<<*junc_id_iter<<'\t';

							ofs_same_score <<endl;
						}
					}

					break;
				}
			}
		}
		else
		{
			vec_sam_rec[0].best = 1;
			ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

			//++count;
		}

		++multipled_mapped;
	}

	//for (tag_iter = sam_rec_map_filter_mismatch.begin(); tag_iter != sam_rec_map_filter_mismatch.end(); ++tag_iter)
	//{
	//	if (tag_iter->second.begin()->second.size() > 1)
	//	{
	//		tag_iter->second.begin()->second[0]->best = 1;
	//		ofs<< tag_iter->second.begin()->second[0]->tostring()<<endl;

	//		++count;
	//					
	//		for (size_t idx = 1; idx < tag_iter->second.begin()->second.size(); ++idx)
	//		{
	//			if (tag_iter->second.begin()->second[idx - 1]->filter_score == tag_iter->second.begin()->second[idx]->filter_score &&
	//				tag_iter->second.begin()->second[idx - 1]->ave_intron_len == tag_iter->second.begin()->second[idx]->ave_intron_len &&
	//				tag_iter->second.begin()->second[idx - 1]->ave_junc_mis == tag_iter->second.begin()->second[idx]->ave_junc_mis)
	//			{
	//				bool matched = true;

	//				//if (tag_iter->second.begin()->second[idx - 1]->chrom_name != tag_iter->second.begin()->second[idx]->chrom_name ||
	//				//	tag_iter->second.begin()->second[idx - 1]->splice_way != tag_iter->second.begin()->second[idx]->splice_way ||
	//				//	tag_iter->second.begin()->second[idx - 1]->start != tag_iter->second.begin()->second[idx]->start)
	//				//	matched = false;

	//				if (matched)
	//				{
	//					tag_iter->second.begin()->second[idx - 1]->best = 1;
	//					ofs<< tag_iter->second.begin()->second[idx - 1]->tostring()<<endl;

	//					++count;
	//				}
	//				else
	//				{
	//					vector<SamRec* >::iterator samrec_ptr_iter;
	//					for (samrec_ptr_iter = tag_iter->second.begin()->second.begin(); samrec_ptr_iter != tag_iter->second.begin()->second.end(); ++samrec_ptr_iter)
	//					{
	//						ofs_same_score<< (*samrec_ptr_iter)->tostring()<<'\t'<<(*samrec_ptr_iter)->filter_score<<'\t' 
	//							<< (*samrec_ptr_iter)->ave_intron_len<<'\t' <<(*samrec_ptr_iter)->ave_junc_mis<<'\t';

	//						vector<string>::iterator junc_id_iter;

	//						for (junc_id_iter = (*samrec_ptr_iter)->junc_id.begin(); junc_id_iter != (*samrec_ptr_iter)->junc_id.end(); ++junc_id_iter)
	//							ofs_same_score<<*junc_id_iter<<'\t';

	//						ofs_same_score <<endl;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				//tag_iter->second.begin()->second.front()->best = 1;
	//				//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

	//				//++count;

	//				break;
	//			}
	//		}
	//	}
	//	else
	//	{
	//		tag_iter->second.begin()->second.front()->best = 1;
	//		ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

	//		++count;
	//	}
	//}

	input_fs.close();

	ofstream stat_ofs(stat_file, ios_base::app);

	stat_ofs << "Filter Multiple Mapped: "<<endl<<"total_mapped: " << count<<endl << "total_mapped_reads: "<< multipled_mapped + exceed_maxhits + unique_mapped << endl
		<<"spliced_multipled_mapped_reads: " <<multipled_mapped<<endl <<"spliced_exceed_maxhits_reads: " <<exceed_maxhits<< endl<<"spliced_unique_mapped_reads: " <<unique_mapped <<endl;

	return count;
}


size_t
SamHandler::FilterMultipleMappedTagByReadFilterByRepeat(const char* sam_file, const char* unique_samfile, const char* filtered_samfile, const char* junc_file, 
										  size_t maxhits, const char* stat_file, double entrpy_weight, double pqlen_weight, double ave_mis_weight)
{
	map<size_t, vector<SamRec> > sam_rec_map_filter_mismatch;

	map<string, vector<JuncForSort> > junc_map;

	cout << "ReadSortedJunction"<<endl; 

	ReadSortedJunction(junc_file, junc_map);

	string filtered_by_samescore = filtered_samfile;filtered_by_samescore.append(".same_score");ofstream ofs_same_score(filtered_by_samescore.c_str());

	string filtered_by_notminmis = filtered_samfile;filtered_by_notminmis.append(".notminmis");ofstream ofs_notminmis (filtered_by_notminmis.c_str());

	string filtered_by_maxhits = filtered_samfile;filtered_by_maxhits.append(".maxhits");ofstream ofs_maxhits (filtered_by_maxhits.c_str());

	ofstream ofs_unique (unique_samfile);

	ofstream ofs_filtered_multiple(filtered_samfile);

	ifstream input_fs(sam_file);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n"); cerr << sam_file<<endl; return 0;
	}

	string line;
	int count=0, tag_count = 0, multipled_mapped = 0, exceed_maxhits = 0, unique_mapped = 0;

	string prev_tagname = "";

	while(getline(input_fs, line))
	{
		if (line == "")
			continue;

		char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
		unsigned short strand, something, mismatch;
		size_t offset;

		char mate_match;

		size_t mate_offest, mate_diff;

		string alterstr = "";

		string qualstr = "I";

		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

		++count;

		if (read_count == 10)
			alterstr = alters;

		qualstr = qual_str;

		size_t tagidx = 0;

		if (prev_tagname.empty() || prev_tagname == tagname)
		{
			sam_rec_map_filter_mismatch[mismatch].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));

			if (prev_tagname.empty())
				prev_tagname = tagname;

			tag_count++;
		}
		else
		{
			//filter

			if (tag_count == 1)
			{
				map<size_t, vector<SamRec> >::iterator mismatch_iter;
				
				for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
				{
					vector<SamRec>::iterator sam_iter;

					for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
					{
						ofs_unique << sam_iter->tostring() << endl;
					}
				}

				++unique_mapped;
			}
			else if (tag_count > maxhits)
			{
				map<size_t, vector<SamRec> >::iterator mismatch_iter;
				
				for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
				{
					vector<SamRec>::iterator sam_iter;

					for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
					{
						ofs_maxhits << sam_iter->tostring() << endl;
					}
				}

				++exceed_maxhits;
			}
			else
			{
				//cout << "FilterByMisMatch"<<endl; 

				FilterByMisMatchByRead(sam_rec_map_filter_mismatch, &ofs_notminmis, &ofs_maxhits);

				//cout << "CalcFilterScore"<<endl; 

				vector<SamRec>& vec_sam_rec = sam_rec_map_filter_mismatch.begin()->second;

				CalcFilterScoreByRead(vec_sam_rec, junc_map, entrpy_weight, pqlen_weight, ave_mis_weight);

				//cout << "start select "<<endl;

				if (vec_sam_rec.size() > 1)
				{
					if (vec_sam_rec[0].filter_score == vec_sam_rec[1].filter_score &&
						vec_sam_rec[0].ave_intron_len == vec_sam_rec[1].ave_intron_len &&
						vec_sam_rec[0].ave_junc_mis == vec_sam_rec[1].ave_junc_mis)
					{
						//vec_sam_rec[idx - 1]->best = 1;
						//ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

						//++count;
					}
					else 
					{
						vec_sam_rec[0].best = 1;
						ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

						//++count;

						for (size_t idx = 1; idx < vec_sam_rec.size(); ++idx)
						{
							if (vec_sam_rec[idx - 1].filter_score == vec_sam_rec[idx].filter_score &&
								vec_sam_rec[idx - 1].ave_intron_len == vec_sam_rec[idx].ave_intron_len &&
								vec_sam_rec[idx - 1].ave_junc_mis == vec_sam_rec[idx].ave_junc_mis)
							{
								//vec_sam_rec[idx - 1]->best = 1;
								ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

								//++count;
							}
							else
							{
								//tag_iter->second.begin()->second.front()->best = 1;
								//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

								//++count;

								for (; idx < vec_sam_rec.size(); ++idx)
								{
									ofs_same_score<< vec_sam_rec[idx].tostring()<<'\t'<<vec_sam_rec[idx].filter_score<<'\t' 
										<< vec_sam_rec[idx].ave_intron_len<<'\t' <<vec_sam_rec[idx].ave_junc_mis<<'\t';

									vector<string>::iterator junc_id_iter;

									for (junc_id_iter = vec_sam_rec[idx].junc_id.begin(); junc_id_iter != vec_sam_rec[idx].junc_id.end(); ++junc_id_iter)
										ofs_same_score<<*junc_id_iter<<'\t';

									ofs_same_score <<endl;
								}

								break;
							}
						}
					}
				}
				else
				{
					vec_sam_rec[0].best = 1;
					ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

					//++count;
				}

				++multipled_mapped;
			}

			sam_rec_map_filter_mismatch.clear();
			sam_rec_map_filter_mismatch[mismatch].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, tagidx, alterstr, qualstr, mate_match, mate_offest, mate_diff));
			prev_tagname = tagname;
			tag_count = 1;
			
		}

		//++count;
	}

	//filter

	if (tag_count == 1)
	{
		map<size_t, vector<SamRec> >::iterator mismatch_iter;

		for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
		{
			vector<SamRec>::iterator sam_iter;

			for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
			{
				ofs_unique << sam_iter->tostring() << endl;
			}
		}

		++unique_mapped;
	}
	else if (tag_count > maxhits)
	{
		map<size_t, vector<SamRec> >::iterator mismatch_iter;

		for (mismatch_iter = sam_rec_map_filter_mismatch.begin(); mismatch_iter != sam_rec_map_filter_mismatch.end(); ++mismatch_iter)
		{
			vector<SamRec>::iterator sam_iter;

			for (sam_iter = mismatch_iter->second.begin(); sam_iter != mismatch_iter->second.end(); ++sam_iter)
			{
				ofs_maxhits << sam_iter->tostring() << endl;
			}
		}

		++exceed_maxhits;
	}
	else
	{
		cout << "FilterByMisMatch"<<endl; 

		FilterByMisMatchByRead(sam_rec_map_filter_mismatch, &ofs_notminmis, &ofs_maxhits);

		cout << "CalcFilterScore"<<endl; 

		vector<SamRec>& vec_sam_rec = sam_rec_map_filter_mismatch.begin()->second;

		CalcFilterScoreByRead(vec_sam_rec, junc_map, entrpy_weight, pqlen_weight, ave_mis_weight);

		cout << "start select "<<endl;

		if (vec_sam_rec.size() > 1)
		{
			if (vec_sam_rec[0].filter_score == vec_sam_rec[1].filter_score &&
				vec_sam_rec[0].ave_intron_len == vec_sam_rec[1].ave_intron_len &&
				vec_sam_rec[0].ave_junc_mis == vec_sam_rec[1].ave_junc_mis)
			{
				//vec_sam_rec[idx - 1]->best = 1;
				//ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

				//++count;
			}
			else 
			{
				vec_sam_rec[0].best = 1;
				ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

				//++count;

				for (size_t idx = 1; idx < vec_sam_rec.size(); ++idx)
				{
					if (vec_sam_rec[idx - 1].filter_score == vec_sam_rec[idx].filter_score &&
						vec_sam_rec[idx - 1].ave_intron_len == vec_sam_rec[idx].ave_intron_len &&
						vec_sam_rec[idx - 1].ave_junc_mis == vec_sam_rec[idx].ave_junc_mis)
					{
						//vec_sam_rec[idx - 1]->best = 1;
						ofs_filtered_multiple<< vec_sam_rec[idx].tostring()<<endl;

						//++count;
					}
					else
					{
						//tag_iter->second.begin()->second.front()->best = 1;
						//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

						//++count;

						for (; idx < vec_sam_rec.size(); ++idx)
						{
							ofs_same_score<< vec_sam_rec[idx].tostring()<<'\t'<<vec_sam_rec[idx].filter_score<<'\t' 
								<< vec_sam_rec[idx].ave_intron_len<<'\t' <<vec_sam_rec[idx].ave_junc_mis<<'\t';

							vector<string>::iterator junc_id_iter;

							for (junc_id_iter = vec_sam_rec[idx].junc_id.begin(); junc_id_iter != vec_sam_rec[idx].junc_id.end(); ++junc_id_iter)
								ofs_same_score<<*junc_id_iter<<'\t';

							ofs_same_score <<endl;
						}

						break;
					}
				}
			}
		}
		else
		{
			vec_sam_rec[0].best = 1;
			ofs_filtered_multiple<< vec_sam_rec[0].tostring()<<endl;

			//++count;
		}

		++multipled_mapped;
	}

	//for (tag_iter = sam_rec_map_filter_mismatch.begin(); tag_iter != sam_rec_map_filter_mismatch.end(); ++tag_iter)
	//{
	//	if (tag_iter->second.begin()->second.size() > 1)
	//	{
	//		tag_iter->second.begin()->second[0]->best = 1;
	//		ofs<< tag_iter->second.begin()->second[0]->tostring()<<endl;

	//		++count;
	//					
	//		for (size_t idx = 1; idx < tag_iter->second.begin()->second.size(); ++idx)
	//		{
	//			if (tag_iter->second.begin()->second[idx - 1]->filter_score == tag_iter->second.begin()->second[idx]->filter_score &&
	//				tag_iter->second.begin()->second[idx - 1]->ave_intron_len == tag_iter->second.begin()->second[idx]->ave_intron_len &&
	//				tag_iter->second.begin()->second[idx - 1]->ave_junc_mis == tag_iter->second.begin()->second[idx]->ave_junc_mis)
	//			{
	//				bool matched = true;

	//				//if (tag_iter->second.begin()->second[idx - 1]->chrom_name != tag_iter->second.begin()->second[idx]->chrom_name ||
	//				//	tag_iter->second.begin()->second[idx - 1]->splice_way != tag_iter->second.begin()->second[idx]->splice_way ||
	//				//	tag_iter->second.begin()->second[idx - 1]->start != tag_iter->second.begin()->second[idx]->start)
	//				//	matched = false;

	//				if (matched)
	//				{
	//					tag_iter->second.begin()->second[idx - 1]->best = 1;
	//					ofs<< tag_iter->second.begin()->second[idx - 1]->tostring()<<endl;

	//					++count;
	//				}
	//				else
	//				{
	//					vector<SamRec* >::iterator samrec_ptr_iter;
	//					for (samrec_ptr_iter = tag_iter->second.begin()->second.begin(); samrec_ptr_iter != tag_iter->second.begin()->second.end(); ++samrec_ptr_iter)
	//					{
	//						ofs_same_score<< (*samrec_ptr_iter)->tostring()<<'\t'<<(*samrec_ptr_iter)->filter_score<<'\t' 
	//							<< (*samrec_ptr_iter)->ave_intron_len<<'\t' <<(*samrec_ptr_iter)->ave_junc_mis<<'\t';

	//						vector<string>::iterator junc_id_iter;

	//						for (junc_id_iter = (*samrec_ptr_iter)->junc_id.begin(); junc_id_iter != (*samrec_ptr_iter)->junc_id.end(); ++junc_id_iter)
	//							ofs_same_score<<*junc_id_iter<<'\t';

	//						ofs_same_score <<endl;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				//tag_iter->second.begin()->second.front()->best = 1;
	//				//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

	//				//++count;

	//				break;
	//			}
	//		}
	//	}
	//	else
	//	{
	//		tag_iter->second.begin()->second.front()->best = 1;
	//		ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

	//		++count;
	//	}
	//}

	input_fs.close();

	ofstream stat_ofs(stat_file, ios_base::app);

	stat_ofs << "Filter Multiple Mapped: "<<endl<<"total_mapped: " << count<<endl << "total_mapped_reads: "<< multipled_mapped + exceed_maxhits + unique_mapped << endl
		<<"spliced_multipled_mapped_reads: " <<multipled_mapped<<endl <<"spliced_exceed_maxhits_reads: " <<exceed_maxhits<< endl<<"spliced_unique_mapped_reads: " <<unique_mapped <<endl;

	return count;
}


size_t
SamHandler::FilterMultipleMappedTag(const char* filtered_samfile, const char* junc_file, size_t maxhits, double entrpy_weight, double pqlen_weight, double ave_mis_weight)
{
	map<string, map<size_t, vector<SamRec* > > > sam_rec_map_filter_mismatch;

	map<string, vector<JuncForSort> > junc_map;

	string filtered_by_samescore = filtered_samfile;filtered_by_samescore.append(".same_score");ofstream ofs_same_score(filtered_by_samescore.c_str());

	cout << "FilterByMisMatch"<<endl; 

	FilterByMisMatch(sam_rec_map_filter_mismatch, maxhits, filtered_samfile);

	cout << "ReadSortedJunction"<<endl; 

	ReadSortedJunction(junc_file, junc_map);

	cout << "CalcFilterScore"<<endl; 

	CalcFilterScore(sam_rec_map_filter_mismatch, junc_map, entrpy_weight, pqlen_weight, ave_mis_weight);

	ofstream ofs(filtered_samfile);

	cout << "start select "<<endl;

	map<string, map<size_t, vector<SamRec* > > >::iterator tag_iter;

	size_t count = 0;

	for (tag_iter = sam_rec_map_filter_mismatch.begin(); tag_iter != sam_rec_map_filter_mismatch.end(); ++tag_iter)
	{
		if (tag_iter->second.begin()->second.size() > 1)
		{
			tag_iter->second.begin()->second[0]->best = 1;
			ofs<< tag_iter->second.begin()->second[0]->tostring()<<endl;

			++count;
						
			for (size_t idx = 1; idx < tag_iter->second.begin()->second.size(); ++idx)
			{
				if (tag_iter->second.begin()->second[idx - 1]->filter_score == tag_iter->second.begin()->second[idx]->filter_score &&
					tag_iter->second.begin()->second[idx - 1]->ave_intron_len == tag_iter->second.begin()->second[idx]->ave_intron_len &&
					tag_iter->second.begin()->second[idx - 1]->ave_junc_mis == tag_iter->second.begin()->second[idx]->ave_junc_mis)
				{
					bool matched = true;

					//if (tag_iter->second.begin()->second[idx - 1]->chrom_name != tag_iter->second.begin()->second[idx]->chrom_name ||
					//	tag_iter->second.begin()->second[idx - 1]->splice_way != tag_iter->second.begin()->second[idx]->splice_way ||
					//	tag_iter->second.begin()->second[idx - 1]->start != tag_iter->second.begin()->second[idx]->start)
					//	matched = false;

					if (matched)
					{
						tag_iter->second.begin()->second[idx - 1]->best = 1;
						ofs<< tag_iter->second.begin()->second[idx - 1]->tostring()<<endl;

						++count;
					}
					else
					{
						vector<SamRec* >::iterator samrec_ptr_iter;
						for (samrec_ptr_iter = tag_iter->second.begin()->second.begin(); samrec_ptr_iter != tag_iter->second.begin()->second.end(); ++samrec_ptr_iter)
						{
							ofs_same_score<< (*samrec_ptr_iter)->tostring()<<'\t'<<(*samrec_ptr_iter)->filter_score<<'\t' 
								<< (*samrec_ptr_iter)->ave_intron_len<<'\t' <<(*samrec_ptr_iter)->ave_junc_mis<<'\t';

							vector<string>::iterator junc_id_iter;

							for (junc_id_iter = (*samrec_ptr_iter)->junc_id.begin(); junc_id_iter != (*samrec_ptr_iter)->junc_id.end(); ++junc_id_iter)
								ofs_same_score<<*junc_id_iter<<'\t';

							ofs_same_score <<endl;
						}
					}
				}
				else
				{
					//tag_iter->second.begin()->second.front()->best = 1;
					//ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

					//++count;

					break;
				}
			}
		}
		else
		{
			tag_iter->second.begin()->second.front()->best = 1;
			ofs<< tag_iter->second.begin()->second.front()->tostring()<<endl;

			++count;
		}
	}

	return count;
}

pair<size_t, size_t>
SamHandler::SepSplicedUnspliced(const char* spliced, const char* unspliced)
{
	ofstream ofs_spliced(spliced);

	ofstream ofs_unspliced(unspliced);

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count_unspliced = 0, count_spliced = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			if (samrec_iter->spliceway_vec.size() == 1)
			{
				ofs_unspliced<<samrec_iter->tostring()<<endl;
				++count_unspliced;
			}
			else
			{
				ofs_spliced<<samrec_iter->tostring()<<endl;
				++count_spliced;
			}
		}

		++lhs_tag_iter;
	}

	return make_pair(count_unspliced, count_spliced);

}

size_t
SamHandler::FilterByAnchor(const char* sam_file, size_t min_anchor)
{
	string sam_file_shortanchor = sam_file; sam_file_shortanchor.append(".shortanchor"); ofstream ofs_shortanchor(sam_file_shortanchor.c_str());

	string sam_file_longanchor = sam_file; sam_file_longanchor.append(".longanchor"); ofstream ofs_longanchor(sam_file_longanchor.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			vector<pair<size_t, int> >::iterator vp_iter;

			bool shortanchor = false;
			
			for (vp_iter = samrec_iter->spliceway_vec.begin(); vp_iter != samrec_iter->spliceway_vec.end(); ++vp_iter)
			{
				if (vp_iter->second < min_anchor)
				{
					shortanchor = true;
					break;
				}
			}

			if (shortanchor)
			{
				ofs_shortanchor << samrec_iter->tostring()<<endl;

				++count;
			}
			else
			{
				ofs_longanchor << samrec_iter->tostring()<<endl;
			}

		}

		++lhs_tag_iter;

	}

	return count;
}

void
GetFastaReads(const char* reads_file, map<string, int>& read_names)
{
	static const size_t input_buffer_size = 1000000;

	ifstream in(reads_file);
	if (!in) 
	{
		cout << "cannot open input file " <<reads_file <<endl;
	}

	string s, name;
	bool first_line = true;
	while (!in.eof())
	{
		string buffer;
		getline(in, buffer);
		if (buffer == "")
			continue;

		if (in.gcount() == static_cast<int>(input_buffer_size))
			cout << "Line in " << name << "\nexceeds max length: " << input_buffer_size;

		// correct for dos carriage returns before newlines
		if (buffer[strlen(buffer.c_str()) - 1] == '\r')
			buffer[strlen(buffer.c_str()) - 1] = '\0';

		if (buffer[0] == '>')
		{
			if (first_line == false && s.length() > 0)
			{
				while (name.find_last_of("~") != string::npos)
				{
					name = name.substr(name.find_last_of("~") + 1);
				}
				read_names[name] = 1;
				//m_reads.push_back(s);
			}
			else first_line = false;
			name = buffer;
			name = name.substr(name.find_first_not_of(">"));
			s = "";
		}
		else s.append(buffer);
	}

	if (!first_line && s.length() > 0)
	{
		while (name.find_last_of("~") != string::npos)
		{
			name = name.substr(name.find_last_of("~") + 1);
		}

		if (name.find_last_of("/") != string::npos)
			name = name.substr(0, name.find_last_of("/"));
		read_names[name] = 1;
		//m_reads.push_back(s);
	}
}

void
GetFastqReads(const char* reads_file, map<string, int>& read_names)
{
	static const size_t input_buffer_size = 1000000;

	ifstream in(reads_file);
	if (!in) 
	{
		cout << "cannot open input file " <<reads_file <<endl;
	}

	string s, name;

	while (!in.eof())
	{
		string buffer, buffer2, buffer3, buffer4;
		getline(in, buffer);

		getline(in, buffer2);

		getline(in, buffer3);

		getline(in, buffer4);

		if (buffer == "" || buffer2 == "" || buffer3 == "" || buffer4 == "")
			continue;

		// correct for dos carriage returns before newlines
		if (buffer[strlen(buffer.c_str()) - 1] == '\r')
			buffer[strlen(buffer.c_str()) - 1] = '\0';

		name = buffer;

		while (name.find_last_of("~") != string::npos)
		{
			name = name.substr(name.find_last_of("~") + 1);
		}

		name = name.substr(name.find_first_not_of("@"));

		if (name.find_last_of("/") != string::npos)
			name = name.substr(0, name.find_last_of("/"));

		read_names[name] = 1;
	}
}

size_t
SamHandler::FilterByRepeated(const char* repeated_reads, size_t fq, const char* sam_file)
{
	map<string, int> read_names;
	if (fq)
		GetFastqReads(repeated_reads, read_names);
	else
		GetFastaReads(repeated_reads, read_names);

	string sam_file_inrepeat = sam_file; sam_file_inrepeat.append(".inrepeat"); ofstream ofs_inrepeat(sam_file_inrepeat.c_str());

	string sam_file_notinrepeat = sam_file; sam_file_notinrepeat.append(".notinrepeat"); ofstream ofs_notinrepeat(sam_file_notinrepeat.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		ofstream* ofs_ptr;

		bool inrepeat = false;

		if (read_names.find(lhs_tag_iter->first) != read_names.end())
		{
			ofs_ptr = &ofs_inrepeat;

			inrepeat = true;
			//count += lhs_tag_iter->second.size();;
		}
		else
			ofs_ptr = &ofs_notinrepeat;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			*ofs_ptr << samrec_iter->tostring()<<endl;
		}

		++lhs_tag_iter;

		if (inrepeat)
			++count;

	}

	return count;
}

size_t
SamHandler::FilterByUnmapped(const char* unmapped_reads, size_t fq, const char* sam_file)
{
	map<string, int> read_names;
	if (fq)
		GetFastqReads(unmapped_reads, read_names);
	else
		GetFastaReads(unmapped_reads, read_names);

	string sam_file_inunmapped = sam_file; sam_file_inunmapped.append(".inunmapped"); ofstream ofs_inunmapped(sam_file_inunmapped.c_str());

	string sam_file_notinunmapped = sam_file; sam_file_notinunmapped.append(".notinunmapped"); ofstream ofs_notinunmapped(sam_file_notinunmapped.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		ofstream* ofs_ptr;

		if (read_names.find(lhs_tag_iter->first) != read_names.end())
		{
			ofs_ptr = &ofs_inunmapped;
		}
		else
		{
			ofs_ptr = &ofs_notinunmapped;
			count += lhs_tag_iter->second.size();;
		}

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			*ofs_ptr << samrec_iter->tostring()<<endl;
		}

		++lhs_tag_iter;

	}

	return count;
}

size_t 
SamHandler::FilterByIntronLenHmer(const char* sam_file, size_t intronlen, size_t seg_len)
{
	string sam_file_inintronlenhmer = sam_file; sam_file_inintronlenhmer.append(".inintronlenhmer"); ofstream ofs_inintronlenhmer(sam_file_inintronlenhmer.c_str());

	string sam_file_notinintronlenhmer = sam_file; sam_file_notinintronlenhmer.append(".notinintronlenhmer"); ofstream ofs_notinintronlenhmer(sam_file_notinintronlenhmer.c_str());

	//string sam_file_smallexon = sam_file; sam_file_smallexon.append(".smallexon"); ofstream ofs_smallexon(sam_file_smallexon.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			vector<pair<size_t, int> >::iterator vp_iter;

			for (vp_iter = samrec_iter->spliceway_vec.begin(); vp_iter != samrec_iter->spliceway_vec.end() - 1; ++vp_iter)
			{
				if (vp_iter->first + vp_iter->second + intronlen < (vp_iter + 1)->first)
				{
					ofs_notinintronlenhmer << samrec_iter->tostring()<<endl;
					++count;

					break;
				}
			}

			if (vp_iter == samrec_iter->spliceway_vec.end() - 1)
				ofs_inintronlenhmer << samrec_iter->tostring()<<endl;

		}

		++lhs_tag_iter;

	}

	return count;
}

size_t
SamHandler::FilterByUnsplicedMapped(SamHandler& unsplice_sam, char* sam_file)
{
	string sam_file_inunsplicedmapped = sam_file; sam_file_inunsplicedmapped.append(".inunsplicedmapped"); ofstream ofs_inunsplicedmapped(sam_file_inunsplicedmapped.c_str());

	string sam_file_notinunsplicedmapped = sam_file; sam_file_notinunsplicedmapped.append(".notinunsplicedmapped"); ofstream ofs_notinunsplicedmapped(sam_file_notinunsplicedmapped.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		if (unsplice_sam.sam_rec_map.find(lhs_tag_iter->first) != unsplice_sam.sam_rec_map.end())
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				ofs_inunsplicedmapped << samrec_iter->tostring()<<endl;
				++count;
			}
		}
		else
		{
			vector<SamRec>::iterator samrec_iter;

			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				ofs_notinunsplicedmapped  << samrec_iter->tostring()<<endl;
			}
		}

		++lhs_tag_iter;
		
	}
	return count;
}

size_t
SamHandler::FilterBySmallExon(const char* sam_file, size_t seg_len)
{
	string sam_file_insmallexon = sam_file; sam_file_insmallexon.append(".insmallexon"); ofstream ofs_insmallexon(sam_file_insmallexon.c_str());

	string sam_file_notinsmallexon = sam_file; sam_file_notinsmallexon.append(".notinsmallexon"); ofstream ofs_notinsmallexon(sam_file_notinsmallexon.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			//vector<pair<size_t, size_t> >::iterator vp_iter;

			//bool anynotinisland = false;
			if (samrec_iter->spliceway_vec.size() > 2)
			{
				bool is_smallexon = false;
				for (size_t i = 1; i < samrec_iter->spliceway_vec.size() - 2; ++i)
				{
					if (samrec_iter->spliceway_vec[i].second < seg_len)
						is_smallexon = true;
				}
	
				if (is_smallexon)
				{
					ofs_insmallexon << samrec_iter->tostring()<<endl;
					++count;
				}
				else
					ofs_notinsmallexon << samrec_iter->tostring()<<endl;
			}
			else
			{
				ofs_notinsmallexon << samrec_iter->tostring()<<endl;
			}
		}

		++lhs_tag_iter;

	}

	return count;
}

bool paircomp ( const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs)
{
	return lhs.first < rhs.first;
}

int
readisland(const char* islands, map<string, vector<pair<size_t, size_t> > >& islandmap)
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

			char chrom[100], group[100], islandname[100];
			size_t start, end;

			sscanf(line.c_str(), "%s\t%s\t%s\t%lu\t%lu", chrom, group, islandname, &start, &end);

			string chromstr = chrom;
			if ((chrom[0] >= '1' && chrom[0] <= '9') || chrom[0] == 'X' || chrom[0] == 'Y')
			{
				string tempchrom = "chr";
				tempchrom.append(chromstr);
				chromstr = tempchrom;
			}

			if (islandmap[chromstr].empty() || islandmap[chromstr].back().second < start)
				islandmap[chromstr].push_back(make_pair(start, end));
			else if(islandmap[chromstr].back().first > start)
			{
				cout << "not ordered!"<<endl;
				exit(0);
			}
			else if(islandmap[chromstr].back().second < end)
				islandmap[chromstr].back().second = end;
			else
			{
				cout << islandmap[chromstr].size() << endl;
				cout << "back " << islandmap[chromstr].back().first << '\t' << islandmap[chromstr].back().second<<endl;
				cout <<"now " <<start  << '\t'<<end<<endl;
				cout << line<<endl;

				getchar();
			}
		}

		ifs.close();
	}

	map<string, vector<pair<size_t, size_t> > >::iterator m_iter;

	for (m_iter = islandmap.begin(); m_iter != islandmap.end(); ++m_iter)
	{
		sort(m_iter->second.begin(), m_iter->second.end(), paircomp);

		cout << m_iter->first << '\t' << m_iter->second.size() <<endl;
	}

	return 0;
}



bool
isinisland(vector<pair<size_t, size_t> > &islandmap, pair<size_t, size_t>& junction)
{
	vector<pair<size_t, size_t> >::iterator lowerbound = lower_bound(islandmap.begin(), islandmap.end(), junction, paircomp);

	if (lowerbound != islandmap.begin())
		--lowerbound;

	//cout << "lowerbound " << lowerbound->first << '\t' <<lowerbound->second<<endl;

	vector<pair<size_t, size_t> >::iterator vpit;

	//cout << "junction "<<junction.first<<'\t'<<junction.second<<endl;

	for (vpit = lowerbound; vpit != islandmap.end(); ++vpit)
	{
		//cout <<"current "<< vpit->first <<'\t'<<vpit->second<<endl;

		if (junction.first >= vpit->first && junction.first <= vpit->second)
			break;
		else if (junction.second < vpit->first)
		{
			//cout << "junction.second < vpit->first"<<endl;
			vpit = islandmap.end();
			break;
		}
	}

	if (vpit == islandmap.end())
	{
		//cout << "first not in"<<endl;
		return false;
	}

	for (vpit = lowerbound; vpit != islandmap.end(); ++vpit)
	{
		//cout <<"current "<< vpit->first <<'\t'<<vpit->second<<endl;
		if (junction.second >= vpit->first && junction.second <= vpit->second)
			break;
		else if (junction.second < vpit->first)
		{
			//cout << "junction.second < vpit->first"<<endl;
			vpit = islandmap.end();
			break;
		}
	}

	if (vpit == islandmap.end())
	{
		//cout << "second not in"<<endl;
		return false;
	}

	return true;
}

size_t
SamHandler::FilterByIsland(const char* islandfile, const char* sam_file)
{
	map<string, vector<pair<size_t, size_t> > > islandmap;
	readisland(islandfile, islandmap);

	string sam_file_inisland = sam_file; sam_file_inisland.append(".inisland"); ofstream ofs_inisland(sam_file_inisland.c_str());

	string sam_file_notinisland = sam_file; sam_file_notinisland.append(".notinisland"); ofstream ofs_notinisland(sam_file_notinisland.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin();

	size_t count = 0;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			vector<pair<size_t, int> >::iterator vp_iter;

			bool anynotinisland = false;

			for (vp_iter = samrec_iter->spliceway_vec.begin(); vp_iter != samrec_iter->spliceway_vec.end() - 1; ++vp_iter)
			{
				pair<size_t, size_t> junction = make_pair(vp_iter->first + vp_iter->second - 1, (vp_iter + 1)->first);

				if (islandmap.find(samrec_iter->chrom_name) == islandmap.end())
				{
					cout <<"warning island not found "<<endl << samrec_iter->tostring()<<endl;
				}
				else
				{
					bool inisland = isinisland(islandmap[samrec_iter->chrom_name], junction);

					if (!inisland)
					{
						//cout << samrec_iter->chrom_name<<'\t' <<junction.first <<'\t' << junction.second<<endl;

						//getchar();
						anynotinisland = true;
						break;
					}
				}
			}

			if (anynotinisland)
			{
				ofs_notinisland << samrec_iter->tostring()<<endl;
				++count;
			}
			else
				ofs_inisland << samrec_iter->tostring()<<endl;
		}

		++lhs_tag_iter;
	}

	return count;
}

size_t
SamHandler::CompareToSyntheticMappedReads(SamHandler& synthetic_sam, const char* output_file, const char* same_offest_output, size_t anchor_width, size_t range)
{
	string output_root = output_file;

	string synthe_vs_cur_1v1 = output_root + ".synthe_vs_cur_1v1";
	ofstream synthe_vs_cur_1v1_ofs(synthe_vs_cur_1v1.c_str());

	string synthe_vs_cur_1vN = output_root + ".synthe_vs_cur_1vN.notbestmatch";
	ofstream synthe_vs_cur_1vN_ofs(synthe_vs_cur_1vN.c_str());

	string synthe_vs_cur_1vN_bestmatch = output_root + ".synthe_vs_cur_1vN.bestmatch";
	ofstream synthe_vs_cur_1vN_bestmatch_ofs(synthe_vs_cur_1vN_bestmatch.c_str());

	string synthe_only = output_root + ".synthe_only";
	ofstream synthe_only_ofs(synthe_only.c_str());

	string cur_only = output_root + ".cur_only";
	ofstream cur_only_ofs(cur_only.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin(), rhs_tag_iter = synthetic_sam.sam_rec_map.begin();

	while (lhs_tag_iter != sam_rec_map.end() &&  rhs_tag_iter != synthetic_sam.sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		if (lhs_tag_iter->first < rhs_tag_iter->first)
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				cur_only_ofs << samrec_iter->tostring()<<endl;
			}

			if (lhs_tag_iter->second.front().isunique)
				++unique_notshared;
			else
				++multiple_notshared;

			if (lhs_tag_iter->second.front().isexonic)
				++exonic_notshared;
			
			if (lhs_tag_iter->second.front().isspliced)
				++spliced_notshared;

			++lhs_tag_iter;
		}
		else if (lhs_tag_iter->first > rhs_tag_iter->first)
		{
			for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
			{
				synthe_only_ofs << samrec_iter->tostring()<<endl;
			}

			if (rhs_tag_iter->second.front().isunique)
				++synthetic_sam.unique_notshared;
			else
				++synthetic_sam.multiple_notshared;

			if (rhs_tag_iter->second.front().isexonic)
				++synthetic_sam.exonic_notshared;
			
			if (rhs_tag_iter->second.front().isspliced)
				++synthetic_sam.spliced_notshared;

			++rhs_tag_iter;
		}
		else
		{
			vector<SamRec>::iterator lhs_samrec_iter = lhs_tag_iter->second.begin(), rhs_samrec_iter = rhs_tag_iter->second.begin();

			if (rhs_tag_iter->second.size() != 1)
				cout <<"warning synthetic reads mapped more than 1 time"<<endl;

			bool matched = false;

			for (; lhs_samrec_iter != lhs_tag_iter->second.end(); ++lhs_samrec_iter)
			{
				if (lhs_samrec_iter->spliceway_vec.size() == rhs_samrec_iter->spliceway_vec.size() && lhs_samrec_iter->chrom_name == rhs_samrec_iter->chrom_name)
				{
					vector<pair<size_t, int> >::iterator lhs_pair_iter = lhs_samrec_iter->spliceway_vec.begin(),  rhs_pair_iter = rhs_samrec_iter->spliceway_vec.begin();

					matched = true;

					while (lhs_pair_iter != lhs_samrec_iter->spliceway_vec.end() && rhs_pair_iter != rhs_samrec_iter->spliceway_vec.end())
					{
						if (lhs_pair_iter->first + range >=  rhs_pair_iter->first && lhs_pair_iter->first <=  rhs_pair_iter->first + range &&
							lhs_pair_iter->second + range >=  rhs_pair_iter->second && lhs_pair_iter->second <=  rhs_pair_iter->second + range)
						{
							//matched = f;
							//break;
						}
						else
						{
							matched = false;
						}

						++lhs_pair_iter;
						++rhs_pair_iter;
					}

					if (matched)
					{
						if (lhs_tag_iter->second.size() == 1)
						{
							synthe_vs_cur_1v1_ofs << lhs_tag_iter->second.front().tostring() <<endl;
							synthe_vs_cur_1v1_ofs << rhs_samrec_iter->tostring() <<endl<<endl;

							if (lhs_tag_iter->second.front().isunique)
								++unique_shared;
							else
								++multiple_shared;

							if (lhs_tag_iter->second.front().isexonic)
								++exonic_shared;

							if (lhs_tag_iter->second.front().isspliced)
								++spliced_shared;

							if (rhs_tag_iter->second.front().isunique)
								++synthetic_sam.unique_shared;
							else
								++synthetic_sam.multiple_shared;

							if (rhs_tag_iter->second.front().isexonic)
								++synthetic_sam.exonic_shared;

							if (rhs_tag_iter->second.front().isspliced)
								++synthetic_sam.spliced_shared;

							if (rhs_tag_iter->second.front().isspliced != lhs_tag_iter->second.front().isspliced)
							{
								cout << "spliced unmatched"<<endl;
								getchar();
							}

							if (synthetic_sam.spliced_shared != spliced_shared)
							{
								cout << "spliced number unmatched"<<endl;
								cout << synthetic_sam.spliced_shared <<":"<<spliced_shared<<endl;
								getchar();
							}
						}
						else
						{
							if ((lhs_samrec_iter)->best == 0/* && rhs_samrec_iter->spliceway_vec.size() > 1*/)
							{
								for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
								{
									synthe_vs_cur_1vN_ofs << samrec_iter->tostring()<<'\t'<<samrec_iter->best <<'\t'
										<< samrec_iter->filter_score<<'\t';

									vector<string>::iterator juncname_iter;

									for (juncname_iter = samrec_iter->junc_id.begin(); juncname_iter != samrec_iter->junc_id.end(); ++juncname_iter)
									{
										synthe_vs_cur_1vN_ofs << *juncname_iter<<'\t';
									}

									synthe_vs_cur_1vN_ofs<<endl;
								}
								synthe_vs_cur_1vN_ofs << rhs_samrec_iter->tostring() <<'\t'<<(lhs_samrec_iter)->best<< endl<<endl;
							}
							else
							{
								for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
								{
									synthe_vs_cur_1vN_bestmatch_ofs << samrec_iter->tostring()<<'\t'<<samrec_iter->best <<'\t'
										<< samrec_iter->filter_score<<'\t';

									vector<string>::iterator juncname_iter;

									for (juncname_iter = samrec_iter->junc_id.begin(); juncname_iter != samrec_iter->junc_id.end(); ++juncname_iter)
									{
										synthe_vs_cur_1vN_bestmatch_ofs << *juncname_iter<<'\t';
									}

									synthe_vs_cur_1vN_bestmatch_ofs<<endl;
								}
								synthe_vs_cur_1vN_bestmatch_ofs << rhs_samrec_iter->tostring() <<'\t'<<(lhs_samrec_iter)->best<< endl<<endl;
							}

														
							if (lhs_tag_iter->second.front().isunique)
								++unique_shared;
							else
								++multiple_shared;

							if (lhs_tag_iter->second.front().isexonic)
								++exonic_shared;

							if (lhs_tag_iter->second.front().isspliced)
								++spliced_shared;

							if (rhs_tag_iter->second.front().isunique)
								++synthetic_sam.unique_shared;
							else
								++synthetic_sam.multiple_shared;

							if (rhs_tag_iter->second.front().isexonic)
								++synthetic_sam.exonic_shared;

							if (rhs_tag_iter->second.front().isspliced)
								++synthetic_sam.spliced_shared;

							if (rhs_tag_iter->second.front().isspliced != lhs_tag_iter->second.front().isspliced)
							{
								cout << "spliced unmatched"<<endl;
								getchar();
							}

							if (synthetic_sam.spliced_shared != spliced_shared)
							{
								cout << "spliced number unmatched"<<endl;
								cout << synthetic_sam.spliced_shared <<":"<<spliced_shared<<endl;
								getchar();
							}
						}

						break;
					}
				}				
			}

			if (!matched)
			{
				for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
				{
					cur_only_ofs << samrec_iter->tostring()<<endl;
				}

				for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
				{
					synthe_only_ofs << samrec_iter->tostring()<<endl;
				}

				if (lhs_tag_iter->second.front().isunique)
					++unique_notshared;
				else
					++multiple_notshared;

				if (lhs_tag_iter->second.front().isexonic)
					++exonic_notshared;

				if (lhs_tag_iter->second.front().isspliced)
					++spliced_notshared;

				if (rhs_tag_iter->second.front().isunique)
					++synthetic_sam.unique_notshared;
				else
					++synthetic_sam.multiple_notshared;

				if (rhs_tag_iter->second.front().isexonic)
					++synthetic_sam.exonic_notshared;

				if (rhs_tag_iter->second.front().isspliced)
					++synthetic_sam.spliced_notshared;
				
			}

			++lhs_tag_iter;
			++rhs_tag_iter;
		}
	}

	vector<SamRec>::iterator samrec_iter;

	while (lhs_tag_iter != sam_rec_map.end())
	{
		for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
		{
			cur_only_ofs << samrec_iter->tostring()<<endl;
		}
		if (lhs_tag_iter->second.front().isunique)
			++unique_notshared;
		else
			++multiple_notshared;

		if (lhs_tag_iter->second.front().isexonic)
			++exonic_notshared;

		if (lhs_tag_iter->second.front().isspliced)
			++spliced_notshared;
		++lhs_tag_iter;
	}

	while (rhs_tag_iter != synthetic_sam.sam_rec_map.end())
	{
		for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
		{
			synthe_only_ofs << samrec_iter->tostring()<<endl;
		}

		if (rhs_tag_iter->second.front().isunique)
			++synthetic_sam.unique_notshared;
		else
			++synthetic_sam.multiple_notshared;

		if (rhs_tag_iter->second.front().isexonic)
			++synthetic_sam.exonic_notshared;

		if (rhs_tag_iter->second.front().isspliced)
			++synthetic_sam.spliced_notshared;
		++rhs_tag_iter;
	}

	CompareTwoSamByOffset(cur_only.c_str(), synthe_only.c_str(), same_offest_output, anchor_width, range);

	return 0;
}

size_t
SamHandler::Compare2Sam(SamHandler& compare_sam, string output_dir)
{
	string output_root = output_dir;

	string mps_vs_comp_1v1_matched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.matched"; ofstream mps_vs_comp_1v1_matched_ofs(mps_vs_comp_1v1_matched.c_str());

	string mps_vs_comp_1v1_notmatched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.notmatched"; ofstream mps_vs_comp_1v1_notmatched_ofs(mps_vs_comp_1v1_notmatched.c_str());

	string mps_vs_comp_1v1_notmatched_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.notmatched.better"; ofstream mps_vs_comp_1v1_notmatched_better_ofs(mps_vs_comp_1v1_notmatched_better.c_str());

	string mps_vs_comp_1v1_notmatched_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.notmatched.equal"; ofstream mps_vs_comp_1v1_notmatched_equal_ofs(mps_vs_comp_1v1_notmatched_equal.c_str());

	string mps_vs_comp_1v1_notmatched_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.notmatched.worse"; ofstream mps_vs_comp_1v1_notmatched_worse_ofs(mps_vs_comp_1v1_notmatched_worse.c_str());

	string mps_vs_comp_1v1_samestartend = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.samestartend"; ofstream mps_vs_comp_1v1_samestartend_ofs(mps_vs_comp_1v1_samestartend.c_str());

	string mps_vs_comp_1v1_samestartend_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.samestartend.better"; ofstream mps_vs_comp_1v1_samestartend_better_ofs(mps_vs_comp_1v1_samestartend_better.c_str());

	string mps_vs_comp_1v1_samestartend_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.samestartend.equal"; ofstream mps_vs_comp_1v1_samestartend_equal_ofs(mps_vs_comp_1v1_samestartend_equal.c_str());

	string mps_vs_comp_1v1_samestartend_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1v1.samestartend.worse"; ofstream mps_vs_comp_1v1_samestartend_worse_ofs(mps_vs_comp_1v1_samestartend_worse.c_str());

	string mps_vs_comp_1vN_matched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.matched"; ofstream mps_vs_comp_1vN_matched_ofs(mps_vs_comp_1vN_matched.c_str());

	string mps_vs_comp_1vN_notmatched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.notmatched"; ofstream mps_vs_comp_1vN_notmatched_ofs(mps_vs_comp_1vN_notmatched.c_str());

	string mps_vs_comp_1vN_notmatched_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.notmatched.better"; ofstream mps_vs_comp_1vN_notmatched_better_ofs(mps_vs_comp_1vN_notmatched_better.c_str());

	string mps_vs_comp_1vN_notmatched_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.notmatched.equal"; ofstream mps_vs_comp_1vN_notmatched_equal_ofs(mps_vs_comp_1vN_notmatched_equal.c_str());

	string mps_vs_comp_1vN_notmatched_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.notmatched.worse"; ofstream mps_vs_comp_1vN_notmatched_worse_ofs(mps_vs_comp_1vN_notmatched_worse.c_str());

	string mps_vs_comp_1vN_samestartend = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.samestartend"; ofstream mps_vs_comp_1vN_samestartend_ofs(mps_vs_comp_1vN_samestartend.c_str());

	string mps_vs_comp_1vN_samestartend_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.samestartend.better"; ofstream mps_vs_comp_1vN_samestartend_better_ofs(mps_vs_comp_1vN_samestartend_better.c_str());

	string mps_vs_comp_1vN_samestartend_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.samestartend.equal"; ofstream mps_vs_comp_1vN_samestartend_equal_ofs(mps_vs_comp_1vN_samestartend_equal.c_str());

	string mps_vs_comp_1vN_samestartend_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".1vN.samestartend.worse"; ofstream mps_vs_comp_1vN_samestartend_worse_ofs(mps_vs_comp_1vN_samestartend_worse.c_str());

	string mps_vs_comp_Nv1_matched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.matched"; ofstream mps_vs_comp_Nv1_matched_ofs(mps_vs_comp_Nv1_matched.c_str());

	string mps_vs_comp_Nv1_notmatched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.notmatched"; ofstream mps_vs_comp_Nv1_notmatched_ofs(mps_vs_comp_Nv1_notmatched.c_str());

	string mps_vs_comp_Nv1_notmatched_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.notmatched.better"; ofstream mps_vs_comp_Nv1_notmatched_better_ofs(mps_vs_comp_Nv1_notmatched_better.c_str());

	string mps_vs_comp_Nv1_notmatched_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.notmatched.equal"; ofstream mps_vs_comp_Nv1_notmatched_equal_ofs(mps_vs_comp_Nv1_notmatched_equal.c_str());

	string mps_vs_comp_Nv1_notmatched_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.notmatched.worse"; ofstream mps_vs_comp_Nv1_notmatched_worse_ofs(mps_vs_comp_Nv1_notmatched_worse.c_str());

	string mps_vs_comp_Nv1_samestartend = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.samestartend"; ofstream mps_vs_comp_Nv1_samestartend_ofs(mps_vs_comp_Nv1_samestartend.c_str());

	string mps_vs_comp_Nv1_samestartend_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.samestartend.better"; ofstream mps_vs_comp_Nv1_samestartend_better_ofs(mps_vs_comp_Nv1_samestartend_better.c_str());

	string mps_vs_comp_Nv1_samestartend_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.samestartend.equal"; ofstream mps_vs_comp_Nv1_samestartend_equal_ofs(mps_vs_comp_Nv1_samestartend_equal.c_str());

	string mps_vs_comp_Nv1_samestartend_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv1.samestartend.worse"; ofstream mps_vs_comp_Nv1_samestartend_worse_ofs(mps_vs_comp_Nv1_samestartend_worse.c_str());

	string mps_vs_comp_NvM_matched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.matched"; ofstream mps_vs_comp_NvM_matched_ofs(mps_vs_comp_NvM_matched.c_str());

	string mps_vs_comp_NvM_notmatched = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.notmatched"; ofstream mps_vs_comp_NvM_notmatched_ofs(mps_vs_comp_NvM_notmatched.c_str());

	string mps_vs_comp_NvM_notmatched_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.notmatched.better"; ofstream mps_vs_comp_NvM_notmatched_better_ofs(mps_vs_comp_NvM_notmatched_better.c_str());

	string mps_vs_comp_NvM_notmatched_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.notmatched.equal"; ofstream mps_vs_comp_NvM_notmatched_equal_ofs(mps_vs_comp_NvM_notmatched_equal.c_str());

	string mps_vs_comp_NvM_notmatched_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.notmatched.worse"; ofstream mps_vs_comp_NvM_notmatched_worse_ofs(mps_vs_comp_NvM_notmatched_worse.c_str());

	string mps_vs_comp_NvM_samestartend = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.samestartend"; ofstream mps_vs_comp_NvM_samestartend_ofs(mps_vs_comp_NvM_samestartend.c_str());

	string mps_vs_comp_NvM_samestartend_better = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.samestartend.better"; ofstream mps_vs_comp_NvM_samestartend_better_ofs(mps_vs_comp_NvM_samestartend_better.c_str());

	string mps_vs_comp_NvM_samestartend_equal = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.samestartend.equal"; ofstream mps_vs_comp_NvM_samestartend_equal_ofs(mps_vs_comp_NvM_samestartend_equal.c_str());

	string mps_vs_comp_NvM_samestartend_worse = output_root + base_file_name + "_" + compare_sam.base_file_name + ".NvM.samestartend.worse"; ofstream mps_vs_comp_NvM_samestartend_worse_ofs(mps_vs_comp_NvM_samestartend_worse.c_str());

	size_t one_v_one_notmatched_better_mis = 0, one_v_one_notmatched_better_qual = 0, one_v_one_notmatched_better_intron = 0, one_v_one_notmatched_equal = 0;
	
	size_t one_v_one_notmatched_worse_mis = 0, one_v_one_notmatched_worse_qual = 0, one_v_one_notmatched_worse_intron = 0;

	size_t one_v_one_samestartend_better_mis = 0, one_v_one_samestartend_better_qual = 0, one_v_one_samestartend_better_intron = 0, one_v_one_samestartend_equal = 0;
	
	size_t one_v_one_samestartend_worse_mis = 0, one_v_one_samestartend_worse_intron = 0, one_v_one_samestartend_worse_qual = 0;

	size_t one_v_N_notmatched_better_mis = 0, one_v_N_notmatched_better_qual = 0, one_v_N_notmatched_better_intron = 0, one_v_N_notmatched_equal = 0;
	
	size_t one_v_N_notmatched_worse_mis = 0, one_v_N_notmatched_worse_qual = 0, one_v_N_notmatched_worse_intron = 0;

	size_t one_v_N_samestartend_better_mis = 0, one_v_N_samestartend_better_qual = 0, one_v_N_samestartend_better_intron = 0, one_v_N_samestartend_equal = 0;
	
	size_t one_v_N_samestartend_worse_mis = 0, one_v_N_samestartend_worse_qual = 0, one_v_N_samestartend_worse_intron = 0;

	size_t N_v_one_notmatched_better_mis = 0, N_v_one_notmatched_better_qual = 0, N_v_one_notmatched_better_intron = 0, N_v_one_notmatched_equal = 0;
	
	size_t N_v_one_notmatched_worse_mis = 0, N_v_one_notmatched_worse_qual = 0, N_v_one_notmatched_worse_intron = 0;

	size_t N_v_one_samestartend_better_mis = 0, N_v_one_samestartend_better_qual = 0, N_v_one_samestartend_better_intron = 0, N_v_one_samestartend_equal = 0;
	
	size_t N_v_one_samestartend_worse_mis = 0, N_v_one_samestartend_worse_qual = 0, N_v_one_samestartend_worse_intron = 0;

	size_t N_v_M_notmatched_better_mis = 0, N_v_M_notmatched_better_qual = 0, N_v_M_notmatched_better_intron = 0, N_v_M_notmatched_equal = 0;
	
	size_t N_v_M_notmatched_worse_mis = 0, N_v_M_notmatched_worse_qual = 0, N_v_M_notmatched_worse_intron = 0;

	size_t N_v_M_samestartend_better_mis = 0, N_v_M_samestartend_better_qual = 0, N_v_M_samestartend_better_intron = 0, N_v_M_samestartend_equal = 0;
	
	size_t N_v_M_samestartend_worse_mis = 0, N_v_M_samestartend_worse_qual = 0, N_v_M_samestartend_worse_intron = 0;

	//string mps_vs_comp_Nv1_bestmatched = output_root + ".mps_vs_comp_Nv1.bestmatched"; ofstream mps_vs_comp_Nv1_bestmatched_ofs(mps_vs_comp_Nv1_bestmatched.c_str());

	//string mps_vs_comp_Nv1_bestnotmatched = output_root + ".mps_vs_comp_Nv1.bestnotmatched"; ofstream mps_vs_comp_Nv1_bestnotmatched_ofs(mps_vs_comp_Nv1_bestnotmatched.c_str());

	//string mps_vs_comp_NvM_bestmatched = output_root + ".mps_vs_comp_NvM.bestmatched"; ofstream mps_vs_comp_NvM_bestmatched_ofs(mps_vs_comp_NvM_bestmatched.c_str());

	//string mps_vs_comp_NvM_bestnotmatched = output_root + ".mps_vs_comp_NvM.bestnotmatched"; ofstream mps_vs_comp_NvM_bestnotmatched_ofs(mps_vs_comp_NvM_bestnotmatched.c_str());

	string mps_vs_comp_0vN = output_root + base_file_name + "_" + compare_sam.base_file_name + ".0vN"; ofstream mps_vs_comp_0vN_ofs(mps_vs_comp_0vN.c_str());

	string mps_vs_comp_Nv0 = output_root + base_file_name + "_" + compare_sam.base_file_name + ".Nv0"; ofstream mps_vs_comp_Nv0_ofs(mps_vs_comp_Nv0.c_str());

	if (sam_rec_map.empty())
		ReadNextTag();

	if (compare_sam.sam_rec_map.empty())
		compare_sam.ReadNextTag();

	while (true)
	{
		if (sam_rec_map.size() > 1 || compare_sam.sam_rec_map.size() > 1)
			cout << "one of sam rec map too large"<<endl;

		if (cur_ifs.eof() && compare_sam.cur_ifs.eof() && sam_rec_map.empty() && compare_sam.sam_rec_map.empty())
			break;

		if (sam_rec_map.empty() && compare_sam.sam_rec_map.empty())
		{
			cout << "both are empty but not finished"<<endl;

			continue;
		}		

		vector<SamRec>::iterator samrec_iter;

		if ((!sam_rec_map.empty() && compare_sam.sam_rec_map.empty()) || 
			(!sam_rec_map.empty() && !compare_sam.sam_rec_map.empty() && strcmp(sam_rec_map.begin()->first.c_str(), compare_sam.sam_rec_map.begin()->first.c_str()) < 0 
			/*comp_str_less_than_linux(sam_rec_map.begin()->first, compare_sam.sam_rec_map.begin()->first)*/ /*sam_rec_map.begin()->first < compare_sam.sam_rec_map.begin()->first*/))
		{
			bool isexonic = false;

			bool isspliced = false;

			for (samrec_iter = sam_rec_map.begin()->second.begin(); samrec_iter != sam_rec_map.begin()->second.end(); ++samrec_iter)
			{
				mps_vs_comp_Nv0_ofs << samrec_iter->cur_line<<endl;

				if (samrec_iter->isexonic)
					isexonic = true;

				if (samrec_iter->isspliced)
					isspliced = true;
			}

			if (sam_rec_map.begin()->second.size() == 1)
				++unique_only;
			else
				++multiple_only;

			if (isexonic)
				++exonic_only;
			
			if (isspliced)
				++spliced_only;

			ReadNextTag();

		}
		else if ((sam_rec_map.empty() && !compare_sam.sam_rec_map.empty()) || 
				 (!sam_rec_map.empty() && !compare_sam.sam_rec_map.empty() && strcmp(sam_rec_map.begin()->first.c_str(), compare_sam.sam_rec_map.begin()->first.c_str()) > 0
				 /*comp_str_less_than_linux(compare_sam.sam_rec_map.begin()->first, sam_rec_map.begin()->first)*/ /*sam_rec_map.begin()->first > compare_sam.sam_rec_map.begin()->first*/))
		{
			bool isexonic = false;

			bool isspliced = false;

			for (samrec_iter = compare_sam.sam_rec_map.begin()->second.begin(); samrec_iter != compare_sam.sam_rec_map.begin()->second.end(); ++samrec_iter)
			{
				mps_vs_comp_0vN_ofs << samrec_iter->cur_line<<endl;

				if (samrec_iter->isexonic)
					isexonic = true;

				if (samrec_iter->isspliced)
					isspliced = true;
			}

			if (compare_sam.sam_rec_map.begin()->second.size() == 1)
				++compare_sam.unique_only;
			else
				++compare_sam.multiple_only;

			if (isexonic)
				++compare_sam.exonic_only;
			
			if (isspliced)
				++compare_sam.spliced_only;

			compare_sam.ReadNextTag();
		}
		else
		{
			map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin(), rhs_tag_iter = compare_sam.sam_rec_map.begin();

			sort(lhs_tag_iter->second.begin(), lhs_tag_iter->second.end(), comp_chromname_offset_spliceway);

			sort(rhs_tag_iter->second.begin(), rhs_tag_iter->second.end(), comp_chromname_offset_spliceway);

			vector<SamRec>::iterator lhs_samrec_iter = lhs_tag_iter->second.begin(), rhs_samrec_iter = rhs_tag_iter->second.begin();

			bool matched = false;

			//bool bestmatched = false;

			int match_id = 1;

			bool same_start_or_end = false;

			bool isexonic = false;

			bool isspliced = false;

			while (lhs_samrec_iter != sam_rec_map.begin()->second.end() && rhs_samrec_iter != compare_sam.sam_rec_map.begin()->second.end())
			{
				bool l_comp_r = comp_chromname_offset_spliceway(*lhs_samrec_iter, *rhs_samrec_iter);
				bool r_comp_l = comp_chromname_offset_spliceway(*rhs_samrec_iter, *lhs_samrec_iter);

				if (lhs_samrec_iter->start == rhs_samrec_iter->start || lhs_samrec_iter->end == rhs_samrec_iter->end)
				{
					same_start_or_end = true;
				}

				if (l_comp_r == true && r_comp_l == false)
				{
					++lhs_samrec_iter;
				}
				else if (l_comp_r == false && r_comp_l == true)
				{
					++rhs_samrec_iter;
				}
				else
				{
					lhs_samrec_iter->matched_id = match_id;

					rhs_samrec_iter->matched_id = match_id;

					if (lhs_samrec_iter->isexonic)
						isexonic = true;

					if (lhs_samrec_iter->isspliced)
						isspliced = true;

					//if (!bestmatched && lhs_samrec_iter->best)
					//	bestmatched = true;

					++match_id;

					++lhs_samrec_iter;
					++rhs_samrec_iter;
					matched = true;
				}
			}

			if (matched)
			{
				ofstream* matched_ofs;

				if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_1v1_matched_ofs;
				}
				else if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() > 1)
				{
					matched_ofs = &mps_vs_comp_1vN_matched_ofs;
				}
				else if (lhs_tag_iter->second.size() > 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_Nv1_matched_ofs;
					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_Nv1_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_Nv1_bestnotmatched_ofs;
				}
				else
				{
					matched_ofs = &mps_vs_comp_NvM_matched_ofs;
					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_NvM_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_NvM_bestnotmatched_ofs;
				}

				//*matched_ofs <<base_file_name<<endl;
				for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
				{
					*matched_ofs<<base_file_name<<'~'<<lhs_samrec_iter->cur_line << endl;

					//vector<string>::iterator juncname_iter;

					//for (juncname_iter = lhs_samrec_iter->junc_id.begin(); juncname_iter != lhs_samrec_iter->junc_id.end(); ++juncname_iter)
					//{
					//	*matched_ofs << *juncname_iter<<'\t';
					//}
					//*matched_ofs <<endl;
				}

				//*matched_ofs<<compare_sam.base_file_name<<endl;

				for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
				{
					*matched_ofs<<compare_sam.base_file_name<<'~'<<rhs_samrec_iter->cur_line<<endl;
				}

				//*matched_ofs<<endl;

				if (lhs_tag_iter->second.size() == 1)
					++unique_shared;
				else
					++multiple_shared;

				if (isexonic)
					++exonic_shared;

				if (isspliced)
					++spliced_shared;

				if (rhs_tag_iter->second.size() == 1)
					++compare_sam.unique_shared;
				else
					++compare_sam.multiple_shared;

				if (isexonic)
					++compare_sam.exonic_shared;

				if (isspliced)
					++compare_sam.spliced_shared;

			}
			else if (same_start_or_end)
			{
				ofstream* matched_ofs;

				ofstream* whichbetter_ofs;

				size_t min_lhs_mis = 255, min_lhs_qual = 255, min_rhs_mis = 255, min_rhs_qual = 255, min_lhs_intron = 10000000000, min_rhs_intron = 100000000000;

				vector<SamRec>::iterator min_lhs_sam_rec_iter, min_rhs_sam_rec_iter;

				for (min_lhs_sam_rec_iter = lhs_tag_iter->second.begin(); min_lhs_sam_rec_iter != lhs_tag_iter->second.end(); ++min_lhs_sam_rec_iter)
				{
					if (min_lhs_mis > min_lhs_sam_rec_iter->mis_match)
						min_lhs_mis = min_lhs_sam_rec_iter->mis_match;

					if (min_lhs_qual > min_lhs_sam_rec_iter->confid)
						min_lhs_qual = min_lhs_sam_rec_iter->confid;

					if (min_lhs_sam_rec_iter->spliceway_vec.size() == 1)
						min_lhs_intron = 0;
					else
					{
						size_t sum_intron = 0;
						for (size_t k = 1; k < min_lhs_sam_rec_iter->spliceway_vec.size(); ++k)
						{
							sum_intron += min_lhs_sam_rec_iter->spliceway_vec[k].first - min_lhs_sam_rec_iter->spliceway_vec[k - 1].first - min_lhs_sam_rec_iter->spliceway_vec[k - 1].second;
						}

						if (min_lhs_intron > sum_intron)
							min_lhs_intron = sum_intron;
					}
				}

				for (min_rhs_sam_rec_iter = rhs_tag_iter->second.begin(); min_rhs_sam_rec_iter != rhs_tag_iter->second.end(); ++min_rhs_sam_rec_iter)
				{
					if (min_rhs_mis > min_rhs_sam_rec_iter->mis_match)
						min_rhs_mis = min_rhs_sam_rec_iter->mis_match;

					if (min_rhs_qual > min_rhs_sam_rec_iter->confid)
						min_rhs_qual = min_rhs_sam_rec_iter->confid;

					if (min_rhs_sam_rec_iter->spliceway_vec.size() == 1)
						min_rhs_intron = 0;
					else
					{
						size_t sum_intron = 0;
						for (size_t k = 1; k < min_rhs_sam_rec_iter->spliceway_vec.size(); ++k)
						{
							sum_intron += min_rhs_sam_rec_iter->spliceway_vec[k].first - min_rhs_sam_rec_iter->spliceway_vec[k - 1].first - min_rhs_sam_rec_iter->spliceway_vec[k - 1].second;
						}

						if (min_rhs_intron > sum_intron)
							min_rhs_intron = sum_intron;
					}
				}

				min_lhs_qual = 255, min_rhs_qual = 255;

				int whichbetter = 0;

				string reason = "";

				if (min_lhs_mis < min_rhs_mis)
				{
					whichbetter = -1;
					reason = "\tXA:Z:bettermis";
				}
				else if (min_lhs_mis > min_rhs_mis)
				{
					whichbetter = 1;
					reason = "\tXA:Z:worsemis";
				}
				else
				{
					if (min_lhs_qual < min_rhs_qual)
					{
						whichbetter = -2;
						reason = "\tXA:Z:betterqual";
					}
					else if (min_lhs_qual > min_rhs_qual)
					{
						whichbetter = 2;
						reason = "\tXA:Z:worsequal";
					}
					else
					{
						if (min_lhs_intron < min_rhs_intron)
						{
							whichbetter = -3;
							reason = "\tXA:Z:betterintron";
						}
						else if (min_lhs_intron > min_rhs_intron)
						{
							whichbetter = 3;
							reason = "\tXA:Z:worseintron";
						}
						else
						{
							whichbetter = 0;
						}
					}
				}
				

				//size_t N_v_M_notmatched_better = 0, N_v_M_notmatched_equal = 0, N_v_M_notmatched_worse = 0;

				if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_1v1_samestartend_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_1v1_samestartend_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_1v1_samestartend_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_1v1_samestartend_equal_ofs;

					if (whichbetter == -3)
						++one_v_one_samestartend_better_intron;
					else if (whichbetter == -2)
						++one_v_one_samestartend_better_qual;
					else if (whichbetter == -1)
						++one_v_one_samestartend_better_mis;
					else if (whichbetter == 0)
						++one_v_one_samestartend_equal;
					else if (whichbetter == 1)
						++one_v_one_samestartend_worse_mis;
					else if (whichbetter == 2)
						++one_v_one_samestartend_worse_qual;
					else if (whichbetter == 3)
						++one_v_one_samestartend_worse_intron;
				}
				else if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() > 1)
				{
					matched_ofs = &mps_vs_comp_1vN_samestartend_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_1vN_samestartend_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_1vN_samestartend_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_1vN_samestartend_equal_ofs;

					if (whichbetter == -3)
						++one_v_N_samestartend_better_intron;
					else if (whichbetter == -2)
						++one_v_N_samestartend_better_qual;
					else if (whichbetter == -1)
						++one_v_N_samestartend_better_mis;
					else if (whichbetter == 0)
						++one_v_N_samestartend_equal;
					else if (whichbetter == 1)
						++one_v_N_samestartend_worse_mis;
					else if (whichbetter == 2)
						++one_v_N_samestartend_worse_qual;
					else if (whichbetter == 3)
						++one_v_N_samestartend_worse_intron;

				}
				else if (lhs_tag_iter->second.size() > 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_Nv1_samestartend_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_Nv1_samestartend_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_Nv1_samestartend_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_Nv1_samestartend_equal_ofs;

					if (whichbetter == -3)
						++N_v_one_samestartend_better_intron;
					else if (whichbetter == -2)
						++N_v_one_samestartend_better_qual;
					else if (whichbetter == -1)
						++N_v_one_samestartend_better_mis;
					else if (whichbetter == 0)
						++N_v_one_samestartend_equal;
					else if (whichbetter == 1)
						++N_v_one_samestartend_worse_mis;
					else if (whichbetter == 2)
						++N_v_one_samestartend_worse_qual;
					else if (whichbetter == 3)
						++N_v_one_samestartend_worse_intron;
					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_Nv1_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_Nv1_bestnotmatched_ofs;
				}
				else
				{
					matched_ofs = &mps_vs_comp_NvM_samestartend_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_NvM_samestartend_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_NvM_samestartend_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_NvM_samestartend_equal_ofs;

					if (whichbetter == -3)
						++N_v_M_samestartend_better_intron;
					else if (whichbetter == -2)
						++N_v_M_samestartend_better_qual;
					else if (whichbetter == -1)
						++N_v_M_samestartend_better_mis;
					else if (whichbetter == 0)
						++N_v_M_samestartend_equal;
					else if (whichbetter == 1)
						++N_v_M_samestartend_worse_mis;
					else if (whichbetter == 2)
						++N_v_M_samestartend_worse_qual;
					else if (whichbetter == 3)
						++N_v_M_samestartend_worse_intron;

					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_NvM_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_NvM_bestnotmatched_ofs;
				}

				bool isexonic = false;

				bool isspliced = false;

				//*matched_ofs <<base_file_name<<endl;
				for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
				{
					*matched_ofs<<base_file_name<<'~'<<lhs_samrec_iter->cur_line <<reason<< endl;

					*whichbetter_ofs<<base_file_name<<'~'<<lhs_samrec_iter->cur_line << reason<< endl;

					if (lhs_samrec_iter->isexonic)
						isexonic = true;

					if (lhs_samrec_iter->isspliced)
						isspliced = true;

					//vector<string>::iterator juncname_iter;

					//for (juncname_iter = lhs_samrec_iter->junc_id.begin(); juncname_iter != lhs_samrec_iter->junc_id.end(); ++juncname_iter)
					//{
					//	*matched_ofs << *juncname_iter<<'\t';
					//}
					//*matched_ofs <<endl;
				}

				if (lhs_tag_iter->second.size() == 1)
					++unique_same_start_or_end;
				else
					++multiple_same_start_or_end;

				if (isexonic)
					++exonic_same_start_or_end;

				if (isspliced)
					++spliced_same_start_or_end;

				isexonic = false;

				isspliced = false;

				//*matched_ofs<<compare_sam.base_file_name<<endl;

				for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
				{
					*matched_ofs<<compare_sam.base_file_name<<'~'<<rhs_samrec_iter->cur_line<<reason<< endl;

					*whichbetter_ofs<<compare_sam.base_file_name<<'~'<<rhs_samrec_iter->cur_line<<reason<< endl;

					if (rhs_samrec_iter->isexonic)
						isexonic = true;

					if (rhs_samrec_iter->isspliced)
						isspliced = true;
				}

				//*matched_ofs<<endl;

				if (rhs_tag_iter->second.size() == 1)
					++compare_sam.unique_same_start_or_end;
				else
					++compare_sam.multiple_same_start_or_end;

				if (isexonic)
					++compare_sam.exonic_same_start_or_end;

				if (isspliced)
					++compare_sam.spliced_same_start_or_end;
			}
			else
			{
				ofstream* matched_ofs;

				ofstream* whichbetter_ofs;

				size_t min_lhs_mis = 255, min_lhs_qual = 255, min_rhs_mis = 255, min_rhs_qual = 255, min_lhs_intron = 10000000000, min_rhs_intron = 100000000000;

				vector<SamRec>::iterator min_lhs_sam_rec_iter, min_rhs_sam_rec_iter;

				for (min_lhs_sam_rec_iter = lhs_tag_iter->second.begin(); min_lhs_sam_rec_iter != lhs_tag_iter->second.end(); ++min_lhs_sam_rec_iter)
				{
					if (min_lhs_mis > min_lhs_sam_rec_iter->mis_match)
						min_lhs_mis = min_lhs_sam_rec_iter->mis_match;

					if (min_lhs_qual > min_lhs_sam_rec_iter->confid)
						min_lhs_qual = min_lhs_sam_rec_iter->confid;

					if (min_lhs_sam_rec_iter->spliceway_vec.size() == 1)
						min_lhs_intron = 0;
					else
					{
						size_t sum_intron = 0;
						for (size_t k = 1; k < min_lhs_sam_rec_iter->spliceway_vec.size(); ++k)
						{
							sum_intron += min_lhs_sam_rec_iter->spliceway_vec[k].first - min_lhs_sam_rec_iter->spliceway_vec[k - 1].first - min_lhs_sam_rec_iter->spliceway_vec[k - 1].second;
						}

						if (min_lhs_intron > sum_intron)
							min_lhs_intron = sum_intron;
					}
				}

				for (min_rhs_sam_rec_iter = rhs_tag_iter->second.begin(); min_rhs_sam_rec_iter != rhs_tag_iter->second.end(); ++min_rhs_sam_rec_iter)
				{
					if (min_rhs_mis > min_rhs_sam_rec_iter->mis_match)
						min_rhs_mis = min_rhs_sam_rec_iter->mis_match;

					if (min_rhs_qual > min_rhs_sam_rec_iter->confid)
						min_rhs_qual = min_rhs_sam_rec_iter->confid;

					if (min_rhs_sam_rec_iter->spliceway_vec.size() == 1)
						min_rhs_intron = 0;
					else
					{
						size_t sum_intron = 0;
						for (size_t k = 1; k < min_rhs_sam_rec_iter->spliceway_vec.size(); ++k)
						{
							sum_intron += min_rhs_sam_rec_iter->spliceway_vec[k].first - min_rhs_sam_rec_iter->spliceway_vec[k - 1].first - min_rhs_sam_rec_iter->spliceway_vec[k - 1].second;
						}

						if (min_rhs_intron > sum_intron)
							min_rhs_intron = sum_intron;
					}
				}

				min_lhs_qual = 255, min_rhs_qual = 255;

				int whichbetter = 0;

				string reason = "";

				if (min_lhs_mis < min_rhs_mis)
				{
					whichbetter = -1;
					reason = "\tXA:Z:bettermis";
				}
				else if (min_lhs_mis > min_rhs_mis)
				{
					whichbetter = 1;
					reason = "\tXA:Z:worsemis";
				}
				else
				{
					if (min_lhs_qual < min_rhs_qual)
					{
						whichbetter = -2;
						reason = "\tXA:Z:betterqual";
					}
					else if (min_lhs_qual > min_rhs_qual)
					{
						whichbetter = 2;
						reason = "\tXA:Z:worsequal";
					}
					else
					{
						if (min_lhs_intron < min_rhs_intron)
						{
							whichbetter = -3;
							reason = "\tXA:Z:betterintron";
						}
						else if (min_lhs_intron > min_rhs_intron)
						{
							whichbetter = 3;
							reason = "\tXA:Z:worseintron";
						}
						else
						{
							whichbetter = 0;
						}
					}
				}

				if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_1v1_notmatched_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_1v1_notmatched_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_1v1_notmatched_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_1v1_notmatched_equal_ofs;

					if (whichbetter == -3)
						++one_v_one_notmatched_better_intron;
					else if (whichbetter == -2)
						++one_v_one_notmatched_better_qual;
					else if (whichbetter == -1)
						++one_v_one_notmatched_better_mis;
					else if (whichbetter == 0)
						++one_v_one_notmatched_equal;
					else if (whichbetter == 1)
						++one_v_one_notmatched_worse_mis;
					else if (whichbetter == 2)
						++one_v_one_notmatched_worse_qual;
					else if (whichbetter == 3)
						++one_v_one_notmatched_worse_intron;
				}
				else if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() > 1)
				{
					matched_ofs = &mps_vs_comp_1vN_notmatched_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_1vN_notmatched_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_1vN_notmatched_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_1vN_notmatched_equal_ofs;

					if (whichbetter == -3)
						++one_v_N_notmatched_better_intron;
					else if (whichbetter == -2)
						++one_v_N_notmatched_better_qual;
					else if (whichbetter == -1)
						++one_v_N_notmatched_better_mis;
					else if (whichbetter == 0)
						++one_v_N_notmatched_equal;
					else if (whichbetter == 1)
						++one_v_N_notmatched_worse_mis;
					else if (whichbetter == 2)
						++one_v_N_notmatched_worse_qual;
					else if (whichbetter == 3)
						++one_v_N_notmatched_worse_intron;
				}
				else if (lhs_tag_iter->second.size() > 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_Nv1_notmatched_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_Nv1_notmatched_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_Nv1_notmatched_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_Nv1_notmatched_equal_ofs;

					if (whichbetter == -3)
						++N_v_one_notmatched_better_intron;
					else if (whichbetter == -2)
						++N_v_one_notmatched_better_qual;
					else if (whichbetter == -1)
						++N_v_one_notmatched_better_mis;
					else if (whichbetter == 0)
						++N_v_one_notmatched_equal;
					else if (whichbetter == 1)
						++N_v_one_notmatched_worse_mis;
					else if (whichbetter == 2)
						++N_v_one_notmatched_worse_qual;
					else if (whichbetter == 3)
						++N_v_one_notmatched_worse_intron;

					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_Nv1_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_Nv1_bestnotmatched_ofs;
				}
				else
				{
					matched_ofs = &mps_vs_comp_NvM_notmatched_ofs;

					if (whichbetter > 0)
						whichbetter_ofs = &mps_vs_comp_NvM_notmatched_worse_ofs;
					else if (whichbetter < 0)
						whichbetter_ofs = &mps_vs_comp_NvM_notmatched_better_ofs;
					else
						whichbetter_ofs = &mps_vs_comp_NvM_notmatched_equal_ofs;

					if (whichbetter == -3)
						++N_v_M_notmatched_better_intron;
					else if (whichbetter == -2)
						++N_v_M_notmatched_better_qual;
					else if (whichbetter == -1)
						++N_v_M_notmatched_better_mis;
					else if (whichbetter == 0)
						++N_v_M_notmatched_equal;
					else if (whichbetter == 1)
						++N_v_M_notmatched_worse_mis;
					else if (whichbetter == 2)
						++N_v_M_notmatched_worse_qual;
					else if (whichbetter == 3)
						++N_v_M_notmatched_worse_intron;
					//if (bestmatched)
					//	matched_ofs = &mps_vs_comp_NvM_bestmatched_ofs;
					//else
					//	matched_ofs = &mps_vs_comp_NvM_bestnotmatched_ofs;
				}

				bool isexonic = false;

				bool isspliced = false;

				//*matched_ofs <<base_file_name<<endl;
				for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
				{
					*matched_ofs<<base_file_name<<'~'<<lhs_samrec_iter->cur_line <<reason<< endl;

					*whichbetter_ofs<<base_file_name<<'~'<<lhs_samrec_iter->cur_line <<reason<< endl;

					if (lhs_samrec_iter->isexonic)
						isexonic = true;

					if (lhs_samrec_iter->isspliced)
						isspliced = true;

					//vector<string>::iterator juncname_iter;

					//for (juncname_iter = lhs_samrec_iter->junc_id.begin(); juncname_iter != lhs_samrec_iter->junc_id.end(); ++juncname_iter)
					//{
					//	*matched_ofs << *juncname_iter<<'\t';
					//}
					//*matched_ofs <<endl;
				}

				if (lhs_tag_iter->second.size() == 1)
					++unique_notshared;
				else
					++multiple_notshared;

				if (isexonic)
					++exonic_notshared;

				if (isspliced)
					++spliced_notshared;

				isexonic = false;

				isspliced = false;

				//*matched_ofs<<compare_sam.base_file_name<<endl;

				for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
				{
					*matched_ofs<<compare_sam.base_file_name<<'~'<<rhs_samrec_iter->cur_line<<reason<<endl;

					*whichbetter_ofs<<compare_sam.base_file_name<<'~'<<rhs_samrec_iter->cur_line<<reason<<endl;

					if (rhs_samrec_iter->isexonic)
						isexonic = true;

					if (rhs_samrec_iter->isspliced)
						isspliced = true;
				}

				if (rhs_tag_iter->second.size() == 1)
					++compare_sam.unique_notshared;
				else
					++compare_sam.multiple_notshared;

				if (isexonic)
					++compare_sam.exonic_notshared;

				if (isspliced)
					++compare_sam.spliced_notshared;

				//*matched_ofs<<endl;

			}

			ReadNextTag();

			compare_sam.ReadNextTag();
		}
	}

	string basecompname = base_file_name + "_" + compare_sam.base_file_name;

	cout <<basecompname<<"\tbetter_mis\tbetter_qual\tbetter_intron\tworse_mis\tworse_qual\tworse_intron\tequal\t" << endl;

	cout << "one_v_one_notmatched\t" << one_v_one_notmatched_better_mis<<'\t'<<one_v_one_notmatched_better_qual<<'\t'<<one_v_one_notmatched_better_intron<<'\t'
		<<one_v_one_notmatched_worse_mis<<'\t'<<one_v_one_notmatched_worse_qual<<'\t'<<one_v_one_notmatched_worse_intron<<'\t'<<one_v_one_notmatched_equal<< endl;
	cout << "one_v_N_notmatched\t" << one_v_N_notmatched_better_mis<<'\t'<<one_v_N_notmatched_better_qual<<'\t'<<one_v_N_notmatched_better_intron<<'\t'
		<<one_v_N_notmatched_worse_mis<<'\t'<<one_v_N_notmatched_worse_qual<<'\t'<<one_v_N_notmatched_worse_intron<<'\t'<<one_v_N_notmatched_equal<<'\t'<< endl;
	cout << "N_v_one_notmatched\t" << N_v_one_notmatched_better_mis<<'\t'<<N_v_one_notmatched_better_qual<<'\t'<<N_v_one_notmatched_better_intron<<'\t'
		<<N_v_one_notmatched_worse_mis<<'\t'<<N_v_one_notmatched_worse_qual<<'\t'<<N_v_one_notmatched_worse_intron<<'\t'<<N_v_one_notmatched_equal<<'\t'<< endl;
	cout << "N_v_M_notmatched\t" << N_v_M_notmatched_better_mis<<'\t'<<N_v_M_notmatched_better_qual<<'\t'<<N_v_M_notmatched_better_intron<<'\t'
		<<N_v_M_notmatched_worse_mis<<'\t'<<N_v_M_notmatched_worse_qual<<'\t'<<N_v_M_notmatched_worse_intron<<'\t'<<N_v_M_notmatched_equal<<'\t'<< endl;
	cout << "one_v_one_samestartend\t" << one_v_one_samestartend_better_mis<<'\t'<<one_v_one_samestartend_better_qual<<'\t'<<one_v_one_samestartend_better_intron<<'\t'
		<<one_v_one_samestartend_worse_mis<<'\t'<<one_v_one_samestartend_worse_qual<<'\t'<<one_v_one_samestartend_worse_intron<<'\t'<<one_v_one_samestartend_equal<<'\t'<< endl;
	cout << "one_v_N_samestartend\t" << one_v_N_samestartend_better_mis<<'\t'<<one_v_N_samestartend_better_qual<<'\t'<<one_v_N_samestartend_better_intron<<'\t'
		<<one_v_N_samestartend_worse_mis<<'\t'<<one_v_N_samestartend_worse_qual<<'\t'<<one_v_N_samestartend_worse_intron<<'\t'<<one_v_N_samestartend_equal<<'\t'<< endl;
	cout << "N_v_one_samestartend\t" << N_v_one_samestartend_better_mis<<'\t'<<N_v_one_samestartend_better_qual<<'\t'<<N_v_one_samestartend_better_intron<<'\t'
		<<N_v_one_samestartend_worse_mis<<'\t'<<N_v_one_samestartend_worse_qual<<'\t'<<N_v_one_samestartend_worse_intron<<'\t'<<N_v_one_samestartend_equal<<'\t'<< endl;
	cout << "N_v_M_samestartend\t" << N_v_M_samestartend_better_mis<<'\t'<<N_v_M_samestartend_better_qual<<'\t'<<N_v_M_samestartend_better_intron<<'\t'
		<<N_v_M_samestartend_worse_mis<<'\t'<<N_v_M_samestartend_worse_qual<<'\t'<<N_v_M_samestartend_worse_intron<<'\t'<<N_v_M_samestartend_equal<<'\t'<< endl;

	return 0;
}

size_t 
SamHandler::Compare2Sam(SamHandler& compare_sam, const char* output_file, const char* same_offest_output, size_t anchor_width, size_t range)
{
	string output_root = output_file;

	string mps_vs_comp_1v1 = output_root + ".mps_vs_comp_1v1"; ofstream mps_vs_comp_1v1_ofs(mps_vs_comp_1v1.c_str());

	string mps_vs_comp_1vN = output_root + ".mps_vs_comp_1vN"; ofstream mps_vs_comp_1vN_ofs(mps_vs_comp_1vN.c_str());

	string mps_vs_comp_Nv1_bestmatched = output_root + ".mps_vs_comp_Nv1.bestmatched"; ofstream mps_vs_comp_Nv1_bestmatched_ofs(mps_vs_comp_Nv1_bestmatched.c_str());

	string mps_vs_comp_Nv1_bestnotmatched = output_root + ".mps_vs_comp_Nv1.bestnotmatched"; ofstream mps_vs_comp_Nv1_bestnotmatched_ofs(mps_vs_comp_Nv1_bestnotmatched.c_str());

	string mps_vs_comp_NvM_bestmatched = output_root + ".mps_vs_comp_NvM.bestmatched"; ofstream mps_vs_comp_NvM_bestmatched_ofs(mps_vs_comp_NvM_bestmatched.c_str());

	string mps_vs_comp_NvM_bestnotmatched = output_root + ".mps_vs_comp_NvM.bestnotmatched"; ofstream mps_vs_comp_NvM_bestnotmatched_ofs(mps_vs_comp_NvM_bestnotmatched.c_str());

	string mps_vs_comp_0vN = output_root + ".mps_vs_comp_0vN"; ofstream mps_vs_comp_0vN_ofs(mps_vs_comp_0vN.c_str());

	string mps_vs_comp_Nv0 = output_root + ".mps_vs_comp_Nv0"; ofstream mps_vs_comp_Nv0_ofs(mps_vs_comp_Nv0.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = sam_rec_map.begin(), rhs_tag_iter = compare_sam.sam_rec_map.begin();

	while (lhs_tag_iter != sam_rec_map.end() &&  rhs_tag_iter != compare_sam.sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		if (lhs_tag_iter->first < rhs_tag_iter->first)
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				mps_vs_comp_Nv0_ofs << samrec_iter->tostring()<<endl;
			}

			if (lhs_tag_iter->second.front().isunique)
				++unique_notshared;
			else
				++multiple_notshared;

			if (lhs_tag_iter->second.front().isexonic)
				++exonic_notshared;
			
			if (lhs_tag_iter->second.front().isspliced)
				++spliced_notshared;

			++lhs_tag_iter;
		}
		else if (lhs_tag_iter->first > rhs_tag_iter->first)
		{
			for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
			{
				mps_vs_comp_0vN_ofs << samrec_iter->tostring()<<endl;
			}

			if (rhs_tag_iter->second.front().isunique)
				++compare_sam.unique_notshared;
			else
				++compare_sam.multiple_notshared;

			if (rhs_tag_iter->second.front().isexonic)
				++compare_sam.exonic_notshared;
			
			if (rhs_tag_iter->second.front().isspliced)
				++compare_sam.spliced_notshared;

			++rhs_tag_iter;
		}
		else
		{
			vector<SamRec>::iterator lhs_samrec_iter = lhs_tag_iter->second.begin(), rhs_samrec_iter = rhs_tag_iter->second.begin();

			bool matched = false;

			bool bestmatched = false;

			int match_id = 1;

			while (lhs_samrec_iter != lhs_tag_iter->second.end() && rhs_samrec_iter != rhs_tag_iter->second.end())
			{
				bool l_comp_r = comp_tagname_offset(*lhs_samrec_iter, *rhs_samrec_iter);
				bool r_comp_l = comp_tagname_offset(*rhs_samrec_iter, *lhs_samrec_iter);
				if (l_comp_r == true && r_comp_l == false)
				{
					++lhs_samrec_iter;
				}
				else if (l_comp_r == false && r_comp_l == true)
				{
					++rhs_samrec_iter;
				}
				else/* if (l_comp_r && r_comp_l)*/
				{
					lhs_samrec_iter->matched_id = match_id;

					rhs_samrec_iter->matched_id = match_id;

					if (!bestmatched && lhs_samrec_iter->best)
						bestmatched = true;

					++match_id;

					++lhs_samrec_iter;
					++rhs_samrec_iter;
					matched = true;
				}
				//else
				//{
				//	cout <<"impossible comp result"<<endl;
				//	cout << lhs_samrec_iter->tostring()<<endl;
				//	cout << rhs_samrec_iter->tostring()<<endl;
				//	getchar();
				//}
			}

			if (matched)
			{
				ofstream* matched_ofs;

				if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() == 1)
				{
					matched_ofs = &mps_vs_comp_1v1_ofs;
				}
				else if (lhs_tag_iter->second.size() == 1 && rhs_tag_iter->second.size() > 1)
				{
					matched_ofs = &mps_vs_comp_1vN_ofs;
				}
				else if (lhs_tag_iter->second.size() > 1 && rhs_tag_iter->second.size() == 1)
				{
					if (bestmatched)
						matched_ofs = &mps_vs_comp_Nv1_bestmatched_ofs;
					else
						matched_ofs = &mps_vs_comp_Nv1_bestnotmatched_ofs;
				}
				else
				{
					if (bestmatched)
						matched_ofs = &mps_vs_comp_NvM_bestmatched_ofs;
					else
						matched_ofs = &mps_vs_comp_NvM_bestnotmatched_ofs;
				}

                *matched_ofs <<"mapsplice"<<endl;
				for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
				{
					*matched_ofs<<lhs_samrec_iter->tostring()<<'\t'<<lhs_samrec_iter->best<<'\t'<<lhs_samrec_iter->matched_id<<'\t'
							<< lhs_samrec_iter->filter_score<<'\t';

					vector<string>::iterator juncname_iter;

					for (juncname_iter = lhs_samrec_iter->junc_id.begin(); juncname_iter != lhs_samrec_iter->junc_id.end(); ++juncname_iter)
					{
						*matched_ofs << *juncname_iter<<'\t';
					}
					*matched_ofs <<endl;
				}
				*matched_ofs<<"tophat"<<endl;

				for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
				{
					*matched_ofs<<rhs_samrec_iter->tostring()<<'\t' <<rhs_samrec_iter->matched_id<<endl;
				}

				*matched_ofs<<endl;

				if (lhs_tag_iter->second.front().isunique)
					++unique_shared;
				else
					++multiple_shared;

				if (lhs_tag_iter->second.front().isexonic)
					++exonic_shared;

				if (lhs_tag_iter->second.front().isspliced)
					++spliced_shared;

				if (rhs_tag_iter->second.front().isunique)
					++compare_sam.unique_shared;
				else
					++compare_sam.multiple_shared;

				if (rhs_tag_iter->second.front().isexonic)
					++compare_sam.exonic_shared;

				if (rhs_tag_iter->second.front().isspliced)
					++compare_sam.spliced_shared;

			}
			else
			{
				for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
				{
					mps_vs_comp_Nv0_ofs<<lhs_samrec_iter->tostring()<<endl;
				}

				for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
				{
					mps_vs_comp_0vN_ofs<<rhs_samrec_iter->tostring()<<endl;
				}

				if (lhs_tag_iter->second.front().isunique)
					++unique_notshared;
				else
					++multiple_notshared;

				if (lhs_tag_iter->second.front().isexonic)
					++exonic_notshared;

				if (lhs_tag_iter->second.front().isspliced)
					++spliced_notshared;

				if (rhs_tag_iter->second.front().isunique)
					++compare_sam.unique_notshared;
				else
					++compare_sam.multiple_notshared;

				if (rhs_tag_iter->second.front().isexonic)
					++compare_sam.exonic_notshared;

				if (rhs_tag_iter->second.front().isspliced)
					++compare_sam.spliced_notshared;
			}

			++rhs_tag_iter;
			++lhs_tag_iter;
		}	  
	}

	for (; lhs_tag_iter != sam_rec_map.end(); ++lhs_tag_iter)
	{
		vector<SamRec>::iterator lhs_samrec_iter;
		for (lhs_samrec_iter = lhs_tag_iter->second.begin(); lhs_samrec_iter !=  lhs_tag_iter->second.end(); ++lhs_samrec_iter)
		{
			mps_vs_comp_Nv0_ofs<<lhs_samrec_iter->tostring()<<endl;
		}

		if (lhs_tag_iter->second.front().isunique)
			++unique_notshared;
		else
			++multiple_notshared;

		if (lhs_tag_iter->second.front().isexonic)
			++exonic_notshared;

		if (lhs_tag_iter->second.front().isspliced)
			++spliced_notshared;

	}

	for (; rhs_tag_iter != compare_sam.sam_rec_map.end(); ++rhs_tag_iter)
	{
		vector<SamRec>::iterator rhs_samrec_iter;
		for (rhs_samrec_iter = rhs_tag_iter->second.begin(); rhs_samrec_iter !=  rhs_tag_iter->second.end(); ++rhs_samrec_iter)
		{
			mps_vs_comp_0vN_ofs<<rhs_samrec_iter->tostring()<<endl;
		}

		if (rhs_tag_iter->second.front().isunique)
			++compare_sam.unique_notshared;
		else
			++compare_sam.multiple_notshared;

		if (rhs_tag_iter->second.front().isexonic)
			++compare_sam.exonic_notshared;

		if (rhs_tag_iter->second.front().isspliced)
			++compare_sam.spliced_notshared;
	}

	return 0;
}



size_t
SamHandler::CompareTwoSamByOffset(const char* cur_file, const char* syn_file, const char* output, size_t anchor_width, size_t range)
{
	SamHandler cur_sam, syn_sam;

	cur_sam.ReadAppendSamFile(cur_file);

	cur_sam.SortSamByTagnameOffset();

	syn_sam.ReadAppendSamFile(syn_file);

	syn_sam.SortSamByTagnameOffset();

	string shared_1v1 = output; shared_1v1.append(".same_offset_1v1"); ofstream ofs_shared_1v1(shared_1v1.c_str());

	string shared_1vN = output; shared_1vN.append(".same_offset_1vN"); ofstream ofs_shared_1vN(shared_1vN.c_str());

	string incur = output; incur.append(".in_"); incur.append(base_name(cur_file)); ofstream ofs_incur(incur.c_str());

	string insyn = output; insyn.append(".in_"); insyn.append(base_name(syn_file)); ofstream ofs_insyn(insyn.c_str());

	map<string, vector<SamRec> >::iterator lhs_tag_iter = cur_sam.sam_rec_map.begin(), rhs_tag_iter = syn_sam.sam_rec_map.begin();

	while (lhs_tag_iter != cur_sam.sam_rec_map.end() &&  rhs_tag_iter != syn_sam.sam_rec_map.end())
	{
		vector<SamRec>::iterator samrec_iter;

		if (lhs_tag_iter->first < rhs_tag_iter->first)
		{
			for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
			{
				ofs_incur << samrec_iter->tostring()<<endl;
			}

			++lhs_tag_iter;
		}
		else if (lhs_tag_iter->first > rhs_tag_iter->first)
		{
			for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
			{
				ofs_insyn << samrec_iter->tostring()<<endl;
			}

			++rhs_tag_iter;
		}
		else
		{
			vector<SamRec>::iterator lhs_samrec_iter = lhs_tag_iter->second.begin(), rhs_samrec_iter = rhs_tag_iter->second.begin();

			if (rhs_tag_iter->second.size() != 1)
				cout <<"warning synthetic reads mapped more than 1 time"<<endl;

			bool matched = false;

			for (; lhs_samrec_iter != lhs_tag_iter->second.end(); ++lhs_samrec_iter)
			{
				if (lhs_samrec_iter->start <= rhs_samrec_iter->start + range && lhs_samrec_iter->start + range >= rhs_samrec_iter->start && lhs_samrec_iter->chrom_name == rhs_samrec_iter->chrom_name)
				{
					matched = true;

					if (lhs_tag_iter->second.size() == 1)
					{
						ofs_shared_1v1 << lhs_tag_iter->second.front().tostring() <<endl;
						ofs_shared_1v1 << rhs_samrec_iter->tostring() <<endl<<endl;
					}
					else
					{
						if ((lhs_samrec_iter)->best == 0/* && rhs_samrec_iter->spliceway_vec.size() > 1*/)
						{
							for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
							{
								ofs_shared_1vN << samrec_iter->tostring()<<'\t'<<samrec_iter->best <<'\t'
									<< samrec_iter->filter_score<<endl;
							}

							ofs_shared_1vN << rhs_samrec_iter->tostring() <<'\t'<<(lhs_samrec_iter)->best<< endl<<endl;
						}
					}

					break;
				}				
			}

			if (!matched)
			{
				for (samrec_iter = lhs_tag_iter->second.begin(); samrec_iter != lhs_tag_iter->second.end(); ++samrec_iter)
				{
					ofs_incur << samrec_iter->tostring()<<endl;
				}

				for (samrec_iter = rhs_tag_iter->second.begin(); samrec_iter != rhs_tag_iter->second.end(); ++samrec_iter)
				{
					ofs_insyn << samrec_iter->tostring()<<endl;
				}
			}

			++lhs_tag_iter;
			++rhs_tag_iter;
		}
	}

	SamHandler syn_only_sam;

	cout << insyn << endl;

	//cout << "anchor width " << anchor_width << endl;

	//cout << "number of records before filtered: "<< syn_only_sam.ReadAppendSamFile(insyn.c_str())<<endl;

	//syn_only_sam.SortSamByTagnameOffset();

	//cout << "number of records filtered: "<< syn_only_sam.FilterByAnchor(insyn.c_str(), anchor_width)<<endl;

	return 0;

}

void
SamHandler::WriteSamByChromoOffset(const char* samfile)
{
	map<string, vector<SamRec> >::iterator tag_iter;

	vector<SamRec* > sam_rec_vec;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		vector<SamRec>::iterator samrec_iter;

		for (samrec_iter = tag_iter->second.begin(); samrec_iter != tag_iter->second.end(); ++samrec_iter)
		{
			sam_rec_vec.push_back(&(*samrec_iter));
		}
	}
	sort(sam_rec_vec.begin(), sam_rec_vec.end(), comp_offset);

	ofstream ofs(samfile);

	vector<SamRec* >::iterator sam_iter;

	for (sam_iter = sam_rec_vec.begin(); sam_iter != sam_rec_vec.end(); ++sam_iter)
	{
		ofs << (*sam_iter)->tostring() << endl;
	}

	ofs.close();
}

void
SamHandler::SortSamByTagnameOffset()
{
	map<string, vector<SamRec> >::iterator tag_iter;

	for (tag_iter = sam_rec_map.begin(); tag_iter != sam_rec_map.end(); ++tag_iter)
	{
		sort(tag_iter->second.begin(), tag_iter->second.end(), comp_tagname_offset);
	}
}

void 
SamHandler::WriteSamByTagnameOffset(const char* samfile)
{
	SortSamByTagnameOffset();

	WriteSamFile(samfile);
}