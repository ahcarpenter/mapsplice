#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdlib.h>

using namespace std;

//#define VS

#ifdef VS
#include <hash_map> //vc only
#include <hash_set>
//#define _CRT_SECURE_NO_WARNINGS 
//#pragma warning(disable:_CRT_SECURE_NO_WARNINGS)
 
#else
#include <ext/hash_map> //g++ only
#include <ext/hash_set>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;

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

#ifdef VS
using namespace stdext;
#endif

#define IS_PAIRED 0x0001
#define IS_PAIRED_MAPPED 0x0002
#define IS_UNMAPPED 0x0004
#define MATE_UNMAPPED 0x0008
#define IS_REVERSE 0x0010
#define IS_MATE_REVERSE 0x0020
#define IS_FIRST_END 0x040
#define IS_SECOND_END 0x0080
#define IS_PRIMARY 0x0100
#define IS_FAILED_QUAL_CHECK 0x0200
#define IS_PCR_DUP 0x0400

struct SamRec {
	string tag_name;
	unsigned short strand_t;
	string chrom_name;
	string qual_str;
	string alters;
	size_t start;
	unsigned short confid;
	string splice_way;
	string mapped_seq;
	unsigned short mis_match;

	vector<pair<size_t, int> > spliceway_vec;

	char paired;

	size_t mate_offset;

	int mate_diff;

	bool wrong_format;

	size_t intron_size;

	size_t mappedlen;

	SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, const string& alt, const string& qualstr) : 
	tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), mis_match(mismatch), qual_str(qualstr), alters(alt), paired('*'), mate_offset(0), mate_diff(0), 
	intron_size(0), mappedlen(0)
	{
		size_t index = 0;

		wrong_format = false;

		string flag_str = " ";

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
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
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

					intron_size += maplen;
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
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));

					intron_size += maplen;
				}
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
	}

	string tostring()
	{
		char sam_rec_char[5000];

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s"/*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
			chrom_name.c_str(), start, confid, splice_way.c_str(), paired, mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, alters.c_str()/*, best, filter_score*/);

		return sam_rec_char;
	}
};

bool comp(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.mis_match < rhs.mis_match;
}

struct PairedSamRec {

	PairedSamRec(int dist, unsigned short tm, size_t is, size_t maplen, SamRec* sam_rec1, SamRec* sam_rec2) : mate_dist(dist), intron_size(is), mappedlen(maplen), total_mismatch(tm) 
	{
		paired_sam_rec = make_pair(sam_rec1, sam_rec2);
	}

	int mate_dist;
	size_t intron_size;
	size_t mappedlen;
	unsigned short total_mismatch;
	pair<SamRec*, SamRec*> paired_sam_rec;
};

bool comp_dist(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (abs(lhs.mate_dist) == abs(rhs.mate_dist))
		{
			if (lhs.total_mismatch == rhs.total_mismatch)
				return lhs.intron_size < rhs.intron_size;
			else
				return lhs.total_mismatch < rhs.total_mismatch;
		}
		else
			return abs(lhs.mate_dist) < abs(rhs.mate_dist);
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

void FilterPairedSamRec(vector<PairedSamRec>& pairedsamrec_vec)
{
	if (pairedsamrec_vec.size() == 0)
		return;

	sort(pairedsamrec_vec.begin(), pairedsamrec_vec.end(), comp_dist);

	vector<PairedSamRec>::iterator vps_iter;

//#define IS_PAIRED 0x0001
//#define IS_PAIRED_MAPPED 0x0002
//#define IS_UNMAPPED 0x0004
//#define MATE_UNMAPPED 0x0008
//#define IS_REVERSE 0x0010
//#define IS_MATE_REVERSE 0x0020
//#define IS_FIRST_END 0x040
//#define IS_SECOND_END 0x0080
//#define IS_PRIMARY 0x0100
//#define IS_FAILED_QUAL_CHECK 0x0200
//#define IS_PCR_DUP 0x0400

	vps_iter = pairedsamrec_vec.begin();

	vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

	vps_iter->paired_sam_rec.first->mate_diff = static_cast<int> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

	vps_iter->paired_sam_rec.first->paired = '=';

	if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

	vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

	vps_iter->paired_sam_rec.second->mate_diff = static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

	vps_iter->paired_sam_rec.second->paired = '=';

	if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

	for (vps_iter = pairedsamrec_vec.begin() + 1; vps_iter != pairedsamrec_vec.end(); ++vps_iter)
	{
		if ((vps_iter - 1)->mappedlen == (vps_iter)->mappedlen && (vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (vps_iter - 1)->total_mismatch == (vps_iter)->total_mismatch && (vps_iter - 1)->intron_size == (vps_iter)->intron_size)
		{
			vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

			vps_iter->paired_sam_rec.first->mate_diff = static_cast<int> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

			vps_iter->paired_sam_rec.first->paired = '=';

			if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
				vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

			vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

			vps_iter->paired_sam_rec.second->mate_diff = static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

			vps_iter->paired_sam_rec.second->paired = '=';

			if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
				vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;
		}
		else
			break;
	}
}



void
RemoveDupMapreads(const char* infile, const char* paired_reads_file, const char* fusion_paired_file, const char* single_reads_file, const char* filtered_paired_reads_file, const char* stat_file, size_t max_pair_dist)
{
	ifstream ifs(infile);

	ofstream paired_reads_ofs(paired_reads_file);
	ofstream fusion_paired_reads_ofs(fusion_paired_file);
	ofstream single_reads_ofs(single_reads_file);
	ofstream filtered_paired_reads_ofs(filtered_paired_reads_file);

	string prev_tagname = "";

	//chrom pairid offset line
	hash_map<string, hash_map<int, map<size_t, vector<SamRec> > > > mapped_reads;

	map<int, vector<pair<char, string> > > paired_reads;

	vector<PairedSamRec> paired_reads_ptr;

	//mapped_reads.

	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	bool end1_mapped = false, end2_mapped = false;

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
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

			if (read_count == 10)
				alterstr = alters;

			string tagnamestr = tagname;

			int pairid = atoi(tagnamestr.substr(tagnamestr.length() - 1, 1).c_str());

			size_t st_idx = tagnamestr.find("~");			

			tagnamestr = tagnamestr.substr(st_idx + 1, tagnamestr.length() - 3 - st_idx);			

			++count;

			if (prev_tagname.empty() || prev_tagname == tagnamestr)
			{
				mapped_reads[chrom][pairid][offset].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, alterstr, qual_str));

				if (pairid == 1)
					end1_mapped = true;

				if (pairid == 2)
					end2_mapped = true;

				if (prev_tagname.empty())
					prev_tagname = tagnamestr;

				//tag_count++;
			}
			else
			{
				//hash_map<string, int>::iterator msi_iter;

				hash_map<string, hash_map<int, map<size_t, vector<SamRec> > > >::iterator chrom_iter;

				for (chrom_iter = mapped_reads.begin(); chrom_iter != mapped_reads.end(); ++chrom_iter)
				{
					hash_map<int, map<size_t, vector<SamRec> > >::iterator pair_iter1, pair_iter2;

					if (chrom_iter->second.size() == 1)
					{
						pair_iter1 = chrom_iter->second.begin();

						map<size_t, vector<SamRec> >::iterator offset_iter;
						for (offset_iter = pair_iter1->second.begin(); offset_iter != pair_iter1->second.end(); ++offset_iter)
						{
							vector<SamRec>::iterator sam_rec_iter;

							for (sam_rec_iter = offset_iter->second.begin(); sam_rec_iter != offset_iter->second.end(); ++sam_rec_iter)
							{
								//paired_reads_ofs << sam_rec_iter->tostring() << endl;

								//paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
							}
						}
					}
					else if (chrom_iter->second.size() == 2)
					{
						pair_iter1 = chrom_iter->second.begin();

						pair_iter2 = pair_iter1;

						++pair_iter2;

						//if (pair_iter1->second.size() == 1 && 

						map<size_t, vector<SamRec> >::iterator offset_iter1 = pair_iter1->second.begin();//, offset_iter2 = pair_iter2->second.begin();

						//map<size_t, vector<SamRec> >& pair1_map_ref = pair_iter1->second;

						map<size_t, vector<SamRec> >& pair2_map_ref = pair_iter2->second;

						while (offset_iter1 != pair_iter1->second.end())//; ++offset_iter1)
						{
							map<size_t, vector<SamRec> >::iterator lower_bound = pair2_map_ref.lower_bound(offset_iter1->first);

							map<size_t, vector<SamRec> >::iterator lower_bound_bef;

							size_t dist_upper = -1, dist_lower = -1;

							if (lower_bound != pair2_map_ref.end())
							{
								dist_upper = lower_bound->first - offset_iter1->first;
							}

							if (lower_bound != pair2_map_ref.begin())
							{
								lower_bound_bef = lower_bound;

								--lower_bound_bef;

								dist_lower = offset_iter1->first - lower_bound_bef->first;
							}

							if (dist_upper < dist_lower && dist_upper < max_pair_dist)
							{
								//if (lower_bound->second.size() > 1)
								//	sort(lower_bound->second.begin(), lower_bound->second.end(), comp);

								//if (offset_iter1->second.size() > 1)
								//	sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

								//lower_bound->second.front().paired = '=';

								//lower_bound->second.front().mate_offset = offset_iter1->second.front().start;

								//lower_bound->second.front().mate_diff = -(int)dist_upper;

								//offset_iter1->second.front().paired = '=';

								//offset_iter1->second.front().mate_offset = lower_bound->second.front().start;

								//offset_iter1->second.front().mate_diff = (int)dist_upper;
							
								vector<SamRec>::iterator sam_rec_iter1, sam_rec_iter2;

								for (sam_rec_iter1 = lower_bound->second.begin(); sam_rec_iter1 != lower_bound->second.end(); ++sam_rec_iter1)
								{
									for (sam_rec_iter2 = offset_iter1->second.begin(); sam_rec_iter2 != offset_iter1->second.end(); ++sam_rec_iter2)
									{
										if ((sam_rec_iter1->strand_t ^ sam_rec_iter2->strand_t) && IS_REVERSE)
										{
											paired_reads_ptr.push_back(PairedSamRec(dist_upper, sam_rec_iter1->mis_match + sam_rec_iter2->mis_match, 
											sam_rec_iter1->intron_size + sam_rec_iter2->intron_size, sam_rec_iter1->mappedlen + sam_rec_iter2->mappedlen, &(*sam_rec_iter1), &(*sam_rec_iter2)));

											paired = true;
										}
									}
								}

								//vector<SamRec>::iterator sam_rec_iter;

								//for (sam_rec_iter = lower_bound->second.begin(); sam_rec_iter != lower_bound->second.end(); ++sam_rec_iter)
								//{
								//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;
								//}

								//for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
								//{
								//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;
								//}

								map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
								//pair_iter1->second.erase(offset_iter1_temp);

								//pair2_map_ref.erase(lower_bound);
							}
							else if (dist_lower < dist_upper && dist_lower < max_pair_dist)
							{
								//if (lower_bound_bef->second.size() > 1)
								//	sort(lower_bound_bef->second.begin(), lower_bound_bef->second.end(), comp);

								//if (offset_iter1->second.size() > 1)
								//	sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

								//lower_bound_bef->second.front().paired = '=';

								//lower_bound_bef->second.front().mate_offset = offset_iter1->second.front().start;

								//lower_bound_bef->second.front().mate_diff = (int)dist_lower;

								//offset_iter1->second.front().paired = '=';

								//offset_iter1->second.front().mate_offset = lower_bound_bef->second.front().start;

								//offset_iter1->second.front().mate_diff = -(int)dist_lower;

								vector<SamRec>::iterator sam_rec_iter1, sam_rec_iter2;

								for (sam_rec_iter1 = lower_bound_bef->second.begin(); sam_rec_iter1 != lower_bound_bef->second.end(); ++sam_rec_iter1)
								{
									for (sam_rec_iter2 = offset_iter1->second.begin(); sam_rec_iter2 != offset_iter1->second.end(); ++sam_rec_iter2)
									{
										if ((sam_rec_iter1->strand_t ^ sam_rec_iter2->strand_t) && IS_REVERSE)
										{
											paired_reads_ptr.push_back(PairedSamRec(dist_lower, sam_rec_iter1->mis_match + sam_rec_iter2->mis_match, 
											sam_rec_iter1->intron_size + sam_rec_iter2->intron_size, sam_rec_iter1->mappedlen + sam_rec_iter2->mappedlen, &(*sam_rec_iter1), &(*sam_rec_iter2)));

											paired = true;
										}
									}
								}

								//vector<SamRec>::iterator sam_rec_iter;

								

								//for (sam_rec_iter = lower_bound_bef->second.begin(); sam_rec_iter != lower_bound_bef->second.end(); ++sam_rec_iter)
								//{
								//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;

								//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));								
								//}

								//for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
								//{
								//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;

								//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								//}

								map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
								//pair_iter1->second.erase(offset_iter1_temp);

								//pair2_map_ref.erase(lower_bound_bef);
							}
							else
								++offset_iter1;
						}

						//for (offset_iter1 = pair_iter1->second.begin(); offset_iter1 != pair_iter1->second.end(); ++offset_iter1)
						//{
						//	vector<SamRec>::iterator sam_rec_iter;

						//	for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
						//	{
						//		paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
						//		//paired_reads_ofs << sam_rec_iter->tostring() << endl;
						//	}
						//}

						//for (offset_iter1 = pair_iter2->second.begin(); offset_iter1 != pair_iter2->second.end(); ++offset_iter1)
						//{
						//	vector<SamRec>::iterator sam_rec_iter;

						//	for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
						//	{
						//		paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
						//		//paired_reads_ofs << sam_rec_iter->tostring() << endl;
						//	}
						//}
						
					}
					else
					{
						cout << "pair id size exceed 2:" << chrom_iter->second.size() << endl;
					}
				}

				FilterPairedSamRec(paired_reads_ptr);

				//hash_map<string, hash_map<int, map<size_t, vector<SamRec> > > >::iterator chrom_iter;

				for (chrom_iter = mapped_reads.begin(); chrom_iter != mapped_reads.end(); ++chrom_iter)
				{
					hash_map<int, map<size_t, vector<SamRec> > >::iterator pair_iter;

					for (pair_iter = chrom_iter->second.begin(); pair_iter != chrom_iter->second.end(); ++pair_iter)
					{
						map<size_t, vector<SamRec> >::iterator offset_iter; 

						for (offset_iter = pair_iter->second.begin(); offset_iter != pair_iter->second.end(); ++offset_iter)
						{
							vector<SamRec>::iterator samrec_iter;

							for (samrec_iter = offset_iter->second.begin(); samrec_iter != offset_iter->second.end(); ++samrec_iter)
							{
								if (!paired)
								{
//#define IS_PAIRED 0x0001
//#define IS_PAIRED_MAPPED 0x0002
//#define IS_UNMAPPED 0x0004
//#define MATE_UNMAPPED 0x0008
//#define IS_REVERSE 0x0010
//#define IS_MATE_REVERSE 0x0020
//#define IS_FIRST_END 0x040
//#define IS_SECOND_END 0x0080
//#define IS_PRIMARY 0x0100
//#define IS_FAILED_QUAL_CHECK 0x0200
//#define IS_PCR_DUP 0x0400

									if (end1_mapped && end2_mapped)
										fusion_paired_reads_ofs << samrec_iter->tostring()<<endl;
									else
									{
										samrec_iter->strand_t |= MATE_UNMAPPED;

										single_reads_ofs << samrec_iter->tostring()<<endl;
									}
								}
								else if (samrec_iter->paired == '=')
									paired_reads_ofs << samrec_iter->tostring()<<endl;
								else
									filtered_paired_reads_ofs << samrec_iter->tostring()<<endl;
							}
						}
					}
				}

				paired_reads.clear();

				paired_reads_ptr.clear();

				mapped_reads.clear();

				end1_mapped = false, end2_mapped = false;

				mapped_reads[chrom][pairid][offset].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, alterstr, qual_str));

				if (pairid == 1)
					end1_mapped = true;

				if (pairid == 2)
					end2_mapped = true;

				paired = false;				

				prev_tagname = tagnamestr;
			}			
		}
		ifs.close();
	}
	else
	{
		cout << "can't open file "<< infile<<endl; exit(1);
	}

	hash_map<string, hash_map<int, map<size_t, vector<SamRec> > > >::iterator chrom_iter;

	for (chrom_iter = mapped_reads.begin(); chrom_iter != mapped_reads.end(); ++chrom_iter)
	{
		hash_map<int, map<size_t, vector<SamRec> > >::iterator pair_iter1, pair_iter2;

		if (chrom_iter->second.size() == 1)
		{
			pair_iter1 = chrom_iter->second.begin();

			map<size_t, vector<SamRec> >::iterator offset_iter;
			for (offset_iter = pair_iter1->second.begin(); offset_iter != pair_iter1->second.end(); ++offset_iter)
			{
				vector<SamRec>::iterator sam_rec_iter;

				for (sam_rec_iter = offset_iter->second.begin(); sam_rec_iter != offset_iter->second.end(); ++sam_rec_iter)
				{
					//paired_reads_ofs << sam_rec_iter->tostring() << endl;

					//paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
				}
			}
		}
		else if (chrom_iter->second.size() == 2)
		{
			pair_iter1 = chrom_iter->second.begin();

			pair_iter2 = pair_iter1;

			++pair_iter2;

			//if (pair_iter1->second.size() == 1 && 

			map<size_t, vector<SamRec> >::iterator offset_iter1 = pair_iter1->second.begin();//, offset_iter2 = pair_iter2->second.begin();

			//map<size_t, vector<SamRec> >& pair1_map_ref = pair_iter1->second;

			map<size_t, vector<SamRec> >& pair2_map_ref = pair_iter2->second;

			while (offset_iter1 != pair_iter1->second.end())//; ++offset_iter1)
			{
				map<size_t, vector<SamRec> >::iterator lower_bound = pair2_map_ref.lower_bound(offset_iter1->first);

				map<size_t, vector<SamRec> >::iterator lower_bound_bef;

				size_t dist_upper = -1, dist_lower = -1;

				if (lower_bound != pair2_map_ref.end())
				{
					dist_upper = lower_bound->first - offset_iter1->first;
				}

				if (lower_bound != pair2_map_ref.begin())
				{
					lower_bound_bef = lower_bound;

					--lower_bound_bef;

					dist_lower = offset_iter1->first - lower_bound_bef->first;
				}

				if (dist_upper < dist_lower && dist_upper < max_pair_dist)
				{
					//if (lower_bound->second.size() > 1)
					//	sort(lower_bound->second.begin(), lower_bound->second.end(), comp);

					//if (offset_iter1->second.size() > 1)
					//	sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

					//lower_bound->second.front().paired = '=';

					//lower_bound->second.front().mate_offset = offset_iter1->second.front().start;

					//lower_bound->second.front().mate_diff = -(int)dist_upper;

					//offset_iter1->second.front().paired = '=';

					//offset_iter1->second.front().mate_offset = lower_bound->second.front().start;

					//offset_iter1->second.front().mate_diff = (int)dist_upper;
					
					vector<SamRec>::iterator sam_rec_iter1, sam_rec_iter2;

					for (sam_rec_iter1 = lower_bound->second.begin(); sam_rec_iter1 != lower_bound->second.end(); ++sam_rec_iter1)
					{
						for (sam_rec_iter2 = offset_iter1->second.begin(); sam_rec_iter2 != offset_iter1->second.end(); ++sam_rec_iter2)
						{
							if ((sam_rec_iter1->strand_t ^ sam_rec_iter2->strand_t) && IS_REVERSE)
							{
								paired_reads_ptr.push_back(PairedSamRec(dist_upper, sam_rec_iter1->mis_match + sam_rec_iter2->mis_match, 
								sam_rec_iter1->intron_size + sam_rec_iter2->intron_size, sam_rec_iter1->mappedlen + sam_rec_iter2->mappedlen, &(*sam_rec_iter1), &(*sam_rec_iter2)));

								paired = true;
							}
						}
					}

					//vector<SamRec>::iterator sam_rec_iter;

					//for (sam_rec_iter = lower_bound->second.begin(); sam_rec_iter != lower_bound->second.end(); ++sam_rec_iter)
					//{
					//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;
					//}

					//for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
					//{
					//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;
					//}

					map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
					//pair_iter1->second.erase(offset_iter1_temp);

					//pair2_map_ref.erase(lower_bound);
				}
				else if (dist_lower < dist_upper && dist_lower < max_pair_dist)
				{
					//if (lower_bound_bef->second.size() > 1)
					//	sort(lower_bound_bef->second.begin(), lower_bound_bef->second.end(), comp);

					//if (offset_iter1->second.size() > 1)
					//	sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

					//lower_bound_bef->second.front().paired = '=';

					//lower_bound_bef->second.front().mate_offset = offset_iter1->second.front().start;

					//lower_bound_bef->second.front().mate_diff = (int)dist_lower;

					//offset_iter1->second.front().paired = '=';

					//offset_iter1->second.front().mate_offset = lower_bound_bef->second.front().start;

					//offset_iter1->second.front().mate_diff = -(int)dist_lower;

					vector<SamRec>::iterator sam_rec_iter1, sam_rec_iter2;

					for (sam_rec_iter1 = lower_bound_bef->second.begin(); sam_rec_iter1 != lower_bound_bef->second.end(); ++sam_rec_iter1)
					{
						for (sam_rec_iter2 = offset_iter1->second.begin(); sam_rec_iter2 != offset_iter1->second.end(); ++sam_rec_iter2)
						{
							if ((sam_rec_iter1->strand_t ^ sam_rec_iter2->strand_t) && IS_REVERSE)
							{
								paired_reads_ptr.push_back(PairedSamRec(dist_lower, sam_rec_iter1->mis_match + sam_rec_iter2->mis_match, 
								sam_rec_iter1->intron_size + sam_rec_iter2->intron_size, sam_rec_iter1->mappedlen + sam_rec_iter2->mappedlen, &(*sam_rec_iter1), &(*sam_rec_iter2)));

								paired = true;
							}
						}
					}

					//vector<SamRec>::iterator sam_rec_iter;

					

					//for (sam_rec_iter = lower_bound_bef->second.begin(); sam_rec_iter != lower_bound_bef->second.end(); ++sam_rec_iter)
					//{
					//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;

					//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));								
					//}

					//for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
					//{
					//	//paired_reads_ofs << sam_rec_iter->tostring() << endl;

					//	paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					//}

					map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
					//pair_iter1->second.erase(offset_iter1_temp);

					//pair2_map_ref.erase(lower_bound_bef);
				}
				else
					++offset_iter1;
			}

			//for (offset_iter1 = pair_iter1->second.begin(); offset_iter1 != pair_iter1->second.end(); ++offset_iter1)
			//{
			//	vector<SamRec>::iterator sam_rec_iter;

			//	for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
			//	{
			//		paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
			//		//paired_reads_ofs << sam_rec_iter->tostring() << endl;
			//	}
			//}

			//for (offset_iter1 = pair_iter2->second.begin(); offset_iter1 != pair_iter2->second.end(); ++offset_iter1)
			//{
			//	vector<SamRec>::iterator sam_rec_iter;

			//	for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
			//	{
			//		paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
			//		//paired_reads_ofs << sam_rec_iter->tostring() << endl;
			//	}
			//}

		}
		else
		{
			cout << "pair id size exceed 2:" << chrom_iter->second.size() << endl;
		}
	}

	FilterPairedSamRec(paired_reads_ptr);

	//hash_map<string, hash_map<int, map<size_t, vector<SamRec> > > >::iterator chrom_iter;

	for (chrom_iter = mapped_reads.begin(); chrom_iter != mapped_reads.end(); ++chrom_iter)
	{
		hash_map<int, map<size_t, vector<SamRec> > >::iterator pair_iter;

		for (pair_iter = chrom_iter->second.begin(); pair_iter != chrom_iter->second.end(); ++pair_iter)
		{
			map<size_t, vector<SamRec> >::iterator offset_iter; 

			for (offset_iter = pair_iter->second.begin(); offset_iter != pair_iter->second.end(); ++offset_iter)
			{
				vector<SamRec>::iterator samrec_iter;

				for (samrec_iter = offset_iter->second.begin(); samrec_iter != offset_iter->second.end(); ++samrec_iter)
				{
					if (!paired)
					{
						if (end1_mapped && end2_mapped)
							fusion_paired_reads_ofs << samrec_iter->tostring()<<endl;
						else
						{
							samrec_iter->strand_t |= MATE_UNMAPPED;

							single_reads_ofs << samrec_iter->tostring()<<endl;
						}
					}
					else if (samrec_iter->paired == '=')
						paired_reads_ofs << samrec_iter->tostring()<<endl;
					else
						filtered_paired_reads_ofs << samrec_iter->tostring()<<endl;
				}
			}
		}
	}

	//ofstream stat_ofs(stat_file, ios_base::app);

	//stat_ofs << "Remove duplication"<<endl << "total_mapped: " << count<<endl<<"unspliced: "<< unspliced <<endl<< 
	//	"unique_spliced_reads: " << unique_spliced<<endl<<"multiple_spliced_reads: "<<  multiple_spliced<<endl;

	//ofs.close();

}

int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << "infile unspliced_reads unique_spliced_reads multiple_spliced_reads" <<endl;
		exit(0);
	}
	const char* infile = argv[1];
	const char* paired_reads = argv[2];
	const char* fusion_paired_reads = argv[3];
	const char* single_reads = argv[4];
	const char* filtered_paired_reads = argv[5];
	const char* stat_file = argv[6];
	int max_pair_dist = atoi(argv[7]);

	RemoveDupMapreads(infile, paired_reads, fusion_paired_reads, single_reads, filtered_paired_reads, stat_file, max_pair_dist);
	//RemoveDupMapreads(infile, unspliced_reads, unique_spliced_reads, multiple_spliced_reads, stat_file);
}