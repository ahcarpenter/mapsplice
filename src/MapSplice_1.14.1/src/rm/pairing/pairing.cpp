/*    
 *    pairing.cpp		
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

	vector<pair<size_t, size_t> > spliceway_vec;

	char paired;

	size_t mate_offset;

	int mate_diff;

	bool wrong_format;

	size_t intron_size;

	SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, const string& alt, const string& qualstr) : 
	tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), mis_match(mismatch), qual_str(qualstr), alters(alt), paired('*'), mate_offset(0), mate_diff(0), intron_size(0)
	{
		size_t index = 0;

		wrong_format = false;

		while (true)
		{
			int maplen, intron;

			if (index == 0)
			{
				sscanf(splice_way.c_str() + index, "%dM", &maplen);
				spliceway_vec.push_back(make_pair(start, maplen));

				if (maplen <=0)
				{
					wrong_format =true;
					break;
				}

				index = splice_way.find("M", index) + 1;

				//cout << "index "<<index << '\t' << maplen << endl;getchar();

			}
			else
			{
				if (index == splice_way.length())
					break;

				sscanf(splice_way.c_str() + index, "%dN%dM", &intron, &maplen);

				if (maplen <=0 || intron <= 0)
				{
					wrong_format = true;
					break;
				}

				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));


				intron_size += intron;

				//if (intron < min_intron || intron > max_intron)
				//	wrong_format = true;

				if (splice_way.find("M", index) == string::npos || splice_way.find("M", index) == splice_way.length() - 1)
					break;

				index = splice_way.find("M", index) + 1;

				//cout << "index "<<index << '\t' << intron  <<'\t' << maplen << endl;getchar();
			}
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
	return lhs.intron_size < rhs.intron_size;
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

								paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
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
								if (lower_bound->second.size() > 1)
									sort(lower_bound->second.begin(), lower_bound->second.end(), comp);

								if (offset_iter1->second.size() > 1)
									sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

								lower_bound->second.front().paired = '=';

								lower_bound->second.front().mate_offset = offset_iter1->second.front().start;

								lower_bound->second.front().mate_diff = -(int)dist_upper;

								offset_iter1->second.front().paired = '=';

								offset_iter1->second.front().mate_offset = lower_bound->second.front().start;

								offset_iter1->second.front().mate_diff = (int)dist_upper;

								paired = true;

								vector<SamRec>::iterator sam_rec_iter;

								for (sam_rec_iter = lower_bound->second.begin(); sam_rec_iter != lower_bound->second.end(); ++sam_rec_iter)
								{
									paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
									//paired_reads_ofs << sam_rec_iter->tostring() << endl;
								}

								for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
								{
									paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
									//paired_reads_ofs << sam_rec_iter->tostring() << endl;
								}

								map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
								pair_iter1->second.erase(offset_iter1_temp);

								pair2_map_ref.erase(lower_bound);
							}
							else if (dist_lower < dist_upper && dist_lower < max_pair_dist)
							{
								if (lower_bound_bef->second.size() > 1)
									sort(lower_bound_bef->second.begin(), lower_bound_bef->second.end(), comp);

								if (offset_iter1->second.size() > 1)
									sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

								lower_bound_bef->second.front().paired = '=';

								lower_bound_bef->second.front().mate_offset = offset_iter1->second.front().start;

								lower_bound_bef->second.front().mate_diff = (int)dist_lower;

								offset_iter1->second.front().paired = '=';

								offset_iter1->second.front().mate_offset = lower_bound_bef->second.front().start;

								offset_iter1->second.front().mate_diff = -(int)dist_lower;

								vector<SamRec>::iterator sam_rec_iter;

								paired = true;

								for (sam_rec_iter = lower_bound_bef->second.begin(); sam_rec_iter != lower_bound_bef->second.end(); ++sam_rec_iter)
								{
									//paired_reads_ofs << sam_rec_iter->tostring() << endl;

									paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));								
								}

								for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
								{
									//paired_reads_ofs << sam_rec_iter->tostring() << endl;

									paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								}

								map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
								pair_iter1->second.erase(offset_iter1_temp);

								pair2_map_ref.erase(lower_bound_bef);
							}
							else
								++offset_iter1;
						}

						for (offset_iter1 = pair_iter1->second.begin(); offset_iter1 != pair_iter1->second.end(); ++offset_iter1)
						{
							vector<SamRec>::iterator sam_rec_iter;

							for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
							{
								paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								//paired_reads_ofs << sam_rec_iter->tostring() << endl;
							}
						}

						for (offset_iter1 = pair_iter2->second.begin(); offset_iter1 != pair_iter2->second.end(); ++offset_iter1)
						{
							vector<SamRec>::iterator sam_rec_iter;

							for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
							{
								paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
								//paired_reads_ofs << sam_rec_iter->tostring() << endl;
							}
						}
						
					}
					else
					{
						cout << "pair id size exceed 2:" << chrom_iter->second.size() << endl;
					}
				}

				map<int, vector<pair<char, string> > >::iterator pairid_iter;

				for (pairid_iter = paired_reads.begin(); pairid_iter != paired_reads.end(); ++pairid_iter)
				{
					vector<pair<char, string> >::iterator mapped_iter;

					for (mapped_iter = pairid_iter->second.begin(); mapped_iter != pairid_iter->second.end(); ++mapped_iter)
					{
						if (!paired)
						{
							if (end1_mapped && end2_mapped)
								fusion_paired_reads_ofs << mapped_iter->second<<endl;
							else
								single_reads_ofs << mapped_iter->second<<endl;
						}
						else if (mapped_iter->first == '=')
							paired_reads_ofs << mapped_iter->second<<endl;
						else
							filtered_paired_reads_ofs << mapped_iter->second<<endl;
					}
				}

				paired_reads.clear();

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

					paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
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
					sort(lower_bound->second.begin(), lower_bound->second.end(), comp);

					sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

					lower_bound->second.front().paired = '=';

					lower_bound->second.front().mate_offset = offset_iter1->second.front().start;

					lower_bound->second.front().mate_diff = -(int)dist_upper;

					offset_iter1->second.front().paired = '=';

					offset_iter1->second.front().mate_offset = lower_bound->second.front().start;

					offset_iter1->second.front().mate_diff = (int)dist_upper;

					paired = true;

					vector<SamRec>::iterator sam_rec_iter;

					for (sam_rec_iter = lower_bound->second.begin(); sam_rec_iter != lower_bound->second.end(); ++sam_rec_iter)
					{
						paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
						//paired_reads_ofs << sam_rec_iter->tostring() << endl;
					}

					for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
					{
						paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
						//paired_reads_ofs << sam_rec_iter->tostring() << endl;
					}

					map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
					pair_iter1->second.erase(offset_iter1_temp);

					pair2_map_ref.erase(lower_bound);
				}
				else if (dist_lower < dist_upper && dist_lower < max_pair_dist)
				{
					sort(lower_bound_bef->second.begin(), lower_bound_bef->second.end(), comp);

					sort(offset_iter1->second.begin(), offset_iter1->second.end(), comp);

					lower_bound_bef->second.front().paired = '=';

					lower_bound_bef->second.front().mate_offset = offset_iter1->second.front().start;

					lower_bound_bef->second.front().mate_diff = (int)dist_lower;

					offset_iter1->second.front().paired = '=';

					offset_iter1->second.front().mate_offset = lower_bound_bef->second.front().start;

					offset_iter1->second.front().mate_diff = -(int)dist_lower;

					paired = true;

					vector<SamRec>::iterator sam_rec_iter;

					for (sam_rec_iter = lower_bound_bef->second.begin(); sam_rec_iter != lower_bound_bef->second.end(); ++sam_rec_iter)
					{
						//paired_reads_ofs << sam_rec_iter->tostring() << endl;

						paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));						
					}

					for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
					{
						//paired_reads_ofs << sam_rec_iter->tostring() << endl;

						paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					}

					map<size_t, vector<SamRec> >::iterator offset_iter1_temp = offset_iter1++;
					pair_iter1->second.erase(offset_iter1_temp);

					pair2_map_ref.erase(lower_bound_bef);
				}
				else
					++offset_iter1;
			}

			for (offset_iter1 = pair_iter1->second.begin(); offset_iter1 != pair_iter1->second.end(); ++offset_iter1)
			{
				vector<SamRec>::iterator sam_rec_iter;

				for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
				{
					paired_reads[pair_iter1->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					//paired_reads_ofs << sam_rec_iter->tostring() << endl;
				}
			}

			for (offset_iter1 = pair_iter2->second.begin(); offset_iter1 != pair_iter2->second.end(); ++offset_iter1)
			{
				vector<SamRec>::iterator sam_rec_iter;

				for (sam_rec_iter = offset_iter1->second.begin(); sam_rec_iter != offset_iter1->second.end(); ++sam_rec_iter)
				{
					paired_reads[pair_iter2->first].push_back(make_pair(sam_rec_iter->paired, sam_rec_iter->tostring()));
					//paired_reads_ofs << sam_rec_iter->tostring() << endl;
				}
			}

		}
		else
		{
			cout << "pair id size exceed 2:" << chrom_iter->second.size() << endl;
		}
	}

	map<int, vector<pair<char, string> > >::iterator pairid_iter;

	for (pairid_iter = paired_reads.begin(); pairid_iter != paired_reads.end(); ++pairid_iter)
	{
		vector<pair<char, string> >::iterator mapped_iter;

		for (mapped_iter = pairid_iter->second.begin(); mapped_iter != pairid_iter->second.end(); ++mapped_iter)
		{
			if (!paired)
			{
				if (end1_mapped && end2_mapped)
					fusion_paired_reads_ofs << mapped_iter->second<<endl;
				else
					single_reads_ofs << mapped_iter->second<<endl;
			}
			else if (mapped_iter->first == '=')
				paired_reads_ofs << mapped_iter->second<<endl;
			else
				filtered_paired_reads_ofs << mapped_iter->second<<endl;
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