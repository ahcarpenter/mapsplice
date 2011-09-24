/*    
 *    ReadRegions.cpp		
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

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string.h>
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

struct PairedRegions
{
	string chrom1;

	size_t strand1;

	size_t start1, end1;

	hash_set<size_t> mapped_reads_index;

	string chrom2;

	size_t strand2;

	size_t start2, end2;

	//vector<size_t> mapped_reads_index2;	

	PairedRegions(const string& chr1, size_t sd1, size_t st1, size_t en1, const string& chr2,  size_t sd2, size_t st2, size_t en2) : chrom1(chr1), strand1(sd1), start1(st1), end1(en1), chrom2(chr2), strand2(sd2), start2(st2), end2(en2)
	{
	}

	string tostring()
	{
		char output_chr[5000];

		sprintf(output_chr, "%s\t%llu\t%llu\t%llu\t%s\t%llu\t%llu\t%llu", chrom1.c_str(), strand1, start1, end1, chrom2.c_str(), strand2, start2, end2);

		return output_chr;		
	}
};

inline char
complement(int i) {
	static const int b2c_size = 20;
	static const char b2c[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
	};
	static const char b2cl[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
	};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

inline string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

typedef vector<PairedRegions> PAIRED_REGIONS_VEC;
typedef vector<PairedRegions>::iterator PAIRED_REGIONS_ITERATOR;
typedef PairedRegions* PAIRED_REGIONS_PTR;
typedef vector<PAIRED_REGIONS_PTR> PAIRED_REGIONS_PTR_VEC;
typedef vector<PAIRED_REGIONS_PTR>::iterator PAIRED_REGIONS_PTR_VEC_ITERATOR;
typedef hash_map<string, PAIRED_REGIONS_PTR_VEC > HASH_CHROM_PAIRED_REGIONS_PTR_VEC;
typedef hash_map<string, PAIRED_REGIONS_PTR_VEC >::iterator HASH_CHROM_PAIRED_REGIONS_PTR_VEC_ITERATOR;

size_t max_dist = 0;

int ReadRegions(const char* region_file, PAIRED_REGIONS_VEC& paired_regions, size_t ext)
{
	ifstream ifs(region_file);

	int count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line.empty())
				continue;

			char chrom1[1000];

			size_t strand1;
			
			size_t start1, end1;

			char chrom2[1000];

			size_t strand2;

			size_t start2, end2;

			sscanf(line.c_str(), "%s\t%llu\t%llu\t%llu\t%s\t%llu\t%llu\t%llu", chrom1, &strand1, &start1, &end1, chrom2, &strand2, &start2, &end2);

			string chrom1str = chrom1;

			string chrom2str = chrom2;

			//if (chrom1str != "chr20" || chrom2str != "chr17")
			//	continue;

			if (start1 > ext)
				start1 -= ext;
			else
				start1 = 1;

			end1 += ext;

			if (start2 > ext)
				start2 -= ext;
			else
				start2 = 1;

			end2 += ext;

			paired_regions.push_back(PairedRegions(chrom1, strand1, start1, end1, chrom2, strand2, start2, end2));

			if (max_dist <  (end1 - start1))
				max_dist = end1 - start1;

			if (max_dist <  (end2 - start2))
				max_dist = end2 - start2;

			++count;
		}
	}
	else
		cout <<"can't open file: " << region_file << endl;

	return count;
}

bool comp1(const PairedRegions* lhs, const PairedRegions* rhs)
{
	if (lhs->start1 == rhs->start1)
		return lhs->end1 < rhs->end1 ;

	return lhs->start1 < rhs->start1;
}

bool comp2(const PairedRegions* lhs, const PairedRegions* rhs)
{
	if (lhs->start2 == rhs->start2)
		return lhs->end2 < rhs->end2;

	return lhs->start2 < rhs->start2;
}

void SortRegion(PAIRED_REGIONS_VEC& paired_regions, HASH_CHROM_PAIRED_REGIONS_PTR_VEC& paired_regions_sorted1, HASH_CHROM_PAIRED_REGIONS_PTR_VEC& paired_regions_sorted2)
{
	PAIRED_REGIONS_ITERATOR paired_regions_iter;

	for (paired_regions_iter = paired_regions.begin(); paired_regions_iter != paired_regions.end(); ++paired_regions_iter)
	{
		paired_regions_sorted1[paired_regions_iter->chrom1].push_back(&(*paired_regions_iter));

		paired_regions_sorted2[paired_regions_iter->chrom2].push_back(&(*paired_regions_iter));
	}

	HASH_CHROM_PAIRED_REGIONS_PTR_VEC_ITERATOR chrom_paired_regions_iter;
	for (chrom_paired_regions_iter = paired_regions_sorted1.begin(); chrom_paired_regions_iter != paired_regions_sorted1.end(); ++chrom_paired_regions_iter)
		sort(chrom_paired_regions_iter->second.begin(), chrom_paired_regions_iter->second.end(), comp1);

	for (chrom_paired_regions_iter = paired_regions_sorted2.begin(); chrom_paired_regions_iter != paired_regions_sorted2.end(); ++chrom_paired_regions_iter)
		sort(chrom_paired_regions_iter->second.begin(), chrom_paired_regions_iter->second.end(), comp2);
}

int MarkMappedReads(const char* mapped_sam_file, hash_map<string, PAIRED_REGIONS_PTR_VEC >& paired_regions_sorted1, hash_map<string, PAIRED_REGIONS_PTR_VEC >& paired_regions_sorted2, hash_set<size_t>& mapped_ids, size_t& max_id)
{
	ifstream ifs(mapped_sam_file);

	int count = 0, matched_count1 = 0, matched_count2 = 0;

	max_id = 0;

	//string prev_chrom1, prev_chrom2;

	//PAIRED_REGIONS_PTR_VEC* paired_regions_ptr1, paired_regions_ptr2;

	//PAIRED_REGIONS_PTR_VEC_ITERATOR st_iter1, st_iter2;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line.empty())
				continue;

			char chromname[1000], readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
			char strand = '+';
			size_t prefixst, strand_t, incorrect, mate_offest, read_index;

			int mate_diff;

			unsigned short mis_match;

			char mate_match;

			sscanf(line.c_str(), "%llu~%s\t%llu\t%s\t%llu\t%llu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu", 
				&read_index, readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, &mate_match, &mate_offest, &mate_diff, chromseq, qualseq, &mis_match);

			//if (prev_chrom1 != chromname)
			//{
			//	prev_chrom1 = chromname;

			//	paired_regions_ptr1 = &(paired_regions_sorted1[chromname]);

			//	st_iter1 = paired_regions_ptr1->begin();
			//}

			//if (prev_chrom2 != chromname)
			//{
			//	prev_chrom2 = chromname;

			//	paired_regions_ptr2 = &(paired_regions_sorted2[chromname]);

			//	st_iter2 = paired_regions_ptr2->begin();
			//}

			//while ((*st_iter1)->end1 < prefixst)
			//	++st_iter1;

			//while ((*st_iter2)->end2 < prefixst)
			//	++st_iter2;

			++count;

			if (count % 100000 == 0)
				cout << matched_count1 <<'\t' <<matched_count2<< endl;

			size_t pair_id = atoi(readname + strlen(readname) - 1);

			size_t comb_offset = read_index;

			if (pair_id == 1)
				++comb_offset;
			else
				--comb_offset;

			//if (max_id < comb_offset)
			//	max_id = comb_offset;

			//mapped_ids.insert(comb_offset);

			//mapped_ids.insert(read_index);

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

					index = splicewaystr.find("M", index) + 1;
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

					if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
						break;

					index = splicewaystr.find("M", index) + 1;
				}
			}

			size_t end = spliceway_vec.back().first + spliceway_vec.back().second;

			PairedRegions cur_region(chromname, strand_t, prefixst, end, chromname, strand_t, prefixst, end);

			PAIRED_REGIONS_PTR_VEC_ITERATOR upper_bound_iter1, upper_bound_iter2, upper_bound_prev_iter1, upper_bound_prev_iter2;

			upper_bound_iter1 = upper_bound(paired_regions_sorted1[chromname].begin(), paired_regions_sorted1[chromname].end(), &cur_region, comp1);

			//upper_bound_prev_iter1 = upper_bound_iter1 - 1;

			while (upper_bound_iter1 != paired_regions_sorted1[chromname].begin())
			{
				upper_bound_prev_iter1 = upper_bound_iter1 - 1;

				if ((*upper_bound_prev_iter1)->end1 + max_dist < cur_region.start1)
					break;

				if ((*upper_bound_prev_iter1)->start1 < cur_region.start1 && (*upper_bound_prev_iter1)->end1 > cur_region.start1) 
				{
					//if ((*upper_bound_prev_iter1)->chrom1 == "chr5" && (*upper_bound_prev_iter1)->start1 == 225649 && (*upper_bound_prev_iter1)->end1 == 226018)
					//{
					//	cout << 1<<endl;
					//	cout << line << endl;
					//	cout << (*upper_bound_prev_iter1)->tostring() << endl;
					//	getchar();
					//}
					if (max_id < comb_offset)
						max_id = comb_offset;

					mapped_ids.insert(comb_offset);

					mapped_ids.insert(read_index);

					(*upper_bound_prev_iter1)->mapped_reads_index.insert(comb_offset);
					++matched_count1;
				}

				upper_bound_iter1 = upper_bound_prev_iter1;
			}

			upper_bound_iter2 = upper_bound(paired_regions_sorted2[chromname].begin(), paired_regions_sorted2[chromname].end(), &cur_region, comp2);

			//upper_bound_prev_iter2 = upper_bound_iter2 - 1;

			while (upper_bound_iter2 != paired_regions_sorted2[chromname].begin())
			{
				upper_bound_prev_iter2 = upper_bound_iter2 - 1;

				if ((*upper_bound_prev_iter2)->end2 + max_dist < cur_region.start1)
					break;

				if ((*upper_bound_prev_iter2)->start2 < cur_region.start1 && (*upper_bound_prev_iter2)->end2 > cur_region.start1)
				{
					//if ((*upper_bound_prev_iter1)->chrom1 == "chr5" && (*upper_bound_prev_iter1)->start1 == 225649 && (*upper_bound_prev_iter1)->end1 == 226018)
					//{
					//	cout << 2<<endl;
					//	cout << line << endl;
					//	cout << (*upper_bound_prev_iter2)->tostring() << endl;
					//	getchar();
					//}

					if (max_id < comb_offset)
						max_id = comb_offset;

					mapped_ids.insert(comb_offset);

					mapped_ids.insert(read_index);

					(*upper_bound_prev_iter2)->mapped_reads_index.insert(comb_offset);
					++matched_count2;
				}

				upper_bound_iter2 = upper_bound_prev_iter2;
			}
		}
	}
	else
		cout <<"can't open file: " << mapped_sam_file << endl;

	return count;
}

//int RemBothMapped(PAIRED_REGIONS_VEC& paired_regions, const char* output_file)
//{
//	ofstream ofs(output_file);
//
//	PAIRED_REGIONS_VEC::iterator paired_regions_iter;
//
//	vector<size_t>::iterator mapped_reads_index_iter1, mapped_reads_index_iter2;
//
//	hash_map<size_t, size_t> mapped_reads_index_map1;
//
//	for (paired_regions_iter = paired_regions.begin(); paired_regions_iter != paired_regions.end(); ++paired_regions_iter)
//	{
//		ofs << paired_regions_iter->chrom1 <<'\t'<< paired_regions_iter->start1 <<'\t'<<  paired_regions_iter->end1 
//			<< '\t'<<  paired_regions_iter->chrom2 << '\t'<< paired_regions_iter->start2 <<'\t'<<  paired_regions_iter->end2 << '\t';
//
//		for (mapped_reads_index_iter1 = paired_regions_iter->mapped_reads_index1.begin(); mapped_reads_index_iter1 != paired_regions_iter->mapped_reads_index1.end(); ++mapped_reads_index_iter1)
//		{
//			ofs << *mapped_reads_index_iter1<<';';			
//		}
//
//		ofs << '\t';
//
//		for (mapped_reads_index_iter2 = paired_regions_iter->mapped_reads_index2.begin(); mapped_reads_index_iter2 != paired_regions_iter->mapped_reads_index2.end(); ++mapped_reads_index_iter2)
//		{
//			ofs << *mapped_reads_index_iter2<<';';
//		}
//
//		ofs << endl;
//	}
//
//	return 0;
//}

int countall=0;

size_t ReadReadsFile(const char* reads_file, vector<string >& read_vec, hash_set<size_t>& mapped_reads_ids, int fastq, string extracted_reads_file)
{
	ifstream input_fs(reads_file);

	ofstream extracted_reads_ofs(extracted_reads_file.c_str());

	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}    

	size_t count = 0;

	string readidline;
	string line;
	string readidline2;
	string qualstr;
	//int countall=0;
	while(!input_fs.eof())
	{
		getline(input_fs, readidline);
		getline(input_fs, line);

		if (fastq)
		{
			getline(input_fs, readidline2);
			getline(input_fs, qualstr);
		}

		if (line == "" || readidline == "")
			continue;

		++count;

		size_t combine_index;

		if (mapped_reads_ids.find(count) != mapped_reads_ids.end())
		{
			if (read_vec.size() < count)
			{
				cout << "read index exceed" << endl;
				cout << readidline << endl;
				return countall;
			}

			string combinedline;
			
			combinedline.append(readidline);combinedline.append("\n");
			combinedline.append(line);combinedline.append("\n");

			if (fastq)
			{
				combinedline.append(readidline2);combinedline.append("\n");
				combinedline.append(qualstr);combinedline.append("\n");
			}

			extracted_reads_ofs << combinedline;

			//read_vec[count - 1] = combinedline;			
		}

		countall++;
	}
	input_fs.close();
	
	cout<<countall<<" reads converted for "<<reads_file<<endl; 

	return countall;
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

			if (line.length() == 0)
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";
}

int ReadRegions(PAIRED_REGIONS_VEC& paired_regions, vector<string>& reserved_reads_space, string out_path, string chrom_dir, string director_file)
{
	int count = 0;

	string chromseq1 = "";
	string chromseq2 = "";
	string prevchrom1 = "";
	string prevchrom2 = "";

	PAIRED_REGIONS_VEC::iterator paired_regions_iter;

	int region_count = 0;

	string sub_dir;

	ofstream director_file_ofs(director_file.c_str());

	ofstream combined_ofs;//(combined_sequence.c_str());

	for (paired_regions_iter = paired_regions.begin(); paired_regions_iter != paired_regions.end(); ++paired_regions_iter)
	{
		++region_count;

		string chrom1 = paired_regions_iter->chrom1;

		string chrom2 = paired_regions_iter->chrom2;

		size_t start1 = paired_regions_iter->start1;

		size_t start2 = paired_regions_iter->start2;

		size_t end1 = paired_regions_iter->end1;

		size_t end2 = paired_regions_iter->end2;

		size_t strand1 = paired_regions_iter->strand1;

		size_t strand2 = paired_regions_iter->strand2;

		if (prevchrom1 != chrom1)
		{
			string chrom_file = chrom_dir;

			chrom_file.append(chrom1);

			chrom_file.append(".fa");

			chromseq1.clear();

			readchrom(chrom_file.c_str(), chromseq1);

			prevchrom1 = chrom1;
		}

		if (prevchrom2 != chrom2)
		{
			string chrom_file = chrom_dir;

			chrom_file.append(chrom2);

			chrom_file.append(".fa");

			chromseq2.clear();

			readchrom(chrom_file.c_str(), chromseq2);

			prevchrom2 = chrom2;
		}

		if (region_count % 200 == 1)
		{
			char sub_dir_chr[1000];

			size_t region_count_10000 = region_count + 199;

			sprintf(sub_dir_chr, "%llu~%llu", region_count, region_count_10000);

			sub_dir = out_path + sub_dir_chr;

			director_file_ofs << sub_dir << '\t'<<out_path<<endl;

			sub_dir.append("/");

			string mkdir_cmd2 = "mkdir ";

			mkdir_cmd2.append(sub_dir);

			system(mkdir_cmd2.c_str());

			string combined_sequence_str = sub_dir + "combined_sequence";

			if (combined_ofs.is_open())
				combined_ofs.close();

			combined_ofs.open(combined_sequence_str.c_str());
		}

		if (end1 > chromseq1.length())
			end1 = chromseq1.length();

		if (end2 > chromseq2.length())
			end2 = chromseq2.length(); 

		char path_chr[1000];

		sprintf(path_chr, "%s~%llu~%llu~%llu~%s~%llu~%llu~%llu", chrom1.c_str(), strand1, start1, end1, chrom2.c_str(), strand2, start2, end2);

		string cur_path = sub_dir;// + path_chr;

		//director_file_ofs << cur_path << '\t'<<path_chr << endl;

		//string mkdir_cmd = "mkdir ";

		//mkdir_cmd.append(cur_path);

		//system(mkdir_cmd.c_str());

		string cur_chrom_file = cur_path + "/";

		cur_chrom_file.append(path_chr); cur_chrom_file.append(".fa");

		ofstream cur_chrom_file_ofs(cur_chrom_file.c_str());

		cur_chrom_file_ofs << '>' <<path_chr<<endl;

		combined_ofs <<'>' <<path_chr<<endl;

		string comb_chrom_seq, chrom1_seq, chrom2_seq;
		
		if (strand1 == 0)
			chrom1_seq = chromseq1.substr(start1, end1 - start1);
		else
			chrom1_seq = revcomp(chromseq1.substr(start1, end1 - start1));
		
		if (strand2 == 0)
			chrom2_seq = chromseq2.substr(start2, end2 - start2);
		else
			chrom2_seq = revcomp(chromseq2.substr(start2, end2 - start2));

		comb_chrom_seq = chrom1_seq + chrom2_seq;

		for (size_t i = 0; i < comb_chrom_seq.length(); i = i + 60)
		{
			cur_chrom_file_ofs << comb_chrom_seq.substr(i, 60)<<endl;
			combined_ofs << comb_chrom_seq.substr(i, 60)<<endl;
		}		

		//string cur_reads_file = cur_path + "/";

		//cur_reads_file.append(path_chr); cur_reads_file.append(".fq");

		//ofstream cur_reads_file_ofs(cur_reads_file.c_str());

		////set<size_t> mapped_reads_ids;

		//hash_set<size_t>::iterator mapped_reads_iter;

		////for (mapped_reads_iter = paired_regions_iter->mapped_reads_index.begin(); mapped_reads_iter != paired_regions_iter->mapped_reads_index.end(); ++mapped_reads_iter)
		////	mapped_reads_ids.insert(*mapped_reads_iter);

		////set<size_t>::iterator read_index_iter;

		//for (mapped_reads_iter =  paired_regions_iter->mapped_reads_index.begin(); mapped_reads_iter !=  paired_regions_iter->mapped_reads_index.end(); ++mapped_reads_iter)
		//{
		//	if (reserved_reads_space[(*mapped_reads_iter) - 1].empty())
		//	{
		//		cout << "read not stored: "<<(*mapped_reads_iter) <<endl;
		//		continue;
		//	}

		//	cur_reads_file_ofs << reserved_reads_space[(*mapped_reads_iter) - 1];
		//}
	}

	return count;
}

int main(int argc, char* argv[])
{
	const char* region_file = argv[1];

	const char* mapped_sam_file = argv[2];

	const char* reads_file = argv[3];

	string out_path = argv[4];

	int fastq = atoi(argv[5]);

	string chrom_dir = argv[6];

	string director_file = argv[7];

	size_t ext = atoi(argv[8]);

	vector<PairedRegions> paired_regions;

	cout << "read regions" << endl;

	ReadRegions(region_file, paired_regions, ext);

	cout << max_dist<<endl;

	HASH_CHROM_PAIRED_REGIONS_PTR_VEC paired_regions_sorted1;

	HASH_CHROM_PAIRED_REGIONS_PTR_VEC paired_regions_sorted2;

	cout << "sort regions" << endl;

	SortRegion(paired_regions, paired_regions_sorted1, paired_regions_sorted2);

	cout << "mark regions" << endl;

	hash_set<size_t> mapped_reads_ids;

	size_t max_id;

	MarkMappedReads(mapped_sam_file, paired_regions_sorted1, paired_regions_sorted2, mapped_reads_ids, max_id);

	vector<string> reserved_reads_space(max_id);

	cout << "read reads files"<<endl;

	ReadReadsFile(reads_file, reserved_reads_space, mapped_reads_ids, fastq, out_path + "extracted_reads_file");

	cout << "read regions"<<endl;

	//python /homes/kai/code/bin/MapSplice_1.5.0/MapSplice_1.13/mapsplice_segments.py --max-hits 10 -o /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/output -X 4 -n 10 -w 50 -S fa -Q fq -c /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/ -u /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/chr2~99735158~99820411~chr2~99989586~99989682.fq -B /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/index -L 25 -I 1 -X 4 --full-running --not-rem-temp --fusion --not-rerun-all 2>/sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/log
	ReadRegions(paired_regions, reserved_reads_space, out_path, chrom_dir, director_file);

	//RemBothMapped(paired_regions, output_file);

}