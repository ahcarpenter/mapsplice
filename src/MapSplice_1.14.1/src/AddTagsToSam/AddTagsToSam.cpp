//#define VS

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>

#ifdef VS
#include <hash_map> //vc only
#include <hash_set> //vc only
#else
#include <ext/hash_map> //g++ only
#include <ext/hash_set> //g++ only
#endif

#ifdef VS
using namespace stdext;
#endif

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;
#endif

using namespace std;

size_t count_del = 0;

struct SamRec {

	string tag_name;

	size_t strand_t;

	string chrom_name;

	size_t start, end;

	unsigned short confid;

	string splice_way;

	string mapped_seq;

	unsigned short mis_match;

	vector<pair<size_t, int> > spliceway_vec;

	//vector<pair<size_t, int> > spliceway_vec;
	//bool wrong_format;

	size_t tagidx;

	//double filter_score;

	//size_t ave_intron_len;

	//double ave_junc_mis;

	//vector<string> junc_id;

	//int best;

	//string alters;

	string qual_str;

	//bool isunique, isexonic, isspliced, issmallins, issmalldel;
	
	//size_t canon_count, noncanon_count;

	//double canon_rate;

	int matched_id;

	char mate_match;

	size_t mate_offest;

	int mate_diff;

	string cur_line;

	size_t ih, hi;

	size_t mappedlen;

	//bool is_insert;

	//size_t mappedlen;

	//bool is_spliced;

	SamRec() {}

	SamRec(const string& tname, size_t strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& qualstr, 
		char matematch, size_t mateoffest, int matediff, int min_del, string line = "") : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
		mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offest(mateoffest), mate_diff(matediff), tagidx(tidx), qual_str(qualstr), cur_line(line), mappedlen(0)
	{
		size_t index = 0;

		string flag_str = " ";

		bool isdel = false;

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

			if (flag == 'N' && maplen <= min_del)
			{
				splice_way[index-1] = 'D';
				isdel = true;
			}

		}

		if (isdel)
			++count_del;
	}

	//SamRec(const string& tname, size_t strand, const string& cname, size_t st, size_t conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
	//	char matematch, size_t mateoffest, int matediff, string line = "") : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
	//	mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offest(mateoffest), mate_diff(matediff), wrong_format(false), tagidx(tidx), filter_score(0), best(0), alters(alt), qual_str(qualstr), 
	//	isunique(false), isexonic(false), isspliced(false), ave_intron_len(0), ave_junc_mis(-1), canon_count(0), noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), cur_line(line),
	//	is_insert(false)/*, is_spliced(false)*/, mappedlen(0)
	//{
	//	//size_t index = 0;

	//	//if ((spliceway.find("I") != string::npos))
	//	//	issmallins = true;
	//	//else if ((spliceway.find("D") != string::npos))
	//	//	issmalldel = true;
	//	//else if (spliceway.find("N") != string::npos)
	//	//	isspliced = true;
	//	//else
	//	//	isexonic = true;		

	//	//string flag_str = " ";

	//	//isexonic = true;

	//	//isspliced = false;

	//	//while (true)
	//	//{
	//	//	if (index >= splice_way.length())
	//	//		break;

	//	//	int maplen;

	//	//	char flag;

	//	//	sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

	//	//	if (flag_str[0] == ' ')
	//	//	{
	//	//		if (flag == 'I')
	//	//			spliceway_vec.push_back(make_pair(start, -maplen));
	//	//		else if (flag == 'M')
	//	//		{
	//	//			spliceway_vec.push_back(make_pair(start, maplen));

	//	//			mappedlen += maplen;
	//	//		}
	//	//		else if (flag == 'N')
	//	//		{
	//	//			cout<<"start with N?"<<endl;
	//	//			spliceway_vec.push_back(make_pair(start + maplen, 0));
	//	//		}
	//	//	}
	//	//	else if (flag_str[0] == 'M')
	//	//	{
	//	//		if (flag == 'I')
	//	//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
	//	//		else if (flag == 'M')
	//	//		{
	//	//			cout << "continue Ms?"<<endl;
	//	//			spliceway_vec.back().second += maplen;

	//	//			mappedlen += maplen;
	//	//		}
	//	//		else if (flag == 'N')
	//	//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
	//	//	}
	//	//	else if (flag_str[0] == 'N')
	//	//	{
	//	//		if (flag == 'I')
	//	//			spliceway_vec.back().second = -maplen;
	//	//		else if (flag == 'M')
	//	//		{
	//	//			spliceway_vec.back().second = maplen;

	//	//			mappedlen += maplen;
	//	//		}
	//	//		else if (flag == 'N')
	//	//		{
	//	//			cout << "continue Ns?"<<endl;
	//	//			spliceway_vec.back().first += maplen;
	//	//		}
	//	//	}
	//	//	else if (flag_str[0] == 'I')
	//	//	{
	//	//		if (flag == 'I')
	//	//		{
	//	//			cout << "continue Is?"<<endl;
	//	//			spliceway_vec.back().second += -maplen;
	//	//		}
	//	//		else if (flag == 'M')
	//	//		{
	//	//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));

	//	//			mappedlen += maplen;
	//	//		}
	//	//		else if (flag == 'N')
	//	//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
	//	//	}

	//	//	flag_str[0] = flag;

	//	//	index = splice_way.find(flag_str, index) + 1;

	//	//}

	//	//end = spliceway_vec.back().first + spliceway_vec.back().second - 1;
	//}

	string tostring()
	{
		if (!cur_line.empty())
			return cur_line;

		char sam_rec_char[5000];

		sprintf(sam_rec_char, "%s\t%llu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu", tag_name.c_str(), strand_t, 
			chrom_name.c_str(), start, confid, splice_way.c_str(), mate_match, mate_offest, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, ih, hi);

		return sam_rec_char;
	}
};

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

void
AddTagsToSam(const char* infile, const char* added_tags_file, const char* stat_file, size_t ispaired, hash_set<size_t>& reads_ids, int min_del)
{
	ifstream infile_ifs(infile);

	string forward = added_tags_file; forward.append(".forward");

	ofstream forward_ofs(forward.c_str());

	string reverse = added_tags_file; reverse.append(".reverse");

	ofstream reverse_ofs(reverse.c_str());

	ofstream added_tags_file_ofs(added_tags_file);

	size_t prev_tag_idx = -1;

	//map<string, int> mapped_reads;

	vector<SamRec> mapped_reads;

	size_t count = 0, multiple_mapped = 0;

	if (infile_ifs.is_open())
	{
		string line;
		while (getline(infile_ifs,line))
		{
			if (line == "")
				continue;

			char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000];
			unsigned short strand, confd, mismatch;
			size_t offset;

			char mate_match;

			size_t mate_offest, mate_diff, tagidx;

			//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
			//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

			int read_count = sscanf(line.c_str(), "%llu~%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu", &tagidx, tagname, &strand, chrom, &offset, &confd, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch);

			//string tagnamestr = tagname;

			//size_t last_idx = tagnamestr.find_first_of("~");

			//size_t tagidx = -1;

			//if (last_idx != string::npos)
			//{
			//	string tagidxstr = tagnamestr.substr(0, last_idx);
			//	tagidx = atoi(tagidxstr.c_str()) - 1;

			//	reads_ids.insert(tagidx);
			//}
			//else
			//{
			//	cout << "read not indexed? "<<endl;
			//	continue;
			//}

			if (prev_tag_idx == -1 || prev_tag_idx == tagidx)
			{
				mapped_reads.push_back(SamRec(tagname, strand, chrom, offset, confd, mapped, seq, mismatch, tagidx, qual_str, mate_match, mate_offest, mate_diff, min_del));

				if (prev_tag_idx == -1)
					prev_tag_idx = tagidx;
			}
			else
			{
				vector<SamRec>::iterator msi_iter;

				++count;

				if (mapped_reads.size() > 1)
					++multiple_mapped;

				for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
				{
					msi_iter->ih = mapped_reads.size();

					msi_iter->hi = msi_iter - mapped_reads.begin() + 1;

					msi_iter->strand_t |= ispaired;

					if (msi_iter->mate_match == '=')
						msi_iter->strand_t |= IS_PAIRED_MAPPED;

					char end_id = msi_iter->tag_name[msi_iter->tag_name.length() - 1];

					if (ispaired)
					{
						if (end_id == '1')
							msi_iter->strand_t |= IS_FIRST_END;
						else if (end_id == '2')
							msi_iter->strand_t |= IS_SECOND_END;
						else
						{
							cout <<"read has no end id ? "<<endl;
						}
					}

					if (msi_iter->strand_t & IS_REVERSE)
						reverse_ofs << msi_iter->tostring() << endl;
					else
						forward_ofs << msi_iter->tostring() << endl;

					//added_tags_file_ofs << msi_iter->tostring() << endl;


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

				}

				reads_ids.insert(tagidx);

				mapped_reads.clear();

				mapped_reads.push_back(SamRec(tagname, strand, chrom, offset, confd, mapped, seq, mismatch, tagidx, qual_str, mate_match, mate_offest, mate_diff, min_del));

				prev_tag_idx = tagidx;
			}			
		}

		infile_ifs.close();
	}
	else
	{
		cout << "can't open file "<< infile<<endl; exit(1);
	}

	vector<SamRec>::iterator msi_iter;

	++count;

	if (mapped_reads.size() > 1)
		++multiple_mapped;

	if (prev_tag_idx != -1)
		reads_ids.insert(prev_tag_idx);

	for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
	{
		msi_iter->ih = mapped_reads.size();

		msi_iter->hi = msi_iter - mapped_reads.begin() + 1;

		msi_iter->strand_t |= ispaired;

		if (msi_iter->mate_match == '=')
			msi_iter->strand_t |= IS_PAIRED_MAPPED;

		char end_id = msi_iter->tag_name[msi_iter->tag_name.length() - 1];

		if (ispaired)
		{
			if (end_id == '1')
				msi_iter->strand_t |= IS_FIRST_END;
			else if (end_id == '2')
				msi_iter->strand_t |= IS_SECOND_END;
			else
			{
				cout <<"read has no end id ? "<<endl;
			}
		}

		if (msi_iter->strand_t & IS_REVERSE)
			reverse_ofs << msi_iter->tostring() << endl;
		else
			forward_ofs << msi_iter->tostring() << endl;

		//added_tags_file_ofs << msi_iter->tostring() << endl;
	}

	ofstream stat_ofs(stat_file, ios_base::app);

	stat_ofs << "total_reads: " << count<<endl<<"multiple mapped: "<< multiple_mapped <<endl;
}

size_t AddUnmappedReadsFQ(const hash_set<size_t>& reads_ids, string reads_file, string unmapped_reads_sam)
{
	ifstream ifs(reads_file.c_str());

	ofstream ofs(unmapped_reads_sam.c_str());

	size_t read_count = 0, count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string reads_name, reads_seq, reads_qual, strand;

			getline(ifs, reads_name);

			if (reads_name.empty())
				continue;

			getline(ifs, reads_seq);

			if (reads_seq.empty())
				continue;

			getline(ifs, reads_qual);

			if (reads_seq.empty())
				continue;

			getline(ifs, strand);

			if (reads_seq.empty())
				continue;

			++read_count;

			if (reads_ids.find(read_count) == reads_ids.end()) // unmapped
			{
				reads_name = reads_name.substr(1);

				size_t last_idx = reads_name.find_first_of("~");

				if (last_idx != string::npos)
				{
					reads_name = reads_name.substr(last_idx + 1);
				}

				ofs << reads_name <<"\t4\t*\t0\t0\t*\t*\t0\t0\t"<<reads_seq <<"\t"<<reads_qual <<"\tIH:i:0\tHI:i:0" <<endl;

				++count;
			}
		}
	}

	return count;
}

size_t AddUnmappedReadsFA(const hash_set<size_t>& reads_ids, string reads_file, string unmapped_reads_sam)
{
	ifstream ifs(reads_file.c_str());

	ofstream ofs(unmapped_reads_sam.c_str());

	size_t read_count = 0, count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string reads_name, reads_seq, reads_qual, strand;

			getline(ifs, reads_name);

			if (reads_name.empty())
				continue;

			getline(ifs, reads_seq);

			if (reads_seq.empty())
				continue;

			if (reads_seq.empty())
				continue;

			++read_count;

			if (reads_ids.find(read_count) == reads_ids.end()) // unmapped
			{
				reads_name = reads_name.substr(1);

				size_t last_idx = reads_name.find_first_of("~");

				if (last_idx != string::npos)
				{
					reads_name = reads_name.substr(last_idx + 1);
				}

				ofs << reads_name <<"\t4\t*\t0\t0\t*\t*\t0\t0\t"<<reads_seq <<"\t"<<string('I', reads_seq.length()) <<"\tIH:i:0\tHI:i:0" <<endl;

				++count;
			}
		}
	}

	return count;
}

int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << "infile added_tags_file stat_file 0|1(is_paired) reads_file unmapped_samfile fq|fa" <<endl;
		exit(0);
	}
	const char* infile = argv[1];

	const char* added_tags_file = argv[2];

	const char* stat_file = argv[3];

	size_t ispaired = atoi(argv[4]);

	string reads_file = argv[5];

	string unmapped_sam_file = argv[6];

	string reads_format = argv[7];

	int min_del = atoi(argv[8]);

	bool isfq = (reads_format == "fq" ? true : false);

	hash_set<size_t> reads_ids;

	cout << "add tag to sam "<<endl;

	AddTagsToSam(infile, added_tags_file, stat_file, ispaired, reads_ids, min_del);

	size_t unmapped_reads;

	cout << "output unmapped sam "<<endl;

	if (isfq)
		unmapped_reads = AddUnmappedReadsFQ(reads_ids, reads_file, unmapped_sam_file);
	else
		unmapped_reads = AddUnmappedReadsFA(reads_ids, reads_file, unmapped_sam_file);

	cout << count_del <<" deletions alignments"<<endl;

}