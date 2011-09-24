#ifndef BWTMAP_INFO_H
#define BWTMAP_INFO_H

#include <string.h>
using namespace std;
#pragma warning(disable:4996)

struct Bwtmap_Info
{
public:
	Bwtmap_Info()
	{ }

	Bwtmap_Info(string line, bool pairend)
	{ 
		if(!pairend)     //non pair end read
			set_value_single(line);
		else   //pair end read
			set_value_pairend(line);
	}

	void set_value_single(string line)
	{
		char read_id_tmp[MAX_READID_LENGTH], seg_no_tmp[MAX_SEGNO_LENGTH], strand_tmp[MAX_STRAND_LENGTH], chrom_tmp[MAX_CHROM_LENGTH], sequence_tmp[MAX_READ_LENGTH], quality_score_tmp[MAX_READ_LENGTH], mismatch_info_tmp[MAX_MISMATCH_LENGTH], unmapped_index_tmp[MAX_UNMAPPED_INDEX];
		int read_count=sscanf(line.c_str(), "%[^~]%*[~]%[^~]%*[~]%s \t    %s\t        %s\t     %d\t       %s\t             %s\t     %d\t     %s    ",
			seg_no_tmp, unmapped_index_tmp, read_id_tmp, strand_tmp, chrom_tmp, &start, sequence_tmp, quality_score_tmp, &other_map   ,mismatch_info_tmp);	
		read_id=read_id_tmp;
		unmapped_index=atoi(unmapped_index_tmp);
		start_seg_no=end_seg_no=atoi(seg_no_tmp);
		strand=strand_tmp;
		chrom=chrom_tmp;
		sequence=sequence_tmp;
		quality_score=quality_score_tmp;
		end=start + (int)sequence.length() -1;
		head=-1;
		tail=-1;
		if(read_count==10)// has mismatch
		{
			char* p;
			p = strtok(mismatch_info_tmp ,mismatch_spliter); 
			while(p != NULL) 
			{ 
				Mismatch new_mismatch;
				new_mismatch.seg_no=start_seg_no;
				new_mismatch.pos=(atoi(p));
				if(strand=="-")
					new_mismatch.pos=(int)(sequence.length()) - new_mismatch.pos - 1;
				p = strtok(NULL , mismatch_spliter); 
				new_mismatch.src_char=(*p);
				new_mismatch.dst_char=*(p+2);
				mismatches.push_back(new_mismatch);
				p = strtok(NULL , mismatch_spliter); 
			} 
		}
	}

	void set_value_pairend(string line)
	{
		char read_id_tmp[MAX_READID_LENGTH], seg_no_tmp[MAX_SEGNO_LENGTH], strand_tmp[MAX_STRAND_LENGTH], chrom_tmp[MAX_CHROM_LENGTH], sequence_tmp[MAX_READ_LENGTH], quality_score_tmp[MAX_READ_LENGTH], mismatch_info_tmp[MAX_MISMATCH_LENGTH], unmapped_index_tmp[MAX_UNMAPPED_INDEX];
		int read_count=sscanf(line.c_str(), "%[^~]%*[~]%[^~]%*[~]%[^/]%*[/]%d \t    %s\t        %s\t     %d\t       %s\t             %s\t     %d\t     %s    ",
			seg_no_tmp, unmapped_index_tmp, read_id_tmp,&pair_no  ,strand_tmp, chrom_tmp, &start, sequence_tmp, quality_score_tmp, &other_map   ,mismatch_info_tmp);	
		read_id=read_id_tmp;
		unmapped_index=atoi(unmapped_index_tmp);
		start_seg_no=end_seg_no=atoi(seg_no_tmp);
		strand=strand_tmp;
		chrom=chrom_tmp;
		sequence=sequence_tmp;
		quality_score=quality_score_tmp;
		end=start + (int)sequence.length() -1;
		head=-1;
		tail=-1;
		if(read_count==11)// has mismatch
		{
			char* p;
			p = strtok(mismatch_info_tmp ,mismatch_spliter); 
			while(p != NULL) 
			{ 
				Mismatch new_mismatch;
				new_mismatch.seg_no=start_seg_no;
				new_mismatch.pos=(atoi(p));
				if(strand=="-")
					new_mismatch.pos=(int)(sequence.length()) - new_mismatch.pos - 1;
				p = strtok(NULL , mismatch_spliter); 
				new_mismatch.src_char=(*p);
				new_mismatch.dst_char=*(p+2);
				mismatches.push_back(new_mismatch);
				p = strtok(NULL , mismatch_spliter); 
			} 
		}	
	}

	~Bwtmap_Info()
	{
	}

	bool check_head(int check_seg_no)
	{
		if(head==-1)
			return true;
		if(strand=="+" && check_seg_no>=head)
			return true;
		else if(strand=="-" && check_seg_no<=head)
			return true;
		return false;
	}

	bool check_tail(int check_seg_no)
	{
		if(tail==-1)
			return true;
		if(strand=="+" && check_seg_no<=tail)
			return true;
		else if(strand=="-" && check_seg_no>=tail)
			return true;
		return false;
	}

	string read_id;
	int unmapped_index;
	int start_seg_no;
	int end_seg_no;
	string strand;
	string chrom;
	int start;
	int end;
	string sequence;
	string quality_score;
	int other_map;
	vector<Mismatch> mismatches;
	int head;
	int tail;
	vector<Splicemap_Info> splice_internal;
	vector<Splicemap_Info> splice_head;
	vector<Splicemap_Info> splice_tail;
	int pair_no;
};

#endif