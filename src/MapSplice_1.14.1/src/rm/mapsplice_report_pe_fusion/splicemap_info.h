#ifndef SPLICEMAP_INFO_H
#define SPLICEMAP_INFO_H

#include <string.h>
using namespace std;
#pragma warning(disable:4996)


struct Splicemap_Info
{
public:
	Splicemap_Info()
	{ }

	Splicemap_Info(string line, bool pairend)
	{
		if(!pairend)
			set_value_single(line);
		else
			set_value_pairend(line);
	}

	void set_value(string line, bool pairend)
	{
		if(!pairend)
			set_value_single(line);
		else
			set_value_pairend(line);
	}

	void set_value_single(string line)
	{
		char read_id_tmp[MAX_READID_LENGTH], chrom_tmp[MAX_CHROM_LENGTH], contig_name_tmp[MAX_READ_LENGTH],jump_code_tmp[MAX_MISMATCH_LENGTH], sequence_tmp[MAX_READ_LENGTH],mismatch_info_tmp[MAX_MISMATCH_LENGTH], unmapped_index_tmp[MAX_UNMAPPED_INDEX];
		int read_count=sscanf(line.c_str(), "%[^~]%*[~]%s\t         %d\t      %s\t         %d\t        %s\t             %d\t          %s\t          %s\t            %d\t          %d\t      %*s\t    %s  ",
			unmapped_index_tmp, read_id_tmp, &strand, chrom_tmp, &start_pos,  contig_name_tmp,  &left_length, jump_code_tmp, sequence_tmp,       &start_contig, &end_contig,  mismatch_info_tmp);	
		read_id=read_id_tmp;
		unmapped_index=atoi(unmapped_index_tmp);
		chrom=chrom_tmp;
		contig_name=contig_name_tmp;
		sequence=sequence_tmp;
		char* p;
		p = strtok(jump_code_tmp ,jumpcode_spliter); 
		while(p != NULL)          //skip offset
		{ 
		jump_code.push_back(atoi(p));
		p = strtok(NULL , jumpcode_spliter); 
		}
		if(read_count==12)// has mismatch
		{
			char* q;
			q = strtok(mismatch_info_tmp ,mismatch_spliter); 
			while(q != NULL) 
			{ 
				Mismatch new_mismatch;
				new_mismatch.seg_no=1;
				new_mismatch.pos=(atoi(q));
				q = strtok(NULL , mismatch_spliter); 
				new_mismatch.src_char=(*q);
				new_mismatch.dst_char=*(q+2);
				mismatches.push_back(new_mismatch);
				q = strtok(NULL , mismatch_spliter); 
			} 
		}
	}

	void set_value_pairend(string line)
	{
		char read_id_tmp[MAX_READID_LENGTH], chrom_tmp[MAX_CHROM_LENGTH], contig_name_tmp[MAX_READ_LENGTH],jump_code_tmp[MAX_MISMATCH_LENGTH], sequence_tmp[MAX_READ_LENGTH],mismatch_info_tmp[MAX_MISMATCH_LENGTH], unmapped_index_tmp[MAX_UNMAPPED_INDEX];
		int read_count=sscanf(line.c_str(), "%[^~]%*[~]%[^/]%*[/]%d         %d\t      %s\t         %d\t        %s\t             %d\t          %s\t          %s\t            %d\t          %d\t      %*s\t    %s  ",
			unmapped_index_tmp, read_id_tmp, &pair_no,  &strand, chrom_tmp, &start_pos,  contig_name_tmp,  &left_length, jump_code_tmp, sequence_tmp,       &start_contig, &end_contig,  mismatch_info_tmp);	
		read_id=read_id_tmp;
		unmapped_index=atoi(unmapped_index_tmp);
		chrom=chrom_tmp;
		contig_name=contig_name_tmp;
		sequence=sequence_tmp;
		char* p;
		p = strtok(jump_code_tmp ,jumpcode_spliter); 
		while(p != NULL)          //skip offset
		{ 
			jump_code.push_back(atoi(p));
			p = strtok(NULL , jumpcode_spliter); 
		}
		if(read_count==13)// has mismatch
		{
			char* q;
			q = strtok(mismatch_info_tmp ,mismatch_spliter); 
			while(q != NULL) 
			{ 
				Mismatch new_mismatch;
				new_mismatch.seg_no=1;
				new_mismatch.pos=(atoi(q));
				q = strtok(NULL , mismatch_spliter); 
				new_mismatch.src_char=(*q);
				new_mismatch.dst_char=*(q+2);
				mismatches.push_back(new_mismatch);
				q = strtok(NULL , mismatch_spliter); 
			} 
		}
	}

	~Splicemap_Info()
	{
	}

	string read_id;
	int unmapped_index;
	int strand;
	string chrom;
	int start_pos;
	string contig_name;
	int left_length;
	vector<int> jump_code;
	vector<Mismatch> mismatches;
	string sequence;
	int start_contig;
	int end_contig;
	int pair_no;
};


#endif