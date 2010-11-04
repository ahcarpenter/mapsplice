#ifndef SAM_INFO_H
#define SAM_INFO_H

#include <string.h>
using namespace std;
#pragma warning(disable:4996)

struct Sam_Info
{
public:
	Sam_Info()
	{
		start_pos=-1;
	}

	~Sam_Info()
	{
	}

	Sam_Info(const Sam_Info& my_sam_info)
	{
		start_pos=my_sam_info.start_pos;
		for(size_t i=0;i<my_sam_info.mismatchs.size();i++)
		{
			mismatchs.push_back(my_sam_info.mismatchs[i]);
		}
		for(size_t i=0;i<my_sam_info.jump_code.size();i++)
		{
			jump_code.push_back(my_sam_info.jump_code[i]);
		}
		start_contig_index=my_sam_info.start_contig_index;
		start_seg_no=my_sam_info.start_seg_no;
		end_seg_no=my_sam_info.end_seg_no;
		chrom=my_sam_info.chrom;
		pair_index=-1;
	}

	Sam_Info(Bwtmap_Info& my_bwt_info, size_t current_index)
	{
		start_pos=my_bwt_info.start + 1;
		for(size_t i=0;i<my_bwt_info.mismatches.size();i++)
		{
			mismatchs.push_back(my_bwt_info.mismatches[i]);
		}
		int match_len=my_bwt_info.end-my_bwt_info.start+1;
		jump_code.push_back(match_len);
		start_contig_index=current_index;
		start_seg_no=my_bwt_info.start_seg_no;
		end_seg_no=my_bwt_info.end_seg_no;
		chrom=my_bwt_info.chrom;
		pair_index=-1;
	}

	string chrom;
	int start_pos;
	vector<Mismatch> mismatchs;
	vector<int> jump_code;
	size_t start_contig_index;
	size_t start_seg_no;
	size_t end_seg_no;
	int pair_index;
};

#endif