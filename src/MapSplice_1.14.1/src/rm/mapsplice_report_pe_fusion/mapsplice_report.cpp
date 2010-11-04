#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include "common.h"
#include "mismatch_info.h"
#include "quality_score.h"
#include "splicemap_info.h"
#include "splicemap_file.h"
#include "bwtmap_info.h"
#include "bwtmap_file.h"
#include "sam_info.h"
#include "reads_file.h"

inline bool bwtmap_compare(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort bwtmap compare function
{
	if(seg1.chrom!=seg2.chrom)
		return seg1.chrom<seg2.chrom;
	if(seg1.strand!=seg2.strand)
		return seg1.strand<seg2.strand;
	if (seg1.start != seg2.start)
		return seg1.start < seg2.start;
	else
		return seg1.end <seg2.end;
}

inline bool bwtmap_compare_stable(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort bwtmap compare function stable version
{
	if(seg1.chrom!=seg2.chrom)
		return seg1.chrom<seg2.chrom;
	if(seg1.strand!=seg2.strand)
		return seg1.strand<seg2.strand;
	if (seg1.start != seg2.start)
		return seg1.start < seg2.start;
	if(seg1.end != seg2.end)
		return seg1.end <seg2.end;
	else
		return seg1.start_seg_no<seg2.start_seg_no;
}

inline bool bwtmap_compare_pe_stable(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort compare function
{
	if(seg1.pair_no!=seg2.pair_no)
		return seg1.pair_no < seg2.pair_no;
	if(seg1.chrom!=seg2.chrom)
		return seg1.chrom<seg2.chrom;
	if(seg1.strand!=seg2.strand)
		return seg1.strand<seg2.strand;
	if (seg1.start != seg2.start)
		return seg1.start < seg2.start;
	if(seg1.end != seg2.end)
		return seg1.end <seg2.end;
	else
		return seg1.start_seg_no<seg2.start_seg_no;
}

inline bool mismatch_compare(const Mismatch& mis1, const Mismatch& mis2)   //sort mismatch compare function
{
	return mis1.pos<mis2.pos;
}


inline bool spliced_sam_compare_stable(const Sam_Info& sam1, const Sam_Info& sam2)
{
	if(sam1.chrom!=sam1.chrom)
		return sam1.chrom<sam2.chrom;
	else if(sam1.start_pos!=sam2.start_pos)
		return sam1.start_pos<sam2.start_pos;
	else
	{
		int intron1=0, intron2=0;
		for(size_t i=0;i<sam1.jump_code.size();i++)
		{
			if(i%2==1)
				intron1+=sam1.jump_code[i];
		}
		for(size_t i=0;i<sam2.jump_code.size();i++)
		{
			if(i%2==1)
				intron2+=sam2.jump_code[i];
		}
		return intron1<intron2;
	}
}

class mapsplice_report
{
public: 
	map<string, size_t> chrom_map;
	vector<Bwtmap_Info> bwt_vector;
	Read_Block r_block;
	vector<Sam_Info> fusion_vec;
	vector<Sam_Info> pair_vec;
	string global_qual_string;
	int segment_length;
	int segment_number;
	int trunc_len;
	int trunc_left_len;
	int extend_len;
	int extend_left_len;
	int max_intron_thresh;
	int seg_thresh;
	int fusion_thresh;
	bool Fusion;
	bool fa;
	bool read_aligned;
	bool read_unspliced;
	bool pair_end;
	bool read_paired;
	bool output_filtered;
	vector<Splicemap_File*> fixed_hole_vec;
	vector<Splicemap_File*> fixed_hmer_vec;
	vector<Splicemap_File*> fixed_head_vec;
	vector<Splicemap_File*> fixed_tail_vec;
	map<string, ofstream*> fusion_fs_vec;
	ofstream sam_fs;
	ofstream sam_pe_fs;
	ofstream sam_fo_fs;


	mapsplice_report(int seg_len, int seg_num, int ext_len, int max_intron_threshold, int seg_threshold, int min_fusion_seg)
	{
		segment_length = seg_len;
		segment_number = seg_num;
		trunc_len=seg_len/2;
		trunc_left_len=seg_len-trunc_len;
		extend_len=ext_len;
		extend_left_len=seg_len-extend_len;
		max_intron_thresh=max_intron_threshold;
		Fusion=false;
		seg_thresh=seg_threshold;
		global_qual_string.assign(seg_num*seg_len,'I');
		fusion_thresh=min_fusion_seg;
		pair_end=false;
		output_filtered=false;
	}

	~mapsplice_report()
	{
	}

	void load_chrom_table(char* chrom_len_file, char* fixhole_path, char* fixhmer_path, char* fixhead_path, char* fixtail_path, char* sam_file, char* fusion_path, char* sam_pe_file, char* sam_fo_file)
	{
		sam_fs.open(sam_file);
		if(sam_pe_file!=NULL)
		{
			pair_end = true;
			sam_pe_fs.open(sam_pe_file);
		}
		if(sam_fo_file!=NULL)
		{
			output_filtered = true;
			sam_fo_fs.open(sam_fo_file);
		}
		ifstream chrom_len_fs(chrom_len_file);
		if( !chrom_len_fs ) 
		{
			fprintf(stderr,"error: open chrom_len_file error\n");
			exit(1);
		}  
		string line;
		int i=0;
		while(getline(chrom_len_fs, line))
		{
			if(line=="")
				continue;
			char chrom_name_tmp[MAX_CHROM_LENGTH];
			size_t chrom_bp_len;
			sscanf(line.c_str(), "%s\t%zd",chrom_name_tmp, &chrom_bp_len);  //zd for linux
			string chrom_name=chrom_name_tmp;
			chrom_map.insert(make_pair(chrom_name,i));

			string hole_filename=fixhole_path + chrom_name + ".txt.fixed";
			string hmer_filename=fixhmer_path + chrom_name + ".txt.fixed";
			string head_filename=fixhead_path + chrom_name + ".txt.fixed";
			string tail_filename=fixtail_path + chrom_name + ".txt.fixed";

			fixed_hole_vec.push_back(new Splicemap_File(hole_filename.c_str(), pair_end));
			fixed_hmer_vec.push_back(new Splicemap_File(hmer_filename.c_str(), pair_end));
			fixed_head_vec.push_back(new Splicemap_File(head_filename.c_str(), pair_end));
			fixed_tail_vec.push_back(new Splicemap_File(tail_filename.c_str(), pair_end));
			i++;
		}
		chrom_len_fs.close();
		if(fusion_path!=NULL)
		{
			Fusion=true;		
			map<string, size_t>::iterator it1;
			for(it1=chrom_map.begin();it1!=chrom_map.end();it1++)
			{
				map<string, size_t>::iterator it2;
				for(it2=chrom_map.begin();it2!=chrom_map.end();it2++)
				{
					string fusion_filename=fusion_path+it1->first+"_"+it2->first+".txt";
					fusion_fs_vec.insert(make_pair(it1->first + it2->first,new ofstream(fusion_filename.c_str())));
				}
			}
		}
	}

	int check_contig(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)  //check 2 bwtmap segment are connected
	{
		if (seg1.chrom!=seg2.chrom || seg1.strand!=seg2.strand || seg2.start - seg1.end > 1)
			return 0;            // case 0: no relationship
		if(pair_end && seg1.pair_no!=seg2.pair_no)
			return 0;
		else if(seg2.start - seg1.end == 1) //position related
		{
			if(seg1.strand=="+"&&seg2.start_seg_no == seg1.end_seg_no + 1)     // strand +
			{
				return 1;        // case 1: contig
			} 
			else if  (seg1.strand=="-"&&seg2.start_seg_no == seg1.end_seg_no - 1)                 // strand -
			{
				return 1;        // case 1: contig
			}
		}
		else if (seg2.start - seg1.end < 1)  //need to continue search
			return -1;
		return 0;                 // case 0: no relationship
	}

	void combine_contig()   //combine segments into contigs
	{
		for (vector<Bwtmap_Info>::iterator it=bwt_vector.begin();it+1!=bwt_vector.end();)
		{
			for(vector<Bwtmap_Info>::iterator it2=it+1;;it2++)
			{
				if(it2==bwt_vector.end())
				{
					it++;
					break;
				}
				int contig_status= check_contig(*it, *(it2));
				if (contig_status == 1)   //contig
				{
					(*it).end_seg_no = (*it2).end_seg_no;
					(*it).end = (*it2).end;
					for(size_t i=0;i<(*it2).mismatches.size();i++)
					{
						(*it).mismatches.push_back((*it2).mismatches[i]);
					}
					bwt_vector.erase(it2);
					break;
				}
				else if(contig_status == 0)
				{
					it++;
					break;
				}
			}
		}
	}

	void generate_sam_head(size_t current_contig_index)
	{
		if(bwt_vector[current_contig_index].splice_head.size()==0)  //no fix head
		{
			Sam_Info new_sam_info;
			new_sam_info.start_contig_index=current_contig_index;
			new_sam_info.start_seg_no=bwt_vector[current_contig_index].start_seg_no;
			new_sam_info.chrom=bwt_vector[current_contig_index].chrom;
			generate_sam_internal(new_sam_info,current_contig_index,false,0);
		}
		else                                                //has fix head
		{
			for(size_t i=0;i<bwt_vector[current_contig_index].splice_head.size();i++)
			{
				Sam_Info new_sam_info;
				new_sam_info.start_contig_index=current_contig_index;
				new_sam_info.chrom=bwt_vector[current_contig_index].chrom;
				if(bwt_vector[current_contig_index].strand=="+")
					new_sam_info.start_seg_no=HEAD_SEG_NO;
				else
					new_sam_info.start_seg_no=segment_number;
				if(bwt_vector[current_contig_index].splice_head[i].jump_code[0]!=0)   // no 0M case
				{
					new_sam_info.start_pos=bwt_vector[current_contig_index].splice_head[i].start_pos+1;
					for(size_t j=0;j<bwt_vector[current_contig_index].splice_head[i].jump_code.size();j++)
					{
						new_sam_info.jump_code.push_back(bwt_vector[current_contig_index].splice_head[i].jump_code[j]);
					}
				}
				else    //0M case
				{
					new_sam_info.start_pos=bwt_vector[current_contig_index].splice_head[i].start_pos+1+bwt_vector[current_contig_index].splice_head[i].jump_code[1];
					for(size_t j=2;j<bwt_vector[current_contig_index].splice_head[i].jump_code.size();j++)
					{
						new_sam_info.jump_code.push_back(bwt_vector[current_contig_index].splice_head[i].jump_code[j]);
					}
				}
				for(size_t j=0;j<bwt_vector[current_contig_index].splice_head[i].mismatches.size();j++)
				{
					new_sam_info.mismatchs.push_back(bwt_vector[current_contig_index].splice_head[i].mismatches[j]);
				}
				generate_sam_internal(new_sam_info,current_contig_index,false,0);
			}
		}
	}

	void generate_sam_internal(Sam_Info &my_sam_info, size_t current_contig_index, bool skip_edge_seg, int overlap_len)
	{
		if(my_sam_info.start_pos==-1) // started without fix head
			my_sam_info.start_pos=bwt_vector[current_contig_index].start+1;
		if(my_sam_info.jump_code.size()==0)
		{
			int match_len=(abs(bwt_vector[current_contig_index].end_seg_no-bwt_vector[current_contig_index].start_seg_no)+1)*segment_length;
			my_sam_info.jump_code.push_back(match_len);
		}
		else
		{
			int match_len=(abs(bwt_vector[current_contig_index].end_seg_no-bwt_vector[current_contig_index].start_seg_no)+1)*segment_length - overlap_len;
			if(skip_edge_seg)
				match_len-=segment_length;
			my_sam_info.jump_code[my_sam_info.jump_code.size()-1]+=match_len;
		}
		if(bwt_vector[current_contig_index].splice_internal.size()==0)   //no more internal splice, generate tail and finish
		{
			for(size_t i=0;i<bwt_vector[current_contig_index].mismatches.size();i++) //add mismatches but skip those in fix hmer or fix hole
			{
				if(skip_edge_seg) //special case
				{
					if(bwt_vector[current_contig_index].mismatches[i].seg_no==bwt_vector[current_contig_index].start_seg_no)
						continue;
					if(bwt_vector[current_contig_index].strand=="+"  &&  bwt_vector[current_contig_index].mismatches[i].seg_no==bwt_vector[current_contig_index].start_seg_no+1 && bwt_vector[current_contig_index].mismatches[i].pos<overlap_len)
						continue;
					if(bwt_vector[current_contig_index].strand=="-"  &&  bwt_vector[current_contig_index].mismatches[i].seg_no==bwt_vector[current_contig_index].start_seg_no-1 && bwt_vector[current_contig_index].mismatches[i].pos<overlap_len)
						continue;
				}
				else //normal case
				{
					if(bwt_vector[current_contig_index].mismatches[i].seg_no==bwt_vector[current_contig_index].start_seg_no && bwt_vector[current_contig_index].mismatches[i].pos<overlap_len)
						continue;
				}
				my_sam_info.mismatchs.push_back(bwt_vector[current_contig_index].mismatches[i]);
				if(bwt_vector[current_contig_index].strand=="+")
					my_sam_info.mismatchs[my_sam_info.mismatchs.size()-1].pos=(bwt_vector[current_contig_index].mismatches[i].seg_no-1)*segment_length + bwt_vector[current_contig_index].mismatches[i].pos;
				else
					my_sam_info.mismatchs[my_sam_info.mismatchs.size()-1].pos=(segment_number - bwt_vector[current_contig_index].mismatches[i].seg_no)*segment_length + bwt_vector[current_contig_index].mismatches[i].pos;
			}
			generate_sam_tail(my_sam_info, current_contig_index); 
		}
		else      //more internal splice
		{
			for(size_t i=0;i<bwt_vector[current_contig_index].splice_internal.size();i++)
			{
				bool next_skip_edge_seg=false;
				int next_overlap_len=0;
				int seg_gap=abs(bwt_vector[current_contig_index].end_seg_no - bwt_vector[bwt_vector[current_contig_index].splice_internal[i].end_contig].start_seg_no);
				if(seg_gap==2)  //gap hole
				{
					next_overlap_len=extend_len;
				}
				else if(seg_gap==1)  // continuous hole
				{
					next_overlap_len=trunc_len;
				}
				else if(seg_gap==0)
				{
					next_overlap_len=extend_len;
					next_skip_edge_seg=true;
				}
				Sam_Info new_sam_info(my_sam_info);
				new_sam_info.jump_code[my_sam_info.jump_code.size()-1]-=next_overlap_len;
				if(next_skip_edge_seg) // special case
					new_sam_info.jump_code[my_sam_info.jump_code.size()-1]-=segment_length;
				new_sam_info.jump_code[my_sam_info.jump_code.size()-1]+=bwt_vector[current_contig_index].splice_internal[i].jump_code[0];
				for(size_t j=1;j<bwt_vector[current_contig_index].splice_internal[i].jump_code.size();j++)
				{
					new_sam_info.jump_code.push_back(bwt_vector[current_contig_index].splice_internal[i].jump_code[j]);
				}
				for(size_t j=0;j<bwt_vector[current_contig_index].mismatches.size();j++)  // add mismatches from bwtmap
				{
					if(skip_edge_seg||next_skip_edge_seg) //special case
					{
						if(skip_edge_seg&&bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].start_seg_no)
							continue;
						if(skip_edge_seg&&bwt_vector[current_contig_index].strand=="+"  &&  bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].start_seg_no+1 && bwt_vector[current_contig_index].mismatches[j].pos<overlap_len)
							continue;
						if(skip_edge_seg&&bwt_vector[current_contig_index].strand=="-"  &&  bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].start_seg_no-1 && bwt_vector[current_contig_index].mismatches[j].pos<overlap_len)
							continue;
						if(!skip_edge_seg&&bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].start_seg_no && bwt_vector[current_contig_index].mismatches[j].pos<overlap_len)
							continue;
						if(next_skip_edge_seg&&bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].end_seg_no)
							continue;
						if(next_skip_edge_seg&&bwt_vector[current_contig_index].strand=="+"  &&  bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].end_seg_no-1 && bwt_vector[current_contig_index].mismatches[j].pos>=segment_length-next_overlap_len)
							continue;
						if(next_skip_edge_seg&&bwt_vector[current_contig_index].strand=="-"  &&  bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].end_seg_no+1 && bwt_vector[current_contig_index].mismatches[j].pos>=segment_length-next_overlap_len)
							continue;
						if(!next_skip_edge_seg&&bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].end_seg_no && bwt_vector[current_contig_index].mismatches[j].pos>=segment_length-next_overlap_len)
							continue;
					}
					else //normal case
					{
						if(bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].start_seg_no && bwt_vector[current_contig_index].mismatches[j].pos<overlap_len)
							continue;
						if(bwt_vector[current_contig_index].mismatches[j].seg_no==bwt_vector[current_contig_index].end_seg_no && bwt_vector[current_contig_index].mismatches[j].pos>=segment_length-next_overlap_len)
							continue;
					}
					new_sam_info.mismatchs.push_back(bwt_vector[current_contig_index].mismatches[j]);
					if(bwt_vector[current_contig_index].strand=="+")
						new_sam_info.mismatchs[new_sam_info.mismatchs.size()-1].pos=(bwt_vector[current_contig_index].mismatches[j].seg_no-1)*segment_length + bwt_vector[current_contig_index].mismatches[j].pos;
					else
						new_sam_info.mismatchs[new_sam_info.mismatchs.size()-1].pos=(segment_number - bwt_vector[current_contig_index].mismatches[j].seg_no)*segment_length + bwt_vector[current_contig_index].mismatches[j].pos;

				}
				for(size_t j=0;j<bwt_vector[current_contig_index].splice_internal[i].mismatches.size();j++)   // add all mismatches in splice region
				{
					new_sam_info.mismatchs.push_back(bwt_vector[current_contig_index].splice_internal[i].mismatches[j]);
				}
				generate_sam_internal(new_sam_info,bwt_vector[current_contig_index].splice_internal[i].end_contig,next_skip_edge_seg,next_overlap_len);
			}
		}
	}

	void generate_sam_tail(Sam_Info &my_sam_info, size_t current_contig_index)
	{
		if(bwt_vector[current_contig_index].splice_tail.size()==0)   // no tail
		{
			my_sam_info.end_seg_no=bwt_vector[current_contig_index].end_seg_no;
			print_saminfo(my_sam_info);
		}
		else              // has tail
		{
			for(size_t i=0;i<bwt_vector[current_contig_index].splice_tail.size();i++)
			{
				Sam_Info new_sam_info(my_sam_info);
				if(bwt_vector[current_contig_index].strand=="+")
					new_sam_info.end_seg_no=segment_number;
				else
					new_sam_info.end_seg_no=HEAD_SEG_NO;
				int match_len=bwt_vector[current_contig_index].splice_tail[i].left_length-1+bwt_vector[current_contig_index].splice_tail[i].jump_code[0]-segment_length;
				new_sam_info.jump_code[my_sam_info.jump_code.size()-1]+=match_len;
				for(size_t j=1;j<bwt_vector[current_contig_index].splice_tail[i].jump_code.size();j++)
				{
					if(j%2==1&&bwt_vector[current_contig_index].splice_tail[i].jump_code[j+1]==0)
						break;
					new_sam_info.jump_code.push_back(bwt_vector[current_contig_index].splice_tail[i].jump_code[j]);
				}
				for(size_t j=0;j<bwt_vector[current_contig_index].splice_tail[i].mismatches.size();j++)
				{
					new_sam_info.mismatchs.push_back(bwt_vector[current_contig_index].splice_tail[i].mismatches[j]);
				}
				print_saminfo(new_sam_info);
			}
		}
	}

	void print_saminfo(Sam_Info &my_sam_info)
	{
		int connected_segment=abs((int)my_sam_info.start_seg_no-(int)my_sam_info.end_seg_no) + 1;
		if(connected_segment>=seg_thresh)   // only output above threshold sam
		{
			read_aligned=true;
			if(pair_end)
			{
				pair_vec.push_back(my_sam_info);
				return;
			}
			/*if(bwt_vector[my_sam_info.start_contig_index].strand=="+" && my_sam_info.start_seg_no!=HEAD_SEG_NO)  //shift mismatch for uncomplete read sam
			{
			int head_lack_len= ((int)my_sam_info.start_seg_no - HEAD_SEG_NO) * segment_length;
			for(size_t i=0;i<my_sam_info.mismatchs.size();i++)
			{
			my_sam_info.mismatchs[i].pos-=head_lack_len;
			}
			}
			else if(bwt_vector[my_sam_info.start_contig_index].strand=="-" && my_sam_info.start_seg_no!=segment_number)
			{
			int head_lack_len= (segment_number - (int)my_sam_info.start_seg_no) * segment_length;
			for(size_t i=0;i<my_sam_info.mismatchs.size();i++)
			{
			my_sam_info.mismatchs[i].pos-=head_lack_len;
			}
			}*/
			string quality_string;
			if(!fa)
			{
				if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //get quality string for score compute
				{
					for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
					{
						quality_string.append(r_block.get_seg_qual((int)i));
					}
				}
				else    
				{
					for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
					{
						quality_string.append(r_block.get_revcom_seg_qual((int)i));
					}
				}
			}
			int quality_score=255;
			if(!fa)
				quality_score = GetQualityScore(my_sam_info.mismatchs,quality_string);
			sam_fs<<bwt_vector[0].unmapped_index<<"~"<<bwt_vector[0].read_id<<"\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")
				sam_fs<<"0"<<"\t";
			else
				sam_fs<<"16"<<"\t";
			sam_fs<<my_sam_info.chrom<<"\t";
			sam_fs<<my_sam_info.start_pos<<"\t";
			sam_fs<<quality_score<<"\t";
			for(size_t i=0;i<my_sam_info.jump_code.size();i++)
			{
				if(i%2==1&&my_sam_info.jump_code[i+1]==0)
				{
					my_sam_info.jump_code[i]+=my_sam_info.jump_code[i+2];
					sam_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
					i+=3;
				}
				else
					sam_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
			}
			sam_fs<<"\t*\t0\t0\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //out put + sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
				{
					sam_fs<<r_block.get_seg_seq((int)i);
				}
				sam_fs<<"\t";
			}
			else                   //out put - sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
				{
					sam_fs<<r_block.get_revcom_seg_seq((int)i);
				}
				sam_fs<<"\t";
			}
			if(!fa)    //fq
				sam_fs<<quality_string<<"\t";
			else       //fa
				sam_fs<<global_qual_string.substr(0,connected_segment*segment_length)<<"\t";
			sam_fs<<"NM:i:"<<my_sam_info.mismatchs.size()<<"\t";	
			/*if(my_sam_info.mismatchs.size()!=0)   //output detailed mismatch
			{
				sort(my_sam_info.mismatchs.begin(),my_sam_info.mismatchs.end(),mismatch_compare);   //sort mismatch by pos
				for(size_t i=0;i<my_sam_info.mismatchs.size();i++)  //output mismatch
				{
					sam_fs<<my_sam_info.mismatchs[i].pos<<":"<<my_sam_info.mismatchs[i].src_char<<">"<<my_sam_info.mismatchs[i].dst_char;
					if(i!=my_sam_info.mismatchs.size()-1)
					{
						sam_fs<<",";
					}
				}
			}*/
			sam_fs<<endl;
		}
		else if(Fusion&&(!read_aligned)&&(!read_unspliced))//add incomplete sam to fusion search candidate
		{
			fusion_vec.push_back(my_sam_info);
		}
	}

	void print_sam_info_pe(Sam_Info &my_sam_info)
	{
		int connected_segment=abs((int)my_sam_info.start_seg_no-(int)my_sam_info.end_seg_no) + 1;
		string quality_string;
		if(!fa)
		{
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //get quality string for score compute
			{
				for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
				{
					quality_string.append(r_block.get_seg_qual((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number));
				}
			}
			else    
			{
				for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
				{
					quality_string.append(r_block.get_revcom_seg_qual((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number));
				}
			}
		}
		int quality_score=255;
		if(!fa)
			quality_score = GetQualityScore(my_sam_info.mismatchs,quality_string);
		if(my_sam_info.pair_index!=-1)  //output paired read
		{
			sam_pe_fs<<bwt_vector[my_sam_info.start_contig_index].unmapped_index<<"~"<<bwt_vector[my_sam_info.start_contig_index].read_id<<"/"<<bwt_vector[my_sam_info.start_contig_index].pair_no<<"\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")
				sam_pe_fs<<"0"<<"\t";
			else
				sam_pe_fs<<"16"<<"\t";
			sam_pe_fs<<my_sam_info.chrom<<"\t";
			sam_pe_fs<<my_sam_info.start_pos<<"\t";
			sam_pe_fs<<quality_score<<"\t";
			for(size_t i=0;i<my_sam_info.jump_code.size();i++)
			{
				if(i%2==1&&my_sam_info.jump_code[i+1]==0)
				{
					my_sam_info.jump_code[i]+=my_sam_info.jump_code[i+2];
					sam_pe_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
					i+=3;
				}
				else
					sam_pe_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
			}
			sam_pe_fs<<"\t=\t"<<pair_vec[my_sam_info.pair_index].start_pos<<"\t"<<pair_vec[my_sam_info.pair_index].start_pos- my_sam_info.start_pos;                ;
			sam_pe_fs<<"\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //out put + sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
				{
					sam_pe_fs<<r_block.get_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_pe_fs<<"\t";
			}
			else                   //out put - sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
				{
					sam_pe_fs<<r_block.get_revcom_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_pe_fs<<"\t";
			}
			if(!fa)    //fq
				sam_pe_fs<<quality_string<<"\t";
			else       //fa
				sam_pe_fs<<global_qual_string.substr(0,connected_segment*segment_length)<<"\t";
			sam_pe_fs<<"NM:i:"<<my_sam_info.mismatchs.size()<<"\t";
			sam_pe_fs<<endl;
		}
		else if(!read_paired)    // single read alignment
		{
			sam_fs<<bwt_vector[my_sam_info.start_contig_index].unmapped_index<<"~"<<bwt_vector[my_sam_info.start_contig_index].read_id<<"/"<<bwt_vector[my_sam_info.start_contig_index].pair_no<<"\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")
				sam_fs<<"0"<<"\t";
			else
				sam_fs<<"16"<<"\t";
			sam_fs<<my_sam_info.chrom<<"\t";
			sam_fs<<my_sam_info.start_pos<<"\t";
			sam_fs<<quality_score<<"\t";
			for(size_t i=0;i<my_sam_info.jump_code.size();i++)
			{
				if(i%2==1&&my_sam_info.jump_code[i+1]==0)
				{
					my_sam_info.jump_code[i]+=my_sam_info.jump_code[i+2];
					sam_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
					i+=3;
				}
				else
					sam_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
			}
			sam_fs<<"\t*\t0\t0\t";                ;
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //out put + sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
				{
					sam_fs<<r_block.get_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_fs<<"\t";
			}
			else                   //out put - sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
				{
					sam_fs<<r_block.get_revcom_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_fs<<"\t";
			}
			if(!fa)    //fq
				sam_fs<<quality_string<<"\t";
			else       //fa
				sam_fs<<global_qual_string.substr(0,connected_segment*segment_length)<<"\t";
			sam_fs<<"NM:i:"<<my_sam_info.mismatchs.size()<<"\t";	
			sam_fs<<endl;
		}
		else if(output_filtered)
		{
			sam_fo_fs<<bwt_vector[my_sam_info.start_contig_index].unmapped_index<<"~"<<bwt_vector[my_sam_info.start_contig_index].read_id<<"/"<<bwt_vector[my_sam_info.start_contig_index].pair_no<<"\t";
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")
				sam_fo_fs<<"0"<<"\t";
			else
				sam_fo_fs<<"16"<<"\t";
			sam_fo_fs<<my_sam_info.chrom<<"\t";
			sam_fo_fs<<my_sam_info.start_pos<<"\t";
			sam_fo_fs<<quality_score<<"\t";
			for(size_t i=0;i<my_sam_info.jump_code.size();i++)
			{
				if(i%2==1&&my_sam_info.jump_code[i+1]==0)
				{
					my_sam_info.jump_code[i]+=my_sam_info.jump_code[i+2];
					sam_fo_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
					i+=3;
				}
				else
					sam_fo_fs<<my_sam_info.jump_code[i]<< jumpcode_type[i%2];
			}
			sam_fo_fs<<"\t*\t0\t0\t";                ;
			if(bwt_vector[my_sam_info.start_contig_index].strand=="+")  //out put + sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i<=my_sam_info.end_seg_no;i++)
				{
					sam_fo_fs<<r_block.get_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_fo_fs<<"\t";
			}
			else                   //out put - sequence
			{
				for(size_t i=my_sam_info.start_seg_no;i>=my_sam_info.end_seg_no;i--)
				{
					sam_fo_fs<<r_block.get_revcom_seg_seq((int)i + (bwt_vector[my_sam_info.start_contig_index].pair_no-1)*segment_number);
				}
				sam_fo_fs<<"\t";
			}
			if(!fa)    //fq
				sam_fo_fs<<quality_string<<"\t";
			else       //fa
				sam_fo_fs<<global_qual_string.substr(0,connected_segment*segment_length)<<"\t";
			sam_fo_fs<<"NM:i:"<<my_sam_info.mismatchs.size()<<"\t";	
			sam_fo_fs<<endl;
		}
	}


	void load_splice_info(Splicemap_File& splice_file, string read_id, int type)  
	{
		if(!splice_file.current_splice.empty()&&splice_file.current_splice[0].read_id==read_id)
		{
			for(size_t i=0;i<splice_file.current_splice.size();i++)
			{
				if(type==HOLE||type==HMER)
				{
					bwt_vector[splice_file.current_splice[i].start_contig].splice_internal.push_back(splice_file.current_splice[i]);
					bwt_vector[splice_file.current_splice[i].end_contig].head = 0;
				}
				else if(type==HEAD)
				{
					bwt_vector[splice_file.current_splice[i].end_contig].splice_head.push_back(splice_file.current_splice[i]);
				}
				else if(type==TAIL)
				{
					bwt_vector[splice_file.current_splice[i].start_contig].splice_tail.push_back(splice_file.current_splice[i]);
				}
			}
			splice_file.load_next_block();
		}
	}

	void search_fusion()
	{
		/*for (size_t i = 0; i < bwt_vector.size(); i++)  // add none spliced contig to fusion search candidate
		{
			if(bwt_vector[i].splice_internal.size()==0&&bwt_vector[i].splice_head.size()==0&&bwt_vector[i].splice_tail.size()==0)
			{
				Sam_Info new_sam_info(bwt_vector[i],i);
				for(size_t j=0;j<new_sam_info.mismatchs.size();j++)
				{
					if(bwt_vector[i].strand=="+")
						new_sam_info.mismatchs[j].pos=(new_sam_info.mismatchs[j].seg_no-1)*segment_length + new_sam_info.mismatchs[j].pos;
					else
						new_sam_info.mismatchs[j].pos=(segment_number - new_sam_info.mismatchs[j].seg_no)*segment_length + new_sam_info.mismatchs[j].pos;
				}
				fusion_vec.push_back(new_sam_info);
			}
		}*/
		for(size_t i=0;i<fusion_vec.size();i++)       //compute mismatch for - strand
		{
			if(bwt_vector[fusion_vec[i].start_contig_index].strand=="-")
			{
				for(size_t j=0;j<fusion_vec[i].mismatchs.size();j++)
				{
					fusion_vec[i].mismatchs[j].pos=segment_number*segment_length - fusion_vec[i].mismatchs[j].pos -1 ; 
				}
			}
		}
		for(size_t i=0; i<fusion_vec.size();i++)   //convert to 0 index
		{
			fusion_vec[i].start_pos -=1;
		}
		for(size_t i=0; i<fusion_vec.size();i++)
		{
			for(size_t j=i+1;j<fusion_vec.size();j++)
			{
				if(pair_end && bwt_vector[fusion_vec[i].start_contig_index].pair_no!= bwt_vector[fusion_vec[j].start_contig_index].pair_no)
					continue;
				int end[2]={fusion_vec[i].start_pos-1,fusion_vec[j].start_pos-1};
				for(size_t k=0;k<fusion_vec[i].jump_code.size();k++)
				{
					end[0]+=fusion_vec[i].jump_code[k];
				}
				for(size_t k=0;k<fusion_vec[j].jump_code.size();k++)
				{
					end[1]+=fusion_vec[j].jump_code[k];
				}
				if(bwt_vector[fusion_vec[i].start_contig_index].strand==bwt_vector[fusion_vec[i].start_contig_index].strand&&bwt_vector[fusion_vec[i].start_contig_index].chrom==bwt_vector[fusion_vec[i].start_contig_index].chrom)
				{
					if(abs(fusion_vec[j].start_pos-end[0])<max_intron_thresh||abs(fusion_vec[i].start_pos-end[1])<max_intron_thresh)
						continue;
				}
				size_t first_index=0;
				size_t second_index=0;
				size_t first_end=0;
				size_t second_end=0;
				if(max(fusion_vec[i].start_seg_no,fusion_vec[i].end_seg_no)<min(fusion_vec[j].start_seg_no,fusion_vec[j].end_seg_no))
				{
					first_index=i;
					second_index=j;
					first_end=0;
					second_end=1;
				}
				else if(max(fusion_vec[j].start_seg_no,fusion_vec[j].end_seg_no)<min(fusion_vec[i].start_seg_no,fusion_vec[i].end_seg_no))
				{
					first_index=j;
					second_index=i;
					first_end=1;
					second_end=0;
				}
				else
					continue;
				int seg_sum=abs((int)fusion_vec[i].end_seg_no- (int)fusion_vec[i].start_seg_no)+ 1 +abs((int)fusion_vec[j].end_seg_no- (int)fusion_vec[j].start_seg_no) + 1;
				if(seg_sum<fusion_thresh)
					continue;
				ofstream* fusion_fs = fusion_fs_vec[bwt_vector[fusion_vec[first_index].start_contig_index].chrom + bwt_vector[fusion_vec[second_index].start_contig_index].chrom];
				size_t seg_gap=min(fusion_vec[second_index].start_seg_no,fusion_vec[second_index].end_seg_no)- max(fusion_vec[first_index].start_seg_no,fusion_vec[first_index].end_seg_no);
				if(seg_gap==1)   // continuous fusion
				{
					if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="+"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="+")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"++"<<"\t"<<end[first_end] - trunc_len +1<<"\t"<<fusion_vec[second_index].start_pos + trunc_len-1<<"\t"<<fusion_vec[first_index].start_seg_no<<"-"<<fusion_vec[first_index].end_seg_no<<"\t"<<fusion_vec[second_index].start_seg_no<<"-"<<fusion_vec[second_index].end_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].start_seg_no;k<fusion_vec[first_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(0,trunc_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(trunc_len);
						for(size_t k=fusion_vec[second_index].start_seg_no+1;k<=fusion_vec[second_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(trunc_left_len)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(0,trunc_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="+"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="-")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"+-"<<"\t"<<end[first_end] - trunc_len +1<<"\t"<<end[second_end]  - trunc_len +1<<"\t"<<fusion_vec[first_index].start_seg_no<<"-"<<fusion_vec[first_index].end_seg_no<<"\t"<<fusion_vec[second_index].end_seg_no<<"-"<<fusion_vec[second_index].start_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].start_seg_no;k<fusion_vec[first_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(0,trunc_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(trunc_len);
						for(size_t k=fusion_vec[second_index].end_seg_no+1;k<=fusion_vec[second_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(trunc_left_len)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(0,trunc_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="-"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="-")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"--"<<"\t"<<fusion_vec[first_index].start_pos + trunc_len -1<<"\t"<<end[second_end] - trunc_len + 1<<"\t"<<fusion_vec[first_index].end_seg_no<<"-"<<fusion_vec[first_index].start_seg_no<<"\t"<<fusion_vec[second_index].end_seg_no<<"-"<<fusion_vec[second_index].start_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].end_seg_no;k<fusion_vec[first_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(0,trunc_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(trunc_len);
						for(size_t k=fusion_vec[second_index].end_seg_no+1;k<=fusion_vec[second_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(trunc_left_len)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(0,trunc_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="-"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="+")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"-+"<<"\t"<<fusion_vec[first_index].start_pos + trunc_len -1<<"\t"<<fusion_vec[second_index].start_pos+ trunc_len -1<<"\t"<<fusion_vec[first_index].end_seg_no<<"-"<<fusion_vec[first_index].start_seg_no<<"\t"<<fusion_vec[second_index].start_seg_no<<"-"<<fusion_vec[second_index].end_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].end_seg_no;k<fusion_vec[first_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(0,trunc_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(trunc_len);
						for(size_t k=fusion_vec[second_index].start_seg_no+1;k<=fusion_vec[second_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(trunc_left_len)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(0,trunc_len);
					}
					(*fusion_fs)<<"\t";
					if(!fa)    //fq
					{
						for(int k=1; k<=segment_number;k++)   // output quality string
						{
							(*fusion_fs)<<r_block.get_seg_qual(k);
						}
					}
					else       //fa
						(*fusion_fs)<<global_qual_string;
					(*fusion_fs)<<"\t";
					for(size_t k=0;k<fusion_vec[first_index].jump_code.size();k++)   //output jump code
					{
						(*fusion_fs)<<fusion_vec[first_index].jump_code[k]<< jumpcode_type[k%2];
					}
					(*fusion_fs)<<"\t";
					for(size_t k=0;k<fusion_vec[second_index].jump_code.size();k++)   //output jump code
					{
						(*fusion_fs)<<fusion_vec[second_index].jump_code[k]<< jumpcode_type[k%2];
					}
					(*fusion_fs)<<"\t";
					int left_cut_pos= (int)max(fusion_vec[first_index].start_seg_no,fusion_vec[first_index].end_seg_no) * segment_length - trunc_len;
					int right_cut_pos= (int)(min(fusion_vec[second_index].start_seg_no,fusion_vec[second_index].end_seg_no)-1) * segment_length + trunc_len;
					for(size_t k=0;k<fusion_vec[first_index].mismatchs.size();k++)  //output mismatches for first part
					{
						if(fusion_vec[first_index].mismatchs[k].pos<left_cut_pos||fusion_vec[first_index].mismatchs[k].pos>=right_cut_pos)
							(*fusion_fs)<<fusion_vec[first_index].mismatchs[k].pos<<":"<<fusion_vec[first_index].mismatchs[k].src_char<<">"<<fusion_vec[first_index].mismatchs[k].dst_char<<",";
					}
					for(size_t k=0;k<fusion_vec[second_index].mismatchs.size();k++) //output mismatches for second part
					{
						if(fusion_vec[second_index].mismatchs[k].pos<left_cut_pos||fusion_vec[second_index].mismatchs[k].pos>=right_cut_pos)
							(*fusion_fs)<<fusion_vec[second_index].mismatchs[k].pos<<":"<<fusion_vec[second_index].mismatchs[k].src_char<<">"<<fusion_vec[second_index].mismatchs[k].dst_char<<",";
					}
					(*fusion_fs)<<endl;
				}
				if(seg_gap==2)   // continuous fusion
				{
					if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="+"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="+")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"++"<<"\t"<<end[first_end] - extend_len +1<<"\t"<<fusion_vec[second_index].start_pos + extend_len-1<<"\t"<<fusion_vec[first_index].start_seg_no<<"-"<<fusion_vec[first_index].end_seg_no<<"\t"<<fusion_vec[second_index].start_seg_no<<"-"<<fusion_vec[second_index].end_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].start_seg_no;k<fusion_vec[first_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(0,extend_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(extend_len);
						for(size_t k=fusion_vec[second_index].start_seg_no+1;k<=fusion_vec[second_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(extend_left_len)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no+1)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(0,extend_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="+"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="-")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"+-"<<"\t"<<end[first_end] - extend_len +1<<"\t"<<end[second_end]  - extend_len +1<<"\t"<<fusion_vec[first_index].start_seg_no<<"-"<<fusion_vec[first_index].end_seg_no<<"\t"<<fusion_vec[second_index].end_seg_no<<"-"<<fusion_vec[second_index].start_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].start_seg_no;k<fusion_vec[first_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(0,extend_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(extend_len);
						for(size_t k=fusion_vec[second_index].end_seg_no+1;k<=fusion_vec[second_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no).substr(extend_left_len)<<r_block.get_seg_seq(fusion_vec[first_index].end_seg_no+1)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(0,extend_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="-"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="-")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"--"<<"\t"<<fusion_vec[first_index].start_pos + extend_len -1<<"\t"<<end[second_end] - extend_len + 1<<"\t"<<fusion_vec[first_index].end_seg_no<<"-"<<fusion_vec[first_index].start_seg_no<<"\t"<<fusion_vec[second_index].end_seg_no<<"-"<<fusion_vec[second_index].start_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].end_seg_no;k<fusion_vec[first_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(0,extend_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(extend_len);
						for(size_t k=fusion_vec[second_index].end_seg_no+1;k<=fusion_vec[second_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(extend_left_len)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no+1)<<r_block.get_seg_seq(fusion_vec[second_index].end_seg_no).substr(0,extend_len);
					}
					else if(bwt_vector[fusion_vec[first_index].start_contig_index].strand=="-"&&bwt_vector[fusion_vec[second_index].start_contig_index].strand=="+")
					{
						(*fusion_fs)<<bwt_vector[fusion_vec[i].start_contig_index].unmapped_index<<"~"<<bwt_vector[fusion_vec[i].start_contig_index].read_id;
						if(pair_end)
							(*fusion_fs)<<"/"<<bwt_vector[fusion_vec[i].start_contig_index].pair_no;
						(*fusion_fs)<<"\t"<<bwt_vector[fusion_vec[first_index].start_contig_index].chrom<<"\t"<<bwt_vector[fusion_vec[second_index].start_contig_index].chrom<<"\t"<<"-+"<<"\t"<<fusion_vec[first_index].start_pos + extend_len -1<<"\t"<<fusion_vec[second_index].start_pos+ extend_len -1<<"\t"<<fusion_vec[first_index].end_seg_no<<"-"<<fusion_vec[first_index].start_seg_no<<"\t"<<fusion_vec[second_index].start_seg_no<<"-"<<fusion_vec[second_index].end_seg_no<<"\t";	
						for(size_t k=fusion_vec[first_index].end_seg_no;k<fusion_vec[first_index].start_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(0,extend_left_len);
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(extend_len);
						for(size_t k=fusion_vec[second_index].start_seg_no+1;k<=fusion_vec[second_index].end_seg_no;k++)
						{
							(*fusion_fs)<<r_block.get_seg_seq(k);
						}
						(*fusion_fs)<<"\t";
						(*fusion_fs)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no).substr(extend_left_len)<<r_block.get_seg_seq(fusion_vec[first_index].start_seg_no+1)<<r_block.get_seg_seq(fusion_vec[second_index].start_seg_no).substr(0,extend_len);
					}
					(*fusion_fs)<<"\t";
					if(!fa)    //fq
					{
						for(int k=1; k<=segment_number;k++)   // output quality string
						{
							(*fusion_fs)<<r_block.get_seg_qual(k);
						}
					}
					else       //fa
						(*fusion_fs)<<global_qual_string;
					(*fusion_fs)<<"\t";
					for(size_t k=0;k<fusion_vec[first_index].jump_code.size();k++)   //output jump code
					{
						(*fusion_fs)<<fusion_vec[first_index].jump_code[k]<< jumpcode_type[k%2];
					}
					(*fusion_fs)<<"\t";
					for(size_t k=0;k<fusion_vec[second_index].jump_code.size();k++)   //output jump code
					{
						(*fusion_fs)<<fusion_vec[second_index].jump_code[k]<< jumpcode_type[k%2];
					}
					(*fusion_fs)<<"\t";
					int left_cut_pos= (int)max(fusion_vec[first_index].start_seg_no,fusion_vec[first_index].end_seg_no) * segment_length - extend_len;
					int right_cut_pos= (int)(min(fusion_vec[second_index].start_seg_no,fusion_vec[second_index].end_seg_no)-1) * segment_length + extend_len;
					for(size_t k=0;k<fusion_vec[first_index].mismatchs.size();k++)  //output mismatches for first part
					{
						if(fusion_vec[first_index].mismatchs[k].pos<left_cut_pos||fusion_vec[first_index].mismatchs[k].pos>=right_cut_pos)
							(*fusion_fs)<<fusion_vec[first_index].mismatchs[k].pos<<":"<<fusion_vec[first_index].mismatchs[k].src_char<<">"<<fusion_vec[first_index].mismatchs[k].dst_char<<",";
					}
					for(size_t k=0;k<fusion_vec[second_index].mismatchs.size();k++) //output mismatches for second part
					{
						if(fusion_vec[second_index].mismatchs[k].pos<left_cut_pos||fusion_vec[second_index].mismatchs[k].pos>=right_cut_pos)
							(*fusion_fs)<<fusion_vec[second_index].mismatchs[k].pos<<":"<<fusion_vec[second_index].mismatchs[k].src_char<<">"<<fusion_vec[second_index].mismatchs[k].dst_char<<",";
					}
					(*fusion_fs)<<endl;
				}
			}
		}
	}

	/*void print_unspliced_read(size_t i)
	{
			string quality_string;
			if(!fa)
			{
				if(bwt_vector[i].strand=="+")  //get quality string for score compute
				{
					for(int i=1;i<=segment_number;i++)
					{
						quality_string.append(r_block.get_seg_qual((int)i));
					}
				}
				else    
				{
					for(int i=segment_number;i>=1;i--)
					{
						quality_string.append(r_block.get_revcom_seg_qual((int)i));
					}
				}
			}
			double quality_score=255;//fa , fq GetQualityScore(my_sam_info.mismatchs,quality_string);
			sam_fs<<bwt_vector[0].unmapped_index<<"~"<<bwt_vector[0].read_id<<"\t";
			if(bwt_vector[i].strand=="+")
				sam_fs<<"0"<<"\t";
			else
				sam_fs<<"16"<<"\t";
			sam_fs<<bwt_vector[i].chrom<<"\t";
			sam_fs<<bwt_vector[i].start+1<<"\t";
			sam_fs<<quality_score<<"\t";
			sam_fs<<segment_length*segment_number<<"M";
			sam_fs<<"\t*\t0\t0\t";
			if(bwt_vector[i].strand=="+")  //out put + sequence
			{
				for(int i=1;i<=segment_number;i++)
				{
					sam_fs<<r_block.get_seg_seq((int)i);
				}
				sam_fs<<"\t";
			}
			else                   //out put - sequence
			{
				for(int i=segment_number;i>=1;i--)
				{
					sam_fs<<r_block.get_revcom_seg_seq((int)i);
				}
				sam_fs<<"\t";
			}
			if(!fa)    //fq
				sam_fs<<quality_string<<"\t";
			else       //fa
				sam_fs<<global_qual_string<<"\t";
			sam_fs<<"NM:i:"<<bwt_vector[i].mismatches.size()<<"\t";	
			sam_fs<<endl;
	}*/

	void output_unspliced_read()
	{
		for(size_t i=0;i<bwt_vector.size();i++)
		{
			if(bwt_vector[i].splice_head.size()==0 && bwt_vector[i].splice_tail.size()==0 && bwt_vector[i].splice_internal.size()==0)
			{
				if(abs(bwt_vector[i].end_seg_no- bwt_vector[i].start_seg_no)+1 ==segment_number)
					read_unspliced=true;
				Sam_Info new_sam_info(bwt_vector[i],i);	
				for(size_t j=0;j<new_sam_info.mismatchs.size();j++)
				{
					if(bwt_vector[i].strand=="+")
						new_sam_info.mismatchs[j].pos=(new_sam_info.mismatchs[j].seg_no-1)*segment_length + new_sam_info.mismatchs[j].pos;
					else
						new_sam_info.mismatchs[j].pos=(segment_number - new_sam_info.mismatchs[j].seg_no)*segment_length + new_sam_info.mismatchs[j].pos;
				}
				print_saminfo(new_sam_info);	
			}
		}
	}

	void report(char* bwtmap_filename, char* reads_filename, bool fasta)
	{
		fa = fasta;
		Bwtmap_File bwt_file(bwtmap_filename);
		Reads_File reads_file(reads_filename, fasta);
		cout<<"reporting sam result"<<endl;
		while(bwt_file.load_next_block(bwt_vector)>0)
		{
			while(reads_file.load_next_read(r_block)!=0)   // load read sequence
			{
				if(bwt_vector[0].read_id==r_block.read_id)
					break;
			}
			for(size_t i=0;i<bwt_vector.size();i++)
			{
				if(bwt_vector[i].strand == "-")
				{
					r_block.get_revcom_info(fasta);
					break;
				}
			}
			read_aligned = false;
			read_unspliced = false;
			sort(bwt_vector.begin(),bwt_vector.end(),bwtmap_compare_stable);
			combine_contig();
			for(size_t i=0;i<chrom_map.size();i++)
			{
				load_splice_info((*fixed_hole_vec[i]), bwt_vector[0].read_id, HOLE);
				load_splice_info((*fixed_hmer_vec[i]), bwt_vector[0].read_id, HMER);
				load_splice_info((*fixed_head_vec[i]), bwt_vector[0].read_id, HEAD);
				load_splice_info((*fixed_tail_vec[i]), bwt_vector[0].read_id, TAIL);
			}
			for(size_t i=0;i<bwt_vector.size();i++)
			{
				if(bwt_vector[i].head==-1&&(bwt_vector[i].splice_internal.size()!=0||bwt_vector[i].splice_head.size()!=0||bwt_vector[i].splice_tail.size()!=0))
					generate_sam_head(i);
			}
			output_unspliced_read();
			if(Fusion&&(!read_aligned)&&(!read_unspliced))
				search_fusion();
			fusion_vec.clear();
		}
		bwt_file.close();
		reads_file.close();
		cout<<"report finished"<<endl;
		close_all();
		if(Fusion)
			close_fusion();
	}

	void pair_2_read(int i, int j)
	{
		pair_vec[i].pair_index=j;
		pair_vec[j].pair_index=i;
		read_paired=true;
	}

	int can_be_paired_rightside(size_t i)
	{
		for(size_t j=i+1;j<pair_vec.size();j++)
		{
			if(pair_vec[j].pair_index!=-1)
				continue;
			if(pair_vec[i].start_pos==pair_vec[j].start_pos)
				continue;
			if(bwt_vector[pair_vec[i].start_contig_index].chrom!=bwt_vector[pair_vec[j].start_contig_index].chrom)
				return -1;
			if(bwt_vector[pair_vec[i].start_contig_index].pair_no==bwt_vector[pair_vec[j].start_contig_index].pair_no&& bwt_vector[pair_vec[i].start_contig_index].strand==bwt_vector[pair_vec[j].start_contig_index].strand)
				return -1;
			int i_j_gap= pair_vec[j].start_pos-pair_vec[i].start_pos;
			if(i_j_gap > max_intron_thresh)
				return -1;
			if(bwt_vector[pair_vec[i].start_contig_index].pair_no==bwt_vector[pair_vec[j].start_contig_index].pair_no|| bwt_vector[pair_vec[i].start_contig_index].strand==bwt_vector[pair_vec[j].start_contig_index].strand)
				continue;
			int j_pair_index=can_be_paired_rightside(j);
			if(j_pair_index ==-1) 
				return (int)j;
			else if(i_j_gap< pair_vec[j_pair_index].start_pos- pair_vec[j].start_pos)
				return (int)j;
			else
			{
				pair_2_read((int)j, j_pair_index);
				return -1;
			}
		}
		return -1;
	}


	void pair_read()
	{
		//cout<<"pairing alignment"<<endl;
		for(size_t i=0;i<pair_vec.size();i++)
		{
			if(pair_vec[i].pair_index==-1)
			{
				int i_pair_index=can_be_paired_rightside(i);
				if(i_pair_index!=-1)
				{
					pair_2_read((int)i, i_pair_index);
				}
			}
		}
		//cout<<"pairing alignment complete"<<endl;
	}

	void report_pairend(char* bwtmap_filename, char* reads_filename, bool fasta)
	{
		fa = fasta;
		Bwtmap_File bwt_file(bwtmap_filename, true);
		Reads_File reads_file(reads_filename, fasta, true);
		cout<<"reporting sam result"<<endl;
		while(bwt_file.load_next_block(bwt_vector)>0)
		{
			//cout<<bwt_vector[0].unmapped_index<<endl;
			while(reads_file.load_next_read(r_block)!=0)   // load read sequence
			{
				if(bwt_vector[0].read_id==r_block.read_id)
					break;
			}
			for(size_t i=0;i<bwt_vector.size();i++)
			{
				if(bwt_vector[i].strand == "-")
				{
					r_block.get_revcom_info(fasta);
					break;
				}
			}
			read_aligned = false;
			read_unspliced = false;
			read_paired =false;
			sort(bwt_vector.begin(),bwt_vector.end(),bwtmap_compare_pe_stable);
			combine_contig();
			for(size_t i=0;i<chrom_map.size();i++)
			{
				load_splice_info((*fixed_hole_vec[i]), bwt_vector[0].read_id, HOLE);
				load_splice_info((*fixed_hmer_vec[i]), bwt_vector[0].read_id, HMER);
				load_splice_info((*fixed_head_vec[i]), bwt_vector[0].read_id, HEAD);
				load_splice_info((*fixed_tail_vec[i]), bwt_vector[0].read_id, TAIL);
			}
			for(size_t i=0;i<bwt_vector.size();i++)
			{
				if(bwt_vector[i].head==-1&&(bwt_vector[i].splice_internal.size()!=0||bwt_vector[i].splice_head.size()!=0||bwt_vector[i].splice_tail.size()!=0))
					generate_sam_head(i);
			}
			output_unspliced_read();
			sort(pair_vec.begin(), pair_vec.end(), spliced_sam_compare_stable);
			pair_read();
			for(size_t i=0;i< pair_vec.size();i++)
			{
				print_sam_info_pe(pair_vec[i]);
			}
			pair_vec.clear();
			if(Fusion&&(!read_aligned)&&(!read_unspliced))
				search_fusion();
			fusion_vec.clear();
		}
		bwt_file.close();
		reads_file.close();
		cout<<"report finished"<<endl;
		close_all();
		if(Fusion)
			close_fusion();
	}

	void close_all()
	{
		for(size_t i=0;i<bwt_vector.size();i++)   // close all splice info file
		{
			(*fixed_hole_vec[i]).close();
			delete fixed_hole_vec[i];
			(*fixed_hmer_vec[i]).close();
			delete fixed_hmer_vec[i];
			(*fixed_head_vec[i]).close();
			delete fixed_head_vec[i];
			(*fixed_tail_vec[i]).close();
			delete fixed_tail_vec[i];
		}
		sam_fs.close();             // close sam output file
		if(pair_end)
			sam_pe_fs.close();
		if(output_filtered)
			sam_fo_fs.close();
	}

	void close_fusion()
	{
		map<string, size_t>::iterator it1, it2;                        // close all fusion output file
		for(it1=chrom_map.begin();it1!=chrom_map.end();it1++)
		{
			for(it2=chrom_map.begin();it2!=chrom_map.end();it2++)
			{
				(*fusion_fs_vec[it1->first + it2->first]).close();
				delete fusion_fs_vec[it1->first + it2->first];
			}
		}
	}
};

void print_usage()
{
	fprintf(stderr,"-l:        segment length, must specify\n");
	fprintf(stderr,"-n:        totoal segment number, must specify\n");
	fprintf(stderr,"-e:        fix hole extend length, default [0]\n");
	fprintf(stderr,"-c:        chromosome length table file name\n");
	fprintf(stderr,"-i1:       bwtmap file\n");
	fprintf(stderr,"-i2:       hole fixed path\n");
	fprintf(stderr,"-i3:       hmer fixed path\n");
	fprintf(stderr,"-i4:       head fixed path\n");
	fprintf(stderr,"-i5:       tail fixed path\n");
	fprintf(stderr,"-r:        reads file name\n");
	fprintf(stderr,"-fa/-fq:   reads file format, default [fa]\n");
	fprintf(stderr,"-o:        sam output file / for pair end read, this specifies the single alignment output\n");
	fprintf(stderr,"-o2:       for pair end read, this specifies the paired alignment output\n");
	fprintf(stderr,"-o3:       for pair end read, this specifies the filterd out single alignment output\n");
	fprintf(stderr,"-I:        max intron size threshold\n");
	fprintf(stderr,"-s:        only output reads have no less than [] number of segment connected\n");
	fprintf(stderr,"-f:        fusion search output file\n");
	fprintf(stderr,"-v:        fusion search minimum mapped segment threshold\n");
	exit(1);
}


int main(int argc, char** argv)
{
	if(argc==1)
		print_usage();
	int seg_len=0;
	int seg_num=0;
	int extend_len=0;
	int trunc_len=0;
	int max_intron_thresh=50000;
	char* chrom_file=NULL;
	char* bwtmap_file=NULL;
	char* hole_fixed_path=NULL;
	char* hmer_fixed_path=NULL;
	char* head_fixed_path=NULL;
	char* tail_fixed_path=NULL;
	char* sam_file=NULL;
	char* sam_pe_file=NULL;
	char* sam_fo_file=NULL;
	char* reads_file=NULL;
	bool fa=true;
	char* fusion_file=NULL;
	int seg_thresh=0;
	int min_fusion_seg=0;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-c")==0)
			chrom_file=argv[++i];
		else if(strcmp(argv[i],"-fa")==0)
			fa=true;
		else if(strcmp(argv[i],"-fq")==0)
			fa=false;
		else if(strcmp(argv[i],"-i1")==0)
			bwtmap_file=argv[++i];
		else if(strcmp(argv[i],"-i2")==0)
			hole_fixed_path=argv[++i];
		else if(strcmp(argv[i],"-i3")==0)
			hmer_fixed_path=argv[++i];
		else if(strcmp(argv[i],"-i4")==0)
			head_fixed_path=argv[++i];
		else if(strcmp(argv[i],"-i5")==0)
			tail_fixed_path=argv[++i];
		else if(strcmp(argv[i],"-o")==0)
			sam_file=argv[++i];
		else if(strcmp(argv[i],"-o2")==0)
			sam_pe_file=argv[++i];
		else if(strcmp(argv[i],"-o3")==0)
			sam_fo_file=argv[++i];
		else if(strcmp(argv[i],"-l")==0)
			seg_len=atoi(argv[++i]);
		else if(strcmp(argv[i],"-n")==0)
			seg_num=atoi(argv[++i]);
		else if(strcmp(argv[i],"-e")==0)
			extend_len=atoi(argv[++i]);
		else if(strcmp(argv[i],"-v")==0)
			min_fusion_seg=atoi(argv[++i]);
		else if(strcmp(argv[i],"-I")==0)
			max_intron_thresh=atoi(argv[++i]);
		else if(strcmp(argv[i],"-s")==0)
			seg_thresh=atoi(argv[++i]);
		else if(strcmp(argv[i],"-r")==0)
			reads_file=argv[++i];
		else if(strcmp(argv[i],"-f")==0)
			fusion_file=argv[++i];
		else 
			continue;
	}
	if(chrom_file==NULL)
	{fprintf(stderr,"error: chromosome length table file not specified\n");exit(-1);}
	if(bwtmap_file==NULL)
	{fprintf(stderr,"error: bwtmap input file not specified\n");		exit(-1);}
	if(hole_fixed_path==NULL)
	{fprintf(stderr,"error: fixed hole input path not specified\n");	exit(-1);}
	if(hmer_fixed_path==NULL)
	{fprintf(stderr,"error: fixed hmer input path not specified\n");	exit(-1);}
	if(head_fixed_path==NULL)
	{fprintf(stderr,"error: fixed head input path not specified\n");	exit(-1);}
	if(tail_fixed_path==NULL)
	{fprintf(stderr,"error: fixed tail input path not specified\n");	exit(-1);}
	if(sam_file==NULL)
	{fprintf(stderr,"error: sam output file not specified\n");	exit(-1);}
	if(seg_len==0)
	{fprintf(stderr,"error: segment length not specified\n");	exit(-1);}
	if(seg_num==0)
	{fprintf(stderr,"error: number of segment not specified\n");	exit(-1);}
	if(reads_file==NULL)
	{fprintf(stderr,"error: reads file not specified\n");	exit(-1);}
	if(seg_thresh==0)
		seg_thresh=seg_num;
	if(min_fusion_seg==0)
		min_fusion_seg=seg_num-1;


	mapsplice_report mr(seg_len,seg_num,extend_len, max_intron_thresh, seg_thresh, min_fusion_seg);
	mr.load_chrom_table(chrom_file, hole_fixed_path, hmer_fixed_path, head_fixed_path, tail_fixed_path, sam_file, fusion_file, sam_pe_file, sam_fo_file);
	if(sam_pe_file==NULL)
		mr.report(bwtmap_file,reads_file, fa);
	else
		mr.report_pairend(bwtmap_file,reads_file, fa);

	/*mapsplice_report mr(25,2,3,50000,4,3);
	mr.load_chrom_table("G://pe_report//chrom_sizes", "G://pe_report//hole//","G://pe_report//hmer//","G://pe_report//head//","G://pe_report//tail//","G://pe_report//accepted_hits.sam",NULL,NULL);
	mr.report("G://pe_report//unspliced_mapped_segments.sorted","G://pe_report//divided_reads.fa" ,false);
	*/

	/*mapsplice_report mr(25,2,3,50000,2,2);
	mr.load_chrom_table("H://pe_report//chrom_sizes", "H://pe_report//hole//","H://pe_report//hmer//","H://pe_report//head//","H://pe_report//tail//","H://pe_report//single.sam",NULL, "H://pe_report//paired.sam", "H://pe_report//filtered.sam");
	mr.report_pairend("H://pe_report//unspliced_mapped_segments.sorted","H://pe_report//divided_reads.fa" ,false);*/
	
}



