#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include "bwtmap_info.h"
#include "bwtmap_file.h"
#include "reads_file.h"

const int MAX_JUMP=3;

inline bool bwtmap_compare(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort compare function
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

inline bool bwtmap_compare_stable(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort compare function
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

inline bool bwtmap_compare_pe_stable(const Bwtmap_Info_Pe& seg1, const Bwtmap_Info_Pe& seg2)   //sort compare function
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

class search_hole
{
public: 
	vector<Bwtmap_Info> bwt_vector;
	vector<Bwtmap_Info_Pe> bwt_vector_pe;
	Read_Block r_block;
	int segment_length;
	int segment_number;
	int pos_threshold;
	int min_threshold;
	int trunc_len;
	int trunc_left_len;
	int extend_len;
	int extend_left_len;
	map<string, size_t> chrom_map;
	map<string, size_t> chrom_len_map;
	vector< vector<bool> > hits_vec;
	vector<ofstream*> hole_fs_vec;
	vector<ofstream*> hmer_fs_vec;
	vector<ofstream*> head_fs_vec;
	vector<ofstream*> tail_fs_vec;
	map< pair<int, string>, pair<size_t, size_t> > pe_index_map;
	size_t index_shift;

	search_hole(int seg_len, int seg_num, int min_thr, int pos_thr, int ext_len)
	{
		segment_length = seg_len;
		trunc_len=seg_len / 2;
		trunc_left_len=seg_len-trunc_len;
		segment_number = seg_num;
		pos_threshold = pos_thr;
		min_threshold=min_thr;
		extend_len=ext_len;
		extend_left_len=seg_len-extend_len;
	}

	~search_hole()
	{}

	void load_chrom_table(char* chrom_len_file, char* hole_path, char* hmer_path, char* head_path, char* tail_path)
	{
		cout<<"allocating memory"<<endl;
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
			char chrom_name_tmp[1000];
			size_t chrom_bp_len;
			sscanf(line.c_str(), "%s\t%zd",chrom_name_tmp, &chrom_bp_len);  //zd for linux
			string chrom_name=chrom_name_tmp;
			chrom_len_map.insert(make_pair(chrom_name,chrom_bp_len));
			chrom_map.insert(make_pair(chrom_name,i));
			vector<bool> new_hits_vec(chrom_bp_len,false);
			hits_vec.push_back(new_hits_vec);
			string hole_file_name=hole_path+chrom_name+".txt";
			string hmer_file_name=hmer_path+chrom_name+".txt";
			string head_file_name=head_path+chrom_name+".txt";
			string tail_file_name=tail_path+chrom_name+".txt";

			hole_fs_vec.push_back(new ofstream(hole_file_name.c_str()));
			hmer_fs_vec.push_back(new ofstream(hmer_file_name.c_str()));
			head_fs_vec.push_back(new ofstream(head_file_name.c_str()));
			tail_fs_vec.push_back(new ofstream(tail_file_name.c_str()));
			i++;
		}
		chrom_len_fs.close();
	}

	int check_contig(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)
	{
		if (seg1.chrom!=seg2.chrom || seg1.strand!=seg2.strand || seg2.start - seg1.end > 1)
			return 0;            // case 0: no relationship
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

	int check_contig(const Bwtmap_Info_Pe& seg1, const Bwtmap_Info_Pe& seg2)
	{
		if (seg1.chrom!=seg2.chrom || seg1.strand!=seg2.strand || seg1.pair_no!=seg2.pair_no ||seg2.start - seg1.end > 1)
			return 0;            // case 0: no relationship
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

	void combine_contig()   
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

	void combine_contig_pe()   
	{
		for (vector<Bwtmap_Info_Pe>::iterator it=bwt_vector_pe.begin();it+1!=bwt_vector_pe.end();)
		{
			for(vector<Bwtmap_Info_Pe>::iterator it2=it+1;;it2++)
			{
				if(it2==bwt_vector_pe.end())
				{
					it++;
					break;
				}
				int contig_status= check_contig(*it, *(it2));
				if (contig_status == 1)   //contig
				{
					(*it).end_seg_no = (*it2).end_seg_no;
					(*it).end = (*it2).end;
					bwt_vector_pe.erase(it2);
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

	void creat_pe_index_map()
	{
		index_shift=0;
		string chrom="";
		string strand="";
		int pair_no=-1;
		bool new_key=false;
		size_t start_index=0;
		size_t end_index=0;
		for(size_t i=0;i<bwt_vector_pe.size();i++)
		{
			if(bwt_vector_pe[i].chrom != chrom || bwt_vector_pe[i].strand != strand || bwt_vector_pe[i].pair_no != pair_no)
			{
				if(bwt_vector_pe[i].pair_no != pair_no)
					index_shift=i;
				if(!new_key)   //first time
				{
					chrom=bwt_vector_pe[i].chrom;
					strand=bwt_vector_pe[i].strand;
					pair_no=bwt_vector_pe[i].pair_no;
					new_key=true;
					start_index=i;
				}
				else
				{
					end_index=i-1;
					pair<int, string> key(pair_no, chrom+strand);
					pair<size_t, size_t> value(start_index, end_index);
					pe_index_map[key]=value;
					chrom=bwt_vector_pe[i].chrom;
					strand=bwt_vector_pe[i].strand;
					pair_no=bwt_vector_pe[i].pair_no;
					start_index=i;
				}
			}
		}
		end_index=bwt_vector_pe.size()-1;
		pair<int, string> key(pair_no, chrom+strand);
		pair<size_t, size_t> value(start_index, end_index);
		pe_index_map[key]=value;
	}

	void set_hits_bit()
	{
		for (size_t i=0; i < bwt_vector.size(); i++)
		{       
			hits_vec[chrom_map[bwt_vector[i].chrom]][bwt_vector[i].start] = true;
		}
	}

	void set_hits_bit_pe()
	{
		for (size_t i=0; i < bwt_vector_pe.size(); i++)
		{       
			hits_vec[chrom_map[bwt_vector_pe[i].chrom]][bwt_vector_pe[i].start] = true;
		}
	}

	void set_whole_mapping_bit(char* whole_bwtmap_file)
	{
		cout<<"setting bit for whole mapping reads"<<endl;
		ifstream wbwt_fs(whole_bwtmap_file);
		string line;
		while(getline(wbwt_fs, line))
		{
			char chrom_tmp[MAX_CHROM_LENGTH];
			int start;
			sscanf(line.c_str(), "%*s\t%*s\t%s\t%d",chrom_tmp,&start);	
			string chrom=chrom_tmp;
			for(int i=0;i<segment_number;i++)
			{
				hits_vec[chrom_map[chrom]][start+i*segment_length] = true;
			}
		}
		wbwt_fs.close();
	}

	void detect_island(int gap_thresh, int ext_len, char* island_file)
	{
		cout<<"detecting island"<<endl;
		ofstream island_fs(island_file);
		island_fs<<"islands"<<endl;
		map<string, size_t>::iterator it;
		for(it=chrom_map.begin();it!=chrom_map.end();it++)
		{
			bool new_island = true;
			int island_start=0;
			int island_end=0;
			int size=(int)hits_vec[it->second].size();
			int total_num = 0;
			for(int i=0;i<size;i++)
			{
				if(hits_vec[it->second][i])
				{
					if(new_island)
					{
						island_start=max(0,i-ext_len);
						island_end=i + segment_length - 1 + ext_len;
						new_island=false;
					}
					else
					{
						if(i>island_end+gap_thresh)
						{
							total_num++;
							island_fs<<it->first<<"\t"<<"MapSplice"<<"\t"<<"island"<<"\t"<<island_start+1<<"\t"<<island_end+1<<"\t"<<"1"<<"\t"<<"."<<"\t"<<"."<<"\t"<<"IL"<<total_num<<endl;
							island_start=max(0,i-ext_len);
							island_end=i + segment_length - 1 + ext_len;
						}
						else
						{
							island_end=i+ segment_length -1 +ext_len;
						}
					}
				}
			}
			if(!new_island)
			{
				total_num++;
				island_fs<<it->first<<"\t"<<"MapSplice"<<"\t"<<"island"<<"\t"<<island_start+1<<"\t"<<island_end+1<<"\t"<<"1"<<"\t"<<"."<<"\t"<<"."<<"\t"<<"IL"<<total_num<<endl;
			}
		}
		island_fs.close();
	}

	int check_relation(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)
	{
		int min_gap=(abs(seg2.start_seg_no - seg1.end_seg_no) - 1) * segment_length + min_threshold;
		if (seg1.chrom==seg2.chrom && seg1.strand==seg2.strand)
		{
			int dist=seg2.start - seg1.end;
			if( dist< min_gap)
				return 1;       // case 1: too close
			else if(dist <= pos_threshold)
				return 2;        // case 2: in range
			else 
				return 3;        // case 3: too far
		}
		else
			return 0;            // case 0: no relationship
	}

	int check_relation(const Bwtmap_Info_Pe& seg1, const Bwtmap_Info_Pe& seg2)
	{
		int min_gap=(abs(seg2.start_seg_no - seg1.end_seg_no) - 1) * segment_length + min_threshold;
		if (seg1.chrom==seg2.chrom && seg1.strand==seg2.strand && seg1.pair_no==seg2.pair_no)
		{
			int dist=seg2.start - seg1.end;
			if( dist< min_gap)
				return 1;       // case 1: too close
			else if(dist <= pos_threshold)
				return 2;        // case 2: in range
			else 
				return 3;        // case 3: too far
		}
		else
			return 0;            // case 0: no relationship
	}

	void search(char* bwtmap_filename, char* reads_filename, bool fasta)
	{
		cout<<"searching hole and hmer"<<endl;
		Bwtmap_File bwt_file(bwtmap_filename);
		Reads_File reads_file(reads_filename, fasta);
		while(bwt_file.load_next_block(bwt_vector)>0)
		{
			while(reads_file.load_next_read(r_block))   // load read sequence
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
			set_hits_bit();
			sort(bwt_vector.begin(),bwt_vector.end(),bwtmap_compare_stable);
			combine_contig();
			for (size_t i = 0; i < bwt_vector.size(); i++)
			{
				for (size_t j = i + 1; j < bwt_vector.size(); j++)
				{
					int relation = check_relation(bwt_vector[i], bwt_vector[j]);
					if(relation == 0)            // not related
						break;
					else if(relation == 1)            // continue search
						continue;
					else  if(relation==2)// in range
					{ 
						if (bwt_vector[i].check_tail(bwt_vector[j].start_seg_no)&&bwt_vector[i].strand =="+")   //strand + 
						{
							int seg_gap=bwt_vector[j].start_seg_no - bwt_vector[i].end_seg_no;
							if(seg_gap==1)	  //continuous hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end - trunc_len + 1<<"\t"<<bwt_vector[j].start+trunc_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no).substr(0,trunc_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[j].start_seg_no).substr(trunc_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no).substr(trunc_left_len)<<r_block.get_seg_seq(bwt_vector[j].start_seg_no).substr(0,trunc_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap== 2)  //gap hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end - extend_len + 1<<"\t"<<bwt_vector[j].start+extend_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no).substr(0,extend_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[j].start_seg_no).substr(extend_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no).substr(extend_left_len)<<r_block.get_seg_seq(bwt_vector[i].end_seg_no+1)<<r_block.get_seg_seq(bwt_vector[j].start_seg_no).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap==0&&bwt_vector[i].start_seg_no!=bwt_vector[i].end_seg_no&&bwt_vector[j].start_seg_no!=bwt_vector[j].end_seg_no)
							{/////// special case
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end -segment_length - extend_len + 1<<"\t"<<bwt_vector[j].start + segment_length + extend_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no-1<<"\t"<<bwt_vector[j].start_seg_no+1<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no-1).substr(0,extend_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[j].start_seg_no+1).substr(extend_len)<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no-1).substr(extend_left_len)<<r_block.get_seg_seq(bwt_vector[i].end_seg_no)<<r_block.get_seg_seq(bwt_vector[j].start_seg_no+1).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap>2&&seg_gap<=MAX_JUMP) //hmer
							{
								(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end + 1<<"\t"<<bwt_vector[j].start -1 <<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no)<<"\t"<<r_block.get_seg_seq(bwt_vector[j].start_seg_no)<<"\t";
								for(int num_jump=1; num_jump<seg_gap;num_jump++)
								{
									(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_seg_seq(bwt_vector[i].end_seg_no+num_jump)<<",";
								}
								(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<i<<"\t"<<j<<endl;
							}
							else
								continue;
							if(bwt_vector[i].end_seg_no!=bwt_vector[j].start_seg_no)  //not special case
							{
								bwt_vector[i].tail=bwt_vector[j].start_seg_no;
								bwt_vector[j].head=bwt_vector[i].end_seg_no;
							}
							else     //special case
							{
								bwt_vector[i].tail=bwt_vector[j].start_seg_no + 1;
								bwt_vector[j].head=bwt_vector[i].end_seg_no - 1;
							}
						}
						else if(bwt_vector[i].check_tail(bwt_vector[j].start_seg_no)&&bwt_vector[i].strand =="-")  //strand - hole
						{
							int seg_gap=bwt_vector[i].end_seg_no - bwt_vector[j].start_seg_no;
							if( seg_gap== 1)  //continuous hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end - trunc_len + 1<<"\t"<<bwt_vector[j].start+trunc_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no).substr(0,trunc_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no).substr(trunc_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no).substr(trunc_left_len)<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no).substr(0,trunc_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap== 2)  //gap hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end - extend_len + 1<<"\t"<<bwt_vector[j].start+extend_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no).substr(0,extend_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no).substr(extend_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no).substr(extend_left_len)<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no-1)<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap==0&&bwt_vector[i].start_seg_no!=bwt_vector[i].end_seg_no&&bwt_vector[j].start_seg_no!=bwt_vector[j].end_seg_no)
							{/////// special case
								(*hole_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end -segment_length - extend_len + 1<<"\t"<<bwt_vector[j].start + segment_length + extend_len-1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no+1<<"\t"<<bwt_vector[j].start_seg_no-1<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no+1).substr(0,extend_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no-1).substr(extend_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no+1).substr(extend_left_len)<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no)<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no-1).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap>2&&seg_gap<=MAX_JUMP) //hmer
							{
								(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end + 1<<"\t"<<bwt_vector[j].start -1<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<bwt_vector[j].start_seg_no<<"-"<<bwt_vector[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[j].start_seg_no)<<"\t";
								for(int num_jump=1; num_jump<seg_gap;num_jump++)
								{
									(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no-num_jump)<<",";
								}
								(*hmer_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<i<<"\t"<<j<<endl;
							}
							else 
								continue;
							if(bwt_vector[i].end_seg_no!=bwt_vector[j].start_seg_no)   //not special case
							{
								bwt_vector[i].tail=bwt_vector[j].start_seg_no;
								bwt_vector[j].head=bwt_vector[i].end_seg_no;
							}
							else     //special case
							{
								bwt_vector[i].tail=bwt_vector[j].start_seg_no - 1;
								bwt_vector[j].head=bwt_vector[i].end_seg_no + 1;
							}
						}	
					}
					/*else  if(relation==3)// for possible fusion junction, do not fix head and tail
					{ 
						int seg_gap=0;
						if (bwt_vector[i].check_tail(bwt_vector[j].start_seg_no)&&bwt_vector[i].strand =="+")   //strand + 
						{
							seg_gap=bwt_vector[j].start_seg_no - bwt_vector[i].end_seg_no;
						}
						else if(bwt_vector[i].check_tail(bwt_vector[j].start_seg_no)&&bwt_vector[i].strand =="-")  //strand - hole
						{
							seg_gap=bwt_vector[i].end_seg_no - bwt_vector[j].start_seg_no;
						}
						if( seg_gap== 1||seg_gap==2)
						{
							bwt_vector[i].tail=bwt_vector[j].start_seg_no;
							bwt_vector[j].head=bwt_vector[i].end_seg_no;
						}
					}*/
				}
			}
			for (unsigned int i = 0; i < bwt_vector.size(); i++)  //jump hmer
			{
				if(bwt_vector[i].strand=="+")
				{
					if (bwt_vector[i].head==-1 && bwt_vector[i].start_seg_no!=1 && (bwt_vector[i].start_seg_no-1 <= MAX_JUMP))
					{//head fix hmer	
						(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<"0"<<"\t"<<bwt_vector[i].start -1<<"\t"<<"0-0"<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<"-"<<"\t"<<r_block.get_seg_seq(bwt_vector[i].start_seg_no)<<"\t";
						for(int num_jump=1; num_jump<bwt_vector[i].start_seg_no;num_jump++)
						{
							(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_seg_seq(bwt_vector[i].start_seg_no - num_jump)<<",";
						}
						(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<"-1"<<"\t"<<i<<endl;
					}
					if (bwt_vector[i].tail==-1 && bwt_vector[i].end_seg_no!= segment_number&& (segment_number - bwt_vector[i].end_seg_no <= MAX_JUMP))
					{//tail fix hmer
						(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end +1 <<"\t"<<"0"<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<"0-0"<<"\t"<<r_block.get_seg_seq(bwt_vector[i].end_seg_no)<<"\t"<<"-"<<"\t";
						for(int num_jump=1; num_jump<=segment_number-bwt_vector[i].end_seg_no;num_jump++)
						{
							(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_seg_seq(bwt_vector[i].end_seg_no + num_jump)<<",";
						}
						(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<i<<"\t"<<"-1"<<endl;
					}
				}
				else if(bwt_vector[i].strand=="-")
				{
					if (bwt_vector[i].head==-1 && bwt_vector[i].start_seg_no!=segment_number&& (segment_number - bwt_vector[i].start_seg_no <= MAX_JUMP))
					{//head fix hmer
						(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<"0"<<"\t"<<bwt_vector[i].start -1 <<"\t"<<"0-0"<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<"-"<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].start_seg_no)<<"\t";
						for(int num_jump=1; num_jump<=segment_number-bwt_vector[i].start_seg_no;num_jump++)
						{
							(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector[i].start_seg_no + num_jump)<<",";
						}
						(*head_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<"-1"<<"\t"<<i<<endl;
					}
					if (bwt_vector[i].tail==-1 && bwt_vector[i].end_seg_no!=1 && (bwt_vector[i].end_seg_no-1 <= MAX_JUMP)) 
					{//tail fix hmer
						(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<bwt_vector[i].unmapped_index<<"~"<<bwt_vector[i].read_id<<"\t"<<bwt_vector[i].chrom<<"\t"<<bwt_vector[i].strand<<"\t"<<bwt_vector[i].end + 1<<"\t"<<"0"<<"\t"<<bwt_vector[i].start_seg_no<<"-"<<bwt_vector[i].end_seg_no<<"\t"<<"0-0"<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no)<<"\t"<<"-"<<"\t";
						for(int num_jump=1; num_jump<bwt_vector[i].end_seg_no;num_jump++)
						{
							(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector[i].end_seg_no - num_jump)<<",";
						}
						(*tail_fs_vec[chrom_map[bwt_vector[i].chrom]])<<"\t"<<i<<"\t"<<"-1"<<endl;
					}
				}
			}
		}
		bwt_file.close();
		reads_file.close();
	}

	int check_pe_boundary(size_t i, int type)  //type 1 == fix head, type 2 == fix tail
	{
		int boundary_ret= -1;
		int search_index;
		int end_index;
		int des_pair_no;
		if(bwt_vector_pe[i].pair_no==1)
			des_pair_no=2;
		else
			des_pair_no=1;
		if(bwt_vector_pe[i].strand=="+"&& type ==1)
		{
			pair<int, string> key(des_pair_no,bwt_vector_pe[i].chrom+"-");
			if(pe_index_map.find(key)!=pe_index_map.end())
			{
				search_index=(int)pe_index_map[key].second;
				end_index=(int)pe_index_map[key].first;
				for(; search_index>= end_index; search_index--)
				{
					if(bwt_vector_pe[search_index].end<bwt_vector_pe[i].start)
					{
						boundary_ret=bwt_vector_pe[search_index].end+1;
						break;
					}
				}
			}
			if(boundary_ret==-1 || bwt_vector_pe[i].start - boundary_ret>pos_threshold)
				boundary_ret=max(bwt_vector_pe[i].start - pos_threshold, 0);
		}
		else if(bwt_vector_pe[i].strand=="+"&& type ==2)
		{
			pair<int, string> key(des_pair_no,bwt_vector_pe[i].chrom+"-");
			if(pe_index_map.find(key)!=pe_index_map.end())
			{
				search_index=(int)pe_index_map[key].first;
				end_index=(int)pe_index_map[key].second;
				for(; search_index<=end_index ; search_index++)
				{
					if(bwt_vector_pe[search_index].start > bwt_vector_pe[i].end)
					{
						boundary_ret=bwt_vector_pe[search_index].start-1;
						break;
					}
				}
			}
			if(boundary_ret==-1 || boundary_ret - bwt_vector_pe[i].end > pos_threshold)
			boundary_ret=min(bwt_vector_pe[i].end + pos_threshold, (int)(chrom_len_map[bwt_vector_pe[i].chrom]));
		}
		else if(bwt_vector_pe[i].strand=="-"&& type ==1)
		{
			pair<int, string> key(des_pair_no,bwt_vector_pe[i].chrom+"+");
			if(pe_index_map.find(key)!=pe_index_map.end())
			{
				search_index=(int)pe_index_map[key].second;
				end_index=(int)pe_index_map[key].first;
				for(; search_index>=end_index ; search_index--)
				{
					if(bwt_vector_pe[search_index].end<bwt_vector_pe[i].start)
					{
						boundary_ret=bwt_vector_pe[search_index].end+1;
						break;
					}
				}
			}
			if(boundary_ret==-1 || bwt_vector_pe[i].start - boundary_ret>pos_threshold)
				boundary_ret=max(bwt_vector_pe[i].start - pos_threshold, 0);
		}
		else if(bwt_vector_pe[i].strand=="-"&& type ==2)
		{
			pair<int, string> key(des_pair_no,bwt_vector_pe[i].chrom+"+");
			if(pe_index_map.find(key)!=pe_index_map.end())
			{
				search_index=(int)pe_index_map[key].first;
				end_index=(int)pe_index_map[key].second;
				for(; search_index<=end_index; search_index++)
				{
					if(bwt_vector_pe[search_index].start > bwt_vector_pe[i].end)
					{
						boundary_ret=bwt_vector_pe[search_index].start-1;
						break;
					}
				}
			}
			if(boundary_ret==-1 || boundary_ret - bwt_vector_pe[i].end > pos_threshold)
			boundary_ret=min(bwt_vector_pe[i].end + pos_threshold, (int)(chrom_len_map[bwt_vector_pe[i].chrom]));
		}
		return boundary_ret;
	}

	void search_pairend(char* bwtmap_filename, char* reads_filename, bool fasta)
	{
		cout<<"searching hole and hmer"<<endl;
		Bwtmap_File bwt_file(bwtmap_filename);
		Reads_File reads_file(reads_filename, fasta);
		while(bwt_file.load_next_block(bwt_vector_pe)>0)
		{
			while(reads_file.load_next_read_pe(r_block))   // load read sequence
			{
				if(bwt_vector_pe[0].read_id==r_block.read_id)
					break;
			}
			for(size_t i=0;i<bwt_vector_pe.size();i++)
			{
				if(bwt_vector_pe[i].strand == "-")
				{
					r_block.get_revcom_info(fasta);
					break;
				}
			}
			set_hits_bit_pe();
			sort(bwt_vector_pe.begin(),bwt_vector_pe.end(),bwtmap_compare_pe_stable);
			combine_contig_pe();
			pe_index_map.clear();
			creat_pe_index_map();
			for (size_t i = 0; i < bwt_vector_pe.size(); i++)
			{
				for (size_t j = i + 1; j < bwt_vector_pe.size(); j++)
				{
					int relation = check_relation(bwt_vector_pe[i], bwt_vector_pe[j]);
					if(relation == 0)            // not related
						break;
					else if(relation == 1)            // continue search
						continue;
					else  if(relation==2)// in range
					{ 
						if (bwt_vector_pe[i].check_tail(bwt_vector_pe[j].start_seg_no)&&bwt_vector_pe[i].strand =="+")   //strand + 
						{
							int seg_gap=bwt_vector_pe[j].start_seg_no - bwt_vector_pe[i].end_seg_no;
							if(seg_gap==1)	  //continuous hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end - trunc_len + 1<<"\t"<<bwt_vector_pe[j].start+trunc_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,trunc_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(trunc_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(trunc_left_len)<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,trunc_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap== 2)  //gap hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end - extend_len + 1<<"\t"<<bwt_vector_pe[j].start+extend_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,extend_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(extend_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(extend_left_len)<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no+1 + (bwt_vector_pe[i].pair_no-1)*segment_number)<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap==0&&bwt_vector_pe[i].start_seg_no!=bwt_vector_pe[i].end_seg_no&&bwt_vector_pe[j].start_seg_no!=bwt_vector_pe[j].end_seg_no)
							{/////// special case
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end -segment_length - extend_len + 1<<"\t"<<bwt_vector_pe[j].start + segment_length + extend_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no-1<<"\t"<<bwt_vector_pe[j].start_seg_no+1<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no-1 + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,extend_left_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no+1 + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(extend_len)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no-1 + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(extend_left_len)<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no+1 + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap>2&&seg_gap<=MAX_JUMP) //hmer
							{
								(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end + 1<<"\t"<<bwt_vector_pe[j].start -1 <<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number)<<"\t";
								for(int num_jump=1; num_jump<seg_gap;num_jump++)
								{
									(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no+num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
								}
								(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<i<<"\t"<<j<<endl;
							}
							else
								continue;
							if(bwt_vector_pe[i].end_seg_no!=bwt_vector_pe[j].start_seg_no)  //not special case
							{
								bwt_vector_pe[i].tail=bwt_vector_pe[j].start_seg_no;
								bwt_vector_pe[j].head=bwt_vector_pe[i].end_seg_no;
							}
							else     //special case
							{
								bwt_vector_pe[i].tail=bwt_vector_pe[j].start_seg_no + 1;
								bwt_vector_pe[j].head=bwt_vector_pe[i].end_seg_no - 1;
							}
						}
						else if(bwt_vector_pe[i].check_tail(bwt_vector_pe[j].start_seg_no)&&bwt_vector_pe[i].strand =="-")  //strand - hole
						{
							int seg_gap=bwt_vector_pe[i].end_seg_no - bwt_vector_pe[j].start_seg_no;
							if( seg_gap== 1)  //continuous hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end - trunc_len + 1<<"\t"<<bwt_vector_pe[j].start+trunc_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,trunc_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(trunc_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(trunc_left_len)<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,trunc_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap== 2)  //gap hole
							{
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end - extend_len + 1<<"\t"<<bwt_vector_pe[j].start+extend_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,extend_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(extend_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(extend_left_len)<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no-1 + (bwt_vector_pe[i].pair_no-1)*segment_number)<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap==0&&bwt_vector_pe[i].start_seg_no!=bwt_vector_pe[i].end_seg_no&&bwt_vector_pe[j].start_seg_no!=bwt_vector_pe[j].end_seg_no)
							{/////// special case
								(*hole_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end -segment_length - extend_len + 1<<"\t"<<bwt_vector_pe[j].start + segment_length + extend_len-1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no+1<<"\t"<<bwt_vector_pe[j].start_seg_no-1<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no+1 + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(0,extend_left_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no-1 + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(extend_len)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no+1 + (bwt_vector_pe[i].pair_no-1)*segment_number).substr(extend_left_len)<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no-1 + (bwt_vector_pe[j].pair_no-1)*segment_number).substr(0,extend_len)<<"\t"<<i<<"\t"<<j<<endl;
							}
							else if(seg_gap>2&&seg_gap<=MAX_JUMP) //hmer
							{
								(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end + 1<<"\t"<<bwt_vector_pe[j].start -1<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<bwt_vector_pe[j].start_seg_no<<"-"<<bwt_vector_pe[j].end_seg_no<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[j].start_seg_no + (bwt_vector_pe[j].pair_no-1)*segment_number)<<"\t";
								for(int num_jump=1; num_jump<seg_gap;num_jump++)
								{
									(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no-num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
								}
								(*hmer_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<i<<"\t"<<j<<endl;
							}
							else 
								continue;
							if(bwt_vector_pe[i].end_seg_no!=bwt_vector_pe[j].start_seg_no)   //not special case
							{
								bwt_vector_pe[i].tail=bwt_vector_pe[j].start_seg_no;
								bwt_vector_pe[j].head=bwt_vector_pe[i].end_seg_no;
							}
							else     //special case
							{
								bwt_vector_pe[i].tail=bwt_vector_pe[j].start_seg_no - 1;
								bwt_vector_pe[j].head=bwt_vector_pe[i].end_seg_no + 1;
							}
						}	
					}
					/*else  if(relation==3)// for possible fusion junction, do not fix head and tail
					{ 
						int seg_gap=0;
						if (bwt_vector_pe[i].check_tail(bwt_vector_pe[j].start_seg_no)&&bwt_vector_pe[i].strand =="+")   //strand + 
						{
							seg_gap=bwt_vector_pe[j].start_seg_no - bwt_vector_pe[i].end_seg_no;
						}
						else if(bwt_vector_pe[i].check_tail(bwt_vector_pe[j].start_seg_no)&&bwt_vector_pe[i].strand =="-")  //strand - hole
						{
							seg_gap=bwt_vector_pe[i].end_seg_no - bwt_vector_pe[j].start_seg_no;
						}
						if( seg_gap== 1||seg_gap==2)
						{
							bwt_vector_pe[i].tail=bwt_vector_pe[j].start_seg_no;
							bwt_vector_pe[j].head=bwt_vector_pe[i].end_seg_no;
						}
					}*/
				}
			}
			for (size_t i = 0; i < bwt_vector_pe.size(); i++)  //single anchored
			{
				if(bwt_vector_pe[i].strand=="+")
				{
					if (bwt_vector_pe[i].head==-1 && bwt_vector_pe[i].start_seg_no!=1 && (bwt_vector_pe[i].start_seg_no-1 <= MAX_JUMP))
					{//head fix hmer	
						int boundary=check_pe_boundary(i, 1);
						(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<boundary<<"\t"<<bwt_vector_pe[i].start -1<<"\t"<<"0-0"<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<"-"<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].start_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t";
						for(int num_jump=1; num_jump<bwt_vector_pe[i].start_seg_no;num_jump++)
						{
							(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_seg_seq(bwt_vector_pe[i].start_seg_no - num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
						}
						(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<"-1"<<"\t"<<i<<endl;
					}
					if (bwt_vector_pe[i].tail==-1 && bwt_vector_pe[i].end_seg_no!= segment_number&& (segment_number - bwt_vector_pe[i].end_seg_no <= MAX_JUMP))
					{//tail fix hmer
						int boundary=check_pe_boundary(i, 2);
						(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end +1 <<"\t"<<boundary<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<"0-0"<<"\t"<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t"<<"-"<<"\t";
						for(int num_jump=1; num_jump<=segment_number-bwt_vector_pe[i].end_seg_no;num_jump++)
						{
							(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_seg_seq(bwt_vector_pe[i].end_seg_no + num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
						}
						(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<i<<"\t"<<"-1"<<endl;
					}
				}
				else if(bwt_vector_pe[i].strand=="-")
				{
					if (bwt_vector_pe[i].head==-1 && bwt_vector_pe[i].start_seg_no!=segment_number&& (segment_number - bwt_vector_pe[i].start_seg_no <= MAX_JUMP))
					{//head fix hmer
						int boundary=check_pe_boundary(i, 1);
						(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<boundary<<"\t"<<bwt_vector_pe[i].start -1 <<"\t"<<"0-0"<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<"-"<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].start_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t";
						for(int num_jump=1; num_jump<=segment_number-bwt_vector_pe[i].start_seg_no;num_jump++)
						{
							(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].start_seg_no + num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
						}
						(*head_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<"-1"<<"\t"<<i<<endl;
					}
					if (bwt_vector_pe[i].tail==-1 && bwt_vector_pe[i].end_seg_no!=1 && (bwt_vector_pe[i].end_seg_no-1 <= MAX_JUMP)) 
					{//tail fix hmer
						int boundary=check_pe_boundary(i, 2);
						(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<bwt_vector_pe[i].unmapped_index<<"~"<<bwt_vector_pe[i].read_id<<"/"<<bwt_vector_pe[i].pair_no<<"\t"<<bwt_vector_pe[i].chrom<<"\t"<<bwt_vector_pe[i].strand<<"\t"<<bwt_vector_pe[i].end + 1<<"\t"<<boundary<<"\t"<<bwt_vector_pe[i].start_seg_no<<"-"<<bwt_vector_pe[i].end_seg_no<<"\t"<<"0-0"<<"\t"<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no + (bwt_vector_pe[i].pair_no-1)*segment_number)<<"\t"<<"-"<<"\t";
						for(int num_jump=1; num_jump<bwt_vector_pe[i].end_seg_no;num_jump++)
						{
							(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<r_block.get_revcom_seg_seq(bwt_vector_pe[i].end_seg_no - num_jump + (bwt_vector_pe[i].pair_no-1)*segment_number)<<",";
						}
						(*tail_fs_vec[chrom_map[bwt_vector_pe[i].chrom]])<<"\t"<<i<<"\t"<<"-1"<<endl;
					}
				}
			}
		}
		bwt_file.close();
		reads_file.close();
	}

	void close_all()
	{
		for(size_t i=0;i<hole_fs_vec.size();i++)
		{
			(*hole_fs_vec[i]).close();
			delete hole_fs_vec[i];
			(*hmer_fs_vec[i]).close();
			delete hmer_fs_vec[i];
			(*head_fs_vec[i]).close();
			delete head_fs_vec[i];
			(*tail_fs_vec[i]).close();
			delete tail_fs_vec[i];
		}
	}
};


void print_usage()
{
	fprintf(stderr,"[argument] [segment mapping file]\n");
	fprintf(stderr,"-i:        minimum gap threshold, default [0]\n");
	fprintf(stderr,"-I:        maximum gap threshold, default [50000]\n");
	fprintf(stderr,"-l:        segment length, must specify\n");
	fprintf(stderr,"-n:        totoal segment number, must specify\n");
	fprintf(stderr,"-e:        fix hole extend length, default [0]\n");
	fprintf(stderr,"-r:		   reads file name\n");
	fprintf(stderr,"-fa/-fq:   reads file format, default [fa]\n");
	fprintf(stderr,"-c:		   chromosome length table file name\n");
	fprintf(stderr,"-m:		   whole reads bowtie mapping file\n");
	fprintf(stderr,"-ie:       island extension length, default [45]\n");
	fprintf(stderr,"-g:        island gap threshold, default [50]\n");
	fprintf(stderr,"-o1:	   fix hole output path\n");
	fprintf(stderr,"-o2:	   fix hmer output path\n");
	fprintf(stderr,"-o3:	   fix head output path\n");
	fprintf(stderr,"-o4:	   fix tail output path\n");
	fprintf(stderr,"-o5:	   island output file name\n");
	fprintf(stderr,"-pe:	   input data is pair end reads\n");
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc==1)
		print_usage();
	int min_thresh=0;
	int max_thresh=50000;
	int seg_len=0;
	int seg_num=0;
	int extend_len=0;
	int island_ext_len=45;
	int island_gap_thresh=50;
	char* reads_file=NULL;
	char* chrom_len_file=NULL;
	char* hole_path=NULL;
	char* hmer_path=NULL;
	char* head_path=NULL;
	char* tail_path=NULL;
	char* island_file=NULL;
	char* whole_bwtmap_file=NULL;
	bool fa=true;
	bool pair_end=false;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-i")==0)
			min_thresh=atoi(argv[++i]);
		else if(strcmp(argv[i],"-I")==0)
			max_thresh=atoi(argv[++i]);
		else if(strcmp(argv[i],"-l")==0)
			seg_len=atoi(argv[++i]);
		else if(strcmp(argv[i],"-n")==0)
			seg_num=atoi(argv[++i]);
		else if(strcmp(argv[i],"-e")==0)
			extend_len=atoi(argv[++i]);
		else if(strcmp(argv[i],"-ie")==0)
			island_ext_len=atoi(argv[++i]);
		else if(strcmp(argv[i],"-g")==0)
			island_gap_thresh=atoi(argv[++i]);
		else if(strcmp(argv[i],"-r")==0)
			reads_file=argv[++i];
		else if(strcmp(argv[i],"-c")==0)
			chrom_len_file=argv[++i];
		else if(strcmp(argv[i],"-o1")==0)
			hole_path=argv[++i];
		else if(strcmp(argv[i],"-o2")==0)
			hmer_path=argv[++i];
		else if(strcmp(argv[i],"-o3")==0)
			head_path=argv[++i];
		else if(strcmp(argv[i],"-o4")==0)
			tail_path=argv[++i];
		else if(strcmp(argv[i],"-o5")==0)
			island_file=argv[++i];
		else if(strcmp(argv[i],"-m")==0)
			whole_bwtmap_file=argv[++i];
		else if(strcmp(argv[i],"-fa")==0)
			fa=true;
		else if(strcmp(argv[i],"-fq")==0)
			fa=false;
		else if(strcmp(argv[i],"-pe")==0)
			pair_end=true;
		else 
			continue;
	}

	if(reads_file==NULL)
	{fprintf(stderr,"error: reads file not specified\n");	exit(-1);}
	if(chrom_len_file==NULL)
	{fprintf(stderr,"error: chromosome length table file not specified\n");	exit(-1);}
	if(hole_path==NULL)
	{fprintf(stderr,"error: fix hole output path not specified\n");	exit(-1);}
	if(hmer_path==NULL)
	{fprintf(stderr,"error: fix hmer output path not specified\n");	exit(-1);}
	if(head_path==NULL)
	{fprintf(stderr,"error: fix head output path not specified\n");	exit(-1);}
	if(tail_path==NULL)
	{fprintf(stderr,"error: fix tail output path not specified\n");	exit(-1);}
	if(island_file==NULL)
	{fprintf(stderr,"error: island output file not specified\n");	exit(-1);}

	search_hole sh(seg_len, seg_num, min_thresh, max_thresh, extend_len);
	sh.load_chrom_table(chrom_len_file, hole_path, hmer_path, head_path, tail_path);
	if(!pair_end)
		sh.search(argv[argc-1], reads_file, fa);
	else
		sh.search_pairend(argv[argc-1], reads_file, fa);
	sh.close_all();
	if(whole_bwtmap_file!=NULL)
	{
		sh.set_whole_mapping_bit(whole_bwtmap_file);
	}
	sh.detect_island(island_gap_thresh, island_ext_len, island_file);

	//search_hole sh(25, 2, 1, 200000, 3);
	//sh.load_chrom_table("D://pair//chrom_sizes", "D://pair//hole//", "D://pair//hmer//","D://pair//head//","D://pair//tail//");
	//sh.set_whole_mapping_bit("D://new_mapsplice//unspliced_map.bwtout");
	//sh.search_pairend("D://pair//test.txt", "D://pair//divided_reads.fa",false);
	//sh.detect_island(50, 45, "D://new_mapsplice//island.gff");
}