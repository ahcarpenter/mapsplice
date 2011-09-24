/*    
 *    junction_seq_construction.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng
 *
 *    Authors: Zeng Zheng
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
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include "junction_data_struct.h"
#include "junction_file.h"
#include "junction_sequence.h"
#include "sequence_file.h"

const int NUM_BP_PER_LINE=70;

class Junction_Seq_Construction
{
public:
	vector<Junction> junction_set;
	vector<Junc_Seq> junction_head_set;
	vector<Junc_Seq> junction_tail_set;
	size_t junc_count;
	size_t junc_seq_count;
	int min_anchor;
	int max_anchor;
	Seq_File sf;
	map<string, size_t> chrom_start_map;
	map<string, size_t> chrom_end_map;
	ofstream output_fs;
	int curr_chrom_len;
	size_t max_seq_thresh;

	Junction_Seq_Construction(int min_anchor_len, int max_anchor_len, int max_sequence_threshold)
	{
		min_anchor=min_anchor_len;
		max_anchor=max_anchor_len;
		junc_count=0;
	    junc_seq_count=0;
		max_seq_thresh= (size_t)max_sequence_threshold;
	}

	int check_link(int pos1, int pos2, int min_dis, int max_dis)
	{
		if(pos2-pos1 <= min_dis)
			return -1;
		else if(pos2-pos1 >= max_dis -1)
			return 1;
		else
			return 0;
	}

	void generate_junc_head(Junc_Seq& my_junc_seq, size_t curr_junc, size_t min_index)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		for(int i=(int)curr_junc - 1;i>=(int)min_index;i--)
		{
			int link_status=check_link(junction_set[i].end, junction_set[curr_junc].start, min_anchor, len_left);
			if(link_status == -1)
				continue;
			else if(link_status == 1)
				break;
			else
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				int inter_seq_len=junction_set[curr_junc].start - junction_set[i].end + 1;
				string inter_sequence=sf.get_sequence(junction_set[i].end, junction_set[curr_junc].start);
				new_junc_seq.sequence.insert(0,inter_sequence);
				new_junc_seq.jump_code.insert(0, int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.insert(0, int_to_str(junction_set[i].end - junction_set[i].start - 1) + "N" );
				generate_junc_head(new_junc_seq, i, min_index);
			}
		}	
		//normal path
		my_junc_seq.start= max(junction_set[curr_junc].start - len_left + 1, 1);
		string start_sequence=sf.get_sequence(my_junc_seq.start, junction_set[curr_junc].start);
		my_junc_seq.sequence.insert(0,start_sequence);
		my_junc_seq.jump_code.insert(0,int_to_str((int)start_sequence.length())+"M");
		junction_head_set.push_back(my_junc_seq);
	}

	void generate_junc_tail(Junc_Seq& my_junc_seq, size_t curr_junc, size_t max_index)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		for(int i=(int)curr_junc+1;i<=(int)max_index;i++)
		{
			int link_status=check_link(junction_set[curr_junc].end, junction_set[i].start, min_anchor, len_left);
			if(link_status == -1)       // too close
				continue;
			else if(link_status == 1)   // too far
				break;
			else                        // fit
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				int inter_seq_len = junction_set[i].start - junction_set[curr_junc].end + 1;
				string inter_sequence=sf.get_sequence(junction_set[curr_junc].end, junction_set[i].start);
				new_junc_seq.sequence.append(inter_sequence);
				new_junc_seq.jump_code.append( int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.append(int_to_str(junction_set[i].end - junction_set[i].start - 1) + "N" );
				generate_junc_tail(new_junc_seq, i, max_index);
			}
		}	
		//normal path
		int seq_end_pos=min(junction_set[curr_junc].end + len_left - 1, curr_chrom_len);
		string end_sequence=sf.get_sequence(junction_set[curr_junc].end, seq_end_pos);
		my_junc_seq.sequence.append(end_sequence);
		my_junc_seq.jump_code.append( int_to_str((int)end_sequence.length())+"M");
		junction_tail_set.push_back(my_junc_seq);
	}

	string int_to_str(int numerical)
	{
		char c[100];
		sprintf(c,"%d",numerical);
		string str(c);
		return str;
	}

	void output_junc_seq(size_t curr_junc, string chrom)
	{
		string inter_jump_code=int_to_str(junction_set[curr_junc].end-junction_set[curr_junc].start - 1)+"N";
		for(size_t i=0;i<junction_head_set.size();i++)
		{
			for(size_t j=0;j<junction_tail_set.size();j++)
			{
				junc_seq_count++;
				output_fs<<">"<<chrom<<"_"<<junction_head_set[i].start<<":"<<junction_head_set[i].jump_code<<inter_jump_code<<junction_tail_set[j].jump_code<<endl;	
				string whole_sequence=junction_head_set[i].sequence+junction_tail_set[j].sequence;
				if((int)whole_sequence.length() > 2*max_anchor)
					cout<<"error: exceed 2 * max_anchor length"<<endl;
				size_t start = 0;
				size_t end= start + NUM_BP_PER_LINE;
				while(true)
				{
					if(end<whole_sequence.length())
					{
						output_fs<<whole_sequence.substr(start,NUM_BP_PER_LINE)<<endl;
					}
					else if(end<whole_sequence.length())
					{
						output_fs<<whole_sequence.substr(start,NUM_BP_PER_LINE)<<endl;
						break;
					}
					else
					{
						output_fs<<whole_sequence.substr(start,whole_sequence.length() - start)<<endl;
						break;
					}
					start=end;
					end=start + NUM_BP_PER_LINE;
				}
			}
		}
	}

	void construct_junc_seq(char* jun_file, char* refseq_path, char* output_file)
	{
		output_fs.open(output_file);
		if( !output_fs ) 
		{
			fprintf(stderr,"error: write junction sequence file error\n");exit(1);
		} 		
		Junction_File junc_file(jun_file);
		junc_count=junc_file.load_all_junctions(junction_set,chrom_start_map, chrom_end_map); //load all junctions
		junc_file.close();
		cout<<junc_count<<" junctions load complete"<<endl;
		map<string, size_t>::iterator it_start = chrom_start_map.begin();
		map<string, size_t>::iterator it_end = chrom_end_map.begin();
		while(it_start!=chrom_start_map.end()&&it_end!=chrom_end_map.end())
		{
			sf.init(refseq_path + it_start->first + ".fa");
			sf.load_next_chrom_seq();
			curr_chrom_len=sf.get_seq_len();
			for(size_t i=it_start->second;i<=it_end->second;i++)
			{
				junction_head_set.clear();
				junction_tail_set.clear();

				Junc_Seq	head_junc_seq;   
				generate_junc_head(head_junc_seq, i,it_start->second);

				Junc_Seq	tail_junc_seq;   
				tail_junc_seq.start=junction_set[i].end;
				generate_junc_tail(tail_junc_seq, i,it_end->second);

				if(junction_tail_set.size()*junction_head_set.size()>max_seq_thresh)
					continue;
					//cout<<"warning: junction:"<<junction_set[i].chrom<<" "<<junction_set[i].start<<" "<<junction_set[i].end<<" has "<<junction_tail_set.size()*junction_head_set.size()<<" junction sequence"<<endl;
				output_junc_seq(i,it_start->first);
			}
			cout<<it_start->first<<" junction sequence generate complete"<<endl;
			it_start++;
			it_end++;
		}
		output_fs.close();
		cout<<"total "<<junc_seq_count<<" junction sequence generated"<<endl;
		sf.close();
	}
};

void print_usage()
{
	fprintf(stderr,"junc_db [min_anchor] [max_anchor] [max_sequence_threshold] [junction_file] [sequence_path] [output_file]\n");
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc ==1 )
		print_usage();
	else if (argc < 7)
	{
		fprintf(stderr,"error: too few arguments\n");
		print_usage();
		exit(1);
	} 
	else if (argc > 7)
	{
		fprintf(stderr,"error: too many arguments\n");
		print_usage();
		exit(1);
	}

	Junction_Seq_Construction my_construct(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	my_construct.construct_junc_seq(argv[4],argv[5],argv[6]);
}
