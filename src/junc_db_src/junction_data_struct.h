/*    
 *    junction_data_struct.h		
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

#ifndef JUNCTION_H
#define JUNCTION_H

#include <string>

using namespace std;

#pragma warning(disable:4996)

struct Junction
{
public:
	Junction() 
	{}

	~Junction()
	{}

	void setValue(string line)
	{
		clear();
		char chrom_tmp[1000], name_tmp[1000];/*, strand_tmp[10], flank_seq_tmp[1000]*/;
		//sscanf(line.c_str(), "%s    \t   %d  \t %d  \t  %s   \t   %d    \t   %s    \t %*d \t %*d \t %*s \t %*d  \t             %d,%d,                   \t    %*d,%d,   \t %*f \t   %d    \t      %s     \t    %f    \t      %f",
		//	                chrom_tmp, &start, &end, name_tmp, &support, strand_tmp,                            &left_flank_length,&right_flank_length,      &gap_length,     &flank_case,  flank_seq_tmp,  &score_pq, &score_intron_size);        
		
		sscanf(line.c_str(), "%s    \t   %d  \t %d",chrom_tmp, &start, &end);        
		
		
		name = name_tmp;
		chrom=chrom_tmp;
		//strand=strand_tmp;
		//flank_seq=flank_seq_tmp;  
	}

	void clear() 
	{ 
		name.clear();
		chrom.clear();
		strand.clear();
		flank_seq.clear();
		start=0;
		end=0;
		support=0;
		left_flank_length=0;
		right_flank_length=0;
		gap_length=0;
		flank_case=0;
		score_pq=0;
		score_intron_size=0;
	}

	string chrom;   //column 1
	int start;      //column 2
	int end;        //column 3
	string name;    //column 4
	int support;    //column 5
	string strand;  //column 6
	int left_flank_length;  //column 11 left part
	int right_flank_length; //column 11 right part
	int gap_length;         //column 12 right part
	int flank_case;         //column 14
	string flank_seq;       //column 15
	float score_pq;         //column 16
	float score_intron_size;//column 17
	//vector<size_t> next_junc;
	//vector<size_t> pre_junc;
};



#endif