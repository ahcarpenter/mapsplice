/*    
 *    sequence_file.h		
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

#ifndef SEQ_FILE_H
#define SEQ_FILE_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

class Seq_File
{
public:
	ifstream seq_fs;
	string old_line;
	string cur_chrom;
	string cur_seq;

	Seq_File(char* filename)
	{
		seq_fs.open(filename);
		if( !seq_fs ) 
		{
			fprintf(stderr,"error: open sequence file error\n");
		}  
	}

	Seq_File()
	{
	}

	~Seq_File()
	{

	}

	void init(char* filename)
	{
		close();
		seq_fs.open(filename);
		if( !seq_fs.is_open()) 
		{
			fprintf(stderr,"error: open sequence file error\n");
			exit(1);
		}  
	}

	void init(string filename)
	{
		close();
		const char *filename_c_str = filename.c_str();
		seq_fs.open(filename_c_str);
		if( !seq_fs.is_open()) 
		{
			fprintf(stderr,"error: open sequence file error\n");
			exit(1);
		}  
	}

	int get_seq_len()
	{
		return (int)(cur_seq.length());
	}

	void load_next_chrom_seq()
	{
		cout<<"loading sequence"<<endl;
		string new_line;
		getline(seq_fs, new_line);
		cur_chrom=new_line.substr(1,new_line.length()-1);
		while(getline(seq_fs, new_line))
		{
			if(new_line=="")
				continue;
			cur_seq.append(new_line);
		}
		cout<<"sequence load complete"<<endl;
	}

	string get_sequence(int start, int end)
	{
		return cur_seq.substr(start-1, end - start + 1 );
	}

	void close()
	{
		cur_seq.clear();
		cur_chrom.clear();
		old_line.clear();
		if(seq_fs.is_open())
		{
			seq_fs.close();
			seq_fs.clear(); 
		}
	}
};

#endif