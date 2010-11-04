/*    
 *    junction_file.h		
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

#ifndef JUNCTION_FILE_H
#define JUNCTION_FILE_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include "junction_data_struct.h"

using namespace std;

class Junction_File
{
private:
	ifstream junc_fs;
	string line;

public:
	Junction_File(char* filename)
	{
		junc_fs.open(filename);
		if( !junc_fs ) 
		{
			fprintf(stderr,"error: open junction file error\n");
			exit(1);
		}  
		getline(junc_fs, line);// skip header
	}

	~Junction_File()
	{

	}

	size_t load_next_junc(Junction& new_junc) // load a single junction
	{
		if(getline(junc_fs, line))
		{
			if(line=="")
				return 0;
			new_junc.setValue(line);
			return 1;
		}
		else
			return 0;
	}

	size_t load_next_block(vector<Junction>& new_junction_set, int number_of_junctions) // load a block of junctions
	{
		new_junction_set.clear();
		while(true)
		{
			Junction new_junc;   
			size_t count=load_next_junc(new_junc);
			if(count)
				new_junction_set.push_back(new_junc);
			else
				break;
			if(new_junction_set.size()==number_of_junctions)
				break;
		}
		return new_junction_set.size();
	}

	size_t load_all_junctions(vector<Junction>& new_junction_set)    // load all the rest of junctions
	{
		//new_junction_set.clear();
		while(true)
		{
			Junction new_junc;   
			size_t count=load_next_junc(new_junc);
			if(count)
				new_junction_set.push_back(new_junc);
			else
				break;
		}
		return new_junction_set.size();
	}

	size_t load_all_junctions(vector<Junction>& new_junction_set, map<string, size_t>& chrom_start_map, map<string, size_t>& chrom_end_map)    // load all the rest of junctions
	{
		new_junction_set.clear();
		string cur_chrom="";
		size_t total_index=0;
		while(true)
		{
			Junction new_junc;   
			size_t count=load_next_junc(new_junc);
			if(count)
			{
				new_junction_set.push_back(new_junc);
				if(new_junc.chrom!=cur_chrom)
				{
					if(cur_chrom!="")
					{
						chrom_end_map.insert(make_pair(cur_chrom,total_index-1));
					}
					chrom_start_map.insert(make_pair(new_junc.chrom,total_index));
					cur_chrom=new_junc.chrom;
				}
				total_index++;
			}
			else
			{
				chrom_end_map.insert(make_pair(cur_chrom,total_index-1));
				break;
			}
		}
		return new_junction_set.size();
	}

	void close()
	{
		junc_fs.close();
	}
};

#endif