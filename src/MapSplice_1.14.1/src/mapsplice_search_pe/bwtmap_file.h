#ifndef BWTMAP_FILE_H
#define BWTMAP_FILE_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>

using namespace std;

class Bwtmap_File
{
public:
	Bwtmap_File(char* filename)
	{
		bwt_fs.open(filename);
		if( !bwt_fs ) 
		{
			fprintf(stderr,"error: open bwtmap file error\n");
			exit(1);
		}  
		endoffile=false;
		pairend=false;
	}

	Bwtmap_File(const char* filename, bool pair_end)
	{
		bwt_fs.open(filename);
		if( !bwt_fs ) 
		{
			fprintf(stderr,"error: open bwtmap file %s error\n",filename);
			exit(1);
		}  
		endoffile=false;
		pairend=pair_end;
	}

	~Bwtmap_File()
	{

	}

	bool load_next_line()
	{
		if(line.empty())
		{
			if( getline(bwt_fs, line))
				return true;
			else
			{
				endoffile=true;			
				return false;
			}
		}
		else
			return true;
	}

	int load_next_block(vector<Bwtmap_Info>& read_mapping_set) // load a block of junctions
	{
		read_mapping_set.clear();
		if(!load_next_line())
		{		
			return -1;
		}
		Bwtmap_Info first_bwt_info(line, pairend);
		read_mapping_set.push_back(first_bwt_info);
		line.clear();
		while(load_next_line())
		{
			Bwtmap_Info new_bwt_info(line, pairend);
			if(new_bwt_info.read_id==first_bwt_info.read_id)
			{
				read_mapping_set.push_back(new_bwt_info);
				line.clear();
			}
			else
			{
				break;
			}
		}
		return (int)read_mapping_set.size();
	}

	void close()
	{
		bwt_fs.close();
	}

	ifstream bwt_fs;
	string line;
	bool endoffile;
	bool pairend;
};

#endif