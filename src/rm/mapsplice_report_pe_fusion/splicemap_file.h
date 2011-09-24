#ifndef SPLICEMAP_FILE_H
#define SPLICEMAP_FILE_H

#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>

using namespace std;

class Splicemap_File
{
public:
	Splicemap_File(const char* filename)
	{
		bwt_fs.open(filename);
		if( !bwt_fs ) 
		{
			fprintf(stderr,"error: open splice file %s error\n",filename);
			exit(1);
		}  
		endoffile=false;
		pairend=false;
		load_next_block();
	}

	Splicemap_File(const char* filename, bool pair_end)
	{
		bwt_fs.open(filename);
		if( !bwt_fs ) 
		{
			fprintf(stderr,"error: open splice file %s error\n",filename);
			exit(1);
		}  
		endoffile=false;
		pairend=pair_end;
		load_next_block();
	}

	~Splicemap_File()
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

	void load_next_block()
	{
		current_splice.clear();
		if(!load_next_line())
			return;
		else
		{
			Splicemap_Info new_splice_info;
			new_splice_info.set_value(line, pairend);
			current_splice.push_back(new_splice_info);
			line.clear();
		}
		while(load_next_line())
		{
			Splicemap_Info new_splice_info;
			new_splice_info.set_value(line, pairend);
			if(new_splice_info.read_id==current_splice[0].read_id)
			{
				current_splice.push_back(new_splice_info);
				line.clear();
			}
			else
				break;
		}
	}

	void close()
	{
		current_splice.clear();
		bwt_fs.close();
	}

	ifstream bwt_fs;
	string line;
	bool endoffile;
	bool pairend;
	vector<Splicemap_Info> current_splice;
};

#endif