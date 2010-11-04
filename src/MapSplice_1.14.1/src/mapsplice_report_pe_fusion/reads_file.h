#ifndef READS_FILE_H
#define READS_FILE_H

#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>

inline char complement(int i) 
{
	static const int b2c_size = 20;
	static const char b2c[] = {'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'};
	static const char b2cl[] = {'t','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

inline string revcomp(const string& s) 
{
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

struct Read_Block
{
public:
	Read_Block()
	{ 
	}

	~Read_Block()
	{
		clear();
	}

	string get_seg_seq(int seg_no)
	{
		return read_seq[seg_no-1];
	}

	string get_revcom_seg_seq(int seg_no)
	{
		return read_seq_revcom[seg_no-1];
	}

	string get_seg_qual(int seg_no)
	{
		return read_quality[seg_no-1];
	}

	string get_revcom_seg_qual(int seg_no)
	{
		return read_quality_revcom[seg_no-1];
	}

   /////////size_t version
	string get_seg_seq(size_t seg_no)
	{
		return read_seq[seg_no-1];
	}

	string get_revcom_seg_seq(size_t seg_no)
	{
		return read_seq_revcom[seg_no-1];
	}

	string get_seg_qual(size_t seg_no)
	{
		return read_quality[seg_no-1];
	}

	string get_revcom_seg_qual(size_t seg_no)
	{
		return read_quality_revcom[seg_no-1];
	}

	void get_revcom_info(bool fasta)
	{
		for(size_t i=0;i<read_seq.size();i++)	//fa
		{
			string new_revcom_sequence=revcomp(read_seq[i]);
			read_seq_revcom.push_back(new_revcom_sequence);
		}
		if(!fasta)   //fq
		{
			for(size_t i=0;i<read_quality.size();i++)	
			{
				string new_quality_revcom=read_quality[i];
				reverse(new_quality_revcom.begin(),new_quality_revcom.end());
				read_quality_revcom.push_back(new_quality_revcom);
			}
		}
	}

	void clear()
	{
		read_id.clear();
		read_strand.clear();
		read_seq.clear();
		read_seq_revcom.clear();
		read_quality.clear();
		read_quality_revcom.clear();
	}

	string read_id;
	string read_strand;
	vector<string> read_seq;
	vector<string> read_seq_revcom;
	vector<string> read_quality;
	vector<string> read_quality_revcom;
};


class Reads_File
{
public:
	ifstream rd_fs;
	string line;
	bool endoffile;
	bool fasta;
	bool pairend;

	Reads_File(char* filename, bool _fasta)
	{
		rd_fs.open(filename);
		if( !rd_fs ) 
		{
			fprintf(stderr,"error: open reads file error\n");
			exit(1);
		}  
		endoffile=false;
		fasta=_fasta;
		pairend=false;
	}

	Reads_File(char* filename, bool _fasta, bool pair_end)
	{
		rd_fs.open(filename);
		if( !rd_fs ) 
		{
			fprintf(stderr,"error: open reads file error\n");
			exit(1);
		}  
		endoffile=false;
		fasta=_fasta;
		pairend=pair_end;
	}

	~Reads_File()
	{

	}

	bool load_next_line()
	{
		if(line.empty())
		{
			if( getline(rd_fs, line))
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
	size_t load_next_read(Read_Block& r_block)
	{
		if(!pairend)
			return load_next_read_single(r_block);
		else
			return load_next_read_pairend(r_block);
	}

	size_t load_next_read_single(Read_Block& r_block)
	{
		r_block.clear();
		while(load_next_line())
		{
			char id_tmp[MAX_READID_LENGTH];
			sscanf(line.c_str(), "%*c%*[^~]%*[~]%*[^~]%*[~]%s",id_tmp);	
			if(r_block.read_id.empty())
			{
				r_block.read_id=id_tmp;
			}
			else if(r_block.read_id!=id_tmp)
				break;
			line.clear();
			load_next_line();
			r_block.read_seq.push_back(line);
			line.clear();
			if(!fasta)    //fq
			{
				load_next_line();
				r_block.read_strand=line;
				line.clear();
				load_next_line();
				r_block.read_quality.push_back(line);
				line.clear();
			}
		}
		return r_block.read_seq.size();
	}

	size_t load_next_read_pairend(Read_Block& r_block)
	{
		r_block.clear();
		while(load_next_line())
		{
			char id_tmp[MAX_READID_LENGTH];
			sscanf(line.c_str(), "%*c%*[^~]%*[~]%*[^~]%*[~]%[^/]",id_tmp);	
			if(r_block.read_id.empty())
			{
				r_block.read_id=id_tmp;
			}
			else if(r_block.read_id!=id_tmp)
				break;
			line.clear();
			load_next_line();
			r_block.read_seq.push_back(line);
			line.clear();
			if(!fasta)    //fq
			{
				load_next_line();
				r_block.read_strand=line;
				line.clear();
				load_next_line();
				r_block.read_quality.push_back(line);
				line.clear();
			}
		}
		return r_block.read_seq.size();
	}

	void close()
	{
		rd_fs.close();
	}
};

#endif