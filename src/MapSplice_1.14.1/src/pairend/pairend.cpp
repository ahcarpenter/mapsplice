#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;

bool
comp(const pair<int, string>& lhs, const pair<int, string>& rhs)
{
	return lhs.first < rhs.first;
}
void
sort_sam(char* fastaFile)
{
	map<string, vector<pair<int, string> > >allsams;

	ifstream input_fs(fastaFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string line;
	while(getline(input_fs, line))
	{
		char chromname[1000], tagname[1000];
		int strand, offset;
		sscanf(line.c_str(), "%s\t%d\t%s\t%d", tagname, &strand, chromname, &offset);

		allsams[chromname].push_back(make_pair(offset, line));
	}
	input_fs.close();

	ofstream ofs(fastaFile);

	map<string, vector<pair<int, string> > >::iterator m_iter;
	for (m_iter = allsams.begin(); m_iter != allsams.end(); ++m_iter)
	{
		sort(m_iter->second.begin(), m_iter->second.end(), comp);
		vector<pair<int, string> >::iterator m_iter2;
		for (m_iter2 = m_iter->second.begin(); m_iter2 != m_iter->second.end(); ++m_iter2)
		{
			ofs<<m_iter2->second<<endl;
		}
	}
	ofs.close();
}

void read_sam(char* extendFile, map<string, map<string, map<string, vector<pair<int, string> > > > >& sam_map)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}    
 
	string extendFilestr = extendFile;

	if (extendFilestr.find(".txt") != string::npos)
	{
		string skipline;
		getline(input_fs, skipline);
	}
	string line;
	int count=0;
	while(getline(input_fs, line))
	{
		if (line == "")
			continue;
		count++;

		char tagname[1000], chrom[100];
		int strand, offset;
		sscanf(line.c_str(), "%s\t%d\t%s\t%d", tagname, &strand, chrom, &offset);

		string tagnamestr = tagname;

		if (tagnamestr.find_last_of("_") != string::npos)
		{
			tagnamestr = tagnamestr.substr(0, tagnamestr.find_last_of("_"));
		}
		string tidstr = tagnamestr.substr(tagnamestr.length() - 1, 1);

		tagnamestr = tagnamestr.substr(0, tagnamestr.length() - 2);

		sam_map[tagnamestr][tidstr][chrom].push_back(make_pair(offset, line));
	}
	input_fs.close();
	//output_fs.close();
	//cout<<count<<" reads converted for "<<extendFile<<endl; 
}

void read_sam_by_read(char* extendFile, int maxmapped, ofstream& ofs_single, ofstream& ofs_multiple)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}    
 
	string extendFilestr = extendFile;

	//if (extendFilestr.find(".txt") != string::npos)
	//{
	//	string skipline;
	//	getline(input_fs, skipline);
	//}

	size_t unique = 0, multiple = 0, multiple_more_than_1 = 0, multiple_1 = 0;
	string prev_tagname = "";

	map<int, vector<string> >mapped_reads;

	string line;
	int count=0;
	while(getline(input_fs, line))
	{
		if (line == "")
			continue;
		count++;

		char tagname[1000], chrom[100], mapped[100];
		unsigned short strand, something;
		size_t offset;

		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped);

		string mappedstr = mapped;

		if (mappedstr.find("I") != string::npos)
			continue;

		string tagnamestr = tagname;

		string tidstr = tagnamestr.substr(tagnamestr.length() - 1, 1);

		int tid = atoi(tidstr.c_str());

		if (tagnamestr.find_first_of("~") != string::npos)
		{
			size_t index = tagnamestr.find_first_of("~");
			tagnamestr = tagnamestr.substr(index + 1, tagnamestr.length() - index - 3);
		}
		else
		{
			tagnamestr = tagnamestr.substr(0, tagnamestr.length() - 2);
		}

		if (prev_tagname.empty() || prev_tagname == tagnamestr)
		{
			mapped_reads[tid].push_back(line);

			if (prev_tagname.empty())
				prev_tagname = tagnamestr;
		}
		else
		{
			map<int, vector<string> >::iterator msi_iter;

			if (mapped_reads.size() == 1)
			{
				++unique;

				for (msi_iter = mapped_reads.begin(); msi_iter !=mapped_reads.end(); ++msi_iter)
				{
					vector<string>::iterator v_iter;
					for (v_iter = msi_iter->second.begin(); v_iter != msi_iter->second.end() && (v_iter -  msi_iter->second.begin())  < maxmapped; ++v_iter)
					{
						ofs_single<<*v_iter<<endl;
					}
				}
			}
			else if (mapped_reads.size() == 2)
			{
				++multiple;

				bool more_than_1 = false;

				for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
				{
					if (msi_iter->second.size() > 1)
						more_than_1 = true;

					vector<string>::iterator v_iter;

					for (v_iter = msi_iter->second.begin(); v_iter != msi_iter->second.end() && (v_iter -  msi_iter->second.begin()) < maxmapped; ++v_iter)
					{
						ofs_multiple<<*v_iter<<endl;
					}
				}

				if (more_than_1)
					++multiple_more_than_1;
				else
					++multiple_1;

			}
			else
			{
				cout << "should not be here"<<endl;

			}

			mapped_reads.clear();

			mapped_reads[tid].push_back(line);

			prev_tagname = tagnamestr;
		}
	}
	input_fs.close();

	cout<<"single reads "<<unique<<endl; 
	cout<<"paired reads "<<multiple<<endl;
	cout << "paired reads more than 1 "<<multiple_more_than_1<<endl;
	cout << "paired reads 1 "<<multiple_1<<endl;
}

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	int maxmapped = atoi(argv[2]);

	ofstream ofs(argv[3]);

	ofstream singleofs(argv[4]);

	read_sam_by_read(argv[1], maxmapped, singleofs, ofs);
	
	return 0;
}