/*    
 *    merge_sam.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
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
void merge_sam(char* extendFile, ofstream& output_fs)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	if( !output_fs ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	}    
	string line;
	int count=0;
	while(getline(input_fs, line))
	{
		if (line == "")
			continue;
		count++;

		output_fs<<line<<endl;
	}
	input_fs.close();
	
	cout<<count<<" reads converted for "<<extendFile<<endl; 
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	ofstream ofs(argv[argc-1]);

	for(int i=1;i<argc-1;i++)
	{
		merge_sam(argv[i], ofs);
	}

	ofs.close();

	//sort_sam(argv[argc-1]);
	return 0;
}