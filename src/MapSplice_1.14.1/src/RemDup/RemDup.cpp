/*    
 *    RemDup.cpp		
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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>

using namespace std;

void
RemoveDupMapreads(const char* infile, /*const char* outfile, *//*const char* unspliced_reads,*/ const char* unique_spliced_reads, const char* multiple_spliced_reads, const char* stat_file)
{
	ifstream ifs(infile);

	//ofstream ofs(outfile);

	//ofstream unspliced_ofs(unspliced_reads);

	ofstream unique_spliced_ofs(unique_spliced_reads);

	ofstream multiple_spliced_ofs(multiple_spliced_reads);

	string dup_reads = infile; dup_reads.append(".dup"); ofstream ofs_dup(dup_reads.c_str());

	string prev_tagname = "";

	map<string, int> mapped_reads;

	size_t count = 0, unspliced = 0, unique_spliced = 0, multiple_spliced = 0;

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "")
				continue;

			char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];
			unsigned short strand, something, mismatch;
			size_t offset;

			char mate_match;

			size_t mate_offest, mate_diff;

			string alterstr = "";

			string qualstr = "I";

			//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
			//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

			size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

			string mappedstr = mapped;

			++count;

			//if (mappedstr.find("N") == string::npos)
			//{
			//	++unspliced;
			//	unspliced_ofs << line << endl;
			//	continue;
			//}

			if (prev_tagname.empty() || prev_tagname == tagname)
			{
				map<string, int>::iterator mapped_reads_iter = mapped_reads.find(line);
				if (mapped_reads_iter == mapped_reads.end())
					mapped_reads[line] = 1;
				else
					mapped_reads_iter->second++;

				if (prev_tagname.empty())
					prev_tagname = tagname;

				//tag_count++;
			}
			else
			{
				map<string, int>::iterator msi_iter;

				ofstream* ofs_ptr;
				if (mapped_reads.size() == 1)
				{
					++unique_spliced;

					ofs_ptr = &unique_spliced_ofs;
				}
				else
				{
					++multiple_spliced;

					ofs_ptr = &multiple_spliced_ofs;
				}

				for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
				{
					(*ofs_ptr)<<msi_iter->first<<endl;

					if (msi_iter->second > 1)
						ofs_dup<<msi_iter->first<<endl;
				}

				mapped_reads.clear();

				mapped_reads[line] = 1;

				prev_tagname = tagname;
			}			
		}
		ifs.close();

	}
	else
	{
		cout << "can't open file "<< infile<<endl; exit(1);
	}

	map<string, int>::iterator msi_iter;

	ofstream* ofs_ptr;
	if (mapped_reads.size() == 1)
	{
		++unique_spliced;

		ofs_ptr = &unique_spliced_ofs;
	}
	else
	{
		++multiple_spliced;

		ofs_ptr = &multiple_spliced_ofs;
	}

	for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
	{
		(*ofs_ptr)<<msi_iter->first<<endl;

		if (msi_iter->second > 1)
			ofs_dup<<msi_iter->first<<endl;
	}

	ofstream stat_ofs(stat_file, ios_base::app);

	stat_ofs << "Remove duplication"<<endl << "total_mapped: " << count<<endl<<"unspliced: "<< unspliced <<endl<< 
		"unique_spliced_reads: " << unique_spliced<<endl<<"multiple_spliced_reads: "<<  multiple_spliced<<endl;

	//ofs.close();

}

int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cout << "infile unspliced_reads unique_spliced_reads multiple_spliced_reads" <<endl;
		exit(0);
	}
	const char* infile = argv[1];
	//const char* unspliced_reads = argv[2];
	const char* unique_spliced_reads = argv[2];
	const char* multiple_spliced_reads = argv[3];
	const char* stat_file = argv[4];
	RemoveDupMapreads(infile, /*unspliced_reads, */unique_spliced_reads, multiple_spliced_reads, stat_file);
}