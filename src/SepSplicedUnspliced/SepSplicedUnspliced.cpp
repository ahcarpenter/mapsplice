#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>

using namespace std;

void
RemoveDupMapreads(const char* infile, /*const char* outfile, *//*const char* unspliced_reads,*/ const char* spliced, const char* unspliced, const char* stat_file)
{
	ifstream ifs(infile);

	//ofstream ofs(outfile);

	//ofstream unspliced_ofs(unspliced_reads);

	ofstream spliced_ofs(spliced);

	ofstream unspliced_ofs(unspliced);

	ofstream stat_ofs(stat_file);

	size_t count = 0, unspliced_count = 0, spliced_count = 0;

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

			if (mappedstr.find("N") == string::npos || mappedstr.find("I") != string::npos)
			{
				++unspliced_count;
				unspliced_ofs << line << endl;
			}
			else
			{
				++spliced_count;
				spliced_ofs << line << endl;			
			}
		}
		ifs.close();
	}

	stat_ofs << "Separate spliced unspliced"<<endl << "total_mapped: " << count<<endl<<"unspliced: "<< unspliced_count <<endl<< 
		"spliced: " << spliced_count<<endl;
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
	const char* spliced = argv[2];
	const char* unspliced = argv[3];
	const char* stat_file = argv[4];
	RemoveDupMapreads(infile, spliced, unspliced, stat_file);
}