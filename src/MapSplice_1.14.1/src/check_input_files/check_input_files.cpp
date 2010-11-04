#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
using namespace std;

string
basename2(string filename) {

	//cout << "bef: "<<filename<<endl;
	const string s(filename.substr(0, filename.find_last_of(".")));
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
	{
		//cout << "aft 1: "<<s.substr(final_slash + 1)<<endl;
		return s.substr(final_slash + 1);
	}
	else
	{
		//cout << "aft 2: "<<s<<endl;
		return s;
	}
}

void convert_ex_fasta(char* extendFile)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string firstline;

	getline(input_fs, firstline);

	if (firstline.find(" ") != string::npos || firstline.find("\t") != string::npos)
	{
		cerr <<"chromosome name contain blank space or tab"<<endl;
		cerr << firstline << endl;
		exit(2);
	}

	string base_chr_file_name = basename2(extendFile);

	string base_chr_name = firstline.substr(1);

	if (base_chr_name != base_chr_file_name)
	{
		cerr << "base name of chromosome file not consistent with chromosome name" << endl;
		cerr << "base chromosoem file name:" << base_chr_file_name << endl;
		cerr << "chromosome name:" <<base_chr_name << endl;
		exit(3);
	}

	input_fs.close();
}

int main(int argc, char* argv[])
{
	if (argc < 2 )
	{
		fprintf(stderr,"./check_input_files chromosome_file1 chromosome_file2 ...\n");exit(1);
	}

	for (int i = 1; i < argc; ++i)
		convert_ex_fasta(argv[i]);

	return 0;
}