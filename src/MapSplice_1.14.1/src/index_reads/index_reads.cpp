#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
using namespace std;

size_t minlen = 10000, maxlen = 0;

int countall=0;

void convert_ex(char* extendFile, ofstream& output_fs, int fasta)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string readidline;
	string line, line2, line3, line4;
	
	while(!input_fs.eof())
	{
		getline(input_fs, line);	
		getline(input_fs, line2);

		countall++;

		output_fs << line[0]<<countall<<'~'<<line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			output_fs << line3 << endl << line4 << endl;
		}
	}
	input_fs.close();
	
	cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}
	ofstream output_fs(argv[argc - 1]);

	int fasta = atoi(argv[argc - 2]);

	if( !output_fs ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	}

	for (size_t i = 1; i < argc - 2; ++i)
		convert_ex(argv[i], output_fs, fasta);

	cout <<"maxlen "<<maxlen<<" minlen "<<minlen<<endl;

	output_fs.close();
	return 0;
}