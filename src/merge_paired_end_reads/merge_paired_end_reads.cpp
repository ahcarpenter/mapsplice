#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
using namespace std;

size_t minlen = 10000, maxlen = 0;

int countall=0;

void convert_ex(char* extendFile, char* extendFile2, ofstream& output_fs, int fasta)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	ifstream input_fs2(extendFile2);
	if( !input_fs2 ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}   

	string readidline;
	string line, line2, line3, line4;
	string line5, line6, line7, line8;
	
	while(!input_fs.eof() && !input_fs2.eof())
	{
		++countall;

		getline(input_fs, line);
		getline(input_fs, line2);

		if (line.empty() || line2.empty())
			continue;

		output_fs <<line[0] <<countall<<'~'<<  line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			if (line3.empty() || line4.empty())
				continue;

			output_fs << line3 << endl << line4 << endl;
		}

		++countall;

		getline(input_fs2, line5);	
		getline(input_fs2, line6);

		if (line5.empty() || line6.empty())
			continue;

		output_fs <<line[0]  <<countall<<'~'<< line5.c_str() + 1 << endl << line6 << endl;

		if (!fasta)
		{
			getline(input_fs2, line7);	
			getline(input_fs2, line8);

			if (line7.empty() || line8.empty())
				continue;

			output_fs << line7 << endl << line8 << endl;
		}
	}

	input_fs.close();
	
	cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

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
		++countall;

		getline(input_fs, line);
		getline(input_fs, line2);

		if (line.empty() || line2.empty())
			continue;

		output_fs <<line[0] <<countall<<'~'<<  line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			if (line3.empty() || line4.empty())
				continue;

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

	int pairend = atoi(argv[argc - 3]);

	if( !output_fs ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	}

	if (pairend && argc % 2 == 0)
	{
		for (size_t i = 1; i < argc - 3; i += 2)
		{
			if (i + 1 < argc - 3)
				convert_ex(argv[i], argv[i + 1], output_fs, fasta);
			else
				convert_ex(argv[i], output_fs, fasta);
		}
	}
	else
	{
		for (size_t i = 1; i < argc - 3; ++i)
			convert_ex(argv[i], output_fs, fasta);
	}

	cout <<"maxlen "<<maxlen<<" minlen "<<minlen<<endl;

	output_fs.close();
	return 0;
}