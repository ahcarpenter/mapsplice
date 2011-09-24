#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
using namespace std;

size_t minlen = 10000, maxlen = 0;

int countall=0;

void convert_ex_fasta(char* extendFile, ofstream& output_fs, int num, int len)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}    

	string readidline;
	string line;
	//int countall=0;
	while(!input_fs.eof())
	{
		getline(input_fs, readidline);
		getline(input_fs, line);
		if (line == "" || readidline == "")
			continue;
		countall++;

		//if (readidline[0] != '>')
		//	readidline[0] = '>';
//		>HWUSI-EAS524:2:1:1:82#0/1
//ANAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACT

		//string firsthalf = line.substr(0, line.length()/2);

		if (line.length() > maxlen)
			maxlen = line.length();

		if (line.length() < minlen)
			minlen = line.length();
			

		if ((len * num) > line.length())
		{
			cerr << "some read length < "<<len * num << endl <<readidline << endl<< line<<endl;
			exit(2);
		}
		//int len = line.length() / num;

		for (int i = 1; i < num; ++i)
		{
			output_fs<< ">" <<i<<"~"<<countall<<"~"<<readidline.c_str() + 1<<endl;
			output_fs<<line.substr((i-1)*len, len)<<endl;
		}

		//if (line.length() - (len*(num-1)) >= len)
		//{
		output_fs<<">" <<num<<"~"<<countall<<"~"<<readidline.c_str() + 1<<endl;
		output_fs<<line.substr((num - 1)*len, len)<<endl;
		//}
	}
	input_fs.close();
	
	cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

void convert_ex_fastq(char* extendFile, ofstream& output_fs, int num, int len)
{
	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}    

	string readidline;
	string line;
	string readidline2;
	string qualstr;
	//int countall=0;
	while(!input_fs.eof())
	{
		getline(input_fs, readidline);
		getline(input_fs, line);
		getline(input_fs, readidline2);
		getline(input_fs, qualstr);

		if (line == "" || readidline == "" || readidline2 == "" || qualstr == "")
			continue;
		countall++;

		if (line.length() > maxlen)
			maxlen = line.length();

		if (line.length() < minlen)
			minlen = line.length();

		//if (readidline[0] != '>')
		//	readidline[0] = '>';
//		>HWUSI-EAS524:2:1:1:82#0/1
//ANAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACT

		//string firsthalf = line.substr(0, line.length()/2);

		//int len = line.length() / num;
		if ((len * num) > line.length() || (len * num) > qualstr.length())
		{
			cerr << "some read or quality length < "<<len * num << endl <<readidline << endl<< line<<endl << readidline2 <<endl<<qualstr << endl;
			exit(2);
		}

		for (int i = 1; i < num; ++i)
		{
			output_fs<<"@" <<i<<"~"<<countall<<"~"<<readidline.c_str() + 1<<endl;
			output_fs<<line.substr((i-1)*len, len)<<endl;
			output_fs<<readidline2<<endl;
			output_fs<<qualstr.substr((i-1)*len, len)<<endl;
		}

		//if (line.length() - (len*(num-1)) >= len)
		//{
		output_fs<<"@" <<num<<"~"<<countall<<"~"<<readidline.c_str() + 1<<endl;
		output_fs<<line.substr((num - 1)*len, len)<<endl;
		output_fs<<readidline2<<endl;
		output_fs<<qualstr.substr((num - 1)*len, len)<<endl;
		//}
	}
	input_fs.close();
	
	cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

int main(int argc, char* argv[])
{
	if (argc < 6)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	ofstream output_fs(argv[argc-4]);

	int fasta = atoi(argv[argc - 3]);

	if( !output_fs ) 
	{
		fprintf(stderr,"error: write fasta file error\n");exit(1);
	}

	if (fasta)
	{
		for(int i=1;i<argc-4;i++)
		{
			convert_ex_fasta(argv[i], output_fs, atoi(argv[argc-2]), atoi(argv[argc-1]));
		}
	}
	else
	{
		for(int i=1;i<argc-4;i++)
		{
			convert_ex_fastq(argv[i], output_fs, atoi(argv[argc-2]), atoi(argv[argc-1]));
		}
	}

	cout << "minlen "<< minlen<<", maxlen "<<maxlen<<endl;

	output_fs.close();
	return 0;
}