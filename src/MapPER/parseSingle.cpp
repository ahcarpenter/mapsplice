/*    
 *    PERfusion_final.h		
 *    MapPER
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
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








#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>


using namespace std;


char targetChromosome[500][50];
int targetChrNum;

class reads
{
public:
	char name[50];
	bool fusionRead;
	int chrNum1; //chromosome number
	int chrNum2;
	char chromosome1[10];
	char chromosome2[10];
	char strand1;
	char strand2;
	long startpoint;
	long endpoint;
	long startpoint2;
	long endpoint2;
	char end[100];
	long junctionStart;
	long junctionEnd;

	reads();
	void compute_endpoint();
	void compute_chrNum();
};

long readAlignIndex[500];
long totalReadNum;
//reads *readAlign[500][1000000];

ofstream outputfile[500];

char lastTimeName[100];
char dirPrefix[400];

ofstream abnormal_reads;
ofstream dropped_reads;

reads::reads()
{
	fusionRead = false;
	for (int i = 0; i < 50; i++)
	{
		name[i] = '\0';
	}
	for (int i = 0; i < 10; i++)
	{
		chromosome1[i] = '\0';
		chromosome2[i] = '\0';
	}
	for (int i = 0; i < 100; i++)
	{
		end[i] = '\0';
	}
	chrNum1 = 0;
	chrNum2 = 0;
	startpoint = 0;
	startpoint2 = 0;
	endpoint = 0;
	endpoint2 = 0;
	junctionStart = 0;
	junctionEnd = 0;
}

void reads::compute_endpoint()
{
	int i, tmp = 0;

	endpoint = startpoint - 1;

	for (i = 0; end[i] != '\0'; i++)
	{
		if (end[i] == 'M')
		{
			endpoint = endpoint + tmp;
			tmp = 0;
		} 
		else if (end[i] == 'N')
		{
			endpoint = endpoint + tmp;
			tmp = 0;
		}
		else
		{
			tmp = tmp * 10 + end[i] - 48;
		}
	}

	return;
}

void reads::compute_chrNum()
{
	int tmp = 0;

	for (int i = 3; chromosome1[i] != '\0'; i++)
	{
		tmp = tmp * 10 + chromosome1[i] - 48;
	}

	chrNum1 = tmp;

	tmp = 0;

	for (int i = 3; chromosome2[i] != '\0'; i++)
	{
		tmp = tmp * 10 + chromosome2[i] - 48;
	}

	chrNum2 = tmp;

	return;
}


void process(reads* curRead)
{
	//output a single read
	char chromosome[50], outputfilename[500];
	int tmp, i;
	long junctionStart = -1, junctionEnd = -1;

	if (curRead->fusionRead == false)
	{
		strcpy(chromosome, curRead->chromosome1);
	} 
	else
	{
		if (curRead->chrNum1 < curRead->chrNum2)
		{
			
// 			if (curRead->strand1 == '+' && curRead->strand2 == '+')
// 			{
// 				junctionStart = curRead->endpoint;
// 				junctionEnd = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '+' && curRead->strand2 == '-')
// 			{
// 				junctionStart = curRead->endpoint;
// 				junctionEnd = curRead->endpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '+')
// 			{
// 				junctionStart = curRead->startpoint;
// 				junctionEnd = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '-')
// 			{
// 				junctionStart = curRead->startpoint;
// 				junctionEnd = curRead->endpoint2;
// 			}

			strcpy(chromosome, curRead->chromosome1);
			strcat(chromosome, curRead->chromosome2);
		} 
		else if (curRead->chrNum1 == curRead->chrNum2) // && curRead->endpoint < curRead->startpoint2)
		{
			strcpy(chromosome, curRead->chromosome1);

// 			if (curRead->strand1 == '+' && curRead->strand2 == '+')
// 			{
// 				junctionStart = curRead->endpoint;
// 				junctionEnd = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '+' && curRead->strand2 == '-')
// 			{
// 				junctionStart = curRead->endpoint;
// 				junctionEnd = curRead->endpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '+')
// 			{
// 				junctionStart = curRead->startpoint;
// 				junctionEnd = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '-')
// 			{
// 				junctionStart = curRead->startpoint;
// 				junctionEnd = curRead->endpoint2;
// 			}
		}
		else if (curRead->chrNum1 > curRead->chrNum2)
		{
// 			if (curRead->strand1 == '+' && curRead->strand2 == '+')
// 			{
// 				junctionEnd = curRead->endpoint;
// 				junctionStart = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '+' && curRead->strand2 == '-')
// 			{
// 				junctionEnd = curRead->endpoint;
// 				junctionStart = curRead->endpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '+')
// 			{
// 				junctionEnd = curRead->startpoint;
// 				junctionStart = curRead->startpoint2;
// 			}
// 			else if (curRead->strand1 == '-' && curRead->strand2 == '-')
// 			{
// 				junctionEnd = curRead->startpoint;
// 				junctionStart = curRead->endpoint2;
// 			}

			strcpy(chromosome, curRead->chromosome2);
			strcat(chromosome, curRead->chromosome1);
		}
		else
		{
			abnormal_reads << curRead->name << '\t' << curRead->chromosome1 << '\t' << curRead->fusionRead << '\t' << curRead->startpoint << '\t' << curRead->end << endl;  
			return;
		}
	}

	if (curRead->strand1 == '+' && curRead->strand2 == '+')
	{
		junctionStart = curRead->endpoint;
		junctionEnd = curRead->startpoint2;
	}
	else if (curRead->strand1 == '+' && curRead->strand2 == '-')
	{
		junctionStart = curRead->endpoint;
		junctionEnd = curRead->endpoint2;
	}
	else if (curRead->strand1 == '-' && curRead->strand2 == '+')
	{
		junctionStart = curRead->startpoint;
		junctionEnd = curRead->startpoint2;
	}
	else if (curRead->strand1 == '-' && curRead->strand2 == '-')
	{
		junctionStart = curRead->startpoint;
		junctionEnd = curRead->endpoint2;
	}

	curRead->junctionStart = junctionStart;
	curRead->junctionEnd = junctionEnd;

	for (tmp = 1; tmp <= targetChrNum; tmp++)
	{
		if (strcmp(targetChromosome[tmp], chromosome) == 0)
		{
			break;
		}
	}

	if (tmp > targetChrNum)
	{
		//chromosome not found
		sprintf(outputfilename, "%sdata/single/%s.txt", dirPrefix, chromosome);
		outputfile[tmp].open(outputfilename);

		targetChrNum++;
		strcpy(targetChromosome[targetChrNum], chromosome);		
	}

	readAlignIndex[tmp]++;
//	readAlign[tmp][readAlignIndex[tmp]] = curRead;

	outputfile[tmp] << curRead->name << "\t"; 

	if (curRead->fusionRead == false)
	{
		outputfile[tmp] << "0\t" << curRead->chromosome1 << "\t" << curRead->strand1 << "\t" << curRead->startpoint << "\t" << curRead->end << endl;
	} 
	else
	{
		outputfile[tmp] << "1\t" << curRead->chromosome1 << "\t" << curRead->startpoint << "\t" << curRead->endpoint << '\t' << curRead->chromosome2 << '\t' << curRead->startpoint2 << "\t" << curRead->endpoint2 << "\t" << curRead->strand1 << curRead->strand2 << "\t" << curRead->junctionStart << '\t' << curRead->junctionEnd << endl;
	}


// 	for (i = readAlignIndex[tmp] - 1; i > 0; i--)
// 	{
// 		if (strcmp(readAlign[tmp][i]->name, curRead->name) == 0)
// 		{
// 			return;
// 		}
// 	}

	if (strcmp(lastTimeName, curRead->name) == 0)
	{
		return;
	}

	strcpy(lastTimeName, curRead->name);
	totalReadNum++;

	delete curRead;
	return;
}


void outputRead()
{
	int i, j;
	char outputfilename[500];
	reads *outputread;
	ofstream outputChrname;
	long totalAlignNum = 0;

	sprintf(outputfilename, "%sdata/single/ChromosomeName.txt", dirPrefix);
	outputChrname.open(outputfilename);

	for (i = 1; i <= targetChrNum; i++)
	{
		outputChrname << targetChromosome[i] << '\t' << readAlignIndex[i] << endl; 
		totalAlignNum += readAlignIndex[i];

		outputfile[i].close();
	}

	outputChrname << totalReadNum << '\t' << totalAlignNum << endl;

	outputChrname.close();
	dropped_reads.close();
	abnormal_reads.close();

	return;
}


void parse(char* inputfilename)
{
	reads *ptr;
	ptr = NULL;

	ifstream inputfile;
	inputfile.open(inputfilename);

	char name[100], curName[100];
	int strand, tmp, i;
	string info;

	char field1[50], field2[50], field3[50], field4[50], field5[50], field6[50], field7[50], field8[50];



	for (tmp = 0; tmp < 100; tmp++)
	{
		name[tmp] = '\0';
	}

	inputfile >> name;

	while (name[1] != '\0')
	{
		for (tmp = 0; tmp < 100; tmp++)
		{
			if (name[tmp] == '/')
			{
				break;
			}
		}
		name[tmp] = '\0'; 

		inputfile >> field1;
		inputfile >> field2;
		inputfile >> field3;
		inputfile >> field4;
		inputfile >> field5;
		inputfile >> field6;
		inputfile >> field7;
		inputfile >> field8;
		getline(inputfile, info);

		ptr = new reads;

		strcpy(ptr->name, name);

		if (strcmp(field6, "*") == 0)
		{
			//not fusion read
			ptr->fusionRead = false;
			if (strcmp(field1, "8") == 0)
			{
				ptr->strand1 = '+';
			} 
			else
			{
				ptr->strand1 = '-';
			}

			strcpy(ptr->chromosome1, field2);
			ptr->startpoint = strtol(field3, NULL, 10);
			strcpy(ptr->end, field5);

			ptr->compute_endpoint();
		} 
		else
		{
			//fusion read
			ptr->fusionRead = true;
			strcpy(ptr->chromosome1, field1);
			ptr->strand1 = field2[0];
			ptr->startpoint = strtol(field3, NULL, 10);
			ptr->endpoint = strtol(field4, NULL, 10);

			strcpy(ptr->chromosome2, field5);
			ptr->strand2 = field6[0];
			ptr->startpoint2 = strtol(field7, NULL, 10);
			ptr->endpoint2 = strtol(field8, NULL, 10);
		}

		ptr->compute_chrNum();

		process(ptr);

		name[1] = '\0';

		inputfile >> name;		
	}


	inputfile.close();

	return;
}


void initialization()
{
	int i;

// 	for (i = 0; i < 500; i++)
// 	{
// 		readAlignIndex[i] = 0;
// 	}

	targetChrNum = 0;
	totalReadNum = 0;
	sprintf(lastTimeName, "XXXX");

	return;
}


int main(int argc, char* argv[])
{
//	sprintf(dirPrefix, ".");

	if (argc != 3)
	{
		cout << argv[0] << "\t<filename>\t<target_path>" << endl;
		return 1;
	}

	char outpurfilename[500];
	sprintf(outpurfilename, "%sdata/single/DroppedReads.txt", dirPrefix);
	dropped_reads.open(outpurfilename);
	sprintf(outpurfilename, "%sdata/single/AbnormalReads.txt", dirPrefix);
	abnormal_reads.open(outpurfilename);

	char inputfilename[500];
//	sprintf(inputfilename, "%s/data/merged_single.sam", dirPrefix);
//	sprintf(inputfilename, "/sdb1/mapsplice/hmy_111/MapSplice_1.11/50bp_0520/merged_single.sam");
	strcpy(inputfilename, argv[1]);

	strcpy(dirPrefix, argv[2]);

	initialization();

	parse(inputfilename);

	outputRead();

	return 0;
}