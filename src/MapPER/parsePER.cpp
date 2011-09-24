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










/************************************************************************/
/* Feb 22nd. ALL PER from start read to end read                        */
/************************************************************************/




#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib> 

using namespace std;

const long MAX = 2000000000;
const int MAX_CHOICE_NUM = 100000;
const int MAX_COMB_LIMIT = 1000;

long PERCnt_same;
long PERCnt_diff;

char targetChromosome[1000][50];
int targetChrNum;
ofstream outputfile[1000];
ofstream outputfile_chrname;
ofstream dropped_reads;
ofstream abnormal_reads;
ofstream overlap_reads;

char dirPrefix[500];

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

	reads();
	void compute_endpoint();
	void compute_chrNum();
};

reads* startlist[MAX_CHOICE_NUM];
reads* endlist[MAX_CHOICE_NUM];
int startlist_count;
int endlist_count;

ofstream allDistanceToFile;


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

long compute_distance(reads* start_read, reads* end_read)
{
	//compute the min distance between the start_read and the end_read

//	return MAX;
	return end_read->startpoint - start_read->endpoint;
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




void outputPEreads(reads* start_read, reads* end_read)
{
	//output a PER
	char outpurfilename[500];
	char chromosome[50];
	int targetChrIndex, tmp;
	long junctionStart = -1, junctionEnd = -1;
	

	if (start_read->fusionRead == false && end_read->fusionRead == false)
	{
		if (start_read->chrNum1 < end_read->chrNum1)
		{
			strcpy(chromosome, start_read->chromosome1);
			strcat(chromosome, end_read->chromosome1);
		} 
		else if (start_read->chrNum1 > end_read->chrNum1)
		{
			strcpy(chromosome, end_read->chromosome1);
			strcat(chromosome, start_read->chromosome1);
		} 
		else
		{
			if (start_read->endpoint < end_read->startpoint || start_read->startpoint > end_read->endpoint)
			{
				strcpy(chromosome, start_read->chromosome1);
				strcat(chromosome, end_read->chromosome1);
			}
			else
			{
				overlap_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl; 
				return;
			}
		}
	} 
	else
	{
		if (start_read->fusionRead == true && end_read->fusionRead == false)
		{
			if (start_read->chrNum1 < start_read->chrNum2)
			{
				if (end_read->chrNum1 == start_read->chrNum1 && end_read->endpoint < start_read->startpoint || end_read->chrNum1 == start_read->chrNum2 && end_read->endpoint < start_read->startpoint2)
				{
				}
				else if (end_read->chrNum1 == start_read->chrNum2 && end_read->startpoint > start_read->endpoint2 || end_read->chrNum1 == start_read->chrNum1 && end_read->startpoint > start_read->endpoint)
				{
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (start_read->strand1 == '+' && start_read->strand2 == '+')
// 				{
// 					junctionStart = start_read->endpoint;
// 					junctionEnd = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '+' && start_read->strand2 == '-')
// 				{
// 					junctionStart = start_read->endpoint;
// 					junctionEnd = start_read->endpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '+')
// 				{
// 					junctionStart = start_read->startpoint;
// 					junctionEnd = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '-')
// 				{
// 					junctionStart = start_read->startpoint;
// 					junctionEnd = start_read->endpoint2;
// 				}

				strcpy(chromosome, start_read->chromosome1);
				strcat(chromosome, start_read->chromosome2);
			} 
			else if (start_read->chrNum1 == start_read->chrNum2)// && start_read->endpoint < start_read->startpoint2)
			{
				if (end_read->chrNum1 <= start_read->chrNum1) // || (end_read->chrNum1 == start_read->chrNum1 && end_read->endpoint < start_read->startpoint))
				{
					strcpy(chromosome, end_read->chromosome1);
					strcat(chromosome, start_read->chromosome1);
				}
				else if (end_read->chrNum1 > start_read->chrNum2) // || (end_read->chrNum1 == start_read->chrNum2 && end_read->startpoint > start_read->endpoint2))
				{
					strcpy(chromosome, start_read->chromosome1);
					strcat(chromosome, end_read->chromosome1);
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (start_read->strand1 == '+' && start_read->strand2 == '+')
// 				{
// 					junctionStart = start_read->endpoint;
// 					junctionEnd = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '+' && start_read->strand2 == '-')
// 				{
// 					junctionStart = start_read->endpoint;
// 					junctionEnd = start_read->endpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '+')
// 				{
// 					junctionStart = start_read->startpoint;
// 					junctionEnd = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '-')
// 				{
// 					junctionStart = start_read->startpoint;
// 					junctionEnd = start_read->endpoint2;
// 				}
			}
			else if (start_read->chrNum1 > start_read->chrNum2)
			{
				if (end_read->chrNum1 == start_read->chrNum2 && end_read->endpoint < start_read->startpoint2 || end_read->chrNum1 == start_read->chrNum1 && end_read->endpoint < start_read->startpoint)
				{
				}
				else if (end_read->chrNum1 == start_read->chrNum1 && end_read->startpoint > start_read->endpoint || end_read->chrNum1 == start_read->chrNum2 && end_read->startpoint > start_read->endpoint2)
				{
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (start_read->strand1 == '+' && start_read->strand2 == '+')
// 				{
// 					junctionEnd = start_read->endpoint;
// 					junctionStart = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '+' && start_read->strand2 == '-')
// 				{
// 					junctionEnd = start_read->endpoint;
// 					junctionStart = start_read->endpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '+')
// 				{
// 					junctionEnd = start_read->startpoint;
// 					junctionStart = start_read->startpoint2;
// 				}
// 				else if (start_read->strand1 == '-' && start_read->strand2 == '-')
// 				{
// 					junctionEnd = start_read->startpoint;
// 					junctionStart = start_read->endpoint2;
// 				}

				strcpy(chromosome, start_read->chromosome2);
				strcat(chromosome, start_read->chromosome1);
			}
			else
			{
				abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
				return;
			}

			if (start_read->strand1 == '+' && start_read->strand2 == '+')
			{
				junctionStart = start_read->endpoint;
				junctionEnd = start_read->startpoint2;
			}
			else if (start_read->strand1 == '+' && start_read->strand2 == '-')
			{
				junctionStart = start_read->endpoint;
				junctionEnd = start_read->endpoint2;
			}
			else if (start_read->strand1 == '-' && start_read->strand2 == '+')
			{
				junctionStart = start_read->startpoint;
				junctionEnd = start_read->startpoint2;
			}
			else if (start_read->strand1 == '-' && start_read->strand2 == '-')
			{
				junctionStart = start_read->startpoint;
				junctionEnd = start_read->endpoint2;
			}
		}
		else if (start_read->fusionRead == false && end_read->fusionRead == true)
		{
			if (end_read->chrNum1 < end_read->chrNum2)
			{
				if (start_read->chrNum1 == end_read->chrNum1 && start_read->endpoint < end_read->startpoint || start_read->chrNum1 == end_read->chrNum2 && start_read->endpoint < end_read->startpoint2)
				{
				}
				else if (start_read->chrNum1 == end_read->chrNum2 && start_read->startpoint > end_read->endpoint2 || start_read->chrNum1 == end_read->chrNum1 && start_read->startpoint > end_read->endpoint)
				{
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (end_read->strand1 == '+' && end_read->strand2 == '+')
// 				{
// 					junctionStart = end_read->endpoint;
// 					junctionEnd = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '+' && end_read->strand2 == '-')
// 				{
// 					junctionStart = end_read->endpoint;
// 					junctionEnd = end_read->endpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '+')
// 				{
// 					junctionStart = end_read->startpoint;
// 					junctionEnd = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '-')
// 				{
// 					junctionStart = end_read->startpoint;
// 					junctionEnd = end_read->endpoint2;
// 				}

				strcpy(chromosome, end_read->chromosome1);
				strcat(chromosome, end_read->chromosome2);
			} 
			else if (end_read->chrNum1 == end_read->chrNum2)// && end_read->endpoint < end_read->startpoint2)
			{
				if (start_read->chrNum1 <= end_read->chrNum1) // || (start_read->chrNum1 == end_read->chrNum1 && start_read->endpoint < end_read->startpoint))
				{
					strcpy(chromosome, start_read->chromosome1);
					strcat(chromosome, end_read->chromosome1);
				}
				else if (start_read->chrNum1 > end_read->chrNum2) // || (start_read->chrNum1 == end_read->chrNum2 && start_read->startpoint > end_read->endpoint2))
				{
					strcpy(chromosome, end_read->chromosome1);
					strcat(chromosome, start_read->chromosome1);
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (end_read->strand1 == '+' && end_read->strand2 == '+')
// 				{
// 					junctionStart = end_read->endpoint;
// 					junctionEnd = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '+' && end_read->strand2 == '-')
// 				{
// 					junctionStart = end_read->endpoint;
// 					junctionEnd = end_read->endpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '+')
// 				{
// 					junctionStart = end_read->startpoint;
// 					junctionEnd = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '-')
// 				{
// 					junctionStart = end_read->startpoint;
// 					junctionEnd = end_read->endpoint2;
// 				}
			}
			else if (end_read->chrNum1 > end_read->chrNum2)
			{
				if (start_read->chrNum1 == end_read->chrNum2 && start_read->endpoint < end_read->startpoint2 || start_read->chrNum1 == end_read->chrNum1 && start_read->endpoint < end_read->startpoint)
				{
				}
				else if (start_read->chrNum1 == end_read->chrNum1 && start_read->startpoint > end_read->endpoint || start_read->chrNum1 == end_read->chrNum2 && start_read->startpoint > end_read->endpoint2)
				{
				} 
				else
				{
					abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
					return;
				}

// 				if (end_read->strand1 == '+' && end_read->strand2 == '+')
// 				{
// 					junctionEnd = end_read->endpoint;
// 					junctionStart = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '+' && end_read->strand2 == '-')
// 				{
// 					junctionEnd = end_read->endpoint;
// 					junctionStart = end_read->endpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '+')
// 				{
// 					junctionEnd = end_read->startpoint;
// 					junctionStart = end_read->startpoint2;
// 				}
// 				else if (end_read->strand1 == '-' && end_read->strand2 == '-')
// 				{
// 					junctionEnd = end_read->startpoint;
// 					junctionStart = end_read->endpoint2;
// 				}

				strcpy(chromosome, end_read->chromosome2);
				strcat(chromosome, end_read->chromosome1);
			}
			else
			{
				abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
				return;
			}

			if (end_read->strand1 == '+' && end_read->strand2 == '+')
			{
				junctionStart = end_read->endpoint;
				junctionEnd = end_read->startpoint2;
			}
			else if (end_read->strand1 == '+' && end_read->strand2 == '-')
			{
				junctionStart = end_read->endpoint;
				junctionEnd = end_read->endpoint2;
			}
			else if (end_read->strand1 == '-' && end_read->strand2 == '+')
			{
				junctionStart = end_read->startpoint;
				junctionEnd = end_read->startpoint2;
			}
			else if (end_read->strand1 == '-' && end_read->strand2 == '-')
			{
				junctionStart = end_read->startpoint;
				junctionEnd = end_read->endpoint2;
			}
		} 
		else
		{
			abnormal_reads << start_read->name << '\t' << start_read->chromosome1 << '\t' << start_read->fusionRead << '\t' << start_read->startpoint << '\t' << start_read->end << '\t' << end_read->chromosome1 << '\t' << end_read->fusionRead << '\t' << end_read->startpoint << '\t' << end_read->end << endl;  
			return;
		}
	}


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
		outputfile_chrname << chromosome << endl;

		sprintf(outpurfilename, "%sdata/PER/%s.txt", dirPrefix, chromosome);
		targetChrNum++;
		strcpy(targetChromosome[targetChrNum], chromosome);
		outputfile[targetChrNum].open(outpurfilename);
	}

	targetChrIndex = tmp;

	outputfile[targetChrIndex] << start_read->name << "\t"; 

	if (start_read->fusionRead == false && end_read->fusionRead == false)
	{
		outputfile[targetChrIndex] << "0\t" << start_read->chromosome1 << "\t" << start_read->strand1 << "\t" << start_read->startpoint << "\t" << start_read->end << "\t";
		outputfile[targetChrIndex] << "0\t" << end_read->chromosome1 << "\t" << end_read->strand1 << "\t" << end_read->startpoint << "\t" << end_read->end << endl;
	} 
	else
	{
		if (start_read->fusionRead == true && end_read->fusionRead == false)
		{
			outputfile[targetChrIndex] << "1\t" << start_read->chromosome1 << "\t" << start_read->startpoint << "\t" << start_read->endpoint << '\t' << start_read->chromosome2 << '\t' << start_read->startpoint2 << "\t" << start_read->endpoint2 << "\t" << start_read->strand1 << start_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << '\t';

			outputfile[targetChrIndex] << "0\t" << end_read->chromosome1 << "\t" << end_read->strand1 << "\t" << end_read->startpoint << "\t" << end_read->end << endl;
		}
		else if (start_read->fusionRead == false && end_read->fusionRead == true)
		{
			outputfile[targetChrIndex] << "0\t" << start_read->chromosome1 << "\t" << start_read->strand1 << "\t" << start_read->startpoint << "\t" << start_read->end << "\t";

			outputfile[targetChrIndex] << "1\t" << end_read->chromosome1 << "\t" << end_read->startpoint << "\t" << end_read->endpoint << '\t' << end_read->chromosome2 << '\t' << end_read->startpoint2 << "\t" << end_read->endpoint2 << "\t" << end_read->strand1 << end_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << endl;
		}
		else if (start_read->fusionRead == true && end_read->fusionRead == true)
		{
			outputfile[targetChrIndex] << "1\t" << start_read->chromosome1 << "\t" << start_read->startpoint << "\t" << start_read->endpoint << '\t' << start_read->chromosome2 << '\t' << start_read->startpoint2 << "\t" << start_read->endpoint2 << "\t" << start_read->strand1 << start_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << '\t';

			outputfile[targetChrIndex] << "1\t" << end_read->chromosome1 << "\t" << end_read->startpoint << "\t" << end_read->endpoint << '\t' << end_read->chromosome2 << '\t' << end_read->startpoint2 << "\t" << end_read->endpoint2 << "\t" << end_read->strand1 << end_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << endl;
		}
	}



// 	long dist;
// 	dist = end_read->startpoint - start_read->endpoint - 1;
// 
// 	if (dist >= 0 && dist < 50000)
// 	{
// 		allDistanceToFile << dist << endl;
// 	}
	return;
}


void outputPEfile()
{
	//extract all possible alignments

	int i, j, ii, jj;
	bool isSame = false;

	if (startlist_count * endlist_count < MAX_COMB_LIMIT)
	{
		//extract all
		for (i = 1; i <= startlist_count; i++)
		{
			for (j = 1; j <= endlist_count; j++)
			{
				if ((startlist[i]->endpoint < endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint < 400000)|| (startlist[i]->startpoint > endlist[j]->endpoint && startlist[i]->startpoint - endlist[j]->endpoint < 400000))
					outputPEreads(startlist[i], endlist[j]);
				if (startlist[i]->fusionRead == false && endlist[j]->fusionRead == false && strcmp(startlist[i]->chromosome1, endlist[j]->chromosome1) == 0)
					isSame = true;
			}
		}
		if (isSame == true)
			PERCnt_same++;
		else
			PERCnt_diff++;

		return;
	} 
	else //if (startlist_count * endlist_count < MAX_COMB_LIMIT * 10)
	{
		for (i = 1; i <= startlist_count; i++)
		{
			for (j = 1; j <= endlist_count; j++)
			{
//				if (startlist[i]->fusionRead == true || endlist[j]->fusionRead == true || (startlist[i]->fusionRead == false && endlist[j]->fusionRead == false && strcmp(startlist[i]->chromosome1, endlist[j]->chromosome1) == 0))
				if (startlist[i]->fusionRead == false && endlist[j]->fusionRead == false && strcmp(startlist[i]->chromosome1, endlist[j]->chromosome1) == 0)
				{
					if ((startlist[i]->endpoint < endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint < 400000)|| (startlist[i]->startpoint > endlist[j]->endpoint && startlist[i]->startpoint - endlist[j]->endpoint < 400000))
						outputPEreads(startlist[i], endlist[j]);
					isSame = true;
				} 
				/*else
				{
					for (ii = 1; ii <= startlist_count; ii++)
					{
						dropped_reads << startlist[ii]->name << "\t1\t" << startlist[ii]->chromosome1 << '\t' << startlist[ii]->startpoint << '\t' << startlist[ii]->end << endl;
					}
					for (jj = 1; jj <= endlist_count; jj++)
					{
						dropped_reads << endlist[jj]->name << "\t2\t" << endlist[jj]->chromosome2 << '\t' << endlist[jj]->startpoint << '\t' << endlist[jj]->end << endl;
					}
				}*/
			}
		}
		if (isSame == true)
			PERCnt_same++;

/*		for (ii = 1; ii <= startlist_count; ii++)
		{
			dropped_reads << startlist[ii]->name << "\t1\t" << startlist[ii]->chromosome1 << '\t' << startlist[ii]->startpoint << '\t' << startlist[ii]->end << endl;
		}
		for (jj = 1; jj <= endlist_count; jj++)
		{
			dropped_reads << endlist[jj]->name << "\t2\t" << endlist[jj]->chromosome1 << '\t' << endlist[jj]->startpoint << '\t' << endlist[jj]->end << endl;
		}*/
	}

	return;
}

void cleanup()
{
	reads* ptr;

	for (int i = 1; i <= startlist_count; i++)
	{
		ptr = startlist[i];
		delete ptr;
	}
	for (int j = 1; j <= endlist_count; j++)
	{
		ptr = endlist[j];
		delete ptr;
	}

	startlist_count = 0;
	endlist_count = 0;

	return;
}

void parse(char* inputfilename)
{
	char outpurfilename[500];
	sprintf(outpurfilename, "%sdata/PER/ChromosomeName.txt", dirPrefix);
	outputfile_chrname.open(outpurfilename);
	sprintf(outpurfilename, "%sdata/PER/DroppedReads.txt", dirPrefix);
	dropped_reads.open(outpurfilename);
	sprintf(outpurfilename, "%sdata/PER/AbnormalReads.txt", dirPrefix);
	abnormal_reads.open(outpurfilename);
	sprintf(outpurfilename, "%sdata/PER/OverlapReads.txt", dirPrefix);
	overlap_reads.open(outpurfilename);

	reads *ptr;
	ptr = NULL;

	startlist_count = 0;
	endlist_count = 0;

	ifstream inputfile;
	inputfile.open(inputfilename);
		
	char name[100];
	int start_end_switch, tmp, strand;
	bool writeswitch;
	writeswitch = false;
	string info, curLine;

	char field1[100], field2[100], field3[100], field4[100], field5[200], field6[100], field7[100], field8[100];

	for (tmp = 0; tmp < 100; tmp++)
	{
		name[tmp] = '\0';
	}

	inputfile >> name;


	while (name[1] != '\0')
	{
		//find start_end_switch
		for (tmp = 0; tmp < 100; tmp++)
		{
			if (name[tmp] == '/')
			{
				break;
			}
		}
		if (name[tmp+1] == '1')
		{
			start_end_switch = 1;
		} 
		else
		{
			start_end_switch = 2;
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
		
		
	
		if (start_end_switch == 1)
		{
			if (writeswitch == true)
			{
				outputPEfile();
				cleanup();
				writeswitch = false;
			}

			//start point
			startlist_count++;
			ptr = new reads;

			strcpy(ptr->name, name);

			if (field6[0] == '=' || field6[0] == '*')
			{
				//not fusion read
				ptr->fusionRead = false;
				if (strcmp(field1, "16") != 0)
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

			startlist[startlist_count] = ptr;
		} 
		else if (start_end_switch == 2)
		{
			//end point
			endlist_count++;
			ptr = new reads;

			strcpy(ptr->name, name);

			if (field6[0] == '=' || field6[0] == '*')
			{
				//not fusion read
				ptr->fusionRead = false;
				if (strcmp(field1, "16") != 0)
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
				//ptr->strand1 = field2[0];
				if (field2[0] == '+')
				{
					ptr->strand1 = '+';
				} 
				else
				{
					ptr->strand1 = '-';
				}
				ptr->startpoint = strtol(field3, NULL, 10);
				ptr->endpoint = strtol(field4, NULL, 10);

				strcpy(ptr->chromosome2, field5);
				//ptr->strand2 = field6[0];
				if (field6[0] == '+')
				{
					ptr->strand2 = '+';
				}
				else
				{
					ptr->strand2 = '-';
				}
				ptr->startpoint2 = strtol(field7, NULL, 10);
				ptr->endpoint2 = strtol(field8, NULL, 10);
			}

			ptr->compute_chrNum();

			endlist[endlist_count] = ptr;

			writeswitch = true;
		}
	


// 		for (tmp = 0; tmp < 100; tmp++)
// 		{
// 			name[tmp] = '\0';
// 		}
		name[1] = '\0';

		inputfile >> name;		

	}

	outputfile_chrname << "PER Count on same chromosome\t" <<  PERCnt_same << endl << "PER Count on different chromosome\t" << PERCnt_diff << endl;

	inputfile.close();
 	outputfile_chrname.close();
 	dropped_reads.close();
	abnormal_reads.close();

	for (tmp = 1; tmp <= targetChrNum; tmp++)
	{
		outputfile[tmp].close();
	}

	return;
}



int main(int argc, char* argv[])
{
	targetChrNum = 0;
	PERCnt_same = 0;
	PERCnt_diff = 0;

	if (argc != 3)
	{
		cout << argv[0] << "\t<filename>\t<target_path>" << endl;
		return 1;
	}

	char inputfilename[500];
//	sprintf(inputfilename, "./data/MCF_2x35.paired");

	strcpy(inputfilename, argv[1]);
//	allDistanceToFile.open("C:\\Users\\yin\\Desktop\\PEmatch\\human35\\allDistance.txt");

	strcpy(dirPrefix, argv[2]);

	parse(inputfilename);

//	cout << PERCnt_same << endl << PERCnt_diff << endl;

//	allDistanceToFile.close();

	return 0;
}



// 
// 
// void outputPEreads(reads* read1, reads* read2)
// {
// 	//output a PER
// 	char outpurfilename[100];
// 	char chromosome[50];
// 	int targetChrIndex, tmp;
// 	reads *start_read, *end_read;
// 	long junctionStart = -1, junctionEnd = -1;
// 
// 
// 	if (read1->fusionRead == false && read2->fusionRead == false)
// 	{
// 		if (read1->chrNum1 < read2->chrNum1)
// 		{
// 			start_read = read1;
// 			end_read = read2;
// 		} 
// 		else if (read1->chrNum1 > read2->chrNum1)
// 		{
// 			start_read = read2;
// 			end_read = read1;
// 		} 
// 		else
// 		{
// 			if (read1->endpoint < read2->startpoint)
// 			{
// 				start_read = read1;
// 				end_read = read2;
// 			} 
// 			else if (read1->startpoint > read2->endpoint)
// 			{
// 				start_read = read2;
// 				end_read = read1;
// 			}
// 			else
// 			{
// 				abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl; 
// 				return;
// 			}
// 		}
// 
// 		strcpy(chromosome, start_read->chromosome1);
// 		strcat(chromosome, end_read->chromosome1);
// 	} 
// 	else
// 	{
// 		if (read1->fusionRead == true && read2->fusionRead == false)
// 		{
// 			if (read1->chrNum1 < read1->chrNum2)
// 			{
// 				if (read2->chrNum1 == read1->chrNum1 && read2->endpoint < read1->startpoint || read2->chrNum1 == read1->chrNum2 && read2->endpoint < read1->startpoint2)
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				}
// 				else if (read2->chrNum1 == read1->chrNum2 && read2->startpoint > read1->endpoint2 || read2->chrNum1 == read1->chrNum1 && read2->startpoint > read1->endpoint)
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read1->strand1 == '+' && read1->strand2 == '+')
// 				{
// 					junctionStart = read1->endpoint;
// 					junctionEnd = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '+' && read1->strand2 == '-')
// 				{
// 					junctionStart = read1->endpoint;
// 					junctionEnd = read1->endpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '+')
// 				{
// 					junctionStart = read1->startpoint;
// 					junctionEnd = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '-')
// 				{
// 					junctionStart = read1->startpoint;
// 					junctionEnd = read1->endpoint2;
// 				}
// 
// 				strcpy(chromosome, read1->chromosome1);
// 				strcat(chromosome, read1->chromosome2);
// 			} 
// 			else if (read1->chrNum1 == read1->chrNum2 && read1->endpoint < read1->startpoint2)
// 			{
// 				if (read2->chrNum1 < read1->chrNum1 || (read2->chrNum1 == read1->chrNum1 && read2->endpoint < read1->startpoint))
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				}
// 				else if (read2->chrNum1 > read1->chrNum2 || (read2->chrNum1 == read1->chrNum2 && read2->startpoint > read1->endpoint2))
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read1->strand1 == '+' && read1->strand2 == '+')
// 				{
// 					junctionStart = read1->endpoint;
// 					junctionEnd = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '+' && read1->strand2 == '-')
// 				{
// 					junctionStart = read1->endpoint;
// 					junctionEnd = read1->endpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '+')
// 				{
// 					junctionStart = read1->startpoint;
// 					junctionEnd = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '-')
// 				{
// 					junctionStart = read1->startpoint;
// 					junctionEnd = read1->endpoint2;
// 				}
// 
// 				strcpy(chromosome, start_read->chromosome1);
// 				strcat(chromosome, end_read->chromosome1);
// 			}
// 			else if (read1->chrNum1 > read1->chrNum2)
// 			{
// 				if (read2->chrNum1 == read1->chrNum2 && read2->endpoint < read1->startpoint2 || read2->chrNum1 == read1->chrNum1 && read2->endpoint < read1->startpoint)
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				}
// 				else if (read2->chrNum1 == read1->chrNum1 && read2->startpoint > read1->endpoint || read2->chrNum1 == read1->chrNum2 && read2->startpoint > read1->endpoint2)
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read1->strand1 == '+' && read1->strand2 == '+')
// 				{
// 					junctionEnd = read1->endpoint;
// 					junctionStart = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '+' && read1->strand2 == '-')
// 				{
// 					junctionEnd = read1->endpoint;
// 					junctionStart = read1->endpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '+')
// 				{
// 					junctionEnd = read1->startpoint;
// 					junctionStart = read1->startpoint2;
// 				}
// 				else if (read1->strand1 == '-' && read1->strand2 == '-')
// 				{
// 					junctionEnd = read1->startpoint;
// 					junctionStart = read1->endpoint2;
// 				}
// 
// 				strcpy(chromosome, read1->chromosome2);
// 				strcat(chromosome, read1->chromosome1);
// 			}
// 			else
// 			{
// 				abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 				return;
// 			}
// 		}
// 		else if (read1->fusionRead == false && read2->fusionRead == true)
// 		{
// 			if (read2->chrNum1 < read2->chrNum2)
// 			{
// 				if (read1->chrNum1 == read2->chrNum1 && read1->endpoint < read2->startpoint || read1->chrNum1 == read2->chrNum2 && read1->endpoint < read2->startpoint2)
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				}
// 				else if (read1->chrNum1 == read2->chrNum2 && read1->startpoint > read2->endpoint2 || read1->chrNum1 == read2->chrNum1 && read1->startpoint > read2->endpoint)
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read2->strand1 == '+' && read2->strand2 == '+')
// 				{
// 					junctionStart = read2->endpoint;
// 					junctionEnd = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '+' && read2->strand2 == '-')
// 				{
// 					junctionStart = read2->endpoint;
// 					junctionEnd = read2->endpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '+')
// 				{
// 					junctionStart = read2->startpoint;
// 					junctionEnd = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '-')
// 				{
// 					junctionStart = read2->startpoint;
// 					junctionEnd = read2->endpoint2;
// 				}
// 
// 				strcpy(chromosome, read2->chromosome1);
// 				strcat(chromosome, read2->chromosome2);
// 			} 
// 			else if (read2->chrNum1 == read2->chrNum2 && read2->endpoint < read2->startpoint2)
// 			{
// 				if (read1->chrNum1 < read2->chrNum1 || (read1->chrNum1 == read2->chrNum1 && read1->endpoint < read2->startpoint))
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				}
// 				else if (read1->chrNum1 > read2->chrNum2 || (read1->chrNum1 == read2->chrNum2 && read1->startpoint > read2->endpoint2))
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read2->strand1 == '+' && read2->strand2 == '+')
// 				{
// 					junctionStart = read2->endpoint;
// 					junctionEnd = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '+' && read2->strand2 == '-')
// 				{
// 					junctionStart = read2->endpoint;
// 					junctionEnd = read2->endpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '+')
// 				{
// 					junctionStart = read2->startpoint;
// 					junctionEnd = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '-')
// 				{
// 					junctionStart = read2->startpoint;
// 					junctionEnd = read2->endpoint2;
// 				}
// 
// 				strcpy(chromosome, start_read->chromosome1);
// 				strcat(chromosome, end_read->chromosome1);
// 			}
// 			else if (read2->chrNum1 > read2->chrNum2)
// 			{
// 				if (read1->chrNum1 == read2->chrNum2 && read1->endpoint < read2->startpoint2 || read1->chrNum1 == read2->chrNum1 && read1->endpoint < read2->startpoint)
// 				{
// 					start_read = read1;
// 					end_read = read2;
// 				}
// 				else if (read1->chrNum1 == read2->chrNum1 && read1->startpoint > read2->endpoint || read1->chrNum1 == read2->chrNum2 && read1->startpoint > read2->endpoint2)
// 				{
// 					start_read = read2;
// 					end_read = read1;
// 				} 
// 				else
// 				{
// 					abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 					return;
// 				}
// 
// 				if (read2->strand1 == '+' && read2->strand2 == '+')
// 				{
// 					junctionEnd = read2->endpoint;
// 					junctionStart = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '+' && read2->strand2 == '-')
// 				{
// 					junctionEnd = read2->endpoint;
// 					junctionStart = read2->endpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '+')
// 				{
// 					junctionEnd = read2->startpoint;
// 					junctionStart = read2->startpoint2;
// 				}
// 				else if (read2->strand1 == '-' && read2->strand2 == '-')
// 				{
// 					junctionEnd = read2->startpoint;
// 					junctionStart = read2->endpoint2;
// 				}
// 
// 				strcpy(chromosome, read2->chromosome2);
// 				strcat(chromosome, read2->chromosome1);
// 			}
// 			else
// 			{
// 				abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 				return;
// 			}
// 		} 
// 		else
// 		{
// 			abnormal_reads << read1->name << '\t' << read1->chromosome1 << '\t' << read1->fusionRead << '\t' << read1->startpoint << '\t' << read1->end << '\t' << read2->chromosome1 << '\t' << read2->fusionRead << '\t' << read2->startpoint << '\t' << read2->end << endl;  
// 			return;
// 		}
// 	}
// 
// 
// 	for (tmp = 1; tmp <= targetChrNum; tmp++)
// 	{
// 		if (strcmp(targetChromosome[tmp], chromosome) == 0)
// 		{
// 			break;
// 		}
// 	}
// 
// 	if (tmp > targetChrNum)
// 	{
// 		//chromosome not found
// 		outputfile_chrname << chromosome << endl;
// 
// 		sprintf(outpurfilename, "./data/PER/%s.txt", chromosome);
// 		targetChrNum++;
// 		strcpy(targetChromosome[targetChrNum], chromosome);
// 		outputfile[targetChrNum].open(outpurfilename);
// 	}
// 
// 	targetChrIndex = tmp;
// 
// 	outputfile[targetChrIndex] << start_read->name << "\t"; 
// 
// 	if (start_read->fusionRead == false && end_read->fusionRead == false)
// 	{
// 		outputfile[targetChrIndex] << "0\t" << start_read->chromosome1 << "\t" << start_read->strand1 << "\t" << start_read->startpoint << "\t" << start_read->end << "\t";
// 		outputfile[targetChrIndex] << "0\t" << end_read->chromosome1 << "\t" << end_read->strand1 << "\t" << end_read->startpoint << "\t" << end_read->end << endl;
// 	} 
// 	else
// 	{
// 		if (start_read->fusionRead == true && end_read->fusionRead == false)
// 		{
// 			outputfile[targetChrIndex] << "1\t";
// 
// 			if (start_read->chrNum1 < start_read->chrNum2)
// 			{
// 				outputfile[targetChrIndex] << start_read->chromosome1 << "\t" << start_read->startpoint << "\t" << start_read->endpoint << '\t' << start_read->chromosome2 << '\t' << start_read->startpoint2 << "\t" << start_read->endpoint2 << "\t" << start_read->strand1 << start_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << '\t';
// 			} 
// 			else if (start_read->chrNum1 == start_read->chrNum2 && start_read->endpoint < start_read->startpoint2)
// 			{
// 				outputfile[targetChrIndex] << start_read->chromosome1 << "\t" << start_read->startpoint << "\t" << start_read->endpoint << '\t' << start_read->chromosome2 << '\t' << start_read->startpoint2 << "\t" << start_read->endpoint2 << "\t" << start_read->strand1 << start_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << '\t';
// 			}
// 			else if (start_read->chrNum1 > start_read->chrNum2)
// 			{
// 				outputfile[targetChrIndex] << start_read->chromosome2 << '\t' << start_read->startpoint2 << "\t" << start_read->endpoint2 << "\t" << start_read->chromosome1 << "\t" << start_read->startpoint << "\t" << start_read->endpoint << '\t' << start_read->strand2 << start_read->strand1 << "\t" << junctionStart << '\t' << junctionEnd << '\t';
// 			}
// 
// 			outputfile[targetChrIndex] << "0\t" << end_read->chromosome1 << "\t" << end_read->strand1 << "\t" << end_read->startpoint << "\t" << end_read->end << endl;
// 		}
// 		else if (start_read->fusionRead == false && end_read->fusionRead == true)
// 		{
// 			outputfile[targetChrIndex] << "0\t" << start_read->chromosome1 << "\t" << start_read->strand1 << "\t" << start_read->startpoint << "\t" << start_read->end << "\t";
// 
// 			outputfile[targetChrIndex] << "1\t";
// 
// 			if (end_read->chrNum1 < end_read->chrNum2)
// 			{
// 				outputfile[targetChrIndex] << end_read->chromosome1 << "\t" << end_read->startpoint << "\t" << end_read->endpoint << '\t' << end_read->chromosome2 << '\t' << end_read->startpoint2 << "\t" << end_read->endpoint2 << "\t" << end_read->strand1 << end_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << endl;
// 			} 
// 			else if (end_read->chrNum1 == end_read->chrNum2 && end_read->endpoint < end_read->startpoint2)
// 			{
// 				outputfile[targetChrIndex] << end_read->chromosome1 << "\t" << end_read->startpoint << "\t" << end_read->endpoint << '\t' << end_read->chromosome2 << '\t' << end_read->startpoint2 << "\t" << end_read->endpoint2 << "\t" << end_read->strand1 << end_read->strand2 << "\t" << junctionStart << '\t' << junctionEnd << endl;
// 			}
// 			else if (end_read->chrNum1 > end_read->chrNum2)
// 			{
// 				outputfile[targetChrIndex] << end_read->chromosome2 << '\t' << end_read->startpoint2 << "\t" << end_read->endpoint2 << "\t" << end_read->chromosome1 << "\t" << end_read->startpoint << "\t" << end_read->endpoint << '\t' << end_read->strand2 << end_read->strand1 << "\t" << junctionStart << '\t' << junctionEnd << endl;
// 			}
// 		}
// 	}
// 
// 
// 
// 	// 	long dist;
// 	// 	dist = end_read->startpoint - start_read->endpoint - 1;
// 	// 
// 	// 	if (dist >= 0 && dist < 50000)
// 	// 	{
// 	// 		allDistanceToFile << dist << endl;
// 	// 	}
// 	return;
// }
