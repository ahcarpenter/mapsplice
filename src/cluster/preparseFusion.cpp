#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib> 

using namespace std;

const long MAX = 2000000000;
const int MAX_CHOICE_NUM = 1000000;
const int MAX_COMB_LIMIT = 10;


char targetChromosome[1000][50];
int targetChrNum;
ofstream outputfile[1000];
ofstream outputfile_chrname;
ofstream dropped_reads;

char dirPrefix[1000];

class reads
{
public:
	char name[200];
	int chrNum; //chromosome number
	char chromosome[50];
	char strand[5];
	long startpoint;
	long endpoint;
	char end[100];

	void compute_endpoint();
	void compute_chrNum();
};

reads* startlist[MAX_CHOICE_NUM];
reads* endlist[MAX_CHOICE_NUM];
int startlist_count;
int endlist_count;

ofstream allDistanceToFile;


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

	for (int i = 3; chromosome[i] != '\0'; i++)
	{
		tmp = tmp * 10 + chromosome[i] - 48;
	}

	chrNum = tmp;

	return;
}




void outputPEreads(reads* start_read, reads* end_read)
{
	//output a single read
	char outpurfilename[1000];
	char chromosome[50];
	int targetChrIndex, tmp;


	strcpy(chromosome, start_read->chromosome);
	strcat(chromosome, end_read->chromosome);

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

		sprintf(outpurfilename, "%sdata/parsedPER/%s.txt", dirPrefix, chromosome);
		targetChrNum++;
		strcpy(targetChromosome[targetChrNum], chromosome);
		outputfile[targetChrNum].open(outpurfilename);
	}

	targetChrIndex = tmp;

	outputfile[targetChrIndex] << start_read->name << "\t" << start_read->chromosome << "\t" << start_read->strand << "\t" << start_read->startpoint << "\t" << start_read->end << "\t" << end_read->chromosome << "\t" << end_read->strand << "\t" << end_read->startpoint << "\t" << end_read->end << endl; 

// 	long dist;
// 	dist = end_read->startpoint - start_read->endpoint - 1;
// 
// 	if (dist >= 0 && dist < 50000)
// 	{
// 		allDistanceToFile << dist << endl;
// 	}
	return;
}

// void outputPEfile()
// {
// 	//unique alignment
// 
// 	int i, j, best_i, best_j;
// 	long distance_i_j, best_distance;
// 
// 	best_distance = MAX;
// 	best_i = 1;
// 	best_j = 1;
// 	
// 	if (startlist_count * endlist_count > MAX_COMB_LIMIT)
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			dropped_reads << startlist[i]->name << "\t1\t" << startlist[i]->chromosome << '\t' << startlist[i]->startpoint << '\t' << startlist[i]->end << endl;
// 		}
// 		for (j = 1; j <= endlist_count; j++)
// 		{
// 			dropped_reads << endlist[j]->name << "\t2\t" << endlist[j]->chromosome << '\t' << endlist[j]->startpoint << '\t' << endlist[j]->end << endl;
// 		}
// 
// 		return;
// 	}
// 	
// 	for (i = 1; i <= startlist_count; i++)
// 	{
// 		for (j = 1; j <= endlist_count; j++)
// 		{
// // 			if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// // 			{
// // 				//not on the same chromosome
// // 				distance_i_j = MAX;
// // 			} 
// // 			else
// // 			{
// // 				//on the same chromosome
// // 				if (startlist[i]->endpoint < endlist[j]->startpoint)
// // 				{
// // 					//i is start and j is end
// // 					distance_i_j = compute_distance(startlist[i], endlist[j]);
// // 				}
// // 				else if (startlist[i]->startpoint >endlist[j]->endpoint)
// // 				{
// // 					//i is end and j is start
// // 					distance_i_j = compute_distance(endlist[j], startlist[i]);
// // 				}
// // 				else
// // 				{
// // 					//overlap
// // 					distance_i_j = 0;
// // 				}
// // 			}
// // 
// // 			if (distance_i_j < best_distance)
// // 			{
// // 				best_i = i;
// // 				best_j = j;
// // 			}
// 
// 			if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 			{
// 				//not on the same chromosome
// 				outputPEreads(startlist[i], endlist[j]);
// 			} 
// 			else
// 			{
// 				//on the same chromosome
// 				if (startlist[i]->startpoint <= endlist[j]->startpoint)
// 				{
// 					//i is start and j is end
// 					outputPEreads(startlist[i], endlist[j]);
// 				}
// 				else
// 				{
// 					//i is end and j is start
// 					outputPEreads(endlist[j], startlist[i]);
// 				}
// 			}
// 		}		
// 	}
// 
// 	//output reads best_i and best_j
// 	//outputPEreads(startlist[best_i], endlist[best_j])
// 
// 	return;
// }

// void outputPEfile()
// {
// 	//multiple alignment, keep alignments on the same chr
// 
// 	int i, j, ii, jj;
// 
// 	if (startlist_count * endlist_count <= MAX_COMB_LIMIT)
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			for (j = 1; j <= endlist_count; j++)
// 			{
// 
// 				if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 				{
// 					//not on the same chromosome
// 					//outputPEreads(startlist[i], endlist[j]);
// 					if (startlist[i]->chrNum < endlist[j]->chrNum)
// 					{
// 						outputPEreads(startlist[i], endlist[j]);
// 					} 
// 					else
// 					{
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				} 
// 				else
// 				{
// 					//on the same chromosome
// 					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
// 					{
// 						//i is start and j is end
// 						outputPEreads(startlist[i], endlist[j]);
// 					}
// 					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
// 					{
// 						//i is end and j is start
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				}
// 			}
// 		}
// 
// 		return;
// 	} 
// 	else
// 	{
// 		for (i = 1; i <= startlist_count; i++)
// 		{
// 			for (j = 1; j <= endlist_count; j++)
// 			{
// 				if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
// 				{
// 					//not on the same chromosome
// // 					for (ii = 1; ii <= startlist_count; ii++)
// // 					{
// // 						dropped_reads << startlist[ii]->name << "\t1\t" << startlist[ii]->chromosome << '\t' << startlist[ii]->startpoint << '\t' << startlist[ii]->end << endl;
// // 					}
// // 					for (jj = 1; jj <= endlist_count; jj++)
// // 					{
// // 						dropped_reads << endlist[jj]->name << "\t2\t" << endlist[jj]->chromosome << '\t' << endlist[jj]->startpoint << '\t' << endlist[jj]->end << endl;
// // 					}
// 				} 
// 				else
// 				{
// 					//on the same chromosome
// 					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
// 					{
// 						//i is start and j is end
// 						outputPEreads(startlist[i], endlist[j]);
// 					}
// 					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
// 					{
// 						//i is end and j is start
// 						outputPEreads(endlist[j], startlist[i]);
// 					}
// 				}
// 			}
// 		}
// 	}
// 
// 	return;
// }

void outputPEfile()
{
	/************************************************************************/
	/* Specifically for fusion clusters!!!                                  */
	/************************************************************************/
	int i, j, ii, jj;
	bool interesting = true;

	for (i = 1; i <= startlist_count; i++)
	{
		for (j = 1; j <= endlist_count; j++)
		{
			if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) == 0 && (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint < 50000 || endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint < 50000))
			{
				interesting = false;
			} 
		}
	}

	if (interesting == true && startlist_count * endlist_count == 1)
	{
		for (i = 1; i <= startlist_count; i++)
		{
			for (j = 1; j <= endlist_count; j++)
			{
				if (strcmp(startlist[i]->chromosome, endlist[j]->chromosome) != 0)
				{
					//not on the same chromosome
					//outputPEreads(startlist[i], endlist[j]);
//					if (startlist[i]->chrNum < endlist[j]->chrNum)
//					{
						outputPEreads(startlist[i], endlist[j]);
//					} 
//					else
//					{
//						outputPEreads(endlist[j], startlist[i]);
//					}
				} 
				else
				{
					//on the same chromosome
//					if (startlist[i]->endpoint <= endlist[j]->startpoint && endlist[j]->startpoint - startlist[i]->endpoint >= 50000)
//					{
						//i is start and j is end
						outputPEreads(startlist[i], endlist[j]);
//					}
//					else if (endlist[j]->endpoint <= startlist[i]->startpoint && startlist[i]->startpoint - endlist[j]->endpoint >= 50000)
//					{
//						//i is end and j is start
//						outputPEreads(endlist[j], startlist[i]);
//					}
				}
			}
		}

		return;
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
	char outpurfilename[1000];
	sprintf(outpurfilename, "%sdata/parsedPER/ChromosomeName.txt", dirPrefix);
	outputfile_chrname.open(outpurfilename);
	sprintf(outpurfilename, "%sdata/parsedPER/DroppedReads.txt", dirPrefix);
//	dropped_reads.open(outpurfilename);

	reads *ptr;
	ptr = NULL;

	startlist_count = 0;
	endlist_count = 0;

	ifstream inputfile;
	inputfile.open(inputfilename);
		
	char tmpChar[100], name[200], prevName[200];
	int start_end_switch, tmp, i;
	bool writeswitch;
	writeswitch = false;
	string info;

	char field1[100], field2[100], field3[100], field4[100], field5[100], field6[100], field7[100], field8[100];


	for (tmp = 0; tmp < 100; tmp++)
	{
		tmpChar[tmp] = '\0';
		name[tmp] = '\0';
		prevName[tmp] = '\0';
	}

	inputfile >> tmpChar;
	sprintf(prevName, "xxx");

	while (tmpChar[1] != '\0')
	{
		//find start_end_switch
		for (tmp = 0; tmp < 100; tmp++)
		{
			if (tmpChar[tmp] == '~')
			{
				break;
			}
		}
		tmp++;
		for (i = 0; tmp < 100; tmp++, i++)
		{
			if (tmpChar[tmp] == '/')
			{
				break;
			}
			name[i] = tmpChar[tmp];
		}
		name[i] = '\0';
		if (tmpChar[tmp+1] == '1')
		{
			start_end_switch = 1;
		} 
		else
		{
			start_end_switch = 2;
		}
		
		inputfile >> field1;
		inputfile >> field2;
		inputfile >> field3;
		inputfile >> field4;
		inputfile >> field5;
		inputfile >> field6;
		inputfile >> field7;
		inputfile >> field8;
		getline(inputfile, info);

		if (strcmp(prevName, name) != 0)
		{
			if (startlist_count * endlist_count > 0)
			{
				outputPEfile();
			}
			strcpy(prevName, name);
			cleanup();
			writeswitch = false;
		}

		if (start_end_switch == 1)
		{
			//start point
			if (field6[0] == '*' || field6[0] == '=')
			{
				//not fusion read				
				startlist_count++;
				ptr = new reads;

				strcpy(ptr->name, name);

				strcpy(ptr->strand, field1);
				strcpy(ptr->chromosome, field2);
				ptr->startpoint = strtol(field3, NULL, 10);
				strcpy(ptr->end, field5);

				ptr->compute_endpoint();
				ptr->compute_chrNum();

				startlist[startlist_count] = ptr;
			} 
			else
			{
				//fusion read, discard
			}
		} 
		else if (start_end_switch == 2)
		{
			//end point
			if (field6[0] == '*' || field6[0] == '=')
			{
				//not fusion read				
				endlist_count++;
				ptr = new reads;

				strcpy(ptr->name, name);
				
				strcpy(ptr->strand, field1);
				strcpy(ptr->chromosome, field2);
				ptr->startpoint = strtol(field3, NULL, 10);
				strcpy(ptr->end, field5);

				ptr->compute_endpoint();
				ptr->compute_chrNum();

				endlist[endlist_count] = ptr;
			} 
			else
			{
				//fusion read, discard
			}

			writeswitch = true;
		}

// 		for (tmp = 0; tmp < 100; tmp++)
// 		{
// 			name[tmp] = '\0';
// 		}
		tmpChar[1] = '\0';

		inputfile >> tmpChar;		
	}


	inputfile.close();
	outputfile_chrname.close();
// 	dropped_reads.close();

	for (tmp = 1; tmp <= targetChrNum; tmp++)
	{
		outputfile[tmp].close();
	}

	return;
}



int main(int argc, char* argv[])
{
	targetChrNum = 0;

	if (argc != 3)
	{
		cout << argv[0] << "\t<filename>\t<targetPath>" << endl;
		return 1;
	}

	char inputfilename[1000];
//	sprintf(inputfilename, "/sdb1/mapsplice/hmy_111/MapSplice_1.11/50bp_0520/merged_paired.sam");

	strcpy(inputfilename, argv[1]);
	strcpy(dirPrefix, argv[2]);
//	allDistanceToFile.open("C:\\Users\\yin\\Desktop\\PEmatch\\human35\\allDistance.txt");

	parse(inputfilename);

//	allDistanceToFile.close();

	return 0;
}