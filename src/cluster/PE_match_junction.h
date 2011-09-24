/////////// Nov 2nd///////////////
//Update
//1. (From last version) Change head-insertion to tail-insertion for range lists. This is for the convenience to output bed file.
//2. Compute the pe_distance
///////////////////////////////////


#include <fstream>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

char dirPrefix[1000];



const long MAX_LENGTH = 10000000;
const long MAX = 2000000000;
const long MAX_JUNCTION_NUM = 1000000;
const long MAX_PEREAD_DISTANCE = 2000000000;
const int MAX_PATH_NUM = 100;

const long SCORE_ADD = 500;
//index starts from 1 for all arrays

class range
{
public:
	long start;
	long end;

	range* next;

	range();
	void range_copy(range *orig);
};

class used_junction
{
public:
	long junctionID;
	used_junction *next;

	used_junction();
};

class pe_distance
{
public:
	long existing_distance;
	long milestone; //the positions till milestone of this reads have been computed 
	used_junction *junc;
	pe_distance *next;

	pe_distance(long start, long end);
};

class paired_end
{
public:
	char chromosome_start[30];
	char chromosome_end[30];
	char strand_start[5];
	char strand_end[5];
	int range_count_start;
	int range_count_end;
	range* start;
	range* end;
	range* tail_start; // tail pointer for start range list
	range* tail_end;//tail pointer for end range list
	long startrange_low;
	long startrange_high;
	long endrange_low;
	long endrange_high;

	pe_distance* dist;

	paired_end();
	void compute_distance();
};


class cluster
{
	//current considering range
	//used in PEmatch

public:
	char chromosome_1[30];
	char chromosome_2[30];
	char strand_start[5];
	char strand_end[5];
	int range_count_1;
	int range_count_2;
	int fulllist_count;
	range* endpoint1;
	range* endpoint2;
	range* tail_end1; 
	range* tail_end2;
	range* fulllist;
	range* tail_fulllist;
	long end1_range_low;
	long end1_range_high;
	long end2_range_low;
	long end2_range_high;

	cluster(paired_end* startPE);
	void merge_range(paired_end* currentPE);//maintain the range of endpoints after adding an overlap read
	bool searchANDmodify(paired_end* currentPE, long *overlap_start, long *overlap_end); //search whether the currentPE is in this range. if so, modify the range. return the overlap num
	bool build_fulllist(); //merge two lists. if two chromosomes are the same, return true (merge to fulllist). otherwise, return false
};

class junction
{
public:
	char chromosome[30];
	long start_start;
	long end_end;
	long start;
	long end;
	char name[30];
	char strand;
	long thickStart;
	long thickEnd;
	int blockCount;
	char blockSizes[50];
	char blockStarts[50];

	junction();
	void compute_range();
};


//used in acknowledge
long current_set[MAX_LENGTH][3];
long current_set_index;
long last_end;

//////////////////////



int deviation_end;
int deviation_distance;






paired_end* PE_list[MAX_LENGTH]; //input from file, all the paired end reads in original order

long PE_list_ordered[MAX_LENGTH]; //index list of PE_list after sorting

long sort_key[MAX_LENGTH];

long stack[MAX_LENGTH];

long PE_num;

long min_point;
long max_point;
ofstream PE_match_outputfile;
ofstream PE_match_region_outputfile;
ofstream PE_match_region_detail_outputfile;
ofstream PE_match_BED_outputfile;
ofstream cluster_distances;
ofstream PE_match_BED_outputfile_1; //when different chromosomes
ofstream PE_match_BED_outputfile_2; //when different chromosomes

char targetChromosome[100];
char inputfilename[1000];
char chromosome[100];

junction* junction_list[MAX_JUNCTION_NUM]; //record all the junctions
long junction_support[MAX_JUNCTION_NUM]; //#reads supporting each junction
long junction_num;

string junction_title;
bool sameChromosome;


ofstream diffChrClusters;

ofstream finalResult;