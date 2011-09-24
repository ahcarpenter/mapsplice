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



///////////////////////////
///Nov. 22nd
//1. cluster based on paths. one path corresponds to one cluster
//2. remove blocks with no reads
//3. if a block has a range larger than delta(e.g. 50bp) without reads, remove this range and disconnect the seperated parts
//junction file and PER file need to be parsed and separated into files for each chromosome
//junction file has been sorted
//
///Dec. 7th
//1. take PER file as only input, extract junctions from PER file
//
///Dec.17th
//1. output path information in sam file
//
///Dec.21st
//1. consider multiple alignments for a PER (!!assume that the alignments for a same PER are located together)
//2. EM algorithm
//
///Jan 2nd
//1. replace quicksort with mergesort for stable sorting
//2. refine & accelerate EM
///Jan 3rd
//1. collect 2 types of junction support
//////////////////////////////////////////////////////////////////////////

/************************************************************************/
/* Latest version for fusion                                            */
/************************************************************************/
///Feb. 19th
//1. strand information is added; direction needs to be considered; junctions have no direction;
//reads have directions, and directions of different reads are independent;
//
//2. Block splitting is updated. All graph edges now have links on both direction.
//
//3. Computing block paths is updated. Path can be not only from low end to high end but also reverse.
//
//4. The input format is changed.
//
//5. Junction ID: when > 0, refer to junction list; when = -1, means getting to the block from
//   the left, i.e. starting from the low end of the block; when = -2, means getting to the 
//   block from the right, i.e. the high end
/************************************************************************/
/* March 1st, fix loop bug                                              */
/************************************************************************/




///*
#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
//*/
/*
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
*/


using namespace std;

/************************************************************************/
/* Configuration file for PERalign                                      */
/************************************************************************/

/************************************************************************/
/* IMPORTANT SETTINGS                                                   */
/************************************************************************/

//distribution of the PER distances
double PERPathDistanceMean; //mean of the distance distribution
double PERPathDistanceSD; //standard deviation of the distance distribution

//size of input data
long MAX_JUNCTION_NUM; //maximum number of input junctions
long MAX_LENGTH; //maximum number of input reads. (2 times input PER alignments)


//setting
int MAX_PATH_NUM; //size of the path queue.
double MIN_PATH_PROBABILITY; //minimum probability of paths. paths with probability less than MIN_PATH_PROBABILITY will be filtered.

//EM hard decision
double MIN_HD_PER_PROB; //minimum read probability to be considered in hard decision -- sequenced segment
double MIN_HD_PATH_PROB; //minimum path probability to be considered in hard decision -- unsequenced segment



/************************************************************************/
/* SETTINGS NOT OFTEN CHANGED                                           */
/************************************************************************/

const long MAX_CHR_LENGTH = 1000000000; //maximal chromosome length. do not change if the chromosome length is less than 10^9 base pair.


//EM
const long MAX_EM_PARA = 500000; //max number of members involved in the EM
const long MAX_EM_MEMBER = 1000; //max number of members using EM at the same time

const int MAX_PATH_CHANGE_DIRECTION = 2; //maximum number of changes on strand of path. Set as 2 to eliminate loops.

//convergence condition for EM
const int MAX_ITER_COUNT = 10; //maximum iteration count
const double MIN_UPDATE_RATE = 0.0001; //minimum updating rate. EM stops iteration if (improvement on log likelihood)/(previous log likelihood) < MIN_UPDATE_RATE. 
const double HILLCLIMBING_STEP = 0.05; //step size for hill climbing.


class block;

class blockGraphEdge
{
public:
	block *linkedBlk;
	long junctionID; //index of the junction used by this edge. -1 indicates no junction used (connected directly)
	blockGraphEdge *nextEdge;
	long arrivedPosition; //the arrived position if passing this edge
	long leavingPosition; //the leaving position for this edge
//	long skippedLength; //sequence length skipped by this edge

	blockGraphEdge();
};

class linkedGraphEdge
{
public:
	blockGraphEdge *linkedEdge; //edge pointing to this block
	linkedGraphEdge *next;

	linkedGraphEdge();
};

class block
{
public:
	bool supported; //true if there are reads supporting this block
	long start_position;
	long end_position;
	char chromosome[6];

	block *next; //physical linkage
	blockGraphEdge *edge; //multiple edges. logical linkage
	linkedGraphEdge *linkEdge; //edges pointing to this block

	block();
	~block();
};


class junction
{
public:
	char chromosome_start[6];
	char chromosome_end[6];
	//	char name[30];
	//	long start_start;
	long start_end;
	long end_start;
	//	long end_end;
	//	char strand;
	//	int blockCount;
	//	char blockSizes[50];
	//	char blockStarts[50];
	double probability;
	double probBackup;
	bool single; //true for junctions from single reads, false for junctions from paired reads
	bool fusion; //true for fusion junctions, false for non-fusion junctions
	bool isTrue; //is true according to hard decision

	bool strand_start; //strand at the start read, true means +
	bool strand_end; //strand at the end read, true means +

	bool exchange;

	int support; //support from the junction file

	junction();
	//	void compute_range();
};

class pathJunction
{
public:
	long junctionID;
	pathJunction *next;

	pathJunction();
};


class path
{
public:
	block *destBlk; //destination block
	//long skippedLength; //skipped length if taking this path
	long pathLength; //length of the path
	long endposition; //end position of the path
	long startposition; //start position of the path
	pathJunction *junctionsInPath;
	path* next;
	int changeOfDirection; //number of changes of the path direction
	bool currentDirection; //current direction of the path, bool for positive strand
	long lastJunctionID; //last time memory of used junction ID

	//for storing path info
	char info[130]; //path info using xxMxxN format
	long lastUnwritten_M; //keep the unwritten M-segment, do this for merging connected M-segs

	//bool strand; //+ strand from left (low to high); - strand from right (high to low)

	path();
	~path();
};

class pathCluster
{
public:
	block *startBlk;
	block *endBlk;
	path *clusterPaths;
	pathCluster *next;

	pathCluster();
};



class range
{
public:
	long start;
	long end;

	range* next;

	range();
	void range_copy(range *orig);
};

class readJunction
{
public:
	//junction of the read region
	long junctionIndex;
	readJunction *next;

	readJunction();
};

class reads
{
public:
	bool startRead; //true if it is a start read, false for a end read
	block* onBlock; //the block where the read is. the position of the read is the end-position for a start read or the start-position for a end read
	char chromosome[6];
	//int strand;
//	int rangeCount;
//	range* readsRange;
//	range* rangeTail; // tail pointer for range list
	long rangeLow;
	long rangeHigh;
	char end[50]; //end information "xxMxxNxxM"

	//two strand strings for fusion reads. For normal reads, the two strand should be the same.
	bool strandLow; //true for +
	bool strandHigh; //true for +
	long crucialPoint; //point representing the 'position', effected by start/end read and strand
	bool crucialStrand; //associated with crucialPoint, the strand at the representing point

	readJunction* readJunc;

	reads();
	~reads();
};

class PERpath
{
	//possible path of a paired-end reads
public:
	path* pathInfo;
	long distance;
	double distanceProb; //probability of the distance according to the empirical distribution
	double probability; //probability that the path is true given the PER, alignment and the path
	double probBackup;
	double theta; //probability that the PER chooses this path on this alignment
	PERpath *next;
//	char infoString[200];

	PERpath();
	void compute_prob(); //compute the probability that the path is true 
};

class PERalignment
{
public:	
	reads *startRead;
	reads *endRead;
	PERpath* possiblePath;
	int pathNum;
	double maxPathProb;
	long perIndex; //index of the PER 
	double probability; //probability of the PER on this alignment
	double probBackup;
	double tau; //probability that the PER is aligned to this position

	PERalignment();
	void compute_prob(); //compute the probability that the PER is true given the alignment
};


class PER
{
public:
	char name[100];
	long firstAlignmentIndex; //index of the first alignment in the PER alignment list
	int numAlignment; //number of alignments for this PER
	double probability; //probability of the PER
	double probBackup;

	PER();
	void compute_prob(); //compute the probability that the PER is true
};


/************************************************************************/
/* GLOBAL VARIABLES                                                     */
/************************************************************************/
char targetChr_start[10];
char targetChr_end[10];
char targetChr[20];
bool CROSSCHROMOSOME; //whether start chromosome and end end chromosome are different. true for different and false for the same
bool writeSameChrSplice; //if cross-chromosome, true for writing splice events on the same chromosome together with cross-chromosomal splice junctions, false for only writing cross-chromosomal splice junctions


//array index starts from 1
string junction_title;

junction** junctionList; //list for all junctions
long* orderedJunctionList; //index list of junctionList after sorting
long* sortKey_Junction;
long junctionNum;
//long* junctionSupportPath;
//double* junctionSupportProbPath;
long* junctionSupportMPPath;
double* junctionSupportProbMPPath;
long* junctionSupportRead;
double* junctionSupportProbRead;

block* chrBlocks; //blocks for the chromosome
block* chrBlocks_tail; //tail pointer to the blocks of the chromosome
block* chrBlocks_2ndStart; //pointer to the start block of the 2nd chromosome
long chrBlockNum; //number of blocks 


pathCluster *pathClusterList;
pathCluster *pathClusterListTail;

//index starts from 1 for readList and perAlignList
reads** readList; //input from file, all the reads in original order
long* orderedReadList; //index list of readList after sorting
long* sortKey_read;
long readNum;

PERalignment** perAlignList; //all the paired-end reads alignments in the original order
long* orderedPERalignlist; //index list of perAlignList after sorting
long* sortKey_PERalign;
long perAlignNum;

PER** perList; //list of all the PERs 
long perNum;


//PER alignments to be processed with EM
PERalignment** perAlignList_EM; 
long perAlignEMlist_Start;
long perAlignEMlist_End;



//long stack[MAX_LENGTH]; //used for quick sort

long* mergeSort_Larray;
long* mergeSort_Rarray;
long* mergeSort_LorderedList;
long* mergeSort_RorderedList;

//path queue for computing all paths
path** pathQueue; 
long pathQueueTail;
long pathQueueHead;
long pathQueueMaxLength;






//pointers to memorize searched results last time
block *lastJunctionSearchBlock;


//distribution of the PER distances
double PERPathDistanceVAR; 
double PERPathDistanceProbScalar;
double MAX_PATH_LENGTH;



/************************************************************************/
/* Output Stream                                                        */
/************************************************************************/
ofstream outputPathSummary;

ofstream outputEMdetail;
ofstream outputEMprobPER;
ofstream outputEMprobJunction;

ofstream outputSummary;
ofstream outputFusion;

clock_t runTimeStart;
double runTime;


void output_PERpath_prob();
void output_junctionProb();
void output_PERpathJunction_RAW();

char dirPrefix[500];


long statNumRead;
long statNumSingleRead;
long statNumPairRead;
long statNumReadHardMap; //PERs with multiple alignments or unique alignment of distance larger than 200bp
long statNumReadPaired;
long statNumReadAbnormalDistance; //distance larger than normal range (e.g. 200bp)
