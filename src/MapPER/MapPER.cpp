#include "MapPER.h"


/************************************************************************/
/* Construction Functions                                               */
/************************************************************************/

blockGraphEdge::blockGraphEdge()
{
	linkedBlk = NULL;
	junctionID = -1;
	nextEdge = NULL;
	arrivedPosition = 0;
	leavingPosition = 0;
	//	skippedLength = 0;

	return;
}

linkedGraphEdge::linkedGraphEdge()
{
	linkedEdge = NULL;
	next = NULL;

	return;
}

block::block()
{
	supported = true; //! true if not eliminating paths with unsupported blocks

	start_position = 0;
	end_position = 0;
	next = NULL;
	edge = NULL;
	linkEdge = NULL;

	return;
}

block::~block()
{
	blockGraphEdge *curEdge, *nextEdge;

	curEdge = edge;
	while (curEdge != NULL)
	{
		nextEdge = curEdge->nextEdge;
		delete curEdge;
		curEdge = nextEdge;
	}

	linkedGraphEdge *curlEdge, *nextlEdge;

	curlEdge = linkEdge;
	while (curlEdge != NULL)
	{
		nextlEdge = curlEdge->next;
		delete curlEdge;
		curlEdge = nextlEdge;
	}
}

junction::junction()
{
	//	start_start = 0;
	start_end = 0;
	end_start = 0;
	//	end_end = 0;
	//	blockCount = 0;
	probability = 1.0;
	probBackup = 1.0;
	single = false;
	fusion = false;
	isTrue = false;

	strand_start = true;
	strand_end = true;

	exchange = false;

	support = 1;
	return;
}


// void junction::compute_range()
// {
// 	long i, blksize = 0;
// 
// 	for (i = 0; blockSizes[i] != '\0'; i++)
// 	{
// 		if (blockSizes[i] == ',')
// 		{
// 			break;
// 		} 
// 		else
// 		{
// 			blksize = blksize * 10 + blockSizes[i] - 48;
// 		}
// 	}
// 
// 	start_start = start_end - blksize;
// 
// 	blksize = 0;
// 	for (i++; blockSizes[i] != '\0'; i++)
// 	{
// 		if (blockSizes[i] == ',')
// 		{
// 			break;
// 		} 
// 		else
// 		{
// 			blksize = blksize * 10 + blockSizes[i] - 48;
// 		}
// 	}
// 
// 	end_end = end_start + blksize;
// 
// 	return;
// }

pathJunction::pathJunction()
{
	junctionID = -1;
	next = NULL;

	return;
}

path::path()
{
	destBlk = NULL;
	//skippedLength = 0;
	pathLength = 0;
	junctionsInPath = NULL;
	next = NULL;
	endposition = 0;
	startposition = 0;
	changeOfDirection = 0;
	lastUnwritten_M = 0;

	return;
}

path::~path()
{
	destBlk = NULL;
	junctionsInPath = NULL;
	next = NULL;
}

pathCluster::pathCluster()
{
	startBlk = NULL;
	endBlk = NULL;
	clusterPaths = NULL;
	next = NULL;

	return;
}

range::range()
{
	start = 0;
	end = 0;

	next = NULL;

	return;
}

void range::range_copy(range *orig)
{
	start = orig->start;
	end = orig->end;

	return;
}

reads::reads()
{
	startRead = true;
	onBlock = NULL;
	//	rangeCount = 0;
	//	readsRange = NULL;
	//	rangeTail = NULL;
	rangeLow = MAX_CHR_LENGTH * 2;
	rangeHigh = 0;
	readJunc = NULL;
	strandLow = true;
	strandHigh = true;

	crucialPoint = 0;
	crucialStrand = true;

	return;
}

reads::~reads()
{
	readJunction *curJunc, *nextJunc;

	curJunc = readJunc;
	while (curJunc != NULL)
	{
		nextJunc = curJunc->next;
		delete curJunc;
		curJunc = nextJunc;
	}
}

PERpath::PERpath()
{
	pathInfo = NULL;
	distance = 0;
	distanceProb = 0.0;
	probability = 0.0;
	theta = 0.0;
	next = NULL;

	return;
}

PERalignment::PERalignment()
{
	startRead = NULL;
	endRead = NULL;
	possiblePath = NULL;
	pathNum = 0;
	perIndex = 0;
	probability = 0.0;
	probBackup = 0.0;
	maxPathProb = 0.0;
	tau = 0.0;

	return;
}

PER::PER()
{
	firstAlignmentIndex = 0;
	numAlignment = 0;
	probability = 1.0;

	return;
}

readJunction::readJunction()
{
	junctionIndex = 0;
	next = NULL;

	return;
}

/************************************************************************/
/*             Compute Probabilities                                    */
/************************************************************************/

void PERpath::compute_prob()
{
	//compute the probability that the path is true given the PER alignment and the path
	//probability contains two parts: distance and junction

	//part 1. distance probability
	//has been stored in the member 'distanceProb' after the computation of possible paths

	//part 2. probability of junctions used by the path
	probability = 1.0;
	if (pathInfo != NULL)
	{
		pathJunction *currentJunc;
		currentJunc = pathInfo->junctionsInPath;

		while (currentJunc != NULL)
		{
			if (currentJunc->junctionID > 0)
			{
				probability = probability * junctionList[currentJunc->junctionID]->probability;
			}

			currentJunc = currentJunc->next;
		}
	}


	probability = probability * distanceProb;

	return;
}


void PERalignment::compute_prob()
{
	//compute the probability that the PER is true given the alignment
	//probability contains two parts: path and read junction

	//part 1. path probability
	double pathProb;
	pathProb = 0.0;
	maxPathProb = 0.0;

	PERpath *currentPERpath;
	currentPERpath = possiblePath;
	while (currentPERpath != NULL)
	{
		if (currentPERpath->probability > maxPathProb)
		{
			maxPathProb = currentPERpath->probability;
		}
		pathProb += currentPERpath->probability * currentPERpath->theta;
		currentPERpath = currentPERpath->next;
	}


	//part 2. probability of read junctions
	double juncProb;
	juncProb = 1.0;

	readJunction *currentJunc;
	currentJunc = startRead->readJunc;
	while (currentJunc != NULL)
	{
		if (currentJunc->junctionIndex > 0)
		{
			juncProb = juncProb * junctionList[currentJunc->junctionIndex]->probability;
		}
		currentJunc = currentJunc->next;
	}

	currentJunc = endRead->readJunc;
	while (currentJunc != NULL)
	{
		if (currentJunc->junctionIndex > 0)
		{
			juncProb = juncProb * junctionList[currentJunc->junctionIndex]->probability;
		}
		currentJunc = currentJunc->next;
	}

	probability = pathProb * juncProb;

	return;
}


void PER::compute_prob()
{
	//compute the probability that the PER is true
	//by summing up the probability on every alignment
	double prob;
	prob = 0.0;

	PERalignment *currentAlignment;
	for (long i = 0; i < numAlignment; i++)
	{
		currentAlignment = perAlignList[firstAlignmentIndex + i];
		prob += currentAlignment->probability * currentAlignment->tau;
	}

	probability = prob;

	return;
}





/************************************************************************/
/* Block Splitting Algorithm                                           */
/************************************************************************/



double probability_pathDistance(long distancePer)
{
	double prob;

	prob = pow(distancePer - PERPathDistanceMean, 2) / (2 * PERPathDistanceVAR);
	prob = 1 / exp(prob);
	//prob = pow(2.71828, -prob);
	//	prob = prob / sqrt(2 * 3.1415926 * PERPathDistanceVAR);

	//prob = exp(- pow(distancePer - PERPathDistanceMean, 2) / (2 * PERPathDistanceVAR)) / sqrt(2 * 3.1415926 * PERPathDistanceVAR);

	//	prob = prob / PERPathDistanceProbScalar;

	return prob;
}


bool block_contains_position(block *checkedBlk, long position)
{
	//check if the position is within the range of the block indexed by blkIndex
	//if the position is to the left of the block, return -1
	//if the position is to the right of the block, return 1
	//if the position is within the range of the block, return 0

	// 	if (position < checkedBlk->start_position)
	// 	{
	// 		//position is to the left of the block
	// 		return -1;
	// 	}
	// 	else if (position >= checkedBlk->end_position)
	// 	{
	// 		//position is to the right of the block
	// 		return 1;
	// 	}
	// 	else
	// 	{
	// 		//block contains the position
	// 		return 0;
	// 	}

	if (position >= checkedBlk->start_position && position <= checkedBlk->end_position)
	{
		return true;
	} 
	else
	{
		return false;
	}
}

block* block_search(long position, char *chromosome, block *startBlock)
{
	//search the block that contains the position on the chromosome
	//return the index of that block
	//return NULL if not found

	block *ptr, *endBlock;

	if (strcmp(startBlock->chromosome, chromosome) == 0)
	{
		ptr = startBlock;
	}
	else
	{
		if (strcmp(chromosome, targetChr_start) == 0)
		{
			ptr = chrBlocks;
		} 
		else
		{
			ptr = chrBlocks_2ndStart;
		}
	}

	if (strcmp(chromosome, targetChr_start) == 0)
	{
		endBlock = chrBlocks_2ndStart;
	} 
	else
	{
		endBlock = NULL;
	}

	while (ptr != endBlock)
	{
		if (block_contains_position(ptr, position) == true && strcmp(ptr->chromosome, chromosome) == 0)
		{
			return ptr;
		}

		ptr = ptr->next;
	}

	return ptr; //not found
}


void block_split(long junctionInsert)
{
	//split the chromosome blocks 

	//find the blocks where the junction stands
	block *blkJuncStart, *blkJuncEnd; //the index of the block which contains the starting/ending point of the junction
	blockGraphEdge *newBlkGraphEdge, *currentMoveEdge, *prevMoveEdge;
	linkedGraphEdge *newLinkedEdge, *currentMoveLinkedEdge, *prevMoveLinkedEdge;

	junction *insertJunction;
	insertJunction = junctionList[junctionInsert];
	long junctionStart, junctionEnd;
	junctionStart = insertJunction->start_end;
	junctionEnd = insertJunction->end_start;

	/*	if (junctionStart > 56800460 && junctionStart < 56800520 && junctionEnd > 1048844800 && junctionEnd < 1048845200)
	{
	cout << "hei" << endl;
	}*/

	blkJuncStart = block_search(junctionStart, insertJunction->chromosome_start, lastJunctionSearchBlock);
	if (blkJuncStart != NULL)
	{
		lastJunctionSearchBlock = blkJuncStart; //since the junctions have been sorted based on their starting position
	} 
	else
	{
		//cannot find the junction position
		cout << "ERROR: junction position out of chromosome range" << endl;
		exit(1);
	}

	blkJuncEnd = block_search(junctionEnd, insertJunction->chromosome_end, lastJunctionSearchBlock);
	if (blkJuncEnd == NULL)
	{
		//cannot find the junction position
		cout << "ERROR: junction position out of chromosome range" << endl;
		exit(1);
	}


	//split the block node
	if (blkJuncStart == blkJuncEnd)
	{
		//start&end on the same block
		block *newBlk1, *newBlk2;
		newBlk1 = new block;
		newBlk2 = new block;

		strcpy(newBlk1->chromosome, blkJuncStart->chromosome);
		strcpy(newBlk2->chromosome, blkJuncStart->chromosome);


		//update block range
		newBlk2->end_position = blkJuncStart->end_position;
		blkJuncStart->end_position = junctionStart;
		newBlk1->start_position = junctionStart + 1;
		newBlk1->end_position = junctionEnd - 1;
		newBlk2->start_position = junctionEnd;



		//update block linkage (logical and physical edges)
		//double direction links for logical edges
		newBlk2->next = blkJuncStart->next;

		//modify edges
		//Previous error: newBlk2->edge = blkJuncStart->edge; 
		prevMoveEdge = NULL;
		currentMoveEdge = blkJuncStart->edge;
		while (currentMoveEdge != NULL)
		{
			if (currentMoveEdge->leavingPosition == blkJuncStart->start_position)
			{
				//do not move
				prevMoveEdge = currentMoveEdge;
				currentMoveEdge = currentMoveEdge->nextEdge;
			} 
			else if (currentMoveEdge->leavingPosition == newBlk2->end_position)
			{
				//move to newBlk2
				if (prevMoveEdge == NULL)
				{
					blkJuncStart->edge = currentMoveEdge->nextEdge;
				} 
				else
				{
					prevMoveEdge->nextEdge = currentMoveEdge->nextEdge;
				}

				currentMoveEdge->nextEdge = newBlk2->edge;
				newBlk2->edge = currentMoveEdge;

				if (prevMoveEdge == NULL)
				{
					currentMoveEdge = blkJuncStart->edge;
				} 
				else
				{
					currentMoveEdge = prevMoveEdge->nextEdge;
				}
			}
			else
			{
				cout << "error!" << endl;
			}
		}

		//modify linked edges
		prevMoveLinkedEdge = NULL;
		currentMoveLinkedEdge = blkJuncStart->linkEdge;
		while (currentMoveLinkedEdge != NULL)
		{
			if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == blkJuncStart->start_position)
			{
				//do not move
				prevMoveLinkedEdge = currentMoveLinkedEdge;
				currentMoveLinkedEdge = currentMoveLinkedEdge->next;
			} 
			else if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == newBlk2->end_position)
			{
				//modify the linked edge
				currentMoveLinkedEdge->linkedEdge->linkedBlk = newBlk2;

				//move to newBlk2
				if (prevMoveLinkedEdge == NULL)
				{
					blkJuncStart->linkEdge = currentMoveLinkedEdge->next;
				} 
				else
				{
					prevMoveLinkedEdge->next = currentMoveLinkedEdge->next;
				}

				currentMoveLinkedEdge->next = newBlk2->linkEdge;
				newBlk2->linkEdge = currentMoveLinkedEdge;

				if (prevMoveLinkedEdge == NULL)
				{
					currentMoveLinkedEdge = blkJuncStart->linkEdge;
				} 
				else
				{
					currentMoveLinkedEdge = prevMoveLinkedEdge->next;
				}
			}
			else
			{
				cout << "error!" << endl;
			}
		}


		blkJuncStart->next = newBlk1;
		newBlk1->next = newBlk2;

		newBlkGraphEdge = new blockGraphEdge;		
		newBlkGraphEdge->junctionID = -1;
		newBlkGraphEdge->nextEdge = blkJuncStart->edge;
		newBlkGraphEdge->linkedBlk = newBlk1;
		newBlkGraphEdge->arrivedPosition = newBlk1->start_position;
		newBlkGraphEdge->leavingPosition = blkJuncStart->end_position;
		//		newBlkGraphEdge->skippedLength = 0;
		blkJuncStart->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = newBlk1->linkEdge;
		newBlk1->linkEdge = newLinkedEdge;

		newBlkGraphEdge = new blockGraphEdge;		
		newBlkGraphEdge->junctionID = -1;
		newBlkGraphEdge->nextEdge = newBlk1->edge;
		newBlkGraphEdge->linkedBlk = newBlk2;
		newBlkGraphEdge->arrivedPosition = newBlk2->start_position;
		newBlkGraphEdge->leavingPosition = newBlk1->end_position;
		//		newBlkGraphEdge->skippedLength = 0;
		newBlk1->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = newBlk2->linkEdge;
		newBlk2->linkEdge = newLinkedEdge;

		newBlkGraphEdge = new blockGraphEdge;
		newBlkGraphEdge->junctionID = junctionInsert;
		newBlkGraphEdge->nextEdge = blkJuncStart->edge;
		newBlkGraphEdge->linkedBlk = newBlk2;
		newBlkGraphEdge->arrivedPosition = newBlk2->start_position;
		newBlkGraphEdge->leavingPosition = blkJuncStart->end_position;
		//		newBlkGraphEdge->skippedLength = junctionEnd - junctionStart - 1;
		blkJuncStart->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = newBlk2->linkEdge;
		newBlk2->linkEdge = newLinkedEdge;


		newBlkGraphEdge = new blockGraphEdge;		
		newBlkGraphEdge->junctionID = -2;
		newBlkGraphEdge->nextEdge = newBlk1->edge;
		newBlkGraphEdge->linkedBlk = blkJuncStart;
		newBlkGraphEdge->arrivedPosition = blkJuncStart->end_position;
		newBlkGraphEdge->leavingPosition = newBlk1->start_position;
		//		newBlkGraphEdge->skippedLength = 0;
		newBlk1->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = blkJuncStart->linkEdge;
		blkJuncStart->linkEdge = newLinkedEdge;

		newBlkGraphEdge = new blockGraphEdge;		
		newBlkGraphEdge->junctionID = -2;
		newBlkGraphEdge->nextEdge = newBlk2->edge;
		newBlkGraphEdge->linkedBlk = newBlk1;
		newBlkGraphEdge->arrivedPosition = newBlk1->end_position;
		newBlkGraphEdge->leavingPosition = newBlk2->start_position;
		//		newBlkGraphEdge->skippedLength = 0;
		newBlk2->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = newBlk1->linkEdge;
		newBlk1->linkEdge = newLinkedEdge;

		newBlkGraphEdge = new blockGraphEdge;
		newBlkGraphEdge->junctionID = junctionInsert;
		newBlkGraphEdge->nextEdge = newBlk2->edge;
		newBlkGraphEdge->linkedBlk = blkJuncStart;
		newBlkGraphEdge->arrivedPosition = blkJuncStart->end_position;
		newBlkGraphEdge->leavingPosition = newBlk2->start_position;
		//		newBlkGraphEdge->skippedLength = -(junctionEnd - junctionStart - 1);
		newBlk2->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = blkJuncStart->linkEdge;
		blkJuncStart->linkEdge = newLinkedEdge;


		if (blkJuncStart == chrBlocks_tail)
		{
			chrBlocks_tail = newBlk2;
		}
		chrBlockNum = chrBlockNum + 2;
	} 
	else
	{
		//start&end on the different block
		block *newBlk;
		block *newEndBlk;

		//split start block		
		if (junctionStart < blkJuncStart->end_position)
		{
			newBlk = new block;

			//update block range
			newBlk->end_position = blkJuncStart->end_position;
			blkJuncStart->end_position = junctionStart;
			newBlk->start_position = junctionStart + 1;

			strcpy(newBlk->chromosome, blkJuncStart->chromosome);

			//update block linkage (logical and physical edges)
			newBlk->next = blkJuncStart->next;


			//modify edges
			//Previous error: newBlk->edge = blkJuncStart->edge;
			prevMoveEdge = NULL;
			currentMoveEdge = blkJuncStart->edge;
			while (currentMoveEdge != NULL)
			{
				if (currentMoveEdge->leavingPosition == blkJuncStart->start_position)
				{
					//do not move
					prevMoveEdge = currentMoveEdge;
					currentMoveEdge = currentMoveEdge->nextEdge;
				} 
				else if (currentMoveEdge->leavingPosition == newBlk->end_position)
				{
					//move to newBlk
					if (prevMoveEdge == NULL)
					{
						blkJuncStart->edge = currentMoveEdge->nextEdge;
					} 
					else
					{
						prevMoveEdge->nextEdge = currentMoveEdge->nextEdge;
					}

					currentMoveEdge->nextEdge = newBlk->edge;
					newBlk->edge = currentMoveEdge;

					if (prevMoveEdge == NULL)
					{
						currentMoveEdge = blkJuncStart->edge;
					} 
					else
					{
						currentMoveEdge = prevMoveEdge->nextEdge;
					}
				}
				else
				{
					cout << "error!" << endl;
				}
			}

			//modify linked edges
			prevMoveLinkedEdge = NULL;
			currentMoveLinkedEdge = blkJuncStart->linkEdge;
			while (currentMoveLinkedEdge != NULL)
			{
				if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == blkJuncStart->start_position)
				{
					//do not move
					prevMoveLinkedEdge = currentMoveLinkedEdge;
					currentMoveLinkedEdge = currentMoveLinkedEdge->next;
				} 
				else if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == newBlk->end_position)
				{
					//modify the linked edge
					currentMoveLinkedEdge->linkedEdge->linkedBlk = newBlk;

					//move to newBlk
					if (prevMoveLinkedEdge == NULL)
					{
						blkJuncStart->linkEdge = currentMoveLinkedEdge->next;
					} 
					else
					{
						prevMoveLinkedEdge->next = currentMoveLinkedEdge->next;
					}

					currentMoveLinkedEdge->next = newBlk->linkEdge;
					newBlk->linkEdge = currentMoveLinkedEdge;

					if (prevMoveLinkedEdge == NULL)
					{
						currentMoveLinkedEdge = blkJuncStart->linkEdge;
					} 
					else
					{
						currentMoveLinkedEdge = prevMoveLinkedEdge->next;
					}
				}
				else
				{
					cout << "error!" << endl;
				}
			}


			blkJuncStart->next = newBlk;

			newBlkGraphEdge = new blockGraphEdge;
			newBlkGraphEdge->junctionID = -1;
			newBlkGraphEdge->nextEdge = blkJuncStart->edge;
			newBlkGraphEdge->linkedBlk = newBlk;
			newBlkGraphEdge->arrivedPosition = newBlk->start_position;
			newBlkGraphEdge->leavingPosition = blkJuncStart->end_position;
			//			newBlkGraphEdge->skippedLength = 0;
			blkJuncStart->edge = newBlkGraphEdge;

			newLinkedEdge = new linkedGraphEdge;
			newLinkedEdge->linkedEdge = newBlkGraphEdge;
			newLinkedEdge->next = newBlk->linkEdge;
			newBlk->linkEdge = newLinkedEdge;

			newBlkGraphEdge = new blockGraphEdge;
			newBlkGraphEdge->junctionID = -2;
			newBlkGraphEdge->nextEdge = newBlk->edge;
			newBlkGraphEdge->linkedBlk = blkJuncStart;
			newBlkGraphEdge->arrivedPosition = blkJuncStart->end_position;
			newBlkGraphEdge->leavingPosition = newBlk->start_position;
			//			newBlkGraphEdge->skippedLength = 0;
			newBlk->edge = newBlkGraphEdge;

			newLinkedEdge = new linkedGraphEdge;
			newLinkedEdge->linkedEdge = newBlkGraphEdge;
			newLinkedEdge->next = blkJuncStart->linkEdge;
			blkJuncStart->linkEdge = newLinkedEdge;

			chrBlockNum = chrBlockNum++;
		} 

		//split end block
		newEndBlk = blkJuncEnd;
		if (junctionEnd > blkJuncEnd->start_position)
		{
			newBlk = new block;
			newEndBlk = newBlk;

			//update block range
			newBlk->end_position = blkJuncEnd->end_position;
			blkJuncEnd->end_position = junctionEnd - 1;
			newBlk->start_position = junctionEnd;

			strcpy(newBlk->chromosome, blkJuncEnd->chromosome);

			//update block linkage (logical and physical edges)
			newBlk->next = blkJuncEnd->next;


			//modify edges
			//Previous error: newBlk->edge = blkJuncEnd->edge; 
			prevMoveEdge = NULL;
			currentMoveEdge = blkJuncEnd->edge;
			while (currentMoveEdge != NULL)
			{
				if (currentMoveEdge->leavingPosition == blkJuncEnd->start_position)
				{
					//do not move
					prevMoveEdge = currentMoveEdge;
					currentMoveEdge = currentMoveEdge->nextEdge;
				} 
				else if (currentMoveEdge->leavingPosition == newBlk->end_position)
				{
					//move to newBlk
					if (prevMoveEdge == NULL)
					{
						blkJuncEnd->edge = currentMoveEdge->nextEdge;
					} 
					else
					{
						prevMoveEdge->nextEdge = currentMoveEdge->nextEdge;
					}

					currentMoveEdge->nextEdge = newBlk->edge;
					newBlk->edge = currentMoveEdge;

					if (prevMoveEdge == NULL)
					{
						currentMoveEdge = blkJuncEnd->edge;
					} 
					else
					{
						currentMoveEdge = prevMoveEdge->nextEdge;
					}
				}
				else
				{
					cout << "error!" << endl;
				}
			}

			//modify linked edges
			prevMoveLinkedEdge = NULL;
			currentMoveLinkedEdge = blkJuncEnd->linkEdge;
			while (currentMoveLinkedEdge != NULL)
			{
				if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == blkJuncEnd->start_position)
				{
					//do not move
					prevMoveLinkedEdge = currentMoveLinkedEdge;
					currentMoveLinkedEdge = currentMoveLinkedEdge->next;
				} 
				else if (currentMoveLinkedEdge->linkedEdge->arrivedPosition == newBlk->end_position)
				{
					//modify the linked edge
					currentMoveLinkedEdge->linkedEdge->linkedBlk = newBlk;

					//move to newBlk
					if (prevMoveLinkedEdge == NULL)
					{
						blkJuncEnd->linkEdge = currentMoveLinkedEdge->next;
					} 
					else
					{
						prevMoveLinkedEdge->next = currentMoveLinkedEdge->next;
					}

					currentMoveLinkedEdge->next = newBlk->linkEdge;
					newBlk->linkEdge = currentMoveLinkedEdge;

					if (prevMoveLinkedEdge == NULL)
					{
						currentMoveLinkedEdge = blkJuncEnd->linkEdge;
					} 
					else
					{
						currentMoveLinkedEdge = prevMoveLinkedEdge->next;
					}
				}
				else
				{
					cout << "error!" << endl;
				}
			}


			blkJuncEnd->next = newBlk;

			newBlkGraphEdge = new blockGraphEdge;
			newBlkGraphEdge->junctionID = -1;
			newBlkGraphEdge->nextEdge = blkJuncEnd->edge;
			newBlkGraphEdge->linkedBlk = newBlk;
			newBlkGraphEdge->arrivedPosition = newBlk->start_position;
			newBlkGraphEdge->leavingPosition = blkJuncEnd->end_position;
			//			newBlkGraphEdge->skippedLength = 0;
			blkJuncEnd->edge = newBlkGraphEdge;

			newLinkedEdge = new linkedGraphEdge;
			newLinkedEdge->linkedEdge = newBlkGraphEdge;
			newLinkedEdge->next = newBlk->linkEdge;
			newBlk->linkEdge = newLinkedEdge;

			newBlkGraphEdge = new blockGraphEdge;
			newBlkGraphEdge->junctionID = -2;
			newBlkGraphEdge->nextEdge = newBlk->edge;
			newBlkGraphEdge->linkedBlk = blkJuncEnd;
			newBlkGraphEdge->arrivedPosition = blkJuncEnd->end_position;
			newBlkGraphEdge->leavingPosition = newBlk->start_position;
			//			newBlkGraphEdge->skippedLength = 0;
			newBlk->edge = newBlkGraphEdge;

			newLinkedEdge = new linkedGraphEdge;
			newLinkedEdge->linkedEdge = newBlkGraphEdge;
			newLinkedEdge->next = blkJuncEnd->linkEdge;
			blkJuncEnd->linkEdge = newLinkedEdge;

			chrBlockNum = chrBlockNum++;
		}

		newBlkGraphEdge = new blockGraphEdge;
		newBlkGraphEdge->junctionID = junctionInsert;
		newBlkGraphEdge->nextEdge = blkJuncStart->edge;
		newBlkGraphEdge->linkedBlk = newEndBlk;
		newBlkGraphEdge->arrivedPosition = newEndBlk->start_position;
		newBlkGraphEdge->leavingPosition = blkJuncStart->end_position;
		//		newBlkGraphEdge->skippedLength = junctionEnd - junctionStart - 1;
		blkJuncStart->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = newEndBlk->linkEdge;
		newEndBlk->linkEdge = newLinkedEdge;

		newBlkGraphEdge = new blockGraphEdge;
		newBlkGraphEdge->junctionID = junctionInsert;
		newBlkGraphEdge->nextEdge = newEndBlk->edge;
		newBlkGraphEdge->linkedBlk = blkJuncStart;
		newBlkGraphEdge->arrivedPosition = blkJuncStart->end_position;
		newBlkGraphEdge->leavingPosition = newEndBlk->start_position;
		//		newBlkGraphEdge->skippedLength = -(junctionEnd - junctionStart - 1);
		newEndBlk->edge = newBlkGraphEdge;

		newLinkedEdge = new linkedGraphEdge;
		newLinkedEdge->linkedEdge = newBlkGraphEdge;
		newLinkedEdge->next = blkJuncStart->linkEdge;
		blkJuncStart->linkEdge = newLinkedEdge;


		if (blkJuncEnd == chrBlocks_tail)
		{
			chrBlocks_tail = newEndBlk;
		}		
	}


	return;
}


void pathQueue_initialization()
{
	pathQueueHead = 0;
	pathQueueTail = 0;
	pathQueueMaxLength = MAX_PATH_NUM;

	return;
}

bool pathQueue_enqueue(path *x)
{
	pathQueue[pathQueueTail] = x;

	if (pathQueueTail == pathQueueMaxLength)
	{
		pathQueueTail = 0;
	}
	else 
	{
		pathQueueTail++;
	}

	return true;
}

path* pathQueue_dequeue()
{
	path *x;

	if (pathQueueHead < pathQueueTail)
	{
		x = pathQueue[pathQueueHead];

		if (pathQueueHead == pathQueueMaxLength)
		{
			pathQueueHead = 0;
		}
		else 
		{
			pathQueueHead++;
		}

		return x;
	} 
	else
	{
		return NULL;
	}

}

bool checkRepeatPathJunc(path *checkPath, long checkJuncID)
{
	//check if the junction has been in the path
	//true if the junction is existing
	if (checkJuncID < 1)
	{
		//sequence
		return false;
	}

	pathJunction *currentJunc;
	currentJunc = checkPath->junctionsInPath;
	while (currentJunc != NULL)
	{
		if (currentJunc->junctionID == checkJuncID)
		{
			return true;
		}
		currentJunc = currentJunc->next;
	}

	return false;
}


path* computeBlkPath(block *startBlk, block *endBlk)
{
	//compute possible paths from startBlk to endBlk
	//return the path list

	/*	if (startBlk->start_position <= 56800516 && startBlk->end_position >= 56800551 && endBlk->start_position <= 1048844982 && endBlk->end_position >= 1048845017)
	{
	cout << "ha" << endl;
	}*/

	path *currentPath, *resultPath, *newPath;
	blockGraphEdge *currentEdge;
	pathJunction *newPathJunction;
	char tmpStr_cat[100];
	long newLength;
	int newChangeOfDir = 0;
	bool changedDirection; //new direction after adding the new edge

	resultPath = NULL;

	pathQueue_initialization();
	//build initial queue
	currentEdge = startBlk->edge;
	while (currentEdge != NULL)
	{
		//if (currentEdge->linkedBlk->start_position <= endBlk->end_position && currentEdge->linkedBlk->supported == true)
		//ALL ARE POSSIBLE
		{
			newPath = new path;
			newPath->pathLength = 0;
			newPath->destBlk = currentEdge->linkedBlk;
			newPath->next = NULL;
			newPath->endposition = currentEdge->arrivedPosition;
			newPath->startposition = currentEdge->leavingPosition;
			newPath->lastJunctionID = currentEdge->junctionID;

			newPath->lastUnwritten_M = 0;

			if (currentEdge->junctionID >= 0)
			{
				//through junction
				newPathJunction = new pathJunction;
				newPathJunction->junctionID = currentEdge->junctionID;
				newPath->junctionsInPath = newPathJunction;
				//newPath->skippedLength = currentEdge->skippedLength;

				//sprintf(newPath->info, "%ldN", junctionList[currentEdge->junctionID]->end_start - junctionList[currentEdge->junctionID]->start_end - 1);
//				sprintf(newPath->info, "%ldN(%ld,%ld)", newPath->endposition - newPath->startposition, newPath->endposition, newPath->startposition);

				/************************************************************************/
				/* Change on July 10th                                                  */
				/************************************************************************/
			
				sprintf(newPath->info, "%ldN", newPath->endposition - newPath->startposition + (newPath->endposition > newPath->startposition ? (-1) : 1));
			} 
			else
			{
				//through no junction
				newPath->junctionsInPath = NULL;
				//newPath->skippedLength = 0;

				//sprintf(newPath->info, "0N");
				sprintf(newPath->info, "");
			}

			if (newPath->endposition >= newPath->startposition)
			{
				newPath->currentDirection = true;
			} 
			else
			{
				newPath->currentDirection = false;
			}
			newPath->changeOfDirection = 0; //generate path.

			pathQueue_enqueue(newPath);
		}

		currentEdge = currentEdge->nextEdge;
	}

	currentPath = pathQueue_dequeue();
	while (currentPath != NULL)
	{
		if (currentPath->pathLength <= MAX_PATH_LENGTH)
		{
			if (currentPath->destBlk == endBlk)
			{
				//currentPath achieves the destination

				if (currentPath->lastUnwritten_M != 0)
				{
					sprintf(tmpStr_cat, "%ldM", currentPath->lastUnwritten_M);
					strcat(currentPath->info, tmpStr_cat);
				}

				if (resultPath == NULL)
				{
					resultPath = currentPath;
				} 
				else
				{
					currentPath->next = resultPath;
					resultPath = currentPath;
				}
			} 
			else
			{
				currentEdge = currentPath->destBlk->edge;
				while (currentEdge != NULL)
				{
					//eliminate loops
					if ((currentEdge->junctionID > 0 && currentEdge->junctionID == currentPath->lastJunctionID) || (currentEdge->junctionID == -1 && currentPath->lastJunctionID == -2) || (currentEdge->junctionID == -2 && currentPath->lastJunctionID == -1))
					{
						currentEdge = currentEdge->nextEdge;
						continue;
					}


					//if (currentEdge->linkedBlk->start_position <= endBlk->end_position && currentEdge->linkedBlk->supported == true)
					//ALL ARE POSSIBLE

					//updated path length if adding the edge
					newLength = currentPath->pathLength + abs(currentEdge->leavingPosition - currentPath->endposition) + 1;

					//updated number of changes on direction if adding the edge
					if (currentPath->currentDirection == true)
					{
						//from low to high
						if (currentEdge->arrivedPosition - currentEdge->leavingPosition >= 0)
						{
							//same direction
							newChangeOfDir = currentPath->changeOfDirection;
							changedDirection = true;
						} 
						else
						{
							//opposite direction
							newChangeOfDir = currentPath->changeOfDirection + 1;
							changedDirection = false;
						}
					} 
					else
					{
						//from high to low
						if (currentEdge->arrivedPosition - currentEdge->leavingPosition >= 0)
						{
							//opposite direction
							newChangeOfDir = currentPath->changeOfDirection + 1;
							changedDirection = true;
						} 
						else
						{
							//same direction
							newChangeOfDir = currentPath->changeOfDirection;
							changedDirection = false;
						}
					}
					if (currentPath == NULL)
					{	cout << "error" << endl;
					exit(1);
					}

					//cout << newLength << '\t' << MAX_PATH_LENGTH << '\t' << newChangeOfDir << "\t" << checkRepeatPathJunc(currentPath, currentEdge->junctionID) << endl;
					if (newLength <= MAX_PATH_LENGTH && newChangeOfDir <= MAX_PATH_CHANGE_DIRECTION && checkRepeatPathJunc(currentPath, currentEdge->junctionID) == false)
					{
						newPath = new path;
						//cout << "newPath\t";
						//compute the path length
						//newPath->pathLength = currentPath->pathLength + currentPath->destBlk->end_position - currentPath->destBlk->start_position + 1;
						newPath->pathLength = newLength;


						newPath->endposition = currentEdge->arrivedPosition;
						newPath->startposition = currentPath->startposition;
						newPath->destBlk = currentEdge->linkedBlk;
						newPath->lastJunctionID = currentEdge->junctionID;
						newPath->next = NULL;
						//cout << "anchor\t";
						strcpy(newPath->info, currentPath->info);
						
						//newPath->lastUnwritten_M = currentPath->lastUnwritten_M + currentEdge->leavingPosition + 1 - currentPath->endposition;
						newPath->lastUnwritten_M = currentPath->lastUnwritten_M + currentEdge->leavingPosition - currentPath->endposition
							+ (currentEdge->leavingPosition >= currentPath->endposition? 1 : -1);

						if (currentEdge->junctionID >= 0)
						{
							//through junction
							newPathJunction = new pathJunction;
							newPathJunction->junctionID = currentEdge->junctionID;
							newPathJunction->next = currentPath->junctionsInPath;
							newPath->junctionsInPath = newPathJunction;
							//newPath->skippedLength = currentPath->skippedLength + junctionList[currentEdge->junctionID]->end_start - junctionList[currentEdge->junctionID]->start_end - 1;

							//sprintf(tmpStr_cat, "%ldM%ldN", newPath->lastUnwritten_M, junctionList[currentEdge->junctionID]->end_start - junctionList[currentEdge->junctionID]->start_end - 1);
//							sprintf(tmpStr_cat, "%ldM%ldN", newPath->lastUnwritten_M, currentEdge->arrivedPosition - currentEdge->leavingPosition);

							/************************************************************************/
							/* Change on July 10th                                                  */
							/************************************************************************/

							sprintf(tmpStr_cat, "%ldM%ldN", newPath->lastUnwritten_M, currentEdge->arrivedPosition - currentEdge->leavingPosition + (currentEdge->arrivedPosition > currentEdge->leavingPosition ? (-1) : 1));




							strcat(newPath->info, tmpStr_cat);
							newPath->lastUnwritten_M = 0;
						} 
						else
						{
							//through no junction
							newPath->junctionsInPath = currentPath->junctionsInPath;
							//newPath->skippedLength = currentPath->skippedLength;
						}

						newPath->changeOfDirection = newChangeOfDir;
						newPath->currentDirection = changedDirection;

						pathQueue_enqueue(newPath);
					}

					currentEdge = currentEdge->nextEdge;
				}
				//				delete currentPath;
			}
		}

		currentPath = pathQueue_dequeue();
	}


	return resultPath;
}



// long partition_PERsort(long p,long r)
// {
// 	long x, i, j, tmp;
// 
// 	//randomized partition
// 	i = p + (double)rand()/ (RAND_MAX) * (r - p);
// 
// 	if (sortKey_PERalign[r] != sortKey_PERalign[i])
// 	{
// 		tmp = sortKey_PERalign[r];
// 		sortKey_PERalign[r] = sortKey_PERalign[i];
// 		sortKey_PERalign[i] = tmp;
// 
// 		tmp = orderedPERalignlist[r];
// 		orderedPERalignlist[r] = orderedPERalignlist[i];
// 		orderedPERalignlist[i] = tmp;
// 	}
// 
// 	x = sortKey_PERalign[r];
// 	i = p - 1;
// 
// 	for(j = p; j < r; j++)
// 	{
// 		if (sortKey_PERalign[j] <= x)
// 		{
// 			i++;
// 
// 			if (sortKey_PERalign[i] != sortKey_PERalign[j])
// 			{
// 				tmp = sortKey_PERalign[j];
// 				sortKey_PERalign[j] = sortKey_PERalign[i];
// 				sortKey_PERalign[i] = tmp;
// 
// 				tmp = orderedPERalignlist[j];
// 				orderedPERalignlist[j] = orderedPERalignlist[i];
// 				orderedPERalignlist[i] = tmp;
// 			}
// 		}
// 	}
// 
// 	if (sortKey_PERalign[r] != sortKey_PERalign[i+1])
// 	{
// 		tmp = sortKey_PERalign[r];
// 		sortKey_PERalign[r] = sortKey_PERalign[i+1];
// 		sortKey_PERalign[i+1]=tmp;
// 
// 		tmp = orderedPERalignlist[r];
// 		orderedPERalignlist[r] = orderedPERalignlist[i+1];
// 		orderedPERalignlist[i+1] = tmp;
// 	}
// 
// 	return i+1;
// }
// 
// 
// void quicksort_PERsort(long length)
// {
// 	long top = 0, p, r, q;
// 
// 	stack[top++] = 1;
// 	stack[top++] = length;
// 
// 	while (top != 0)
// 	{
// 		r = stack[--top];
// 		p = stack[--top];
// 
// 		if(p>=r)
// 			continue;
// 
// 		q = partition_PERsort(p, r);
// 
// 		stack[top++] = p;
// 		stack[top++] = q - 1;
// 
// 		stack[top++] = q + 1;
// 		stack[top++] = r;
// 	}
// 
// 	return;
// }


void merge_PERalignSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_PERalign[p + i - 1];
		mergeSort_LorderedList[i] = orderedPERalignlist[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_PERalign[q + j];
		mergeSort_RorderedList[j] = orderedPERalignlist[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_PERalign[k] = mergeSort_Larray[i];
			orderedPERalignlist[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey_PERalign[k] = mergeSort_Rarray[j];
			orderedPERalignlist[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_PERalignSort()
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = perAlignNum;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_PERalignSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}



void compute_PER_distance()
{
	//compute all PER's possible distances

	//sort PERs based on position
	for (long i = 1; i <= perAlignNum; i++)
	{
		orderedPERalignlist[i] = i;
	}

	for (long i = 1; i <= perAlignNum; i++)
	{
		sortKey_PERalign[i] = perAlignList[orderedPERalignlist[i]]->endRead->crucialPoint;
	}
	//quicksort_PERsort(perAlignNum);
	mergeSort_PERalignSort();

	for (long i = 1; i <= perAlignNum; i++)
	{
		sortKey_PERalign[i] = perAlignList[orderedPERalignlist[i]]->startRead->crucialPoint;
	}
	//quicksort_PERsort(perAlignNum);
	mergeSort_PERalignSort();


	//compute possible paths and its distanceProb for every PER
	PERalignment *currentPER;
	pathCluster *lastPathCluster; //memorize last time starting position
	lastPathCluster = pathClusterList;
	block *lastStartBlk; //memorize last time starting block
	lastStartBlk = NULL;
	pathCluster *currentPathCluster;
	path *currentPath, *newPath;
	double currentProbability;
	PERpath *newPERpath, *tmpPERpath;
	long pathCount, pathNumEqual0 = 0, pathNumEqual1 = 0, pathNumLarger2Cnt = 0, maxPathCount = 0, totalPathCount = 0, currentLength = 0;
	bool validpath;
	char newPathInfo[250], tmpInfo[50];
	//double temp_prob; //for normalizing probabilities

	for (long i = 1; i <= perAlignNum; i++)
	{
		currentPER = perAlignList[orderedPERalignlist[i]];
		pathCount = 0;
		//temp_prob = 0.0;

		/*		if (orderedPERalignlist[i] == 8257)
		{
		cout << "8257" << endl;
		}*/

		//find the corresponding path cluster
		currentPathCluster = lastPathCluster;
		//currentPathCluster = pathClusterList;
		while (currentPathCluster != NULL && (currentPathCluster->startBlk != currentPER->startRead->onBlock || currentPathCluster->endBlk != currentPER->endRead->onBlock))
		{
			currentPathCluster = currentPathCluster->next;
		}

		if (currentPathCluster == NULL)
		{
			//path cluster not generated
			currentPathCluster = new pathCluster;
			currentPathCluster->startBlk = currentPER->startRead->onBlock;
			currentPathCluster->endBlk = currentPER->endRead->onBlock;
			currentPathCluster->next = NULL;
			currentPathCluster->clusterPaths = computeBlkPath(currentPathCluster->startBlk, currentPathCluster->endBlk);

			if (pathClusterList == NULL)
			{
				pathClusterList = currentPathCluster;
				pathClusterListTail = currentPathCluster;
			} 
			else
			{
				pathClusterListTail->next = currentPathCluster;
				pathClusterListTail = currentPathCluster;
			}
		}

		//update last time memory
		if (lastPathCluster == NULL || lastStartBlk != currentPathCluster->startBlk)
		{
			lastPathCluster = currentPathCluster;
			lastStartBlk = currentPathCluster->startBlk;
		}


		//compute possible paths for current PER
		currentPath = currentPathCluster->clusterPaths;
		while (currentPath != NULL)
		{
			validpath = true;
			currentLength = currentPath->pathLength - 4;

			if (currentPER->endRead->crucialStrand == false && currentPath->endposition <= currentPER->endRead->crucialPoint)
			{
				currentLength += currentPER->endRead->crucialPoint - currentPath->endposition;
//				sprintf(newPathInfo, "%ldM%s", currentPER->endRead->crucialPoint - currentPath->endposition, currentPath->info);
			}
			else if (currentPER->endRead->crucialStrand == true && currentPath->endposition >= currentPER->endRead->crucialPoint)
			{
				currentLength += currentPath->endposition - currentPER->endRead->crucialPoint;
//				sprintf(newPathInfo, "-%ldM%s", currentPath->endposition - currentPER->endRead->crucialPoint, currentPath->info);
			} 
			else
			{
				validpath = false;
			}

			if (validpath == true)
			{
				if (currentPER->startRead->crucialStrand == true && currentPath->startposition >= currentPER->startRead->crucialPoint)
				{
					currentLength += currentPath->startposition - currentPER->startRead->crucialPoint;
//					sprintf(tmpInfo, "%ldM", currentPath->startposition - currentPER->startRead->crucialPoint);
//					strcat(newPathInfo, tmpInfo);
				}
				else if (currentPER->startRead->crucialStrand == false && currentPath->startposition <= currentPER->startRead->crucialPoint)
				{
					currentLength += currentPER->startRead->crucialPoint - currentPath->startposition;
//					sprintf(tmpInfo, "-%ldM", currentPER->startRead->crucialPoint - currentPath->startposition);
//					strcat(newPathInfo, tmpInfo);
				} 
				else
				{
					validpath = false;
				}
			}

			//constraint on the changes of directions
			if (currentPER->startRead->crucialStrand == true && currentPER->endRead->crucialStrand == false && currentPER->startRead->crucialPoint <= currentPER->endRead->crucialPoint)
			{
				if (currentPath->changeOfDirection > 0)
				{
					validpath = false;
				}
			}
			else if (currentPER->startRead->crucialStrand == false && currentPER->endRead->crucialStrand == true && currentPER->startRead->crucialPoint >= currentPER->endRead->crucialPoint)
			{
				if (currentPath->changeOfDirection > 0)
				{
					validpath = false;
				}
			}
			else if (currentPER->startRead->crucialStrand == true && currentPER->endRead->crucialStrand == true)
			{
				if (currentPath->changeOfDirection > 1)
				{
					validpath = false;
				}
			}
			else if (currentPER->startRead->crucialStrand == false && currentPER->endRead->crucialStrand == false)
			{
				if (currentPath->changeOfDirection > 1)
				{
					validpath = false;
				}
			}

			if (validpath == true)
			{
				currentProbability = probability_pathDistance(currentLength);
				if (currentProbability > MIN_PATH_PROBABILITY)
				{
					newPERpath = new PERpath;
					newPERpath->pathInfo = currentPath;
					newPERpath->distance = currentLength;
					newPERpath->distanceProb = currentProbability;
//					strcpy(newPERpath->infoString, newPathInfo);

					if (currentPER->possiblePath == NULL)
					{
						currentPER->possiblePath = newPERpath;
					} 
					else
					{
						newPERpath->next = currentPER->possiblePath;
						currentPER->possiblePath = newPERpath;
					}

					//temp_prob += currentProbability;

					pathCount++;
				}
			}

			currentPath = currentPath->next;
		}	
		// 		//normalize probabilities
		// 		tmpPERpath = currentPER->possiblePath;
		// 		while (tmpPERpath != NULL)
		// 		{
		// 			tmpPERpath->distanceProb = tmpPERpath->distanceProb / temp_prob;
		// 			tmpPERpath = tmpPERpath->next;
		// 		}

		if (currentPER->possiblePath == NULL && currentPER->startRead->onBlock == currentPER->endRead->onBlock)
		{
			//no possible paths, start and end on the same block

			if (currentPER->startRead->crucialStrand == true && currentPER->endRead->crucialStrand == false && currentPER->startRead->crucialPoint <= currentPER->endRead->crucialPoint)
			{
				currentLength = currentPER->endRead->crucialPoint - currentPER->startRead->crucialPoint - 4 - 1;
			} 
			else if (currentPER->startRead->crucialStrand == false && currentPER->endRead->crucialStrand == true && currentPER->endRead->crucialPoint <= currentPER->startRead->crucialPoint)
			{
				currentLength = currentPER->startRead->crucialPoint - currentPER->endRead->crucialPoint - 4 - 1;
			} 
			else
			{
				currentLength = -1;
			}

			if (currentLength > 0 && probability_pathDistance(currentLength) > MIN_PATH_PROBABILITY)
			{
				// 			newPath = new path;
				// 			newPath->destBlk = currentPER->startRead->onBlock;
				// 			newPath->junctionsInPath = NULL;
				// 			newPath->next = NULL;
				// 			newPath->pathLength = 0;
				// 			newPath->skippedLength = 0;

				newPERpath = new PERpath;
				//			newPERpath->pathInfo = newPath;
				newPERpath->pathInfo = NULL;

				newPERpath->distance = currentLength;
				newPERpath->distanceProb = probability_pathDistance(currentLength); //no path

				currentPER->possiblePath = newPERpath;
				pathCount++;
			}
		}

		currentPER->pathNum = pathCount;

		if (pathCount == 0)
		{
			pathNumEqual0++;
		}
		else if (pathCount == 1)
		{
			pathNumEqual1++;
		}
		else if (pathCount > 1)
		{
			pathNumLarger2Cnt++;
			if (pathCount > maxPathCount)
			{
				maxPathCount = pathCount;
			}
		}

		totalPathCount += pathCount;
	}
	outputPathSummary << targetChr << '\t' << pathNumLarger2Cnt << '\t' << maxPathCount << endl;
	outputSummary << totalPathCount << '\t' << maxPathCount << '\t' << pathNumEqual0 << '\t' << pathNumEqual1 << '\t' << pathNumLarger2Cnt << '\t';

	return;
}


// long partition_JunctionSort(long p,long r)
// {
// 	long x, i, j, tmp;
// 
// 	//randomized partition
// 	i = p + (double)rand()/ (RAND_MAX) * (r - p);
// 
// 	if (sortKey_Junction[r] != sortKey_Junction[i])
// 	{
// 		tmp = sortKey_Junction[r];
// 		sortKey_Junction[r] = sortKey_Junction[i];
// 		sortKey_Junction[i] = tmp;
// 
// 		tmp = orderedJunctionList[r];
// 		orderedJunctionList[r] = orderedJunctionList[i];
// 		orderedJunctionList[i] = tmp;
// 	}
// 
// 	x = sortKey_Junction[r];
// 	i = p - 1;
// 
// 	for(j = p; j < r; j++)
// 	{
// 		if (sortKey_Junction[j] <= x)
// 		{
// 			i++;
// 
// 			if (sortKey_Junction[i] != sortKey_Junction[j])
// 			{
// 				tmp = sortKey_Junction[j];
// 				sortKey_Junction[j] = sortKey_Junction[i];
// 				sortKey_Junction[i] = tmp;
// 
// 				tmp = orderedJunctionList[j];
// 				orderedJunctionList[j] = orderedJunctionList[i];
// 				orderedJunctionList[i] = tmp;
// 			}
// 		}
// 	}
// 
// 	if (sortKey_Junction[r] != sortKey_Junction[i+1])
// 	{
// 		tmp = sortKey_Junction[r];
// 		sortKey_Junction[r] = sortKey_Junction[i+1];
// 		sortKey_Junction[i+1]=tmp;
// 
// 		tmp = orderedJunctionList[r];
// 		orderedJunctionList[r] = orderedJunctionList[i+1];
// 		orderedJunctionList[i+1] = tmp;
// 	}
// 
// 	return i+1;
// }
// 
// 
// void quicksort_JunctionSort(long length)
// {
// 	long top = 0, p, r, q;
// 
// 	stack[top++] = 1;
// 	stack[top++] = length;
// 
// 	while (top != 0)
// 	{
// 		r = stack[--top];
// 		p = stack[--top];
// 
// 		if ( p >= r)
// 			continue;
// 
// 		q = partition_JunctionSort(p, r);
// 
// 		stack[top++] = p;
// 		stack[top++] = q - 1;
// 
// 		stack[top++] = q + 1;
// 		stack[top++] = r;
// 	}
// 
// 	return;
// }


void merge_JunctionSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_Junction[p + i - 1];
		mergeSort_LorderedList[i] = orderedJunctionList[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_Junction[q + j];
		mergeSort_RorderedList[j] = orderedJunctionList[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_Junction[k] = mergeSort_Larray[i];
			orderedJunctionList[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey_Junction[k] = mergeSort_Rarray[j];
			orderedJunctionList[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_JunctionSort()
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = junctionNum;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_JunctionSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}


void input_junction(char *junctionFilename)
{
	ifstream junctionFile;
	junctionFile.open(junctionFilename);

	junction *newJunction;
	char name[30], strand[5];
	bool junctionExist;
	long temp;
	int support;
	for (int i = 0; i < 30; i++)
	{
		name[i] = '\0';
	}

	junctionFile >> name;
	while (name[1] != '\0')
	{
		newJunction = new junction;

		junctionFile >> newJunction->chromosome_start;
		junctionFile >> newJunction->start_end;
		junctionFile >> newJunction->chromosome_end;
		junctionFile >> newJunction->end_start;
		junctionFile >> newJunction->exchange;
		junctionFile >> strand;
		junctionFile >> newJunction->support;

		if (strand[0] == '+')
		{
			newJunction->strand_start = true;
		} 
		else
		{
			newJunction->strand_start = false;
		}		

		if (strand[1] == '+')
		{
			newJunction->strand_end = true;
		} 
		else
		{
			newJunction->strand_end = false;
		}

		if (strcmp(targetChr_start, newJunction->chromosome_start) != 0)
		{
			newJunction->start_end += MAX_CHR_LENGTH;
		}

		if (strcmp(targetChr_start, newJunction->chromosome_end) != 0)
		{
			newJunction->end_start += MAX_CHR_LENGTH;
		}

		newJunction->single = false;

		if (strcmp(newJunction->chromosome_start, newJunction->chromosome_end) == 0 && abs(newJunction->end_start - newJunction->start_end) < 50000)
		{
			newJunction->fusion = false;
		}
		else
		{
			newJunction->fusion = true;
		}

		junctionExist = false;
		for (temp = 1; temp <= junctionNum; temp++)
		{
			if (strcmp(junctionList[temp]->chromosome_start, newJunction->chromosome_start) == 0 && strcmp(junctionList[temp]->chromosome_end, newJunction->chromosome_end) == 0 && junctionList[temp]->start_end == newJunction->start_end && junctionList[temp]->end_start == newJunction->end_start)
			{
				//				(junctionList[temp]->support)++;
				junctionExist = true;
				break;
			}
		}

		if (junctionExist == false)// && newJunction->support > 1)
		{
			junctionList[++junctionNum] = newJunction;
		}
		else
		{
			delete newJunction;
		}


		name[1] = '\0';
		junctionFile >> name;	
	}

	junctionFile.close();

	return;
}


void input_PERsSingle(char* PERfileName)
{
	//derive junctions from single reads file
	ifstream inputfile;
	inputfile.open(PERfileName);

	if (inputfile.is_open() == false)
	{
		return;
	}

	char name[100], chromosome[30], end[50], strand[5], prevChromosome[30], juncStartChr[30], juncEndChr[30];
	long start;
	long tmp, i, chrPositionBase;
	junction *newJunction;
	bool junctionExist, isFusion;
	long tempStart, tempEnd, temp, juncStart, juncEnd, endpoint;


	for (tmp = 0; tmp < 100; tmp++)
	{
		name[tmp] = '\0';
	}


	inputfile >> name;

	while (name[1] != '\0')
	{
		inputfile >> tmp;
		if (tmp == 1)
		{
			isFusion = true;
		}
		else
		{
			isFusion = false;
		}

		if (isFusion == false)
		{
			inputfile >> chromosome;
			inputfile >> strand;
			inputfile >> start;
			inputfile >> end;	

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;

			tmp = 0;
			for (i = 0; end[i] != '\0'; i++)
			{
				if (end[i] == 'M')
				{
					if (tmp == 0)
					{
						continue;
					}

					start = start + tmp;
					tmp = 0;
				} 
				else if (end[i] == 'N')
				{
					//derive a junction
					if (tmp > 50)
					{
						junctionExist = false;
						tempStart = start - 1;
						tempEnd = start + tmp;

						for (temp = 1; temp <= junctionNum; temp++)
						{
							if (strcmp(junctionList[temp]->chromosome_start, chromosome) == 0 && strcmp(junctionList[temp]->chromosome_end, chromosome) == 0 && junctionList[temp]->start_end == tempStart && junctionList[temp]->end_start == tempEnd)
							{
								(junctionList[temp]->support)++;
								junctionExist = true;
								break;
							}
						}

						if (junctionExist == false)
						{
							newJunction = new junction;
							strcpy(newJunction->chromosome_start, chromosome);
							strcpy(newJunction->chromosome_end, chromosome);
							newJunction->start_end = tempStart;
							newJunction->end_start = tempEnd;
							newJunction->single = true;
							newJunction->fusion = false;

							junctionList[++junctionNum] = newJunction;
						}
					}

					start = start + tmp;
					tmp = 0;
				}
				else
				{
					tmp = tmp * 10 + end[i] - 48;
				}
			}
		} 
		else
		{
			//seg 1
			inputfile >> prevChromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(prevChromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;


			//seg 2
			inputfile >> chromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;

			//input strand
			inputfile >> strand;


			/************************************************************************/
			/* input fusion junction                                                */
			/************************************************************************/

			inputfile >> tempStart;
			inputfile >> tempEnd;

			if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) == 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart;
					juncEnd = tempEnd;
				} 
				else
				{
					juncStart = tempEnd;
					juncEnd = tempStart;
				}
			}
			else if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) != 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				juncStart = tempStart;
				juncEnd = tempEnd + MAX_CHR_LENGTH;
			}
			else if (strcmp(prevChromosome, targetChr_start) != 0 && strcmp(chromosome,targetChr_start) == 0)
			{
				strcpy(juncStartChr, chromosome);
				strcpy(juncEndChr, prevChromosome);
				juncStart = tempEnd;
				juncEnd = tempStart + MAX_CHR_LENGTH;
			} 
			else
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart + MAX_CHR_LENGTH;
					juncEnd = tempEnd + MAX_CHR_LENGTH;
				} 
				else
				{
					juncStart = tempEnd + MAX_CHR_LENGTH;
					juncEnd = tempStart + MAX_CHR_LENGTH;
				}
			}

			junctionExist = false;
			for (temp = 1; temp <= junctionNum; temp++)
			{
				if (strcmp(junctionList[temp]->chromosome_start, juncStartChr) == 0 && strcmp(junctionList[temp]->chromosome_end, juncEndChr) == 0 && junctionList[temp]->start_end == juncStart && junctionList[temp]->end_start == juncEnd)
				{
					(junctionList[temp]->support)++;
					junctionExist = true;
					break;
				}
			}

			if (junctionExist == false)
			{
				newJunction = new junction;
				strcpy(newJunction->chromosome_start, juncStartChr);
				strcpy(newJunction->chromosome_end, juncEndChr);
				newJunction->start_end = juncStart;
				newJunction->end_start = juncEnd;
				newJunction->single = false;
				newJunction->fusion = true;

				junctionList[++junctionNum] = newJunction;
			}

		}

		statNumSingleRead++;

		name[1] = '\0';

		inputfile >> name;		
	}

	inputfile.close();

	return;
}

void input_PERs(char* PERfileName, char* SingleFilename1, char* SingleFilename2, char* SingleFilename3)
{
	//input PE reads from input file
	ifstream inputfile;
	inputfile.open(PERfileName);

	if (inputfile.is_open() == false)
	{
		cout << "no file     ";
		exit(1);
	}

	char name[100], chromosome[30], end[50], prevName[100], prevChromosome[30], strand[5], juncStartChr[30], juncEndChr[30];
	long start, PE_list_index, minPosition, maxPosition, endpoint, chrPositionBase;
	PERalignment* newPERalign;
	PER *newPER;
	readJunction *newReadJunc;
	reads* newRead;
	range* newRange;
	long tmp, i;
	junction *newJunction;
	bool junctionExist, isFusion;
	long tempStart, tempEnd, temp, juncStart, juncEnd;
	long numPERsingleAlign, numPERmultiAlign, maxNumPERalign;

	minPosition = MAX_CHR_LENGTH * 2;
	maxPosition = 0;

	PE_list_index = 0;
	numPERsingleAlign = 0;
	numPERmultiAlign = 0;
	maxNumPERalign = 0;

	for (tmp = 0; tmp < 100; tmp++)
	{
		name[tmp] = '\0';
	}

	sprintf(prevName, "  ");
	inputfile >> name;

	while (name[1] != '\0')
	{
		newPERalign = new PERalignment;

		/************************************************************************/
		/* input start read                                                     */
		/************************************************************************/

		newRead = new reads;
		newRead->startRead = true;

		inputfile >> tmp;
		if (tmp == 1)
		{
			isFusion = true;
		} 
		else
		{
			isFusion = false;
		}

		if (isFusion == false)
		{
			inputfile >> chromosome;
			inputfile >> strand;
			inputfile >> start;
			inputfile >> end;

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;

			if (strand[0] == '+')
			{
				newRead->strandLow = true;
				newRead->strandHigh = true;
			} 
			else
			{
				newRead->strandLow = false;
				newRead->strandHigh = false;
			}

			//strcpy(newPERalign->name, name);
			strcpy(newRead->chromosome, chromosome);
			strcpy(newRead->end, end);

			statNumPairRead++;

			tmp = 0;
			for (i = 0; end[i] != '\0'; i++)
			{
				if (end[i] == 'M')
				{
					if (tmp == 0)
					{
						continue;
					}

					// 					newRange = new range;
					// 					newRange->start = start;
					// 					newRange->end = start + tmp - 1;
					// 
					// 					if (newRead->rangeTail == NULL)
					// 					{
					// 						newRead->readsRange = newRange;
					// 						newRead->rangeTail = newRange;
					// 					} 
					// 					else
					// 					{
					// 						newRead->rangeTail->next = newRange;
					// 						newRead->rangeTail = newRange;
					// 					}

					if (start < newRead->rangeLow)
					{
						newRead->rangeLow = start;
					}

					if (start + tmp > newRead->rangeHigh)
					{
						newRead->rangeHigh = start + tmp - 1;
					}

					//					(newRead->rangeCount)++;

					start = start + tmp;
					tmp = 0;
				} 
				else if (end[i] == 'N')
				{
					//derive a junction
					if (tmp > 50)
					{
						junctionExist = false;
						tempStart = start - 1;
						tempEnd = start + tmp;

						for (temp = 1; temp <= junctionNum; temp++)
						{
							if (strcmp(junctionList[temp]->chromosome_start, chromosome) == 0 && strcmp(junctionList[temp]->chromosome_end, chromosome) == 0 && junctionList[temp]->start_end == tempStart && junctionList[temp]->end_start == tempEnd)
							{
								(junctionList[temp]->support)++;
								junctionExist = true;
								break;
							}
						}

						if (junctionExist == false)
						{
							newJunction = new junction;
							strcpy(newJunction->chromosome_start, chromosome);
							strcpy(newJunction->chromosome_end, chromosome);
							newJunction->start_end = tempStart;
							newJunction->end_start = tempEnd;
							newJunction->single = false;
							newJunction->fusion = false;

							junctionList[++junctionNum] = newJunction;
						}

						newReadJunc = new readJunction;
						newReadJunc->junctionIndex = temp;
						newReadJunc->next = newRead->readJunc;
						newRead->readJunc = newReadJunc;

						//block_split(junctionNum-1); //insert a new junction, split the chromosome blocks
					}

					start = start + tmp;
					tmp = 0;
				}
				else
				{
					tmp = tmp * 10 + end[i] - 48;
				}
			}

			if (newRead->strandHigh == true)
			{
				newRead->crucialStrand = true;
				newRead->crucialPoint = newRead->rangeHigh - 2; //modify the boundary for junction error
			} 
			else
			{
				newRead->crucialStrand = false;
				newRead->crucialPoint = newRead->rangeLow + 2;
			}
		} 
		else
		{
			/************************************************************************/
			/* input 2 segments of the read                                         */
			/************************************************************************/
			//seg 1
			inputfile >> prevChromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(prevChromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;
			sprintf(newRead->end, "%ldM", endpoint - start + 1);
			newRead->rangeLow = start;

			// 			newRange = new range;
			// 			newRange->start = start;
			// 			newRange->end = endpoint;
			// 			newRead->readsRange = newRange;
			// 			newRead->rangeTail = newRange;
			// 			(newRead->rangeCount)++;

			//seg 2
			inputfile >> chromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;
			newRead->rangeHigh = endpoint;

			// 			newRange = new range;
			// 			newRange->start = start;
			// 			newRange->end = endpoint;
			// 			newRead->rangeTail->next = newRange;
			// 			newRead->rangeTail = newRange;

			//input strand
			inputfile >> strand;

			if (strand[0] == '+')
			{
				newRead->strandLow = true;
			} 
			else
			{
				newRead->strandLow = false;
			}
			if (strand[1] == '+')
			{
				newRead->strandHigh = true;
			} 
			else
			{
				newRead->strandHigh = false;
			}

			//strcpy(newPERalign->name, name);
			sprintf(newRead->chromosome, "%s_%s", prevChromosome, chromosome);
			sprintf(newRead->end, "%s%ldM", newRead->end, endpoint - start + 1);

			statNumPairRead++;


			if (newRead->strandHigh == true)
			{
				newRead->crucialStrand = true;
				newRead->crucialPoint = endpoint - 2;
			} 
			else
			{
				newRead->crucialStrand = false;
				newRead->crucialPoint = start + 2;
			}


			/************************************************************************/
			/* input fusion junction                                                */
			/************************************************************************/

			inputfile >> tempStart;
			inputfile >> tempEnd;

			if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) == 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart;
					juncEnd = tempEnd;
				} 
				else
				{
					juncStart = tempEnd;
					juncEnd = tempStart;
				}
			}
			else if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) != 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				juncStart = tempStart;
				juncEnd = tempEnd + MAX_CHR_LENGTH;
			}
			else if (strcmp(prevChromosome, targetChr_start) != 0 && strcmp(chromosome,targetChr_start) == 0)
			{
				strcpy(juncStartChr, chromosome);
				strcpy(juncEndChr, prevChromosome);
				juncStart = tempEnd;
				juncEnd = tempStart + MAX_CHR_LENGTH;
			} 
			else
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart + MAX_CHR_LENGTH;
					juncEnd = tempEnd + MAX_CHR_LENGTH;
				} 
				else
				{
					juncStart = tempEnd + MAX_CHR_LENGTH;
					juncEnd = tempStart + MAX_CHR_LENGTH;
				}
			}

			junctionExist = false;
			for (temp = 1; temp <= junctionNum; temp++)
			{
				if (strcmp(junctionList[temp]->chromosome_start, juncStartChr) == 0 && strcmp(junctionList[temp]->chromosome_end, juncEndChr) == 0 && junctionList[temp]->start_end == juncStart && junctionList[temp]->end_start == juncEnd)
				{
					(junctionList[temp]->support)++;
					junctionExist = true;
					break;
				}
			}

			if (junctionExist == false)
			{
				newJunction = new junction;
				strcpy(newJunction->chromosome_start, juncStartChr);
				strcpy(newJunction->chromosome_end, juncEndChr);
				newJunction->start_end = juncStart;
				newJunction->end_start = juncEnd;
				newJunction->single = false;
				newJunction->fusion = true;

				junctionList[++junctionNum] = newJunction;
			}

			newReadJunc = new readJunction;
			newReadJunc->junctionIndex = temp;
			newReadJunc->next = newRead->readJunc;
			newRead->readJunc = newReadJunc;
		}

		newPERalign->startRead = newRead;

		readNum++;
		readList[readNum] = newRead;

		if (newRead->rangeLow < minPosition)
		{
			minPosition = newRead->rangeLow;
		}
		if (newRead->rangeHigh > maxPosition)
		{
			maxPosition = newRead->rangeHigh;
		}



		/************************************************************************/
		/* input end read                                                       */
		/************************************************************************/

		newRead = new reads;
		newRead->startRead = false;

		inputfile >> tmp;
		if (tmp == 1)
		{
			isFusion = true;
		} 
		else
		{
			isFusion = false;
		}

		if (isFusion == false)
		{
			inputfile >> chromosome;
			inputfile >> strand;
			inputfile >> start;
			inputfile >> end;

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;

			if (strand[0] == '+')
			{
				newRead->strandLow = true;
				newRead->strandHigh = true;
			} 
			else
			{
				newRead->strandLow = false;
				newRead->strandHigh = false;
			}


			strcpy(newRead->chromosome, chromosome);
			strcpy(newRead->end, end);

			statNumPairRead++;

			tmp = 0;
			for (i = 0; end[i] != '\0'; i++)
			{
				if (end[i] == 'M')
				{
					if (tmp == 0)
					{
						continue;
					}

					newRange = new range;
					newRange->start = start;
					newRange->end = start + tmp - 1;

					// 					if (newRead->rangeTail == NULL)
					// 					{
					// 						newRead->readsRange = newRange;
					// 						newRead->rangeTail = newRange;
					// 					} 
					// 					else
					// 					{
					// 						newRead->rangeTail->next = newRange;
					// 						newRead->rangeTail = newRange;
					// 					}
					//					(newRead->rangeCount)++;

					if (start < newRead->rangeLow)
					{
						newRead->rangeLow = start;
					}

					if (start + tmp > newRead->rangeHigh)
					{
						newRead->rangeHigh = start + tmp - 1;
					}



					start = start + tmp;
					tmp = 0;
				} 
				else if (end[i] == 'N')
				{
					//derive a junction
					if (tmp > 50)
					{
						junctionExist = false;
						tempStart = start - 1;
						tempEnd = start + tmp;

						for (temp = 1; temp <= junctionNum; temp++)
						{
							if (strcmp(junctionList[temp]->chromosome_start, chromosome) == 0 && strcmp(junctionList[temp]->chromosome_end, chromosome) == 0 && junctionList[temp]->start_end == tempStart && junctionList[temp]->end_start == tempEnd)
							{
								(junctionList[temp]->support)++;
								junctionExist = true;
								break;
							}
						}

						if (junctionExist == false)
						{
							newJunction = new junction;
							strcpy(newJunction->chromosome_start, chromosome);
							strcpy(newJunction->chromosome_end, chromosome);
							newJunction->start_end = tempStart;
							newJunction->end_start = tempEnd;
							newJunction->single = false;
							newJunction->fusion = false;

							junctionList[++junctionNum] = newJunction;
						}

						newReadJunc = new readJunction;
						newReadJunc->junctionIndex = temp;
						newReadJunc->next = newRead->readJunc;
						newRead->readJunc = newReadJunc;

						//block_split(junctionNum-1); //insert a new junction, split the chromosome blocks
					}

					start = start + tmp;
					tmp = 0;
				}
				else
				{
					tmp = tmp * 10 + end[i] - 48;
				}
			}

			if (newRead->strandLow == true)
			{
				newRead->crucialStrand = true;
				newRead->crucialPoint = newRead->rangeHigh - 2;
			} 
			else
			{
				newRead->crucialStrand = false;
				newRead->crucialPoint = newRead->rangeLow + 2;
			}

		} 
		else
		{
			/************************************************************************/
			/* input 2 segments of the read                                         */
			/************************************************************************/
			//seg 1
			inputfile >> prevChromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(prevChromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;
			tempStart = start;
			tempEnd = endpoint;
			sprintf(newRead->end, "%ldM", endpoint - start + 1);
			newRead->rangeLow = start;

			// 			newRange = new range;
			// 			newRange->start = start;
			// 			newRange->end = endpoint;
			// 			newRead->readsRange = newRange;
			// 			newRead->rangeTail = newRange;
			// 			(newRead->rangeCount)++;

			//seg 2
			inputfile >> chromosome;
			inputfile >> start;
			inputfile >> endpoint;

			if (strcmp(chromosome, targetChr_start) == 0)
			{
				chrPositionBase = 0;
			} 
			else
			{
				chrPositionBase = MAX_CHR_LENGTH;
			}
			start += chrPositionBase;
			endpoint += chrPositionBase;
			newRead->rangeHigh = endpoint;

			// 			newRange = new range;
			// 			newRange->start = start;
			// 			newRange->end = endpoint;
			// 			newRead->rangeTail->next = newRange;
			// 			newRead->rangeTail = newRange;

			//input strand
			inputfile >> strand;

			if (strand[0] == '+')
			{
				newRead->strandLow = true;
			} 
			else
			{
				newRead->strandLow = false;
			}
			if (strand[1] == '+')
			{
				newRead->strandHigh = true;
			} 
			else
			{
				newRead->strandHigh = false;
			}

			//strcpy(newPERalign->name, name);
			sprintf(newRead->chromosome, "%s_%s", prevChromosome, chromosome);
			sprintf(newRead->end, "%s%ldM", newRead->end, endpoint - start + 1);

			statNumPairRead++;


			if (newRead->strandLow == true)
			{
				newRead->crucialStrand = true;
				newRead->crucialPoint = tempEnd - 2;
			} 
			else
			{
				newRead->crucialStrand = false;
				newRead->crucialPoint = tempStart + 2;
			}


			/************************************************************************/
			/* input fusion junction                                                */
			/************************************************************************/

			inputfile >> tempStart;
			inputfile >> tempEnd;

			if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) == 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart;
					juncEnd = tempEnd;
				} 
				else
				{
					juncStart = tempEnd;
					juncEnd = tempStart;
				}
			}
			else if (strcmp(prevChromosome, targetChr_start) == 0 && strcmp(chromosome, targetChr_start) != 0)
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				juncStart = tempStart;
				juncEnd = tempEnd + MAX_CHR_LENGTH;
			}
			else if (strcmp(prevChromosome, targetChr_start) != 0 && strcmp(chromosome,targetChr_start) == 0)
			{
				strcpy(juncStartChr, chromosome);
				strcpy(juncEndChr, prevChromosome);
				juncStart = tempEnd;
				juncEnd = tempStart + MAX_CHR_LENGTH;
			} 
			else
			{
				strcpy(juncStartChr, prevChromosome);
				strcpy(juncEndChr, chromosome);
				if (tempStart <= tempEnd)
				{
					juncStart = tempStart + MAX_CHR_LENGTH;
					juncEnd = tempEnd + MAX_CHR_LENGTH;
				} 
				else
				{
					juncStart = tempEnd + MAX_CHR_LENGTH;
					juncEnd = tempStart + MAX_CHR_LENGTH;
				}
			}

			junctionExist = false;
			for (temp = 1; temp <= junctionNum; temp++)
			{
				if (strcmp(junctionList[temp]->chromosome_start, juncStartChr) == 0 && strcmp(junctionList[temp]->chromosome_end, juncEndChr) == 0 && junctionList[temp]->start_end == juncStart && junctionList[temp]->end_start == juncEnd)
				{
					(junctionList[temp]->support)++;
					junctionExist = true;
					break;
				}
			}

			if (junctionExist == false)
			{
				newJunction = new junction;
				strcpy(newJunction->chromosome_start, juncStartChr);
				strcpy(newJunction->chromosome_end, juncEndChr);
				newJunction->start_end = juncStart;
				newJunction->end_start = juncEnd;
				newJunction->single = false;
				newJunction->fusion = true;

				junctionList[++junctionNum] = newJunction;
			}

			newReadJunc = new readJunction;
			newReadJunc->junctionIndex = temp;
			newReadJunc->next = newRead->readJunc;
			newRead->readJunc = newReadJunc;

		}

		newPERalign->endRead = newRead;

		readNum++;
		readList[readNum] = newRead;

		if (newRead->rangeLow < minPosition)
		{
			minPosition = newRead->rangeLow;
		}
		if (newRead->rangeHigh > maxPosition)
		{
			maxPosition = newRead->rangeHigh;
		}

		PE_list_index++;
		perAlignList[PE_list_index] = newPERalign;


		if (strcmp(name, prevName) == 0)
		{
			//multi alignment
			perList[perNum]->numAlignment = perList[perNum]->numAlignment + 1;
		} 
		else
		{
			//new PER
			if (perNum >= 1)
			{
				if (perList[perNum]->numAlignment == 1)
				{
					numPERsingleAlign++;
				} 
				else
				{
					numPERmultiAlign++;
					if (perList[perNum]->numAlignment > maxNumPERalign)
					{
						maxNumPERalign = perList[perNum]->numAlignment;
					}
				}
			}

			newPER = new PER;
			strcpy(newPER->name, name);
			newPER->firstAlignmentIndex = PE_list_index;
			newPER->numAlignment = 1;

			perNum++;
			perList[perNum] = newPER;

			strcpy(prevName, name);
		}

		newPERalign->perIndex = perNum;

		/************************************************************************/
		/* Modify the boundary of the PER alignment for junction mismatch       */
		/************************************************************************/
		//newPERalign->startRead->rangeHigh = newPERalign->startRead->rangeHigh - 2;
		//newPERalign->endRead->rangeLow = newPERalign->endRead->rangeLow + 2;


		name[1] = '\0';

		inputfile >> name;		
	}

	perAlignNum = PE_list_index;

	inputfile.close();

	input_PERsSingle(SingleFilename1);
	if (CROSSCHROMOSOME == true)
	{
		input_PERsSingle(SingleFilename2);
		input_PERsSingle(SingleFilename3);

	}

	//split chromosome blocks

	//sort junctions based on position
	for (long i = 1; i <= junctionNum; i++)
	{
		orderedJunctionList[i] = i;
	}

	for (long i = 1; i <= junctionNum; i++)
	{
		sortKey_Junction[i] = junctionList[orderedJunctionList[i]]->end_start;
	}
	//quicksort_JunctionSort(junctionNum);
	mergeSort_JunctionSort();

	for (long i = 1; i <= junctionNum; i++)
	{
		sortKey_Junction[i] = junctionList[orderedJunctionList[i]]->start_end;
	}
	//quicksort_JunctionSort(junctionNum);
	mergeSort_JunctionSort();

	runTimeStart = clock();
	//split
	for (long i = 1; i <= junctionNum; i++)
	{
		block_split(orderedJunctionList[i]);
	}

	runTime = ((double)clock() - runTimeStart) / CLOCKS_PER_SEC;


	// 	chrBlocks->start_position = minPosition;
	// 	chrBlocks_tail->end_position = maxPosition;


	if (perList[perNum]->numAlignment == 1)
	{
		numPERsingleAlign++;
	} 
	else
	{
		numPERmultiAlign++;
		if (perList[perNum]->numAlignment > maxNumPERalign)
		{
			maxNumPERalign = perList[perNum]->numAlignment;
		}
	}

	statNumRead = statNumSingleRead + statNumPairRead;
	outputSummary << statNumRead << '\t' << statNumSingleRead << '\t' << statNumPairRead << '\t' << perNum << '\t' << numPERsingleAlign << '\t' << numPERmultiAlign << '\t' << perAlignNum << '\t' << maxNumPERalign << '\t' << junctionNum << '\t'<< runTime << '\t' << chrBlockNum << '\t';

	return;
}


// long partition_readSort(long p,long r)
// {
// 	long x, i, j, tmp;
// 
// 	//randomized partition
// 	i = p + (double)rand()/ (RAND_MAX) * (r - p);
// 
// 	tmp = sortKey_read[r];
// 	sortKey_read[r] = sortKey_read[i];
// 	sortKey_read[i] = tmp;
// 
// 	tmp = orderedReadList[r];
// 	orderedReadList[r] = orderedReadList[i];
// 	orderedReadList[i] = tmp;
// 
// 
// 	x = sortKey_read[r];
// 	i = p - 1;
// 
// 	for(j = p; j < r; j++)
// 	{
// 		if (sortKey_read[j] <= x)
// 		{
// 			i++;
// 
// 			tmp = sortKey_read[j];
// 			sortKey_read[j] = sortKey_read[i];
// 			sortKey_read[i] = tmp;
// 
// 			tmp = orderedReadList[j];
// 			orderedReadList[j] = orderedReadList[i];
// 			orderedReadList[i] = tmp;
// 		}
// 	}
// 
// 	tmp = sortKey_read[r];
// 	sortKey_read[r] = sortKey_read[i+1];
// 	sortKey_read[i+1]=tmp;
// 
// 	tmp = orderedReadList[r];
// 	orderedReadList[r] = orderedReadList[i+1];
// 	orderedReadList[i+1] = tmp;
// 
// 	return i+1;
// }
// 
// 
// void quicksort_readSort(long length)
// {
// 	long top = 0, p, r, q;
// 
// 	stack[top++] = 1;
// 	stack[top++] = length;
// 
// 	while (top != 0)
// 	{
// 		r = stack[--top];
// 		p = stack[--top];
// 
// 		if(p>=r)
// 			continue;
// 
// 		q = partition_readSort(p, r);
// 
// 		stack[top++] = p;
// 		stack[top++] = q - 1;
// 
// 		stack[top++] = q + 1;
// 		stack[top++] = r;
// 	}
// 
// 	return;
// }



void merge_ReadSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_read[p + i - 1];
		mergeSort_LorderedList[i] = orderedReadList[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_read[q + j];
		mergeSort_RorderedList[j] = orderedReadList[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_read[k] = mergeSort_Larray[i];
			orderedReadList[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey_read[k] = mergeSort_Rarray[j];
			orderedReadList[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_ReadSort()
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = readNum;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_ReadSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}




void reads_map_block()
{
	//map reads to blocks
	long i;

	//sort the reads based on position
	for (i = 1; i <= readNum; i++)
	{
		orderedReadList[i] = i;
	}

	for (i = 1; i <= readNum; i++)
	{
		sortKey_read[i] = readList[orderedReadList[i]]->crucialPoint;
		// 		if (readList[orderedReadList[i]]->startRead == true)
		// 		{
		// 			sortKey_read[i] = readList[orderedReadList[i]]->crucialPoint;
		// 		} 
		// 		else
		// 		{
		// 			sortKey_read[i] = readList[orderedReadList[i]]->crucialPoint;
		// 		}
	}
	//quicksort_readSort(readNum);
	mergeSort_ReadSort();


	//map reads
	block *currentBlk;
	currentBlk = chrBlocks;

	for (i = 1; i <= readNum; i++)
	{
		//examine read readList[orderedReadList[i]]
		if (currentBlk == NULL)
		{
			cout << sortKey_read[i];
			exit(1);
		}

		//find the block where the read is
		while (currentBlk->end_position < sortKey_read[i])
		{
			currentBlk = currentBlk->next;
		}

		if (orderedReadList[i] > readNum)
		{
			cout << orderedReadList[i];
			exit(1);
		}

		readList[orderedReadList[i]]->onBlock = currentBlk;
		currentBlk->supported = true;
	}

	return;
}









/**********************************************************
***********************************************************
************************* EM ******************************
***********************************************************
**********************************************************/

/************************************************************************/
/* Compute Probabilities for All PER alignments                         */
/************************************************************************/

void compute_prob_path(long startIndex, long endIndex)
{
	//compute probabilities for all the paths
	PERpath *currentPERpath;

	for (long i = startIndex; i <= endIndex; i++)
	{
		currentPERpath = perAlignList[i]->possiblePath;

		while (currentPERpath != NULL)
		{
			currentPERpath->compute_prob();
			currentPERpath = currentPERpath->next;
		}
	}

	return;
}

void compute_prob_PERalign(long startIndex, long endIndex)
{
	//compute probabilities for all the PER alignments
	for (long i = startIndex; i <= endIndex; i++)
	{
		perAlignList[i]->compute_prob();
	}

	return;
}

void compute_prob_PER()
{
	//compute probabilities for all the PERs
	for (long i = 1; i <= perNum; i++)
	{
		perList[i]->compute_prob();
	}

	return;
}


void compute_prob_junction()
{
	//compute probabilities for all the junctions
	double p_pos, p_path;
	PERpath *currentPERpath;
	pathJunction *pathJunc;
	long i, j;

	for (long i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probability = 1.0;
	}

	for (long i = 1; i <= perAlignNum; i++)
	{
		p_pos = perAlignList[i]->probability * perAlignList[i]->tau;
		currentPERpath = perAlignList[i]->possiblePath;
		while (currentPERpath != NULL)
		{
			p_path = currentPERpath->probability * currentPERpath->theta * p_pos;

			if (p_path > 0)
			{
				if (currentPERpath->pathInfo != NULL)
				{
					pathJunc = currentPERpath->pathInfo->junctionsInPath;
				} 
				else
				{
					pathJunc = NULL;
				}

				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						junctionList[j]->probability *= (1 - p_path);
					}
					pathJunc = pathJunc->next;
				}
			}

			currentPERpath = currentPERpath->next;
		}
	}

	for (long i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probability = 1.0 - junctionList[i]->probability;
	}

	return;
}


/************************************************************************/
/* Compute Probabilities for EM set                                     */
/************************************************************************/

void compute_prob_path_EMset()
{
	//compute probabilities for all the paths of PER alignments being estimated
	PERpath *currentPERpath;

	for (long i = perAlignEMlist_Start; i <= perAlignEMlist_End; i++)
	{
		currentPERpath = perAlignList_EM[i]->possiblePath;

		while (currentPERpath != NULL)
		{
			currentPERpath->compute_prob();
			currentPERpath = currentPERpath->next;
		}
	}

	return;
}

void compute_prob_PERalign_EMset()
{
	//compute probabilities for all the PER alignments being estimated
	for (long i = perAlignEMlist_Start; i <= perAlignEMlist_End; i++)
	{
		perAlignList_EM[i]->compute_prob();
	}

	return;
}


void compute_prob_PER_EMset(long startIndex, long endIndex)
{
	//compute probabilities for all the PERs being estimated
	for (long i = startIndex; i <= endIndex; i++)
	{
		perList[perAlignList[orderedPERalignlist[i]]->perIndex]->compute_prob();
	}

	return;
}

void compute_prob_junction_EMset()
{
	//compute probabilities for all the junctions under estimation
	double p_pos, p_path;
	PERpath *currentPERpath;
	pathJunction *pathJunc;
	long i, j;

	for (i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probability = 1.0;
	}

	for (i = perAlignEMlist_Start; i <= perAlignEMlist_End; i++)
	{
		p_pos = perAlignList_EM[i]->probability * perAlignList_EM[i]->tau;
		currentPERpath = perAlignList_EM[i]->possiblePath;
		while (currentPERpath != NULL)
		{
			p_path = currentPERpath->probability * currentPERpath->theta * p_pos;

			if (p_path > 0)
			{
				if (currentPERpath->pathInfo != NULL)
				{
					pathJunc = currentPERpath->pathInfo->junctionsInPath;
				} 
				else
				{
					pathJunc = NULL;
				}

				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						junctionList[j]->probability *= (1 - p_path);
					}
					pathJunc = pathJunc->next;
				}
			}

			currentPERpath = currentPERpath->next;
		}
	}

	for (i = 1; i <= junctionNum; i++)
	{
		if (junctionList[i]->probability < 1.0)
		{
			junctionList[i]->probability = 1.0 - junctionList[i]->probability;
		}
	}

	return;
}



/************************************************************************/
/* EM                                                                   */
/************************************************************************/



long separate_subset(long startIndex)
{
	//separate all the PERs into subset which are independent with each other
	//input the index of starting PER alignment  
	//output the index of ending PER alignment

	long endBoard; //ending position of the subset
	long i, tmp;
	PERalignment *currentAlign;
	perAlignEMlist_Start = perAlignEMlist_End + 1;

	//the 'startIndex'th alignment is initiated in the closure
	i = startIndex;
	currentAlign = perAlignList[orderedPERalignlist[i]];
	while ((currentAlign->probability < 0.3 || currentAlign->maxPathProb > 0.7) && i < perAlignNum)
	{
		i++;
		currentAlign = perAlignList[orderedPERalignlist[i]];
	}
	if (currentAlign->pathNum > 1)
	{
		//parameter of current PER alignment needs to be estimated 
		perAlignEMlist_End++;
		perAlignList_EM[perAlignEMlist_End] = currentAlign;
	}
	endBoard = currentAlign->endRead->crucialPoint;

	for (i++; i <= perAlignNum; i++)
	{
		currentAlign = perAlignList[orderedPERalignlist[i]];

		if (currentAlign->probability < 0.3 || currentAlign->maxPathProb > 0.7)
		{
			continue;
		}

		if (currentAlign->startRead->crucialPoint > endBoard || perAlignEMlist_End - perAlignEMlist_Start >= MAX_EM_MEMBER)
		{
			//alignment i out of current closure 
			break;
		}

		tmp = currentAlign->endRead->crucialPoint;

		if (tmp > endBoard)
		{
			endBoard = tmp;
		}

		//build EM set		
		if (currentAlign->pathNum > 1)
		{
			//parameter of current PER alignment needs to be estimated 
			perAlignEMlist_End++;
			perAlignList_EM[perAlignEMlist_End] = currentAlign;
		}
	}

	//	outputEMdetail << endl << "---" << endBoard << '\t';
	return (i - 1);
}

double compute_Eval(long startIndex, long endIndex)
{
	//compute evaluation function value for given range of PER alignments (given closure)
	double eval;
	long curPERindex;
	eval = 0.0;

	for (long i = startIndex; i <= endIndex; i++)
	{
		curPERindex = perAlignList[orderedPERalignlist[i]]->perIndex;
		if (perList[curPERindex]->probability > 0.0001)
		{
			eval += perAlignList[orderedPERalignlist[i]]->tau * log10(perList[curPERindex]->probability);
		}
		else
		{
			eval += perAlignList[orderedPERalignlist[i]]->tau * -4;
		}
	}

	return eval;
}


double compute_Eval_ALL()
{
	//compute evaluation function value for all PERs
	double eval;
	eval = 0.0;

	PERalignment *currentAlignment;

	for (long i = 1; i <= perNum; i++)
	{
		for (long j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlignment = perAlignList[perList[i]->firstAlignmentIndex + j];
			if (perList[i]->probability > 0.0001)
			{
				eval += currentAlignment->tau * log10(perList[i]->probability);
			}
			else
			{
				eval += currentAlignment->tau * -4;
			}
		}
	}

	return eval;
}

void probability_backup(long startIndex, long endIndex)
{
	long i;

	PERpath *currentPERpath;

	for (i = startIndex; i <= endIndex; i++)
	{
		currentPERpath = perAlignList[orderedPERalignlist[i]]->possiblePath;

		while (currentPERpath != NULL)
		{
			currentPERpath->probBackup = currentPERpath->probability;
			currentPERpath = currentPERpath->next;
		}
	}

	for (i = startIndex; i <= endIndex; i++)
	{
		perAlignList[orderedPERalignlist[i]]->probBackup = perAlignList[orderedPERalignlist[i]]->probability;
	}

	for (i = startIndex; i <= endIndex; i++)
	{
		perList[perAlignList[orderedPERalignlist[i]]->perIndex]->probBackup = perList[perAlignList[orderedPERalignlist[i]]->perIndex]->probability;
	}

	for (i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probBackup = junctionList[i]->probability;
	}

	return;
}


void probability_recover(long startIndex, long endIndex)
{
	long i;

	PERpath *currentPERpath;

	for (i = startIndex; i <= endIndex; i++)
	{
		currentPERpath = perAlignList[orderedPERalignlist[i]]->possiblePath;

		while (currentPERpath != NULL)
		{
			currentPERpath->probability = currentPERpath->probBackup;
			currentPERpath = currentPERpath->next;
		}
	}

	for (i = startIndex; i <= endIndex; i++)
	{
		perAlignList[orderedPERalignlist[i]]->probability = perAlignList[orderedPERalignlist[i]]->probBackup;
	}

	for (i = startIndex; i <= endIndex; i++)
	{
		perList[perAlignList[orderedPERalignlist[i]]->perIndex]->probability = perList[perAlignList[orderedPERalignlist[i]]->perIndex]->probBackup;
	}

	for (i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probability = junctionList[i]->probBackup;
	}

	return;
}

void hill_climbing(long startIndex, long endIndex)
{
	//hill climbing algorithm for parameter estimation
	//on given closure [startIndex, endIndex]
	double eval, tmpEval, lastIterEval, step, bestProb;
	long i;
	PERpath *currentPERpath, *tmpPERpath;
	long currentPERpathNum, iterCnt, max_iterCnt;

	max_iterCnt = 3;
	iterCnt = 1;
	eval = compute_Eval(startIndex, endIndex);
	lastIterEval = eval * 2;

	//cout << endl << eval << endl;
	/*outputEMdetail << startIndex << " -- " << endIndex << endl;*/


	while (abs(lastIterEval - eval) / abs(lastIterEval) > MIN_UPDATE_RATE * 3 && iterCnt <= max_iterCnt)
	{
		lastIterEval = eval;

		for (i = perAlignEMlist_Start; i <= perAlignEMlist_End; i++)
		{
			currentPERpathNum = perAlignList_EM[i]->pathNum;

			tmpPERpath = perAlignList_EM[i]->possiblePath;
			currentPERpath = NULL;
			bestProb = 0.1;
			while (currentPERpath != NULL)
			{
				if (tmpPERpath->probability > bestProb)
				{
					currentPERpath = tmpPERpath;
				}
			}
			if (currentPERpath != NULL)
			{
				if ((double)rand()/ (RAND_MAX) <= currentPERpath->probability && currentPERpath->theta < 1.0)
				{
					//check a move
					if (currentPERpath->theta + HILLCLIMBING_STEP > 1)
					{
						step = 1 - currentPERpath->theta;
					}
					else
					{
						step = HILLCLIMBING_STEP;
					}

					currentPERpath->theta += step + step/(currentPERpathNum - 1);

					tmpPERpath = perAlignList_EM[i]->possiblePath;
					while (tmpPERpath != NULL)
					{
						tmpPERpath->theta = tmpPERpath->theta - step/(currentPERpathNum - 1);
						tmpPERpath = tmpPERpath->next;
					}

					probability_backup(startIndex, endIndex);

					//reflect the effects of the move
					compute_prob_PERalign_EMset();
					compute_prob_junction_EMset();

					//recalculate the probabilities after the move
					compute_prob_path(startIndex, endIndex);
					compute_prob_PERalign(startIndex, endIndex);
					compute_prob_PER_EMset(startIndex, endIndex);

					//recalculate the evaluation function value
					tmpEval = compute_Eval(startIndex, endIndex);

					if (tmpEval > eval)
					{
						//take the move
						//outputEMdetail << perAlignList_EM[i]->perIndex << '\t' << '(' << perAlignList_EM[i]->startRead->rangeLow << ',' << perAlignList_EM[i]->endRead->rangeHigh << ")\tQn: " << eval << " -> " << tmpEval << endl;

						eval = tmpEval;

						//cout << eval << endl;
					} 
					else
					{
						//undo the move
						currentPERpath->theta = currentPERpath->theta - step - step/(currentPERpathNum - 1);

						tmpPERpath = perAlignList_EM[i]->possiblePath;
						while (tmpPERpath != NULL)
						{
							tmpPERpath->theta += step/(currentPERpathNum - 1);
							tmpPERpath = tmpPERpath->next;
						}

						probability_recover(startIndex, endIndex);
					}
				}

				//currentPERpath = currentPERpath->next;
			}
		}

		iterCnt++;
	}

	return;
}


void maximization()
{
	//maximization step
	long startIndex, endIndex;
	startIndex = 1;
	endIndex = 1;

	while (startIndex <= perAlignNum)
	{
		endIndex = separate_subset(startIndex);

		if (perAlignEMlist_End - perAlignEMlist_Start >= 0L)
		{
			hill_climbing(startIndex, endIndex);
		}

		startIndex = endIndex + 1;
	}

	return;
}

void expectation()
{
	//expectation step
	long i, j;
	double tempSum;
	long startIndex;
	PER *currentPER;

	for (i = 1; i <= perAlignEMlist_End; i++)
	{
		tempSum = 0.0;
		currentPER = perList[perAlignList_EM[i]->perIndex];

		startIndex = currentPER->firstAlignmentIndex;

		for (j = 0; j < currentPER->numAlignment; j++)
		{
			tempSum += perAlignList[startIndex + j]->probability;
		}

		if (tempSum > 0)
		{
			for (j = 0; j < currentPER->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = perAlignList[startIndex + j]->probability / tempSum;
			}
		}
		else
		{
			for (j = 0; j < currentPER->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = 1.0 / currentPER->numAlignment;
			}
		}
	}

	return;
}

double compute_loglikelihood()
{
	//compute the value of log likelihood
	double ln, likelihood;
	//likelihood = 1.0;
	ln = 0.0;

	for (long i = 1; i <= perNum; i++)
	{
		if (perList[i]->probability > 0.0001)
		{
			likelihood = (perList[i]->probability);
			ln += log10(likelihood);
		}
		else
		{
			ln += -4;
		}
	}

	return ln;
}



bool em_convergence(long iterCnt, double currentValue, double prevValue)
{
	//test whether the EM algorithm has converged
	//	return true;

	if (iterCnt < 2)
		return false;

	if (iterCnt > MAX_ITER_COUNT)
	{
		return true;
	}

	double updateRate;
	if (prevValue != 0)
	{
		updateRate = double(abs(currentValue - prevValue)) / abs(prevValue);
		if (updateRate < MIN_UPDATE_RATE)
		{
			return true;
		}
	} 
	else
	{
		if (currentValue - prevValue == 0)
		{
			return true;
		}
	}

	return false;
}




void em_initialization()
{
	//initialize parameters for em
	long i, j;
	double tempSum;
	long startIndex;
	PERpath *currentPERpath;

	//re-order PER alignments with starting position
	for (i = 1; i <= perAlignNum; i++)
	{
		sortKey_PERalign[i] = perAlignList[orderedPERalignlist[i]]->endRead->crucialPoint;
	}
	//quicksort_PERsort(perAlignNum);
	mergeSort_PERalignSort();

	for (i = 1; i <= perAlignNum; i++)
	{
		sortKey_PERalign[i] = perAlignList[orderedPERalignlist[i]]->startRead->crucialPoint;
	}
	//quicksort_PERsort(perAlignNum);
	mergeSort_PERalignSort();


	delete [] sortKey_PERalign;
	delete [] sortKey_read;
	delete [] sortKey_Junction;

	delete [] mergeSort_Larray;
	delete [] mergeSort_Rarray;
	delete [] mergeSort_LorderedList;
	delete [] mergeSort_RorderedList;

	delete [] pathQueue;

	junctionSupportMPPath = new long[MAX_JUNCTION_NUM+1];
	junctionSupportProbMPPath = new double[MAX_JUNCTION_NUM+1];
	junctionSupportRead = new long[MAX_JUNCTION_NUM+1];
	junctionSupportProbRead = new double[MAX_JUNCTION_NUM+1];

	perAlignList_EM = new PERalignment* [MAX_EM_PARA +1];

	//initiate the probabilities
	compute_prob_path(1, perAlignNum);

	//assign initial theta
	for (i = 1; i <= perAlignNum; i++)
	{
		tempSum = 0.0;
		currentPERpath = perAlignList[i]->possiblePath;
		while (currentPERpath != NULL)
		{
			tempSum += currentPERpath->probability;
			currentPERpath = currentPERpath->next;
		}

		if (tempSum > 0)
		{
			currentPERpath = perAlignList[i]->possiblePath;
			while (currentPERpath != NULL)
			{
				currentPERpath->theta = currentPERpath->probability / tempSum;
				currentPERpath = currentPERpath->next;
			}	
		}
		else
		{
			currentPERpath = perAlignList[i]->possiblePath;
			while (currentPERpath != NULL)
			{
				currentPERpath->theta = 1.0 / perAlignList[i]->pathNum;
				currentPERpath = currentPERpath->next;
			}	
		}
	}

	compute_prob_PERalign(1, perAlignNum);

	//assign initial tau
	for (i = 1; i <= perNum; i++)
	{
		tempSum = 0;
		startIndex = perList[i]->firstAlignmentIndex;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			tempSum += perAlignList[startIndex + j]->probability;
		}

		if (tempSum > 0)
		{
			for (j = 0; j < perList[i]->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = perAlignList[startIndex + j]->probability / tempSum;
			}
		}
		else
		{
			for (j = 0; j < perList[i]->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = 1.0 / perList[i]->numAlignment;
			}
		}
	}

	//calculate initial probabilities for junctions
	compute_prob_junction();

	compute_prob_path(1, perAlignNum);

	//reassign theta
	for (i = 1; i <= perAlignNum; i++)
	{
		tempSum = 0.0;
		currentPERpath = perAlignList[i]->possiblePath;
		while (currentPERpath != NULL)
		{
			tempSum += currentPERpath->probability;
			currentPERpath = currentPERpath->next;
		}

		if (tempSum > 0)
		{
			currentPERpath = perAlignList[i]->possiblePath;
			while (currentPERpath != NULL)
			{
				currentPERpath->theta = currentPERpath->probability / tempSum;
				currentPERpath = currentPERpath->next;
			}	
		}	
		else
		{
			currentPERpath = perAlignList[i]->possiblePath;
			while (currentPERpath != NULL)
			{
				currentPERpath->theta = 1.0 / perAlignList[i]->pathNum;
				currentPERpath = currentPERpath->next;
			}	
		}
	}

	compute_prob_PERalign(1, perAlignNum);

	//reassign tau
	for (i = 1; i <= perNum; i++)
	{
		tempSum = 0;
		startIndex = perList[i]->firstAlignmentIndex;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			tempSum += perAlignList[startIndex + j]->probability;
		}

		if (tempSum > 0)
		{
			for (j = 0; j < perList[i]->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = perAlignList[startIndex + j]->probability / tempSum;
			}
		}
		else
		{
			for (j = 0; j < perList[i]->numAlignment; j++)
			{
				perAlignList[startIndex + j]->tau = 1.0 / perList[i]->numAlignment;
			}
		}
	}

	compute_prob_PER();

	return;
}


void EMoutput_junctionProb()
{
	for (long i = 1; i <= junctionNum; i++)
	{
		outputEMprobJunction << junctionList[i]->chromosome_start << '\t' << junctionList[i]->start_end << '\t' << junctionList[i]->end_start << '\t' << junctionList[i]->probability << endl;
	}

	return;
}

void EMoutput_PERprob()
{
	long i, startIndex;
	long j;
	PERalignment *currentAlign;

	for (i = 1; i <= perNum; i++)
	{
		startIndex = perList[i]->firstAlignmentIndex;

		outputEMprobPER << perList[i]->name << '\t' << perList[i]->probability << "\t" ;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			outputEMprobPER << currentAlign->probability << '\t' << currentAlign->tau;
		}		 		

		outputEMprobPER << endl;
	}

	return;
}



void compute_junction_MPprob()
{
	//compute the junction probabilities with hard decision

	double p_pos, p_path, bestProb;
	PERalignment *currentAlign, *bestAlign;
	PERpath *currentPERpath, *bestPERpath;
	pathJunction *pathJunc;
	long i, j, startIndex;

	for (i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probBackup = 1.0;
	}

	for (i = 1; i <= perNum; i++)
	{
		startIndex = perList[i]->firstAlignmentIndex;
		bestAlign = NULL;
		bestProb = -1.0;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			if (currentAlign->tau > bestProb)
			{
				bestProb = currentAlign->tau;
				bestAlign = currentAlign;
			}
		}		

		if (bestAlign == NULL)
		{
			continue;
		}

		currentPERpath = bestAlign->possiblePath;
		p_pos = bestAlign->probability;

		bestProb = -1.0;
		bestPERpath = NULL;

		while (currentPERpath != NULL)
		{
			if (currentPERpath->theta > bestProb)
			{
				bestProb = currentPERpath->theta;
				bestPERpath = currentPERpath;
			}

			currentPERpath = currentPERpath->next;
		}

		if (bestPERpath != NULL)
		{
			p_path = bestPERpath->probability * p_pos;

			if (p_path > 0)
			{
				if (bestPERpath->pathInfo != NULL)
				{
					pathJunc = bestPERpath->pathInfo->junctionsInPath;
				} 
				else
				{
					pathJunc = NULL;
				}

				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						junctionList[j]->probBackup *= (1 - p_path);
					}
					pathJunc = pathJunc->next;
				}
			}
		}		
	}

	for (i = 1; i <= junctionNum; i++)
	{
		junctionList[i]->probBackup = 1.0 - junctionList[i]->probBackup;
	}

	return;
}



void check_junction_HD()
{
	//check junctions according to hard decision
	long i, j, startIndex;
	PERalignment *currentAlign, *bestAlign;
	pathJunction *pathJunc;
	readJunction *currentReadJunc;
	PERpath *currentPERpath, *bestPERpath;
	double bestProb;

	char filename[500];
	//	ofstream outputPathJunc;
	//	sprintf(filename, "%sresult/junction_support/path_junction_%s.txt", dirPrefix, targetChr);
	//	outputPathJunc.open(filename);

	for (i = 1; i <= perNum; i++)
	{
		//		outputPathJunc << perList[i]->name << '\t';

		startIndex = perList[i]->firstAlignmentIndex;
		bestAlign = NULL;
		bestProb = -1.0;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			if (currentAlign->tau > bestProb)
			{
				bestProb = currentAlign->tau;
				bestAlign = currentAlign;
			}
		}		

		if (bestAlign == NULL)
		{
			//			outputPathJunc << endl;
			continue;
		}

		//		outputPathJunc << "[" << bestAlign->startRead->rangeLow << "," << bestAlign->startRead->rangeHigh << "]\t[" << bestAlign->endRead->rangeLow << "," << bestAlign->endRead->rangeHigh << "]\t" << bestAlign->startRead->crucialStrand << bestAlign->endRead->crucialStrand << "\t"; 

		//sequenced segment
		if (bestAlign->probability > MIN_HD_PER_PROB)
		{
			currentReadJunc = bestAlign->startRead->readJunc;
			while (currentReadJunc != NULL)
			{
				j = currentReadJunc->junctionIndex;
				if (j > 0)
				{
					//					outputPathJunc << j << "_,";
					junctionList[j]->isTrue = true;
				}
				currentReadJunc = currentReadJunc->next;
			}
		}


		//un-sequenced segment
		currentPERpath = bestAlign->possiblePath;

		bestProb = -1.0;
		bestPERpath = NULL;

		while (currentPERpath != NULL)
		{
			if (currentPERpath->theta > bestProb)
			{
				bestProb = currentPERpath->theta;
				bestPERpath = currentPERpath;
			}

			currentPERpath = currentPERpath->next;
		}

		//		outputPathJunc << bestAlign->startRead->end;

		if (bestPERpath != NULL)
		{
			if (bestPERpath->pathInfo != NULL)
			{
				//				outputPathJunc << bestAlign->startRead->onBlock->end_position - bestAlign->startRead->rangeHigh << "M" << bestPERpath->pathInfo->info << bestAlign->endRead->rangeLow - bestAlign->endRead->onBlock->start_position << "M" << bestAlign->endRead->end << '\t';
				pathJunc = bestPERpath->pathInfo->junctionsInPath;
			} 
			else
			{
				//				outputPathJunc << bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1 << "M" << bestAlign->endRead->end;
				pathJunc = NULL;
			}

			if (bestPERpath->probability >= MIN_HD_PATH_PROB)
			{
				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						//						outputPathJunc << j << ",";
						junctionList[j]->isTrue = true;
					}
					pathJunc = pathJunc->next;
				}
			}
		}	
		else
		{
			//			outputPathJunc << bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1 << "V" << bestAlign->endRead->end;
		}

		//		outputPathJunc << endl;
	}

	//	outputPathJunc.close();

	return;
}


void expectation_maximization()
{
	//apply EM on all the PERs
	double prevLn, curLn, eval, eval2, test;
	int iterCnt;

	//cout << "EM   ";

	char filename[500];
	ofstream output_em_log;
	//	sprintf(filename, "%sresult/EM_log/EMlog_%s.txt", dirPrefix, targetChr);
	//	output_em_log.open(filename);
	//	output_em_log << "iterCnt\tprevLn\teval\tevalM\tLnM\tLnE" << endl;

	//	sprintf(filename, "%sresult/EM_detail/EMdetail_%s.txt", dirPrefix, targetChr);
	//	outputEMdetail.open(filename);

	//cout << "init   ";
	em_initialization();

	//cout << "computLog    ";

	curLn = compute_loglikelihood();
	prevLn = curLn * 2;
	iterCnt = 0;
	sprintf(filename, "%s_RAW", targetChr);
	//	output_PERpath_prob(filename);
	sprintf(filename, "_RAW");
	//	output_junctionProb(filename);
	//	output_PERpathJunction_RAW();

	outputSummary << curLn << '\t';
	runTimeStart = clock();

	while (em_convergence(iterCnt, curLn, prevLn) == false)
	{
		// 		outputEMdetail << endl << endl << "======================================================================" << "Iteration " << iterCnt << endl;
		// 		sprintf(filename, "C:\\Users\\yin\\Desktop\\PER\\PER_36bp\\result\\Prob_Junction_Iter%d.txt", iterCnt);
		// 		outputEMprobJunction.open(filename);
		// 		sprintf(filename, "C:\\Users\\yin\\Desktop\\PER\\PER_36bp\\result\\Prob_PER_Iter%d.txt", iterCnt);
		// 		outputEMprobPER.open(filename);
		// 		EMoutput_PERprob();
		// 		EMoutput_junctionProb();
		//cout << "iter     " << iterCnt << endl;

		perAlignEMlist_End = 0;
		prevLn = curLn;

		eval2 = compute_Eval_ALL();
		maximization();

		for (long i = 1; i <= perAlignEMlist_End; i++)
		{
			perAlignList_EM[i]->compute_prob();
		}

		eval = compute_Eval_ALL();
		test = compute_loglikelihood();

		expectation();

		compute_prob_PER();

		curLn = compute_loglikelihood();

		//		output_em_log << iterCnt << '\t' << prevLn << '\t' << eval2 << '\t' << eval << '\t' << test << '\t' << curLn << endl;


		// 		outputEMprobPER.close();
		// 		outputEMprobJunction.close();

		iterCnt++;
	}

	// 	output_em_log << iterCnt << '\t' << curLn << endl;
	// 	output_em_log.close();

	/*	outputEMdetail.close();*/

	compute_prob_junction();

	runTime = ((double)clock() - runTimeStart)/ CLOCKS_PER_SEC;
	outputSummary << curLn << '\t' << iterCnt << '\t' << runTime << '\t';

	compute_junction_MPprob();
	check_junction_HD();

	return;
}




/**********************************************************
***********************************************************
********************* OUTPUT ******************************
***********************************************************
**********************************************************/


// void pathInfoProcess(char *info)
// {
// 	//process path information
// 	char newInfo[100], tmpInfo[100];
// 	long curM, curN, tmp, i, j;
// 
// 	sprintf(newInfo, "");
// 	tmp = 0;
// 	curM = 0;
// 	curN = 0;
// 
// 	for (i = 0; info[i] != '\0'; i++)
// 	{
// 		if (info[i] == 'M')
// 		{
// 			if (tmp == 0)
// 			{
// 				//continue;
// 			}
// 			else
// 			{
// 				if (curN == 0)
// 				{
// 					curM += tmp;
// 				}
// 				else
// 				{
// 					sprintf(tmpInfo, "%ldN", curN);
// 					strcat(newInfo, tmpInfo);
// 					curN = 0;
// 					curM = tmp;
// 				}
// 			}
// 
// 			tmp = 0;
// 		}
// 		else if (info[i] == 'N')
// 		{
// 			if (tmp == 0)
// 			{
// 				//continue;
// 			}
// 			else
// 			{
// 				if (curM == 0)
// 				{
// 					curN += tmp;
// 				}
// 				else
// 				{
// 					sprintf(tmpInfo, "%ldM", curM);
// 					strcat(newInfo, tmpInfo);
// 					curM = 0;
// 					curN = tmp;
// 				}
// 			}
// 
// 			tmp = 0;
// 		}
// 		else if (info[i] == 'V')
// 		{
// 			if (tmp == 0)
// 			{
// 				//continue;
// 			}
// 			else
// 			{
// 				if (curN != 0)
// 				{
// 					sprintf(tmpInfo, "%ldN%ldV", curN, tmp);
// 					strcat(newInfo, tmpInfo);
// 					curN = 0;
// 				}
// 
// 				if (curM != 0)
// 				{
// 					sprintf(tmpInfo, "%ldM%ldV", curM, tmp);
// 					strcat(newInfo, tmpInfo);
// 					curM = 0;
// 				}
// 			}
// 
// 			tmp = 0;
// 		}
// 		else
// 		{
// 			tmp = tmp * 10 + info[i] - 48;
// 		}
// 	}
// 
// 	if (curN != 0)
// 	{
// 		sprintf(tmpInfo, "%ldN", curN);
// 		strcat(newInfo, tmpInfo);
// 	}
// 
// 	if (curM != 0)
// 	{
// 		sprintf(tmpInfo, "%ldM", curM);
// 		strcat(newInfo, tmpInfo);
// 	}
// 
// 	strcpy(info, newInfo);
// 
// 	return;
// }



void pathInfoProcess(char *info, bool changeDirection, char *startRead, char *endRead)
{
	//process path information
	//changeDirection: true for negative strand
	char newInfo[200], tmpInfo[200];
	long curM, curN, tmp, i;

	sprintf(newInfo, "");
	tmp = 0;
	curM = 0;
	curN = 0;

	if (changeDirection == true)
	{
		for (i = 0; info[i] != '\0'; i++)
		{
			if (info[i] == 'M')
			{
				if (tmp == 0)
				{
					//continue;
				}
				else
				{
					if (curN == 0)
					{
						curM += tmp;
					}
					else
					{
						strcpy(tmpInfo, newInfo);
						sprintf(newInfo, "%ldN", curN);
						strcat(newInfo, tmpInfo);
						curN = 0;
						curM = tmp;
					}
				}

				tmp = 0;
			}
			else if (info[i] == 'N')
			{
				if (tmp == 0)
				{
					//continue;
				}
				else
				{
					if (curM == 0)
					{
						curN += tmp;
					}
					else
					{
						strcpy(tmpInfo, newInfo);
						sprintf(newInfo, "%ldM", curM);
						strcat(newInfo, tmpInfo);
						curM = 0;
						curN = tmp;
					}
				}

				tmp = 0;
			}
			else if(info[i] == '-')
			{
				//ignore
			}
			else
			{
				tmp = tmp * 10 + info[i] - 48;
			}
		}

		if (curN != 0)
		{
			strcpy(tmpInfo, newInfo);
			sprintf(newInfo, "%ldN", curN);
			strcat(newInfo, tmpInfo);
			curN = 0;
		}

		if (curM != 0)
		{
			strcpy(tmpInfo, newInfo);
			sprintf(newInfo, "%ldM", curM);
			strcat(newInfo, tmpInfo);
			curM = 0;
		}

		strcpy(info, endRead);
		strcat(info, newInfo);
		strcat(info, startRead);
	}
	else
	{
		strcpy(newInfo, info);
		strcpy(info, startRead);
		strcat(info, newInfo);
		strcat(info, endRead);
	}

	sprintf(newInfo, "");

	for (i = 0; info[i] != '\0'; i++)
	{
		if (info[i] == 'M')
		{
			if (tmp == 0)
			{
				//continue;
			}
			else
			{
				if (curN == 0)
				{
					curM += tmp;
				}
				else
				{
					sprintf(tmpInfo, "%ldN", curN);
					strcat(newInfo, tmpInfo);
					curN = 0;
					curM = tmp;
				}
			}

			tmp = 0;
		}
		else if (info[i] == 'N')
		{
			if (tmp == 0)
			{
				//continue;
			}
			else
			{
				if (curM == 0)
				{
					curN += tmp;
				}
				else
				{
					sprintf(tmpInfo, "%ldM", curM);
					strcat(newInfo, tmpInfo);
					curM = 0;
					curN = tmp;
				}
			}

			tmp = 0;
		}
		else if(info[i] == '-')
		{
			sprintf(info, "*");
			return;
		}
		else
		{
			tmp = tmp * 10 + info[i] - 48;
		}
	}

	if (curN != 0)
	{
		sprintf(tmpInfo, "%ldN", curN);
		strcat(newInfo, tmpInfo);
		curN = 0;
	}

	if (curM != 0)
	{
		sprintf(tmpInfo, "%ldM", curM);
		strcat(newInfo, tmpInfo);
		curM = 0;
	}

	strcpy(info, newInfo);

	return;
}

void output_PERpathJunction_RAW()
{
	//output path junctions of all PER paths
	long i, j, startIndex;
	PERalignment *currentAlign;
	pathJunction *pathJunc;
	PERpath *currentPERpath;

	char filename[500];
	ofstream outputPathJunc;
	sprintf(filename, "%sresult/junction_support/path_junction_RAW_%s.txt", dirPrefix, targetChr);
	outputPathJunc.open(filename);

	for (i = 1; i <= perNum; i++)
	{
		outputPathJunc << perList[i]->name << endl;

		startIndex = perList[i]->firstAlignmentIndex;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];

			if (currentAlign == NULL)
			{
				continue;
			}

			outputPathJunc << "[" << currentAlign->startRead->rangeLow << ", " << currentAlign->startRead->rangeHigh << "]\t[" << currentAlign->endRead->rangeLow << ", " << currentAlign->endRead->rangeHigh << "]\t"; 

			currentPERpath = currentAlign->possiblePath;

			while (currentPERpath != NULL)
			{
				if (currentPERpath->pathInfo != NULL)
				{
					outputPathJunc << currentPERpath->pathInfo->info << '\t';
					pathJunc = currentPERpath->pathInfo->junctionsInPath;
				} 
				else
				{
					pathJunc = NULL;
				}

				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						outputPathJunc << j << ",";
					}
					pathJunc = pathJunc->next;
				}

				outputPathJunc << '\t';

				currentPERpath = currentPERpath->next;
			}

			outputPathJunc << endl;
		}		
	}

	outputPathJunc.close();

	return;
}

void output_PERpath_prob()
{
	//output the distanceProb of paths of each PER
	ofstream perPathProb;
	char filename[500];
	// 	sprintf(filename, "C:\\Users\\yin\\Desktop\\PER\\PER_36bp\\result\\path_prob\\PER_pathProb_%s.txt", targetChr);
	// 	perPathProb.open(filename);
	// 
	// 	PERalignment *currentPER;
	// 	PERpath *currentPERpath;
	// 	path *currentPath;
	// 
	// 	for (long i = 1; i <= perAlignNum; i++)
	// 	{
	// 		currentPER = perAlignList[orderedPERalignlist[i]];
	// 		perPathProb << currentPER->name << '\t' << currentPER->startRead->chromosome << '\t' << currentPER->startRead->rangeLow << '\t' << currentPER->startRead->end << '\t' << currentPER->endRead->chromosome << '\t' << currentPER->endRead->rangeLow << '\t' << currentPER->endRead->end << '\t';
	// 
	// 		//output PER paths
	// 		currentPERpath = currentPER->possiblePath;
	// 		while (currentPERpath != NULL)
	// 		{
	// 			//output exact path info using xxMxxN format
	// 			if (currentPERpath->pathInfo == NULL)
	// 			{
	// 				//on the same block
	// 				perPathProb << currentPER->endRead->rangeLow - currentPER->startRead->rangeHigh - 1 << 'M' << '\t';
	// 			} 
	// 			else
	// 			{
	// 				//on different block
	// 				perPathProb << currentPER->startRead->onBlock->end_position - currentPER->startRead->rangeHigh << 'M' << currentPERpath->pathInfo->info << currentPER->endRead->rangeLow - currentPER->endRead->onBlock->start_position << 'M' << '\t';
	// 			}
	// 
	// 
	// 			//output path distance and prob
	// 			perPathProb << currentPERpath->distance << '\t' << currentPERpath->distanceProb << '\t';
	// 			currentPERpath = currentPERpath->next;
	// 		}
	// 		perPathProb << endl;
	// 	}
	// 
	// 	perPathProb.close();

	/////////////////////////////////
	//   output in sam format
	/////////////////////////////////
	sprintf(filename, "%sresult/path_prob/PER_pathProb_%s.sam", dirPrefix, targetChr);
	perPathProb.open(filename);

	PERalignment *currentPER;
	PERpath *currentPERpath, *outputPERpath;
	/*	path *currentPath;*/
	double best_prob;
	char pathInfo[300], curInfo[300];

	for (long i = 1; i <= perAlignNum; i++)
	{
		currentPER = perAlignList[orderedPERalignlist[i]];
		perPathProb << perList[currentPER->perIndex]->name << '\t' << '0' << '\t' << currentPER->startRead->chromosome << '\t' << currentPER->startRead->rangeLow << '\t' << '0' << '\t';

		//output PER paths
		sprintf(pathInfo, "%s", currentPER->startRead->end);
		currentPERpath = currentPER->possiblePath;
		if (currentPERpath == NULL)
		{
			//perPathProb << currentPER->endRead->rangeLow - currentPER->startRead->rangeHigh - 1 << 'M';
			sprintf(curInfo, "%ldV", currentPER->endRead->rangeLow - currentPER->startRead->rangeHigh - 1);
			strcat(pathInfo, curInfo);
		}
		else
		{
			best_prob = -1.0;
			outputPERpath = NULL;

			while (currentPERpath != NULL)
			{
				if (currentPERpath->probability > best_prob)
				{
					best_prob = currentPERpath->probability;
					outputPERpath = currentPERpath;
				}

				currentPERpath = currentPERpath->next;
			}

			//output exact path info using xxMxxN format
			if (outputPERpath != NULL)
			{
				if (outputPERpath->pathInfo == NULL)
				{
					//on the same block
					//perPathProb << currentPER->endRead->rangeLow - currentPER->startRead->rangeHigh - 1 << 'M';
					sprintf(curInfo, "%ldM", currentPER->endRead->rangeLow - currentPER->startRead->rangeHigh - 1);
					strcat(pathInfo, curInfo);
				} 
				else
				{
					//on different block
					//perPathProb << currentPER->startRead->onBlock->end_position - currentPER->startRead->rangeHigh << 'M' << outputPERpath->pathInfo->info << currentPER->endRead->rangeLow - currentPER->endRead->onBlock->start_position << 'M';
					sprintf(curInfo, "%ldM%s%ldM", currentPER->startRead->onBlock->end_position - currentPER->startRead->rangeHigh, outputPERpath->pathInfo->info, currentPER->endRead->rangeLow - currentPER->endRead->onBlock->start_position);
					strcat(pathInfo, curInfo);
				}
			}
		}

		strcat(pathInfo, currentPER->endRead->end);

		perPathProb << pathInfo << '\t';
		//pathInfoProcess(pathInfo);
		//perPathProb << pathInfo << '\t';


		perPathProb << "*\t0\t0\tAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIII\t" << currentPER->probability << endl;
	}

	perPathProb.close();



	return;
}


void output_junctionSupport()
{
	long i, j, startIndex;
	double p_pos, p_path;
	PERpath *currentPERpath, *outputPERpath;
	PERalignment *currentAlign, *bestAlign;
	pathJunction *pathJunc;


	//compute junction support and support probability
	//output junction support by path support
	// 	for (i = 0; i <= junctionNum; i++)
	// 	{
	// 		junctionSupportPath[i] = 0;
	// 		junctionSupportProbPath[i] = 0.0;
	// 	}
	// 
	// 	for (i = 1; i <= perAlignNum; i++)
	// 	{
	// 		p_pos = perAlignList[i]->probability * perAlignList[i]->tau;
	// 		currentPERpath = perAlignList[i]->possiblePath;
	// 		while (currentPERpath != NULL)
	// 		{
	// 			p_path = currentPERpath->probability * currentPERpath->theta * p_pos;
	// 
	// 			if (currentPERpath->pathInfo != NULL)
	// 			{
	// 				pathJunc = currentPERpath->pathInfo->junctionsInPath;
	// 			} 
	// 			else
	// 			{
	// 				pathJunc = NULL;
	// 			}
	// 
	// 			while (pathJunc != NULL)
	// 			{
	// 				j = pathJunc->junctionID;
	// 				if (j > 0)
	// 				{
	// 					junctionSupportPath[j]++;
	// 					junctionSupportProbPath[j] += p_path;
	// 				}
	// 				pathJunc = pathJunc->next;
	// 			}
	// 
	// 			currentPERpath = currentPERpath->next;
	// 		}
	// 	}


	//compute junction support and support probability with most probable paths
	//output junction support by path support
	double best_prob;

	for (i = 0; i <= junctionNum; i++)
	{
		junctionSupportMPPath[i] = 0;
		junctionSupportProbMPPath[i] = 0.0;
	}

	//output junction support by read junction
	readJunction *currentJunc;

	//compute junction support and support probability
	for (i = 0; i <= junctionNum; i++)
	{
		junctionSupportRead[i] = 0;
		junctionSupportProbRead[i] = 0.0;
	}

	for (i = 1; i <= perNum; i++)
	{
		startIndex = perList[i]->firstAlignmentIndex;
		bestAlign = NULL;
		best_prob = -1.0;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			if (currentAlign->tau > best_prob)
			{
				best_prob = currentAlign->tau;
				bestAlign = currentAlign;
			}
		}		

		if (bestAlign == NULL)
		{
			continue;
		}

		p_pos = bestAlign->probability;

		if (bestAlign->probability > MIN_HD_PER_PROB)
		{
			currentJunc = bestAlign->startRead->readJunc;
			while (currentJunc != NULL)
			{
				j = currentJunc->junctionIndex;
				if (j > 0)
				{
					junctionSupportRead[j]++;
					junctionSupportProbRead[j] += p_pos;
				}
				currentJunc = currentJunc->next;
			}

			currentJunc = bestAlign->endRead->readJunc;
			while (currentJunc != NULL)
			{
				j = currentJunc->junctionIndex;
				if (j > 0)
				{
					junctionSupportRead[j]++;
					junctionSupportProbRead[j] += p_pos;
				}
				currentJunc = currentJunc->next;
			}
		}


		currentPERpath = bestAlign->possiblePath;

		best_prob = -1.0;
		outputPERpath = NULL;

		while (currentPERpath != NULL)
		{
			if (currentPERpath->theta > best_prob)
			{
				best_prob = currentPERpath->theta;
				outputPERpath = currentPERpath;
			}

			currentPERpath = currentPERpath->next;
		}


		if (outputPERpath != NULL)
		{
			p_path = outputPERpath->probability * p_pos;

			if (outputPERpath->probability > MIN_HD_PATH_PROB)
			{
				if (outputPERpath->pathInfo != NULL)
				{
					pathJunc = outputPERpath->pathInfo->junctionsInPath;
				} 
				else
				{
					pathJunc = NULL;
				}

				while (pathJunc != NULL)
				{
					j = pathJunc->junctionID;
					if (j > 0)
					{
						junctionSupportMPPath[j]++;
						junctionSupportProbMPPath[j] += p_path;
					}
					pathJunc = pathJunc->next;
				}
			}
		}		
	}

	return;
}

void output_junctionProb(char *suffix)
{
	char filename[500];
	long temp;
	ofstream outputJuncSupport;

	sprintf(filename, "%sresult/junction_support/splice%s.txt", dirPrefix, suffix);
	string s,line;

	ifstream prevFile;
	prevFile.open(filename);

	while (getline(prevFile, line))
		s += line + "\n";

	prevFile.close();


	outputJuncSupport.open(filename);
	outputJuncSupport << s;




	long i;
	for (long j = 1; j <= junctionNum; j++)
	{
		i = orderedJunctionList[j];

		if (junctionList[i]->fusion == false)
		{
			if (CROSSCHROMOSOME == true && writeSameChrSplice == false)
			{
			}
			else
			{
				if (junctionList[i]->exchange == false)
					outputJuncSupport << i << '\t' <<  junctionList[i]->chromosome_start << '\t' << junctionList[i]->strand_start << junctionList[i]->strand_end << '\t' << junctionList[i]->start_end << '\t' << junctionList[i]->end_start << '\t' << junctionList[i]->isTrue << '\t' << junctionSupportMPPath[i] << '\t' << junctionSupportRead[i] << endl;
				else
					outputJuncSupport << i << '\t' <<  junctionList[i]->chromosome_start << '\t' << junctionList[i]->strand_end << junctionList[i]->strand_start << '\t' << junctionList[i]->end_start << '\t' << junctionList[i]->start_end << '\t' << junctionList[i]->isTrue << '\t' << junctionSupportMPPath[i] << '\t' << junctionSupportRead[i] << endl;
			}
		}
		else 
		{
			if (CROSSCHROMOSOME == true && writeSameChrSplice == false && strcmp(junctionList[i]->chromosome_start, junctionList[i]->chromosome_end) == 0)
			{
			}
			else
			{
				if (junctionList[i]->exchange == false)
				{
					outputFusion << i << '\t' << junctionList[i]->chromosome_start << '_' << junctionList[i]->chromosome_end << '\t' << junctionList[i]->strand_start << junctionList[i]->strand_end << '\t';
					if (junctionList[i]->start_end > MAX_CHR_LENGTH)
					{
						temp = junctionList[i]->start_end - MAX_CHR_LENGTH;
						outputFusion << temp << '\t';
					} 
					else
					{
						outputFusion << junctionList[i]->start_end << '\t';
					}

					if (junctionList[i]->end_start > MAX_CHR_LENGTH)
					{
						temp = junctionList[i]->end_start - MAX_CHR_LENGTH;
						outputFusion << temp << '\t';
					} 
					else
					{
						outputFusion << junctionList[i]->end_start << '\t';
					}
					outputFusion << junctionList[i]->isTrue << '\t' << junctionSupportMPPath[i] << '\t' << junctionSupportRead[i] << endl;
				}
				else
				{
					outputFusion << i << '\t' << junctionList[i]->chromosome_end << '_' << junctionList[i]->chromosome_start << '\t' << junctionList[i]->strand_end << junctionList[i]->strand_start << '\t';
					if (junctionList[i]->end_start > MAX_CHR_LENGTH)
					{
						temp = junctionList[i]->end_start - MAX_CHR_LENGTH;
						outputFusion << temp << '\t';
					} 
					else
					{
						outputFusion << junctionList[i]->end_start << '\t';
					}

					if (junctionList[i]->start_end > MAX_CHR_LENGTH)
					{
						temp = junctionList[i]->start_end - MAX_CHR_LENGTH;
						outputFusion << temp << '\t';
					} 
					else
					{
						outputFusion << junctionList[i]->start_end << '\t';
					}
					outputFusion << junctionList[i]->isTrue << '\t' << junctionSupportMPPath[i] << '\t' << junctionSupportRead[i] << endl;
				}
			}
		}
	}

	outputJuncSupport.close();

	return;
}

/*
void output_PERprob()
{
	ofstream output_PERprob;
	char tmp[500];
	sprintf(tmp, "%sresult/PER_prob/PERprob_%s.sam", dirPrefix, targetChr);
	output_PERprob.open(tmp);

	long i, startIndex, j;
	PERalignment *currentAlign, *bestAlign;
	PERpath *currentPERpath, *outputPERpath;
	double best_prob;
	char pathInfo[300], curInfo[300];

	for (i = 1; i <= perNum; i++)
	{
		startIndex = perList[i]->firstAlignmentIndex;

		if (perList[i]->numAlignment == 1 && (perAlignList[startIndex]->endRead->rangeLow - perAlignList[startIndex]->startRead->rangeHigh > 200 || perAlignList[startIndex]->startRead->rangeLow - perAlignList[startIndex]->endRead->rangeHigh > 200))
		{
		}
		else 
		{
			statNumReadHardMap++;
		}


		//output_PERprob << perList[i]->name << '\t' << perList[i]->probability << '\t' ;

		bestAlign = NULL;
		best_prob = -1.0;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			if (currentAlign->tau > best_prob)
			{
				best_prob = currentAlign->tau;
				bestAlign = currentAlign;
			}
		}		

		if (bestAlign == NULL)// || bestAlign->probability <= MIN_HD_PER_PROB)
		{
			continue;
		}

		//output PER paths
		currentPERpath = bestAlign->possiblePath;
		sprintf(pathInfo, "%s", bestAlign->startRead->end);
		outputPERpath = NULL;

		if (currentPERpath == NULL)
		{
			//perPathProb << bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1 << 'M';
			sprintf(curInfo, "%ldV", bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1);
			strcat(pathInfo, curInfo);
		}
		else
		{
			best_prob = -1.0;

			while (currentPERpath != NULL)
			{
				if (currentPERpath->probability > best_prob)
				{
					best_prob = currentPERpath->probability;
					outputPERpath = currentPERpath;
				}

				currentPERpath = currentPERpath->next;
			}

			//output exact path info using xxMxxN format
			if (outputPERpath != NULL)
			{
				if (outputPERpath->pathInfo == NULL)
				{
					//on the same block
					//perPathProb << bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1 << 'M';

					sprintf(curInfo, "%ldM", bestAlign->endRead->crucialPoint - bestAlign->startRead->crucialPoint);
					strcat(pathInfo, curInfo);
				} 
				else
				{
					//on different block
					//perPathProb << bestAlign->startRead->onBlock->end_position - bestAlign->startRead->rangeHigh << 'M' << outputPERpath->pathInfo->info << bestAlign->endRead->rangeLow - bestAlign->endRead->onBlock->start_position << 'M';
					if (bestAlign->startRead->crucialStrand == true)
					{
						sprintf(curInfo, "%ldM", bestAlign->startRead->onBlock->end_position - bestAlign->startRead->crucialPoint);
					}
					else 
					{
						sprintf(curInfo, "-%ldM", bestAlign->startRead->crucialPoint - bestAlign->startRead->onBlock->start_position);
					}
					strcat(curInfo, outputPERpath->pathInfo->info);
					if (bestAlign->endRead->crucialStrand == true)
					{
						sprintf(curInfo, "%s-%ldM", curInfo, bestAlign->endRead->onBlock->end_position - bestAlign->endRead->crucialPoint);
					}
					else 
					{
						sprintf(curInfo, "%s%ldM", curInfo, bestAlign->endRead->crucialPoint - bestAlign->endRead->onBlock->start_position);
					}

					strcat(pathInfo, curInfo);
				}

				statNumReadPaired++;
				if (bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh > 200 || bestAlign->startRead->rangeLow - bestAlign->endRead->rangeHigh > 200)
				{
					statNumReadAbnormalDistance++;
				}
			}
			else
			{
				cout << "abnormal" << endl;
				exit(1);
			}
		}

		strcat(pathInfo, bestAlign->endRead->end);



		//		output_PERprob << perList[i]->name << '\t' << bestAlign->probability << '\t' << outputPERpath->probability << '\t' << bestAlign->startRead->crucialStrand << bestAlign->endRead->crucialStrand << '\t' << bestAlign->startRead->chromosome << '_' << bestAlign->endRead->chromosome << '\t' << bestAlign->startRead->rangeLow << '\t' << '0' << '\t';
		output_PERprob << perList[i]->name << '\t' << bestAlign->startRead->crucialStrand << bestAlign->endRead->crucialStrand << '\t' << bestAlign->startRead->chromosome << '_' << bestAlign->endRead->chromosome << '\t' << bestAlign->startRead->rangeLow << '\t' << '0' << '\t';

		output_PERprob << pathInfo;
		//pathInfoProcess(pathInfo);
		//perPathProb << pathInfo << '\t';


		output_PERprob << "\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIII\t" << endl;
	}



	output_PERprob.close();

	return;
}*/

void output_PERprob()
{
	ofstream output_PERprob, output_fusionRead;
	char tmp[500];
	sprintf(tmp, "%sresult/PER_prob/PERprob_%s.sam", dirPrefix, targetChr);
	output_PERprob.open(tmp);
	sprintf(tmp, "%sresult/fusionRead/fusionRead_%s.sam", dirPrefix, targetChr);
	output_fusionRead.open(tmp);

	long i, startIndex, j;
	PERalignment *currentAlign, *bestAlign;
	PERpath *currentPERpath, *outputPERpath;
	double best_prob;
	char pathInfo[300];

	for (i = 1; i <= perNum; i++)
	{
		startIndex = perList[i]->firstAlignmentIndex;

		if (perList[i]->numAlignment == 1 && (perAlignList[startIndex]->endRead->rangeLow - perAlignList[startIndex]->startRead->rangeHigh > 200 || perAlignList[startIndex]->startRead->rangeLow - perAlignList[startIndex]->endRead->rangeHigh > 200))
		{
		}
		else 
		{
			statNumReadHardMap++;
		}


		//output_PERprob << perList[i]->name << '\t' << perList[i]->probability << '\t' ;

		bestAlign = NULL;
		best_prob = -1.0;

		for (j = 0; j < perList[i]->numAlignment; j++)
		{
			currentAlign = perAlignList[startIndex + j];
			if (currentAlign->tau > best_prob)
			{
				best_prob = currentAlign->tau;
				bestAlign = currentAlign;
			}
		}		

		if (bestAlign == NULL)// || bestAlign->probability <= MIN_HD_PER_PROB)
		{
			continue;
		}

		//output PER paths
		currentPERpath = bestAlign->possiblePath;
		outputPERpath = NULL;

		if (currentPERpath == NULL || bestAlign->startRead->crucialStrand == bestAlign->endRead->crucialStrand || bestAlign->probability <= MIN_HD_PER_PROB)
		{
			//perPathProb << bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1 << 'M';
			sprintf(pathInfo, "*");

			output_PERprob << perList[i]->name << '\t' << "0" << '\t' << bestAlign->startRead->chromosome << '\t' << bestAlign->startRead->rangeLow << '\t' << '0' << '\t' << "*";
			output_PERprob << "\t*\t0\t0\tNA\tNA\t" << endl;
		}
		else
		{		
			best_prob = -1.0;

			while (currentPERpath != NULL)
			{
				if (currentPERpath->probability > best_prob)
				{
					best_prob = currentPERpath->probability;
					outputPERpath = currentPERpath;
				}

				currentPERpath = currentPERpath->next;
			}

			//output exact path info using xxMxxN format
			if (outputPERpath != NULL && outputPERpath->probability > MIN_HD_PATH_PROB)
			{
				if (outputPERpath->pathInfo == NULL)
				{
					//on the same block
//					sprintf(pathInfo, "%ldM", bestAlign->endRead->crucialPoint - bestAlign->startRead->crucialPoint);
					if (bestAlign->endRead->rangeLow >= bestAlign->startRead->rangeHigh)
					{
						sprintf(pathInfo, "%ldM", bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh - 1);
					} 
					else if (bestAlign->startRead->rangeLow >= bestAlign->endRead->rangeHigh)
					{
						sprintf(pathInfo, "-%ldM", bestAlign->startRead->rangeLow - bestAlign->endRead->rangeHigh - 1);
					}
					else
					{
						sprintf(pathInfo, "*");
					}
				} 
				else
				{
					//on different block
//output_PERprob << bestAlign->startRead->end << '\t' << outputPERpath->pathInfo->info << '\t' << bestAlign->endRead->end << '\t';
					if (bestAlign->startRead->crucialStrand == true)
					{
						if (bestAlign->startRead->onBlock->end_position >= bestAlign->startRead->rangeHigh)
						{
							sprintf(pathInfo, "%ldM", bestAlign->startRead->onBlock->end_position - bestAlign->startRead->rangeHigh);
						} 
						else
						{
							sprintf(pathInfo, "*");
						}
					}
					else 
					{
						if (bestAlign->startRead->rangeLow >= bestAlign->startRead->onBlock->start_position)
						{
							sprintf(pathInfo, "-%ldM", bestAlign->startRead->rangeLow - bestAlign->startRead->onBlock->start_position);
						} 
						else
						{
							sprintf(pathInfo, "*");
						}
					}
					if (pathInfo[0] != '*')
					{
						strcat(pathInfo, outputPERpath->pathInfo->info);
						if (bestAlign->endRead->crucialStrand == true)
						{
							if (bestAlign->endRead->onBlock->end_position >= bestAlign->endRead->rangeHigh)
							{
								sprintf(pathInfo, "%s-%ldM", pathInfo, bestAlign->endRead->onBlock->end_position - bestAlign->endRead->rangeHigh);
							} 
							else
							{
								sprintf(pathInfo, "*");
							}
						}
						else 
						{
							if (bestAlign->endRead->rangeLow >= bestAlign->endRead->onBlock->start_position)
							{
								sprintf(pathInfo, "%s%ldM", pathInfo, bestAlign->endRead->rangeLow - bestAlign->endRead->onBlock->start_position);
							} 
							else
							{
								sprintf(pathInfo, "*");
							}
						}
					}
				}

				statNumReadPaired++;
				if (bestAlign->endRead->rangeLow - bestAlign->startRead->rangeHigh > 200 || bestAlign->startRead->rangeLow - bestAlign->endRead->rangeHigh > 200)
				{
					statNumReadAbnormalDistance++;
				}
			}
			else
			{
				sprintf(pathInfo, "*");
			}



			if (bestAlign->startRead->crucialStrand == true && bestAlign->endRead->crucialStrand == false)
			{
				if (bestAlign->startRead->rangeHigh <= bestAlign->endRead->rangeLow)
				{
//output_PERprob << bestAlign->startRead->end << '\t' << pathInfo << '\t' << bestAlign->endRead->end << '\t' << "+\t" << bestAlign->startRead->rangeLow << "\t" << bestAlign->endRead->rangeLow << "\t";
					if (pathInfo[0] != '*')
						pathInfoProcess(pathInfo, false, bestAlign->startRead->end, bestAlign->endRead->end);
					else
						sprintf(pathInfo, "*");
					output_PERprob << perList[i]->name << '\t' << "0" << '\t' << bestAlign->startRead->chromosome << '\t' << bestAlign->startRead->rangeLow << '\t' << '0' << '\t' << pathInfo;
					output_PERprob << "\t*\t0\t0\tNA\tNA\t" << endl;
				} 
				else
				{
					output_fusionRead << perList[i]->name << '\t' << "0" << '\t' << bestAlign->startRead->chromosome << '\t' << bestAlign->startRead->rangeLow << '\t' << '0' << '\t';
					output_fusionRead << bestAlign->startRead->end << '\t' << pathInfo << '\t' << bestAlign->endRead->end << endl;
				}				
			} 
			else if (bestAlign->startRead->crucialStrand == false && bestAlign->endRead->crucialStrand == true)
			{
				if (bestAlign->startRead->rangeLow >= bestAlign->endRead->rangeHigh)
				{
//output_PERprob << bestAlign->startRead->end << '\t' << pathInfo << '\t' << bestAlign->endRead->end << '\t' << "-\t" << bestAlign->startRead->rangeLow << "\t" << bestAlign->endRead->rangeLow << "\t";
					if (pathInfo[0] != '*')
						pathInfoProcess(pathInfo, true, bestAlign->startRead->end, bestAlign->endRead->end);
					else
						sprintf(pathInfo, "*");
					output_PERprob << perList[i]->name << '\t' << "16" << '\t' << bestAlign->startRead->chromosome << '\t' << bestAlign->endRead->rangeLow << '\t' << '0' << '\t' << pathInfo;
					output_PERprob << "\t*\t0\t0\tNA\tNA\t" << endl;
				} 
				else
				{
					output_fusionRead << perList[i]->name << '\t' << "16" << '\t' << bestAlign->startRead->chromosome << '\t' << bestAlign->endRead->rangeLow << '\t' << '0' << '\t';
					output_fusionRead << bestAlign->endRead->end << '\t' << pathInfo << '\t' << bestAlign->startRead->end << endl;
				}
			}
		}
	}

	output_PERprob.close();
	output_fusionRead.close();

	return;
}

void output_result()
{
	//	output_PERpath_prob();

	output_junctionSupport();

	char suffix[10];
	sprintf(suffix, "");
	output_junctionProb(suffix);

	output_PERprob();

	return;
}








/**********************************************************
***********************************************************
*********************** MAIN ******************************
***********************************************************
**********************************************************/

void initialization()
{
	junctionNum = 0;
	chrBlockNum = 0;

	//build initial block for the whole chromosome
	block *newBlock;

	if (CROSSCHROMOSOME == true)
	{
		newBlock = new block;
		newBlock->start_position = MAX_CHR_LENGTH; 
		newBlock->end_position = MAX_CHR_LENGTH * 2;
		strcpy(newBlock->chromosome, targetChr_end);
		newBlock->next = NULL;
		newBlock->edge = NULL;
		chrBlocks_2ndStart = newBlock;
		chrBlocks_tail = newBlock;
		chrBlockNum++;
	} 
	else
	{
		chrBlocks_2ndStart = NULL;
		chrBlocks_tail = NULL;
	}


	newBlock = new block;
	newBlock->start_position = 0; 
	newBlock->end_position = MAX_CHR_LENGTH;
	strcpy(newBlock->chromosome, targetChr_start);
	newBlock->next = chrBlocks_2ndStart;
	newBlock->edge = NULL;
	chrBlocks = newBlock;
	if (chrBlocks_tail == NULL)
	{
		chrBlocks_tail = chrBlocks;
	}
	chrBlockNum++;

	lastJunctionSearchBlock = chrBlocks;

	readNum = 0;
	perAlignNum = 0;
	perNum = 0;

	pathClusterList = NULL;
	pathClusterListTail = NULL;

	statNumRead = 0;
	statNumSingleRead = 0;
	statNumPairRead = 0;
	statNumReadPaired = 0;
	statNumReadAbnormalDistance = 0;
	statNumReadHardMap = 0;

	srand ( time(NULL) );

	//initialize arrays
	junctionList = new junction* [MAX_JUNCTION_NUM+1]; //list for all junctions
	orderedJunctionList = new long[MAX_JUNCTION_NUM+1]; //index list of junctionList after sorting
	sortKey_Junction = new long[MAX_JUNCTION_NUM+1];
	// 	junctionSupportPath = new long[MAX_JUNCTION_NUM+1];
	// 	junctionSupportProbPath = new double[MAX_JUNCTION_NUM+1];


	readList = new reads* [MAX_LENGTH+1];
	orderedReadList = new long[MAX_LENGTH+1];
	sortKey_read = new long[MAX_LENGTH+1];

	perAlignList = new PERalignment* [MAX_LENGTH / 2 +1];
	orderedPERalignlist = new long[MAX_LENGTH / 2 +1];
	sortKey_PERalign = new long[MAX_LENGTH / 2 +1];

	perList = new PER* [MAX_LENGTH / 2 +1];



	mergeSort_Larray = new long[MAX_LENGTH +1];
	mergeSort_Rarray = new long[MAX_LENGTH +1];
	mergeSort_LorderedList = new long[MAX_LENGTH +1];
	mergeSort_RorderedList = new long[MAX_LENGTH +1];

	pathQueue = new path* [MAX_PATH_NUM + 1];

	return;	
}



void process_single_chromosome()
{
	//process single chromosome
	char junctionFileName[500];
	char PERfileName[500], singleFilename1[500], singleFilename2[500], SingleFilename3[500];

	sprintf(PERfileName, "%sdata/PER/%s.txt", dirPrefix, targetChr);
	//	sprintf(PERfileName, "%sdata\\PER\\%s.txt", dirPrefix, targetChr);


	sprintf(singleFilename1, "%sdata/single/%s.txt", dirPrefix, targetChr_start);
	sprintf(singleFilename2, "%sdata/single/%s.txt", dirPrefix, targetChr_end);
	sprintf(SingleFilename3, "%sdata/single/%s%s.txt", dirPrefix, targetChr_start, targetChr_end);
	/*	
	sprintf(singleFilename1, "%sdata\\single\\%s.txt", dirPrefix, targetChr_start);
	sprintf(singleFilename2, "%sdata\\single\\%s.txt", dirPrefix, targetChr_end);
	sprintf(SingleFilename3, "%sdata\\single\\%s%s.txt", dirPrefix, targetChr_start, targetChr_end);
	*/


	initialization();

	if (CROSSCHROMOSOME == true)
	{
		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_start);
		input_junction(junctionFileName);
		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_end);
		input_junction(junctionFileName);
		sprintf(junctionFileName, "%sdata/junction/%s%s.txt", dirPrefix, targetChr_start, targetChr_end);
		input_junction(junctionFileName);


		/*		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_start);
		input_junction(junctionFileName);
		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_end);
		input_junction(junctionFileName);
		sprintf(junctionFileName, "%sdata/junction/%s%s.txt", dirPrefix, targetChr_start, targetChr_end);
		input_junction(junctionFileName);*/
	} 
	else
	{
		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_start);
		input_junction(junctionFileName);


		/*		sprintf(junctionFileName, "%sdata/junction/%s.txt", dirPrefix, targetChr_start);
		input_junction(junctionFileName);*/
	}

	input_PERs(PERfileName, singleFilename1, singleFilename2, SingleFilename3); //input PERs and construct junctions

	runTimeStart = clock();

	//	cout << "input complete" << endl;

	reads_map_block();

	//	cout << "reads mapping complete" << endl;

	compute_PER_distance();

	//	cout << "distance computing complete" << endl;

	runTime = ((double)clock() - runTimeStart)/ CLOCKS_PER_SEC;
	outputSummary << runTime << '\t';


	PERalignment *currentAlign;
	for (long i = 1; i <= perAlignNum; i++)
	{
		/************************************************************************/
		/* Modify the boundary of the PER alignment for junction mismatch       */
		/************************************************************************/
		currentAlign = perAlignList[i];
		if (currentAlign->startRead->crucialStrand == true)
			currentAlign->startRead->crucialPoint += 2;
		else
			currentAlign->startRead->crucialPoint -= 2;

		if (currentAlign->endRead->crucialStrand == true)
			currentAlign->endRead->crucialPoint += 2;
		else
			currentAlign->endRead->crucialPoint -= 2;
	}


	//	cout << "EMbegins   " << endl;
	expectation_maximization();

	output_result();
	return;
}

bool input_setting(char *settingfile)
{
	ifstream setting;
	setting.open(settingfile);

	string info;
	size_t found;
	int count = 0;
	double temp, temp2;

	while (getline(setting, info))
	{
		found = info.find("PERPathDistanceMean");
		if (found != string::npos)
		{
			setting >> PERPathDistanceMean;
			count++;
		}

		found = info.find("PERPathDistanceSD");
		if (found != string::npos)
		{
			setting >> PERPathDistanceSD;
			count++;
		}
		found = info.find("MAX_JUNCTION_NUM");
		if (found != string::npos)
		{
			setting >> MAX_JUNCTION_NUM;
			count++;
		}
		found = info.find("MAX_LENGTH");
		if (found != string::npos)
		{
			setting >> MAX_LENGTH;
			count++;
		}
		found = info.find("MAX_PATH_NUM");
		if (found != string::npos)
		{
			setting >> MAX_PATH_NUM;
			count++;
		}
		found = info.find("MIN_PATH_PROBABILITY");
		if (found != string::npos)
		{
			setting >> temp2;
			count++;
		}
		found = info.find("MIN_HD_PER_PROB");
		if (found != string::npos)
		{
			setting >> MIN_HD_PER_PROB;
			count++;
		}
		found = info.find("MIN_HD_PATH_PROB");
		if (found != string::npos)
		{
			setting >> temp; 
			count++;
		}
	}

	PERPathDistanceVAR = PERPathDistanceSD * PERPathDistanceSD; 
	//	PERPathDistanceProbScalar = 1 / sqrt(2 * 3.1415926 * PERPathDistanceVAR); //pdf on the mean
	MAX_PATH_LENGTH = PERPathDistanceMean + PERPathDistanceSD * 5;
	//	MIN_HD_PATH_PROB = 1 / exp(temp * temp / (2 * PERPathDistanceVAR));
	MIN_PATH_PROBABILITY = 1 / exp(temp2 * temp2 / (2 * PERPathDistanceVAR));
	MIN_HD_PATH_PROB = temp;

	if (count == 8)
		return true;
	else
		return false;
}

void clearAll()
{
	long i;

	for (i = 1; i <= junctionNum; i++)
	{
		delete junctionList[i];
	}
	for (i = 1; i <= readNum; i++)
	{
		delete readList[i];
	}
	for (i = 1; i <= perAlignNum; i++)
	{
		delete perAlignList[i];
	}
	for (i = 1; i <= perNum; i++)
	{
		delete perList[i];
	}

	pathCluster *curPathCluster, *nextPathCluster;
	path *curPath, *nextPath;
	pathJunction *curJunc, *nextJunc;
	curPathCluster = pathClusterList;
	// 	while (curPathCluster != NULL)
	// 	{
	// 		curPath = curPathCluster->clusterPaths;
	// 		while (curPath != NULL)
	// 		{
	// 			curJunc = curPath->junctionsInPath;
	// 			while (curJunc != NULL)
	// 			{
	// 				nextJunc = curJunc->next;
	// 				delete curJunc;
	// 				curJunc = nextJunc;
	// 			}
	// 
	// 			nextPath = curPath->next;
	// 			delete curPath;
	// 			curPath = nextPath;
	// 		}
	// 		nextPathCluster = curPathCluster->next;
	// 		delete curPathCluster;
	// 		curPathCluster = nextPathCluster;
	// 	}

	block *curBlk, *nextBlk;
	curBlk = chrBlocks;
	while (curBlk != NULL)
	{
		nextBlk = curBlk->next;
		delete curBlk;
		curBlk = nextBlk;
	}

	delete [] junctionList;
	delete [] orderedJunctionList;

	delete [] readList;
	delete [] orderedReadList;

	delete [] perAlignList;
	delete [] orderedPERalignlist;

	delete [] perList;

	delete [] junctionSupportMPPath;
	delete [] junctionSupportProbMPPath;
	delete [] junctionSupportRead;
	delete [] junctionSupportProbRead;

	delete [] perAlignList_EM;
}

int main(int argc, char* argv[])
{
	char chromosome1[30], chromosome2[30], targetChromosome[100], filename[500];
	///*
	if (argc != 6)
	{
		cout << argv[0] << "\t<chr#>\t<chr#>\t<exepath>\t<targetpath>\t<all/cross>" << endl;
		return 1;
	}
	//sprintf(targetChr_start, "chr1");
	strcpy(targetChr_start, argv[1]);
	strcpy(targetChr_end, argv[2]);
	//*/

	//	strcpy(targetChr_start, "chr17");
	//	strcpy(targetChr_end, "chr20");

	sprintf(targetChr, "%s%s", targetChr_start, targetChr_end);

	if (strcmp(targetChr_start, targetChr_end) == 0)
	{
		CROSSCHROMOSOME = false;
	} 
	else
	{
		CROSSCHROMOSOME = true;
	}


	//	sprintf(dirPrefix, "./");	
	strcpy(dirPrefix, argv[4]);
	//	sprintf(dirPrefix, "C:\\Users\\yin\\Documents\\PER\\");

	sprintf(filename, "%ssetting.txt", argv[3]);
	if (input_setting(filename) == false)
	{
		cout << "Setting file error!" << endl;
		return 1;
	}

	if (strcmp(argv[5], "all") == 0)
		writeSameChrSplice = true;
	else if (strcmp(argv[5], "cross") == 0)
		writeSameChrSplice = false;
	else
	{
		cout << "Unrecognized parameter: 5th parameter should be all or cross." << endl;
		return 1;
	}


	if (CROSSCHROMOSOME == true)
	{
		sprintf(filename, "%sresult/outputPathSummary.txt", dirPrefix);
	} 
	else
	{
		sprintf(filename, "%sresult/outputPathSummarySame.txt", dirPrefix);
	}
	//	outputPathSummary.open(filename);


	sprintf(filename, "%sresult/Summary.txt", dirPrefix);
	//sprintf(filename, "%sresult\\Summary.txt", dirPrefix);

	string s,line;

	ifstream prevFile;
	prevFile.open(filename);

	while (getline(prevFile, line))
		s += line + "\n";

	prevFile.close();


	outputSummary.open(filename);
	outputSummary << s << targetChr_start << '_' << targetChr_end << "\t";



	/************************************************************************/
	/* OUTPUT FUSION EVENTS                                                 */
	/************************************************************************/
	sprintf(filename, "%sresult/junction_support/fusion.txt", dirPrefix);

	string sfusion;

	prevFile.open(filename);

	while (getline(prevFile, line))
		sfusion += line + "\n";

	prevFile.close();

	outputFusion.open(filename);
	outputFusion << sfusion;


	process_single_chromosome();

	outputSummary << perNum << '\t' << statNumReadHardMap << '\t' << statNumReadPaired << '\t' << statNumReadAbnormalDistance;
	outputFusion.close();
	//	outputPathSummary.close();
	outputSummary.close();


	//	clearAll();

	// 	ifstream chrNameFile;
	// 	chrNameFile.open("C:\\Users\\yin\\Desktop\\PER\\PER_36bp\\data\\parsedPER\\ChromosomeName.txt");
	// 
	// 	int tmp, i, j;
	// 
	// 	for (tmp = 0; tmp < 100; tmp++)
	// 	{
	// 		targetChromosome[tmp] = '\0';
	// 	}
	// 
	// 	chrNameFile >> targetChromosome;
	// 
	// 
	// 	while (targetChromosome[1] != '\0')
	// 	{
	// 		
	// 		chromosome1[0] = 'c';
	// 
	// 		for (i = 1; targetChromosome[i] != 'c'; i++)
	// 		{
	// 			chromosome1[i] = targetChromosome[i];
	// 		}
	// 		chromosome1[i] = '\0';
	// 
	// 		for (j = 0; targetChromosome[i] != '\0'; i++, j++)
	// 		{
	// 			chromosome2[j] = targetChromosome[i];
	// 		}
	// 		chromosome2[j] = '\0';
	// 
	// 		if (strcmp(chromosome1, chromosome2) == 0)
	// 		{
	// 			process_single_chromosome(chromosome1);
	// 		}
	// 
	// 		for (tmp = 0; tmp < 100; tmp++)
	// 		{
	// 			targetChromosome[tmp] = '\0';
	// 		}
	// 
	// 		chrNameFile >> targetChromosome;	
	// 	}
	// 
	// 
	// 	chrNameFile.close();

	return 0;
}