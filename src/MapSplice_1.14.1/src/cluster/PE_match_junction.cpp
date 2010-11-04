#include "PE_match_junction.h"

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

junction::junction()
{
	start = 0;
	end = 0;

	return;
}


used_junction::used_junction()
{
	junctionID = -1;
	next = NULL;
	
	return;
}

pe_distance::pe_distance(long start, long end)
{
	existing_distance = end - start + 1;
	milestone = start;
	
	junc = NULL;
	next = NULL;
	
	return;
}

void junction::compute_range()
{
	int i, blksize = 0;

	for (i = 0; blockSizes[i] != '\0'; i++)
	{
		if (blockSizes[i] == ',')
		{
			break;
		} 
		else
		{
			blksize = blksize * 10 + blockSizes[i] - 48;
		}
	}

	start_start = start - blksize;

	blksize = 0;
	for (i++; blockSizes[i] != '\0'; i++)
	{
		if (blockSizes[i] == ',')
		{
			break;
		} 
		else
		{
			blksize = blksize * 10 + blockSizes[i] - 48;
		}
	}

	end_end = end + blksize;
	
	return;
}

paired_end :: paired_end()
{
	for (int i = 0; i < 5; i++)
	{
		strand_start[i] = '\0';
		strand_end[i] = '\0';
	}

	start = NULL;
	end = NULL;
	tail_start = NULL;
	tail_end = NULL;

	range_count_start = 0;
	range_count_end = 0;

	startrange_low = MAX;
	startrange_high = 0;
	endrange_low = MAX;
	endrange_high = 0;

	dist = NULL;

	return;
}

void paired_end::compute_distance()
{
	if (strcmp(chromosome_start, chromosome_end) != 0)
	{
		return;
	}
	//assume the junctions list has been sorted based on starting point
	long junction_index;
	junction_index = 1;

	pe_distance *current_dist, *new_dist;
	current_dist = NULL;
	new_dist = NULL;

	used_junction *tmp_usedjunc;
	tmp_usedjunc = NULL;

	int path_num = 1;

	if (endrange_low - startrange_high >= MAX_PEREAD_DISTANCE)
	{
		//ignore
		return;
	}

	while (junction_index <= junction_num)
	{
		if (junction_list[junction_index]->start > endrange_low || path_num > MAX_PATH_NUM)
		{
			//computation completed
			return;
		}

		if (junction_list[junction_index]->start >= startrange_high && junction_list[junction_index]->end <= endrange_low)
		{
			//read supports this junction
			junction_support[junction_index]++;

			current_dist = dist;
			while (current_dist != NULL)
			{
				if (current_dist->milestone <= junction_list[junction_index]->start)
				{
					//can use this junction
					new_dist = new pe_distance(0, 0);
					new_dist->existing_distance = current_dist->existing_distance - (junction_list[junction_index]->end - junction_list[junction_index]->start - 1);
					new_dist->milestone = junction_list[junction_index]->end;

					tmp_usedjunc = new used_junction;
					tmp_usedjunc->junctionID = junction_index;
					tmp_usedjunc->next = current_dist->junc;

					new_dist->junc = tmp_usedjunc;
					new_dist->next = current_dist->next;
					current_dist->next = new_dist;

					path_num++;

					current_dist = current_dist->next;
				}

				current_dist = current_dist->next;
			}
		}

		junction_index++;
	}

	return;
}


cluster::cluster(paired_end* startPE)
{
	strcpy(chromosome_1, startPE->chromosome_start);
	strcpy(chromosome_2, startPE->chromosome_end);

	strcpy(strand_start, startPE->strand_start);
	strcpy(strand_end, startPE->strand_end);

	endpoint1 = startPE->start;
	endpoint2 = startPE->end;
	tail_end1 = startPE->tail_start;
	tail_end2 = startPE->tail_end;

	fulllist = NULL;
	tail_fulllist = NULL;

	range_count_1 = startPE->range_count_start;
	range_count_2 = startPE->range_count_end;
	fulllist_count = 0;

	end1_range_low = startPE->startrange_low;
	end1_range_high = startPE->startrange_high;
	end2_range_low = startPE->endrange_low;
	end2_range_high = startPE->endrange_high;

	return;
}

void cluster::merge_range(paired_end* currentPE)
{
	range *nextrange, *nextPErange, *merge_current, *merge_next;

	if (strcmp(currentPE->chromosome_start, chromosome_1) == 0)
	{
		if (strcmp(currentPE->chromosome_end, chromosome_2) == 0)
		{
			//merge start point 
			nextrange = endpoint1;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->start;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						//overlapping range

						//aggregate current range
						if (nextPErange->start < nextrange->start)
						{
							nextrange->start = nextPErange->start;
						}

						if (nextPErange->end > nextrange->end)
						{
							nextrange->end = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//merge start point ranges
			merge_current = endpoint1;

			if (merge_current != NULL)
			{
				merge_next = merge_current->next;

				while (merge_next != NULL)
				{
					if (merge_current->end >= merge_next->start && merge_current->start <= merge_next->end)
					{
						//merge
						merge_current->start = merge_current->start < merge_next->start ? merge_current->start : merge_next->start;
						merge_current->end = merge_current->end > merge_next->end ? merge_current->end : merge_next->end;

						merge_current->next = merge_next->next;
						delete merge_next;

						range_count_1--;
					}
					else
					{
						merge_current = merge_current->next;
					}

					merge_next = merge_current->next;
				}
			}

			

			//merge end point
			nextrange = endpoint2;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->end;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						//overlapping range

						//aggregate current range
						if (nextPErange->start < nextrange->start)
						{
							nextrange->start = nextPErange->start;
						}

						if (nextPErange->end > nextrange->end)
						{
							nextrange->end = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//merge end point ranges
			merge_current = endpoint2;

			if (merge_current != NULL)
			{
				merge_next = merge_current->next;

				while (merge_next != NULL)
				{
					if (merge_current->end >= merge_next->start && merge_current->start <= merge_next->end)
					{
						//merge
						merge_current->start = merge_current->start < merge_next->start ? merge_current->start : merge_next->start;
						merge_current->end = merge_current->end > merge_next->end ? merge_current->end : merge_next->end;

						merge_current->next = merge_next->next;
						delete merge_next;

						range_count_2--;
					}
					else
					{
						merge_current = merge_current->next;
					}

					merge_next = merge_current->next;
				}
			}
		}
	} 
	else if (strcmp(currentPE->chromosome_start, chromosome_2) == 0)
	{
		if (strcmp(currentPE->chromosome_end, chromosome_1) == 0)
		{
			//check start point 
			nextrange = endpoint2;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->start;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						//overlapping range

						//aggregate current range
						if (nextPErange->start < nextrange->start)
						{
							nextrange->start = nextPErange->start;
						}

						if (nextPErange->end > nextrange->end)
						{
							nextrange->end = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//merge start point ranges
			merge_current = endpoint2;

			if (merge_current != NULL)
			{
				merge_next = merge_current->next;

				while (merge_next != NULL)
				{
					if (merge_current->end >= merge_next->start && merge_current->start <= merge_next->end)
					{
						//merge
						merge_current->start = merge_current->start < merge_next->start ? merge_current->start : merge_next->start;
						merge_current->end = merge_current->end > merge_next->end ? merge_current->end : merge_next->end;

						merge_current->next = merge_next->next;
						delete merge_next;

						range_count_2--;
					}
					else
					{
						merge_current = merge_current->next;
					}

					merge_next = merge_current->next;
				}
			}
			
			
			
			
			//merge end point
			nextrange = endpoint1;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->end;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						//overlapping range

						//aggregate current range
						if (nextPErange->start < nextrange->start)
						{
							nextrange->start = nextPErange->start;
						}

						if (nextPErange->end > nextrange->end)
						{
							nextrange->end = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//merge end point ranges
			merge_current = endpoint1;

			if (merge_current != NULL)
			{
				merge_next = merge_current->next;

				while (merge_next != NULL)
				{
					if (merge_current->end >= merge_next->start && merge_current->start <= merge_next->end)
					{
						//merge
						merge_current->start = merge_current->start < merge_next->start ? merge_current->start : merge_next->start;
						merge_current->end = merge_current->end > merge_next->end ? merge_current->end : merge_next->end;

						merge_current->next = merge_next->next;
						delete merge_next;

						range_count_1--;
					}
					else
					{
						merge_current = merge_current->next;
					}

					merge_next = merge_current->next;
				}
			}
		}
	}


	return;
}


bool cluster::searchANDmodify(paired_end* currentPE, long *overlap_start, long *overlap_end)
{
	bool overlap = false;
	bool overlap_tag_start, overlap_tag_end;

	long tmp1, tmp2, tmp3, tmp4;
	tmp1 = end1_range_low;
	tmp2 = end1_range_high;
	tmp3 = end2_range_low;
	tmp4 = end2_range_high;

	*overlap_start = 0;
	*overlap_end = 0;

	range *nextrange, *nextPErange;

	if (strcmp(currentPE->chromosome_start, chromosome_1) == 0 && strcmp(currentPE->strand_start, strand_start) == 0)
	{
		if (strcmp(currentPE->chromosome_end, chromosome_2) == 0 && strcmp(currentPE->strand_end, strand_end) == 0)
		{
			//check start point 
			nextrange = endpoint1;

			overlap_tag_start = false;
			overlap_tag_end = false;


			while (nextrange != NULL)
			{
				nextPErange = currentPE->start;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						overlap_tag_start = true;

						if (nextPErange->start >= nextrange->start)
						{
							if (nextPErange->end >= nextrange->end)
							{
								*overlap_start += nextrange->end - nextPErange->start + 1;
							} 
							else
							{
								*overlap_start += nextPErange->end - nextPErange->start + 1;
							}
						} 
						else if (nextPErange->end >= nextrange->end)
						{
							*overlap_start += nextrange->end - nextrange->start + 1;
						}
						else
						{
							*overlap_start += nextPErange->end - nextrange->start + 1; 
						}

						if (nextPErange->start < tmp1)
						{
							tmp1 = nextPErange->start;
						}

						if (nextPErange->end > tmp2)
						{
							tmp2 = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//check end point
			nextrange = endpoint2;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->end;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						overlap_tag_end = true;
						
						if (nextPErange->start >= nextrange->start)
						{
							if (nextPErange->end >= nextrange->end)
							{
								*overlap_end += nextrange->end - nextPErange->start + 1;
							} 
							else
							{
								*overlap_end += nextPErange->end - nextPErange->start + 1;
							}
						} 
						else if (nextPErange->end >= nextrange->end)
						{
							*overlap_end += nextrange->end - nextrange->start + 1;
						}
						else
						{
							*overlap_end += nextPErange->end - nextrange->start + 1; 
						}

						if (nextPErange->start < tmp3)
						{
							tmp3 = nextPErange->start;
						}

						if (nextPErange->end > tmp4)
						{
							tmp4 = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//check whether overlap or not
			if (overlap_tag_start == true && overlap_tag_end == true)
			{
				overlap = true;

				end1_range_low = tmp1;
				end1_range_high = tmp2;
				end2_range_low = tmp3;
				end2_range_high = tmp4;

				merge_range(currentPE);
			}
		}
	} 
	else if (strcmp(currentPE->chromosome_start, chromosome_2) == 0 && strcmp(currentPE->strand_start, strand_end) == 0)
	{
		if (strcmp(currentPE->chromosome_end, chromosome_1) == 0 && strcmp(currentPE->strand_end, strand_start) == 0)
		{
			//check start point 
			nextrange = endpoint2;

			overlap_tag_start = false;
			overlap_tag_end = false;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->start;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						overlap_tag_start = true;
						
						if (nextPErange->start >= nextrange->start)
						{
							if (nextPErange->end >= nextrange->end)
							{
								*overlap_start += nextrange->end - nextPErange->start + 1;
							} 
							else
							{
								*overlap_start += nextPErange->end - nextPErange->start + 1;
							}
						} 
						else if (nextPErange->end >= nextrange->end)
						{
							*overlap_start += nextrange->end - nextrange->start + 1;
						}
						else
						{
							*overlap_start += nextPErange->end - nextrange->start + 1; 
						}

						if (nextPErange->start < tmp3)
						{
							tmp3 = nextPErange->start;
						}

						if (nextPErange->end > tmp4)
						{
							tmp4 = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//check end point
			nextrange = endpoint1;

			while (nextrange != NULL)
			{
				nextPErange = currentPE->end;

				while (nextPErange != NULL)
				{
					if (nextPErange->start <= nextrange->end + deviation_end && nextPErange->end >= nextrange->start - deviation_end)
					{
						overlap_tag_end = true;
						
						if (nextPErange->start >= nextrange->start)
						{
							if (nextPErange->end >= nextrange->end)
							{
								*overlap_end += nextrange->end - nextPErange->start + 1;
							} 
							else
							{
								*overlap_end += nextPErange->end - nextPErange->start + 1;
							}
						} 
						else if (nextPErange->end >= nextrange->end)
						{
							*overlap_end += nextrange->end - nextrange->start + 1;
						}
						else
						{
							*overlap_end += nextPErange->end - nextrange->start + 1; 
						}

						if (nextPErange->start < tmp1)
						{
							tmp1 = nextPErange->start;
						}

						if (nextPErange->end > tmp2)
						{
							tmp2 = nextPErange->end;
						}
					}

					nextPErange = nextPErange->next;
				}

				nextrange = nextrange->next;
			}

			//check whether overlap or not
			if (overlap_tag_start = true && overlap_tag_end == true)
			{
				overlap = true;

				end1_range_low = tmp1;
				end1_range_high = tmp2;
				end2_range_low = tmp3;
				end2_range_high = tmp4;

				merge_range(currentPE);
			}
		}
	}


	return overlap;
}



bool cluster::build_fulllist()
{
	if (strcmp(chromosome_1, chromosome_2) != 0)
	{
		//two clone ends are on different chromosomes
		return false;
	}

	//merge the two lists into one 

	range *ptr1, *ptr2, *new_range;
	ptr1 = endpoint1;
	ptr2 = endpoint2;

	while (ptr1 != NULL && ptr2 != NULL)
	{
		if (ptr1->start >= ptr1->end)
		{
			ptr1 = ptr1->next;
			continue;
		}

		if (ptr2->start >= ptr2->end)
		{
			ptr2 = ptr2->next;
			continue;
		}

		if (ptr1->end < ptr2->start)
		{
			//no overlap and *ptr1 is ahead from *ptr2
			new_range = new range;

			new_range->range_copy(ptr1);

			fulllist_count++;

			if (fulllist == NULL)
			{
				fulllist = new_range;
				tail_fulllist = new_range;
			} 
			else
			{
				tail_fulllist->next = new_range;
				tail_fulllist = new_range;
			}

			ptr1 = ptr1->next;
		} 
		else if (ptr1->start > ptr2->end)
		{
			//no overlap and *ptr2 is ahead from *ptr1
			new_range = new range;

			new_range->range_copy(ptr2);

			fulllist_count++;

			if (fulllist == NULL)
			{
				fulllist = new_range;
				tail_fulllist = new_range;
			} 
			else
			{
				tail_fulllist->next = new_range;
				tail_fulllist = new_range;
			}

			ptr2 = ptr2->next;
		}
		else
		{
			//overlap
			new_range = new range;

			new_range->start = ptr1->start < ptr2->start? ptr1->start : ptr2->start;
			new_range->end = ptr1->end > ptr2->end? ptr1->end : ptr2->end;

			fulllist_count++;

			if (fulllist == NULL)
			{
				fulllist = new_range;
				tail_fulllist = new_range;
			} 
			else
			{
				tail_fulllist->next = new_range;
				tail_fulllist = new_range;
			}

			ptr1 = ptr1->next;
			ptr2 = ptr2->next;
		}
	}

	while (ptr1 != NULL)
	{
		new_range = new range;
		new_range->range_copy(ptr1);

		fulllist_count++;

		if (fulllist == NULL)
		{
			fulllist = new_range;
			tail_fulllist = new_range;
		} 
		else
		{
			tail_fulllist->next = new_range;
			tail_fulllist = new_range;
		}

		ptr1 = ptr1->next;
	}

	while (ptr2 != NULL)
	{
		new_range = new range;
		new_range->range_copy(ptr2);

		fulllist_count++;

		if (fulllist == NULL)
		{
			fulllist = new_range;
			tail_fulllist = new_range;
		} 
		else
		{
			tail_fulllist->next = new_range;
			tail_fulllist = new_range;
		}

		ptr2 = ptr2->next;
	}

	ptr1 = fulllist;
	ptr2 = ptr1->next;

	while (ptr2 != NULL && ptr1 != NULL)
	{
		if (ptr2->end <= ptr2->start)
		{
			if (ptr2 == tail_fulllist)
			{
				tail_fulllist = ptr1;
			}

			ptr1->next = ptr2->next;
			delete ptr2;
			ptr2 = ptr1->next;

			fulllist_count--;
		}
		else
		{
			ptr1 = ptr2;
			ptr2 = ptr1->next;
		}
	}

	return true;
}



void PE_matchset_acknowledge_initialization()
{
	//initialize the variables used in the PE_matchset_acknowledge
	//such as setting starting index to 1 (if store in variables) or initializing a ofstream object (if store in file)

	char PE_match_outputfilename[1000];
	char PE_match_region_outputfilename[1000];
	char PE_match_region_detail_outputfilename[1000];
	char PE_match_BED_outputfilename[1000];

	sprintf(PE_match_outputfilename, "%sresult/PEmatch/%s.txt", dirPrefix, targetChromosome);
	sprintf(PE_match_region_outputfilename, "%sresult/PEmatchregion/%s.txt", dirPrefix, targetChromosome);
	sprintf(PE_match_region_detail_outputfilename, "%sresult/PEmatchregion_detail/%s.txt", dirPrefix, targetChromosome);
	sprintf(PE_match_BED_outputfilename, "%sresult/PEmatch_region_BED/%s.bed", dirPrefix, targetChromosome);

	PE_match_outputfile.open(PE_match_outputfilename);
	PE_match_region_outputfile.open(PE_match_region_outputfilename);
	PE_match_region_detail_outputfile.open(PE_match_region_detail_outputfilename);
	PE_match_BED_outputfile.open(PE_match_BED_outputfilename);

	

	if (sameChromosome == false)
	{
		sprintf(PE_match_BED_outputfilename, "%sresult/PEmatch_region_BED/%s_1.bed", dirPrefix, targetChromosome);
		PE_match_BED_outputfile_1.open(PE_match_BED_outputfilename);
		//PE_match_BED_outputfile_1 << "track name=pairedReads description=\"Clone Paired Reads\" useScore=1" << endl;

		sprintf(PE_match_BED_outputfilename, "%sresult/PEmatch_region_BED/%s_2.bed", dirPrefix, targetChromosome);
		PE_match_BED_outputfile_2.open(PE_match_BED_outputfilename);
		//PE_match_BED_outputfile_2 << "track name=pairedReads description=\"Clone Paired Reads\" useScore=1" << endl;
	}
	else
	{
		PE_match_BED_outputfile << "track name=pairedReads description=\"Clone Paired Reads\" useScore=1" << endl;
	}

	last_end = 0;
	current_set_index = 0;
	/*for (long tmp = 0; tmp < MAX_LENGTH; tmp++)
	{
	current_set[tmp] = 0;
	}*/


	min_point = MAX_LENGTH;
	max_point = 0;

	return;
}

/*
void PE_matchset_acknowledge(int newmember, long overlap_start, long overlap_end, cluster *c_range)
{
	//add newmember to current PE match set
	//if newmember is 0, close the current set as a complete set
	range* output_range;

	if (newmember == 0)
	{
		if (current_set_index > 1 && current_set[current_set_index][0] > last_end)
		{
			//output current set
			 
			for (long tmp = 1; tmp <= current_set_index; tmp++)
			{
				PE_match_outputfile << PE_list_ordered[current_set[tmp][0]] << "(" << current_set[tmp][1] << ", " << current_set[tmp][2] << ")\t";
			}

			PE_match_outputfile << endl;

			
			PE_match_region_outputfile << c_range->chromosome_1 << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high << '\t' << c_range->chromosome_2 << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high << '\t' << current_set_index << endl; 

			
			
			//region detail			
			PE_match_region_detail_outputfile << c_range->chromosome_1 << "\t[" << c_range->end1_range_low << ", " << c_range->end1_range_high << "]\t{";
			
			output_range = c_range->endpoint1;
			while (output_range != NULL)
			{
				PE_match_region_detail_outputfile << "(" << output_range->start << ", " << output_range->end << ") ";

				output_range = output_range->next;
			}

			PE_match_region_detail_outputfile << "}\t" << c_range->chromosome_2 << "\t[" << c_range->end2_range_low << ", " << c_range->end2_range_high << "]\t{";

			output_range = c_range->endpoint2;
			while (output_range != NULL)
			{
				PE_match_region_detail_outputfile << "(" << output_range->start << ", " << output_range->end << ") ";

				output_range = output_range->next;
			}

			PE_match_region_detail_outputfile << "}\t" << current_set_index << endl;


			last_end = current_set[current_set_index][0];
		} 

		current_set_index = 0;

		return;
	}

	current_set_index++;
	current_set[current_set_index][0] = newmember;
	current_set[current_set_index][1] = overlap_start;
	current_set[current_set_index][2] = overlap_end;

	return;
}
*/


/*void outputBlocksSize(range *p)
{
	if (p == NULL)
	{
		return;
	}

	outputBlocksSize(p->next);

	PE_match_BED_outputfile << p->end - p->start << ",";

	return;
}

void outputBlocksStart(range *p, long base)
{
	if (p == NULL)
	{
		return;
	}

	outputBlocksStart(p->next, base);

	PE_match_BED_outputfile << p->start - base << ",";

	return;
}*/

void outputFullList(range *fulllist, range *tail_fulllist, long base)
{
	//output blocks to BED file
	range *current;
	
	//output block size	
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile << current->end - current->start + 1 << ",";

		current = current->next;
	}

	PE_match_BED_outputfile << "\t";

	//output block start
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile << current->start - base << ",";

		current =  current->next;
	}

	//PE_match_BED_outputfile << current->start - base << endl;
	PE_match_BED_outputfile << endl;

	return;
}

void outputFullList_diffChrClusters(range *fulllist, range *tail_fulllist, long base)
{
	//output blocks to BED file
	range *current;

	//output block size	
	current = fulllist;
	while (current != NULL)
	{
		diffChrClusters << current->end - current->start + 1 << ",";

		current = current->next;
	}

	PE_match_BED_outputfile << "\t";

	//output block start
	current = fulllist;
	while (current != NULL)
	{
		diffChrClusters << current->start - base << ",";

		current =  current->next;
	}

	//PE_match_BED_outputfile << current->start - base << endl;
	diffChrClusters << endl;

	return;
}

void outputEndList1(range *fulllist, range *tail_fulllist, long base)
{
	//output blocks to BED file
	range *current;

	//output block size	
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile_1 << current->end - current->start + 1 << ",";

		current = current->next;
	}

	PE_match_BED_outputfile_1 << "\t";

	//output block start
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile_1 << current->start - base << ",";

		current =  current->next;
	}

	//PE_match_BED_outputfile_1 << current->start - base << endl;
	PE_match_BED_outputfile_1 << endl;

	return;
}

void outputEndList2(range *fulllist, range *tail_fulllist, long base)
{
	//output blocks to BED file
	range *current;

	//output block size	
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile_2 << current->end - current->start + 1 << ",";

		current = current->next;
	}

	PE_match_BED_outputfile_2 << "\t";

	//output block start
	current = fulllist;
	while (current != NULL)
	{
		PE_match_BED_outputfile_2 << current->start - base << ",";

		current =  current->next;
	}

	//PE_match_BED_outputfile_2 << current->start - base << endl;
	PE_match_BED_outputfile_2 << endl;

	return;
}


//output in BED file format
void PE_matchset_acknowledge(int newmember, long overlap_start, long overlap_end, cluster *c_range)
{
	//add newmember to current PE match set
	//if newmember is 0, close the current set as a complete set
	range* output_range;
	long chromStart, chromEnd;
	pe_distance *tmp_dist;

	if (newmember == 0)
	{
		if (current_set_index > 1 && current_set[current_set_index][0] > last_end)// && c_range->end2_range_low - c_range->end1_range_high < 100000)
		{
			//output current set

			for (long tmp = 1; tmp <= current_set_index; tmp++)
			{
				PE_match_outputfile << PE_list_ordered[current_set[tmp][0]] << "(" << current_set[tmp][1] << ", " << current_set[tmp][2] << ")\t";

				tmp_dist = PE_list[PE_list_ordered[current_set[tmp][0]]]->dist;
				while (tmp_dist != NULL)
				{
					cluster_distances << tmp_dist->existing_distance << '\t';
					tmp_dist = tmp_dist->next;
				}
				cluster_distances << "|\t";
			}

			PE_match_outputfile << endl;
			cluster_distances << endl;


			PE_match_region_outputfile << c_range->chromosome_1 << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high << '\t' << c_range->chromosome_2 << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high << '\t' << current_set_index << endl; 
	
			if (strcmp(c_range->chromosome_1, c_range->chromosome_2) != 0 || (c_range->end1_range_high < c_range->end2_range_low || c_range->end1_range_low > c_range->end2_range_high))
			{
				finalResult << c_range->chromosome_1 << "\t" << c_range->strand_start << "\t" << c_range->end1_range_low << "\t" << c_range->end1_range_high << "\t" << c_range->chromosome_2 << "\t";
				if (strcmp(c_range->strand_end, "0") == 0)
					finalResult << "16";
				else
					finalResult << "0";
				finalResult << "\t" << c_range->end2_range_low << "\t" << c_range->end2_range_high << endl;
			}

			//region detail			
			PE_match_region_detail_outputfile << c_range->chromosome_1 << "\t[" << c_range->end1_range_low << ", " << c_range->end1_range_high << "]\t{";

			output_range = c_range->endpoint1;
			while (output_range != NULL)
			{
				PE_match_region_detail_outputfile << "(" << output_range->start << ", " << output_range->end << ") ";

				output_range = output_range->next;
			}

			PE_match_region_detail_outputfile << "}\t" << c_range->chromosome_2 << "\t[" << c_range->end2_range_low << ", " << c_range->end2_range_high << "]\t{";

			output_range = c_range->endpoint2;
			while (output_range != NULL)
			{
				PE_match_region_detail_outputfile << "(" << output_range->start << ", " << output_range->end << ") ";

				output_range = output_range->next;
			}

			PE_match_region_detail_outputfile << "}\t" << current_set_index << endl;


			//output BED file
			if (sameChromosome == true)
			{
				chromStart = ((c_range->end1_range_low < c_range->end2_range_low)? c_range->end1_range_low : c_range->end2_range_low);
				chromEnd = ((c_range->end1_range_high > c_range->end2_range_high)? c_range->end1_range_high : c_range->end2_range_high);

				PE_match_BED_outputfile << c_range->chromosome_1 << '\t' << chromStart << '\t' << chromEnd + 1 << '\t';
				PE_match_BED_outputfile << "clone" << '\t';

				if (current_set_index + SCORE_ADD < 1000)
				{
					PE_match_BED_outputfile << current_set_index + SCORE_ADD << "\t";
				}
				else
				{
					PE_match_BED_outputfile << "999\t";
				} 


				PE_match_BED_outputfile << '-' << '\t' << chromStart << '\t' << chromEnd << '\t';

				PE_match_BED_outputfile << "0\t";

				if (c_range->build_fulllist() == true)
				{
					PE_match_BED_outputfile << c_range->fulllist_count << '\t';
					outputFullList(c_range->fulllist, c_range->tail_fulllist, chromStart);
				}
			} 
			else
			{
				//on different chromosomes. translocation
				PE_match_BED_outputfile << c_range->chromosome_1 << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
				PE_match_BED_outputfile << "clone" << '\t';
				if (current_set_index + SCORE_ADD < 1000)
				{
					PE_match_BED_outputfile << current_set_index + SCORE_ADD << "\t";
				}
				else
				{
					PE_match_BED_outputfile << "999\t";
				} 
				PE_match_BED_outputfile << '-' << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
				PE_match_BED_outputfile << "0\t";
				PE_match_BED_outputfile << c_range->range_count_1 << '\t';
				outputFullList(c_range->endpoint1, c_range->tail_end1, c_range->end1_range_low);

				PE_match_BED_outputfile << c_range->chromosome_2 << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << '\t';
				PE_match_BED_outputfile << "clone" << '\t';
				if (current_set_index + SCORE_ADD < 1000)
				{
					PE_match_BED_outputfile << current_set_index + SCORE_ADD << "\t";
				}
				else
				{
					PE_match_BED_outputfile << "999\t";
				} 
				PE_match_BED_outputfile << '-' << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << '\t';
				PE_match_BED_outputfile << "0\t";
				PE_match_BED_outputfile << c_range->range_count_2 << '\t';
				outputFullList(c_range->endpoint2, c_range->tail_end2, c_range->end2_range_low);

				PE_match_BED_outputfile_1 << c_range->chromosome_1 << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
				PE_match_BED_outputfile_1 << "clone" << '\t';
				if (current_set_index + SCORE_ADD < 1000)
				{
					PE_match_BED_outputfile_1 << current_set_index + SCORE_ADD << "\t";
				}
				else
				{
					PE_match_BED_outputfile_1 << "999\t";
				} 
				PE_match_BED_outputfile_1 << '-' << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
				PE_match_BED_outputfile_1 << "0\t";
				PE_match_BED_outputfile_1 << c_range->range_count_1 << '\t';
				outputEndList1(c_range->endpoint1, c_range->tail_end1, c_range->end1_range_low);

				PE_match_BED_outputfile_2 << c_range->chromosome_2 << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << '\t';
				PE_match_BED_outputfile_2 << "clone" << '\t';
				if (current_set_index + SCORE_ADD < 1000)
				{
					PE_match_BED_outputfile_2 << current_set_index + SCORE_ADD << "\t";
				}
				else
				{
					PE_match_BED_outputfile_2 << "999\t";
				} 
				PE_match_BED_outputfile_2 << '-' << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << '\t';
				PE_match_BED_outputfile_2 << "0\t";
				PE_match_BED_outputfile_2 << c_range->range_count_2 << '\t';
				outputEndList2(c_range->endpoint2, c_range->tail_end2, c_range->end2_range_low);

				diffChrClusters << c_range->chromosome_1 << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
				diffChrClusters << current_set_index + SCORE_ADD << "\t";
// 				diffChrClusters << '-' << '\t' << c_range->end1_range_low << '\t' << c_range->end1_range_high + 1 << '\t';
// 				diffChrClusters << "0\t";
// 				diffChrClusters << c_range->range_count_1 << '\t';
// 				outputFullList_diffChrClusters(c_range->endpoint1, c_range->tail_end1, c_range->end1_range_low);

				diffChrClusters << c_range->chromosome_2 << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << endl; //<< '\t';
//				diffChrClusters << current_set_index + SCORE_ADD << "\t";
// 				diffChrClusters << '-' << '\t' << c_range->end2_range_low << '\t' << c_range->end2_range_high + 1 << '\t';
// 				diffChrClusters << "0\t";
// 				diffChrClusters << c_range->range_count_2 << '\t';
// 				outputFullList_diffChrClusters(c_range->endpoint2, c_range->tail_end2, c_range->end2_range_low);
			}
			
			

			last_end = current_set[current_set_index][0];
		} 

		current_set_index = 0;

		return;
	}

	current_set_index++;
	current_set[current_set_index][0] = newmember;
	current_set[current_set_index][1] = overlap_start;
	current_set[current_set_index][2] = overlap_end;

	return;
}





void PE_matchset_acknowledge_end()
{
	//clean up the useless variables in the process of set acknowledge, such as close files

	PE_match_outputfile.close();

	PE_match_region_outputfile.close();

	PE_match_region_detail_outputfile.close();

	PE_match_BED_outputfile.close();

	if (sameChromosome == false)
	{
		PE_match_BED_outputfile_1.close();
		PE_match_BED_outputfile_2.close();
	}

	return;
}



void PE_match()
{
	//find all matched paired end reads for every paired end read
	//brute force
	//clustering

	long i, j, overlap_start, overlap_end;
	

	PE_matchset_acknowledge_initialization();

	for (i = 1; i < PE_num; i++)
	{
		//PE_matchset_acknowledge(PE_list_ordered[i], 50, 50, NULL);
		PE_matchset_acknowledge(i, 50, 50, NULL);

		j = i + 1;

		cluster c_range(PE_list[PE_list_ordered[i]]);
		
		while (j <= PE_num && PE_list[PE_list_ordered[j]]->startrange_low <= c_range.end1_range_high)// && PE_list[PE_list_ordered[j]]->endrange_low <= c_range.end2_range_high)
		{
			if (c_range.searchANDmodify(PE_list[PE_list_ordered[j]], &overlap_start, &overlap_end) == true)
			{
				//PE_matchset_acknowledge(PE_list_ordered[j], overlap_start, overlap_end, NULL);
				PE_matchset_acknowledge(j, overlap_start, overlap_end, NULL);
			}
			
			j++;
		}

		PE_matchset_acknowledge(0, 0, 0, &c_range); //end the current set		
	}

	PE_matchset_acknowledge_end();

	return;
}





void input_from_file()
{
	//input PE reads from input file
	ifstream inputfile;
	inputfile.open(inputfilename);

	char name[100], chromosome1[30], chromosome2[30], end1[50], end2[50], strand1[5], strand2[5];
	long start1, start2, PE_list_index;
	paired_end* new_PE;
	range* new_range;

	int tmp, i, cccccount = 1;

	PE_list_index = 0;

	for (int tmp = 0; tmp < 100; tmp++)
	{
		name[tmp] = '\0';
	}

	inputfile >> name;

	while (name[1] != '\0')
	{
		inputfile >> chromosome1;
		inputfile >> strand1;
		inputfile >> start1;
		inputfile >> end1;

		inputfile >> chromosome2;
		inputfile >> strand2;
		inputfile >> start2;
		inputfile >> end2;

		new_PE = new paired_end;

		strcpy(new_PE->chromosome_start, chromosome1);
		strcpy(new_PE->strand_start, strand1);

		tmp = 0;
		for (i = 0; end1[i] != '\0'; i++)
		{
			if (end1[i] == 'M')
			{
				if (tmp == 0)
				{
					continue;
				}

				new_range = new range;
				new_range->start = start1;
				new_range->end = start1 + tmp - 1;

				if (new_PE->tail_start == NULL)
				{
					new_PE->start = new_range;
					new_PE->tail_start = new_range;
				} 
				else
				{
					new_PE->tail_start->next = new_range;
					new_PE->tail_start = new_range;
				}

				if (start1 < new_PE->startrange_low)
				{
					new_PE->startrange_low = start1;
				}

				if (start1 + tmp > new_PE->startrange_high)
				{
					new_PE->startrange_high = start1 + tmp - 1;
				}

				(new_PE->range_count_start)++;

				start1 = start1 + tmp;
				tmp = 0;
			} 
			else if (end1[i] == 'N')
			{
				start1 = start1 + tmp;
				tmp = 0;
			}
			else
			{
				tmp = tmp * 10 + end1[i] - 48;
			}
		}

		strcpy(new_PE->chromosome_end, chromosome2);
		strcpy(new_PE->strand_end, strand2);

		tmp = 0;
		for (i = 0; end2[i] != '\0'; i++)
		{
			if (end2[i] == 'M')
			{
				new_range = new range;
				new_range->start = start2;
				new_range->end = start2 + tmp - 1;

				if (new_PE->tail_end == NULL)
				{
					new_PE->end = new_range;
					new_PE->tail_end = new_range;
				} 
				else
				{
					new_PE->tail_end->next = new_range;
					new_PE->tail_end = new_range;
				}

				if (start2 < new_PE->endrange_low)
				{
					new_PE->endrange_low = start2;
				}

				if (start2 + tmp > new_PE->endrange_high)
				{
					new_PE->endrange_high = start2 + tmp - 1;
				}

				(new_PE->range_count_end)++;

				start2 = start2 + tmp;
				tmp = 0;
			} 
			else if (end2[i] == 'N')
			{
				start2 = start2 + tmp;
				tmp = 0;
			}
			else
			{
				tmp = tmp * 10 + end2[i] - 48;
			}
		}

		if (strcmp(new_PE->chromosome_start, new_PE->chromosome_end) == 0)
		{
			new_PE->dist = new pe_distance(new_PE->startrange_high, new_PE->endrange_low);
		} 
		else
		{
			new_PE->dist = new pe_distance(1, MAX);
		}
		
		//new_PE->pe_distance = start1 - start2;
		//cout << cccccount++ << endl;

// 	//	if (new_PE->endrange_low - new_PE->startrange_high > MAX_PEREAD_DISTANCE)
// 		if (new_PE->endrange_low - new_PE->startrange_high < MAX_PEREAD_DISTANCE)
// 		{
// 			//ignore
// 		} 
/*		else*/
		{
			//new_PE->compute_distance();

			PE_list_index++;
			PE_list[PE_list_index] = new_PE;

		}
		
		for (int tmp = 0; tmp < 100; tmp++)
		{
			name[tmp] = '\0';
		}

		inputfile >> name;		
	}

	PE_num = PE_list_index;

	inputfile.close();

	return;
}


/*

long quicksort_partition(long p, long r)
{
	long x, i, j, tmp;
	
	x = sort_key[r];
	i = p - 1;

	for (j = p; j <= r - 1; j++)
	{
		if (sort_key[j] <= x)
		{
			i++;

			tmp = sort_key[j];
			sort_key[j] = sort_key[i];
			sort_key[i] = tmp;

			tmp = PE_list_ordered[j];
			PE_list_ordered[j] = PE_list_ordered[i];
			PE_list_ordered[i] = tmp;
		}
	}


	tmp = sort_key[r];
	sort_key[r] = sort_key[i+1];
	sort_key[i+1] = tmp;

	tmp = PE_list_ordered[r];
	PE_list_ordered[r] = PE_list_ordered[i+1];
	PE_list_ordered[i+1] = tmp;

	return i+1;
}

long quicksort_randomized_partition(long p, long r)
{
	long i, tmp;
	i = p + (double)rand()/ (RAND_MAX) * (r - p);


	tmp = sort_key[r];
	sort_key[r] = sort_key[i];
	sort_key[i] = tmp;

	tmp = PE_list_ordered[r];
	PE_list_ordered[r] = PE_list_ordered[i];
	PE_list_ordered[i] = tmp;

	return quicksort_partition(p, r);
}

void quicksort(long p, long r)
{
	long q;

	if (p < r)
	{
		q = quicksort_randomized_partition(p, r);

		quicksort(p, q-1);
		quicksort(q+1, r);
	}

	return;
}

*/



long partition(long p,long r)
{
	long x, i, j, tmp;

	//randomized partition
	i = p + (double)rand()/ (RAND_MAX) * (r - p);

	tmp = sort_key[r];
	sort_key[r] = sort_key[i];
	sort_key[i] = tmp;

	tmp = PE_list_ordered[r];
	PE_list_ordered[r] = PE_list_ordered[i];
	PE_list_ordered[i] = tmp;


	x = sort_key[r];
	i = p - 1;
	
	for(j = p; j < r; j++)
	{
		if (sort_key[j] <= x)
		{
			i++;

			tmp = sort_key[j];
			sort_key[j] = sort_key[i];
			sort_key[i] = tmp;

			tmp = PE_list_ordered[j];
			PE_list_ordered[j] = PE_list_ordered[i];
			PE_list_ordered[i] = tmp;
		}
	}

	tmp = sort_key[r];
	sort_key[r] = sort_key[i+1];
	sort_key[i+1]=tmp;

	tmp = PE_list_ordered[r];
	PE_list_ordered[r] = PE_list_ordered[i+1];
	PE_list_ordered[i+1] = tmp;

	return i+1;
}


void quicksort(long length)
{
	long top = 0, p, r, q;

	stack[top++] = 1;
	stack[top++] = length;
	
	while (top != 0)
	{
		r = stack[--top];
		p = stack[--top];
		
		if(p>=r)
			continue;
		
		q = partition(p, r);
		
		stack[top++] = p;
		stack[top++] = q - 1;
		
		stack[top++] = q + 1;
		stack[top++] = r;
	}

	return;
}


void sort_PE_list()
{
	for (int i = 1; i <= PE_num; i++)
	{
		PE_list_ordered[i] = i;
	}


	//sort based on endrange_high
	for (int i = 1; i <= PE_num; i++)
	{
		sort_key[i] = PE_list[PE_list_ordered[i]]->endrange_high;
	}
	quicksort(PE_num);

	//sort based on startrange_high
	for (int i = 1; i <= PE_num; i++)
	{
		sort_key[i] = PE_list[PE_list_ordered[i]]->startrange_high;
	}
	quicksort(PE_num);

	//sort based on endrange_low
	for (int i = 1; i <= PE_num; i++)
	{
		sort_key[i] = PE_list[PE_list_ordered[i]]->endrange_low;
	}
	quicksort(PE_num);

	//sort based on startrange_low
	for (int i = 1; i <= PE_num; i++)
	{
		sort_key[i] = PE_list[PE_list_ordered[i]]->startrange_low;
	}
	quicksort(PE_num);


/*	//sort based on end_chr
	for (int i = 1; i <= PE_num; i++)
	{
		sort_key[i] = ((PE_list[PE_list_ordered[i]]->chromosome_end)[3] - 48) * 10;

		if ((PE_list[PE_list_ordered[i]]->chromosome_end)[4] != '\0')
		{
			sort_key[i] += (PE_list[PE_list_ordered[i]]->chromosome_end)[4] - 48;
		}
	}
	quicksort(PE_num);
*/


	return;
}


void input_from_file_junction(char *targetChr)
{
	ifstream inputfile;
	inputfile.open(inputfilename);

	junction_num = 0;

	long start, end;
	int temp;
	string otherinfo;
	junction *newjunc;
	char tmpchar[30];

	for (int tmp = 0; tmp < 100; tmp++)
	{
		chromosome[tmp] = '\0';
	}

	getline(inputfile, junction_title);

	inputfile >> chromosome;

	while (chromosome[1] != '\0')
	{
		if (strcmp(chromosome, targetChr) == 0)
		{
			newjunc = new junction;

			strcpy(newjunc->chromosome, chromosome);
			inputfile >> newjunc->start;
			inputfile >> newjunc->end;
			inputfile >> newjunc->name;
			inputfile >> temp;
			inputfile >> newjunc->strand;
			inputfile >> newjunc->thickStart;
			inputfile >> newjunc->thickEnd;
			inputfile >> tmpchar;
			inputfile >> newjunc->blockCount;
			inputfile >> newjunc->blockSizes;
			inputfile >> newjunc->blockStarts;
			getline(inputfile, otherinfo);
			//newjunc->info = otherinfo;
			newjunc->compute_range();

			junction_num++;
			junction_list[junction_num] = newjunc;
		} 
		else
		{
			getline(inputfile, otherinfo);
		}


// 		for (int tmp = 0; tmp < 100; tmp++)
// 		{
// 			chromosome[tmp] = '\0';
// 		}
		chromosome[1] = '\0';

		inputfile >> chromosome;		

	}


	inputfile.close();

	return;
}

void cleanup()
{
	paired_end *ptr;
	

	for (int i = 1; i <= PE_num; i++)
	{
		ptr = PE_list[i];
		delete ptr;
	}

	PE_num = 0;


	return;
}


void processSingle()
{
	int i, j;
	char outputfilename[1000], chromosome1[50], chromosome2[50];
	
	chromosome1[0] = 'c';

	for (i = 1; targetChromosome[i] != 'c'; i++)
	{
		chromosome1[i] = targetChromosome[i];
	}
	chromosome1[i] = '\0';

	for (j = 0; targetChromosome[i] != '\0'; i++, j++)
	{
		chromosome2[j] = targetChromosome[i];
	}
	chromosome2[j] = '\0';


	//sprintf(inputfilename, "C:\\Users\\yin\\Desktop\\PEmatch\\junctions_filtered.txt");
	//sprintf(inputfilename, "C:\\Users\\yin\\Desktop\\PEmatch\\human35\\junction\\junctions_rmapped.txt");
	if (strcmp(chromosome1, chromosome2) == 0)
	{
		//sprintf(inputfilename, "C:\\Users\\yin\\Desktop\\PEmatch\\human35\\junction\\junctions.txt");
		//input_from_file_junction(chromosome1);
		junction_num = 0;
		sameChromosome = true;
	}
	else
	{
		//return;

		junction_num = 0;
		sameChromosome = false;
	}
		
	sprintf(outputfilename, "%sresult/distance/distances_%s.txt", dirPrefix, targetChromosome);
	cluster_distances.open(outputfilename);

	sprintf(inputfilename, "%sdata/parsedPER/%s.txt", dirPrefix, targetChromosome);
	input_from_file();

	sort_PE_list();

	PE_match();

	cluster_distances.close();

// 	junction *ptr_junc;
// 
// 	if (strcmp(chromosome1, chromosome2) == 0)
// 	{
// 		ofstream output_junction_support;
// 		sprintf(outputfilename, "C:\\Users\\yin\\Desktop\\PEmatch\\human35\\junction_support\\JunctionSupport_%s.bed", targetChromosome);
// 		output_junction_support.open(outputfilename);
// 		output_junction_support << junction_title << endl;
// 		for (i = 1; i <= junction_num; i++)
// 		{
// 			output_junction_support << junction_list[i]->chromosome << '\t' <<junction_list[i]->start_start << '\t' << junction_list[i]->end_end << '\t' << junction_list[i]->name << '\t';
// 			if (junction_support[i] > 500)
// 			{
// 				output_junction_support << "999";
// 			} 
// 			else if (junction_support[i] == 0)
// 			{
// 				output_junction_support << "167";
// 			}
// 			else
// 			{
// 				output_junction_support << junction_support[i] + 500;
// 			}
// 			output_junction_support << '\t' << junction_list[i]->strand << '\t' << junction_list[i]->start_start << '\t' << junction_list[i]->end_end << '\t' << '0' << '\t' << junction_list[i]->blockCount << '\t' << junction_list[i]->blockSizes << '\t' << "0," << junction_list[i]->end - junction_list[i]->start_start << endl;
// 
// 
// 			junction_support[i] = 0;
// 
// 			ptr_junc = junction_list[i];
// 			delete ptr_junc;
// 		}
// 		output_junction_support.close();
// 	}

	cleanup();

	return;
}





int main(int argc, char* argv[])
{
	char filename[1000];
	deviation_end = 30;
	deviation_distance = MAX;

	if (argc != 2)
	{
		cout << argv[0] << "\t<targetpath>" << endl;
		return 1;
	}

	strcpy(dirPrefix, argv[1]);
	ifstream chrNameFile;
	sprintf(filename, "%sdata/parsedPER/ChromosomeName.txt", dirPrefix);
	chrNameFile.open(filename);

//	diffChrClusters.open("%sresult/cluster_DiffChr.txt", dirPrefix);
	sprintf(filename, "%sresult/cluster.txt", dirPrefix);
	finalResult.open(filename);

	int tmp;

	for (tmp = 0; tmp < 100; tmp++)
	{
		targetChromosome[tmp] = '\0';
	}

	chrNameFile >> targetChromosome;


	while (targetChromosome[1] != '\0')
	{
		processSingle();

		for (tmp = 0; tmp < 100; tmp++)
		{
			targetChromosome[tmp] = '\0';
		}

		chrNameFile >> targetChromosome;	
	}


	chrNameFile.close();
//	diffChrClusters.close();
	finalResult.close();

	return 0;
}
