/*    
 *    SamHandler.cpp	
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
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


#include "SamHandler.h"


/////merge sort sam
//int
//main(int argc, char** argv)
//{
//	SamHandler allsam;
//
//	for (size_t i = 1; i < argc - 2; ++i)
//	{
//		allsam.ReadAppendSamFile(argv[i]);
//	}
//
//	allsam.WriteSamByTagnameOffset(argv[argc - 2]);
//}


//int
//main(int argc, char** argv)
//{
//	char* mps_sam = argv[1];
//
//	char* synthetic_sam = argv[2];
//
//	char* compared_file = argv[3];
//
//	size_t range = atoi(argv[4]);
//
//	char* junc_file = argv[5];
//
//	char* filtered_multiple = argv[6];
//
//	size_t maxhits = atoi(argv[7]);
//
//	size_t read_width = atoi(argv[8]);
//
//	char* chrom_dir = argv[9];
//
//	size_t min_intron = atoi(argv[10]);
//
//	size_t max_intron = atoi(argv[11]);
//
//	SamHandler mps_samhandler, synthetic_samhandler;
//
//	cout << "read sam " << mps_sam << endl;
//
//	mps_samhandler.ReadAppendSamFile(mps_sam);
//
//	cout << "sort sam " << mps_sam << endl;
//
//	mps_samhandler.SortSamByTagnameOffset();
//
//	cout << "sam2junc " << mps_sam << endl;
//
//	mps_samhandler.Sam2Junc(junc_file, read_width, chrom_dir, min_intron, max_intron);
//
//	cout << "filter multiple mapped " << mps_sam << endl;
//
//	mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits);
//
//	cout << "read sam " << synthetic_sam << endl;
//
//	synthetic_samhandler.ReadAppendSamFile(synthetic_sam);
//
//	cout << "sort sam " << synthetic_sam << endl;
//
//	synthetic_samhandler.SortSamByTagnameOffset();
//
//	cout << "compare sam " << mps_sam <<'\t' <<synthetic_sam << endl;
//
//	mps_samhandler.CompareToSyntheticMappedReads(synthetic_samhandler, compared_file, range);
//}

