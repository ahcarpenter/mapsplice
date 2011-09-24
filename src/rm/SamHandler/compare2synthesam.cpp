/*    
 *    compare2synthesam.cpp		
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

///compare to synthetic sam
int
main(int argc, char** argv)
{
	char* mps_sam = argv[1];

	char* synthetic_sam = argv[2];

	char* compared_file = argv[3];

	size_t range = atoi(argv[4]);

	char* junc_file = argv[5];

	char* filtered_multiple = argv[6];

	size_t maxhits = atoi(argv[7]);

	char* unique_mapped = argv[8];

	char* compared_offset = argv[9];

	char* stat_file = argv[10];

	int anchor_width = atoi(argv[11]);

	

	SamHandler mps_samhandler, synthetic_samhandler;

	cout << "read sam " << mps_sam << endl;

	mps_samhandler.ReadAppendSamFile(mps_sam);

	mps_samhandler.CollectStat();

	cout << "sort sam " << mps_sam << endl;

	mps_samhandler.SortSamByTagnameOffset();

	cout << "filter multiple mapped " << mps_sam << endl;

	if (argc >= 15)
	{
		mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits, atof(argv[12]), atof(argv[13]), atof(argv[14]));
	}
	else
	{
		mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits);
	}

	cout << "write unique mapped " << mps_sam << endl;

	mps_samhandler.WriteUniqueMappedTagSam(unique_mapped);

	cout << "read sam " << synthetic_sam << endl;

	synthetic_samhandler.ReadAppendSamFile(synthetic_sam);

	synthetic_samhandler.CollectStat();

	cout << "sort sam " << synthetic_sam << endl;

	synthetic_samhandler.SortSamByTagnameOffset();

	cout << "compare sam " << mps_sam <<'\t' <<synthetic_sam << endl;

	mps_samhandler.CompareToSyntheticMappedReads(synthetic_samhandler, compared_file, compared_offset, anchor_width, range);

	mps_samhandler.WriteCompStat(synthetic_samhandler, stat_file);
}