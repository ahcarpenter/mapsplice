/*    
 *    FilterMultipleMapped.cpp		
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

	char* junc_file = argv[2];

	char* filtered_multiple = argv[3];

	size_t maxhits = atoi(argv[4]);

	char* unique_mapped = argv[5];

	SamHandler mps_samhandler;

	cout << "read sam " << mps_sam << endl;

	mps_samhandler.ReadAppendSamFile(mps_sam);

	mps_samhandler.CollectStat();

	cout << "sort sam " << mps_sam << endl;

	mps_samhandler.SortSamByTagnameOffset();

	cout << "filter multiple mapped " << mps_sam << endl;

	if (argc >= 9)
	{
		mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits, atof(argv[6]), atof(argv[7]), atof(argv[8]));
	}
	else
	{
		mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits);
	}

	//mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits);

	cout << "write unique mapped " << mps_sam << endl;

	mps_samhandler.WriteUniqueMappedTagSam(unique_mapped);

	
}