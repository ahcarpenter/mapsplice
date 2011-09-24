/*    
 *    FilterReadsByCanonNoncanon.cpp		
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

	char* filtered_canon = argv[3];

	char* filtered_noncanon = argv[4];

	char* fitlered_noncanon_canon = argv[5];

	//size_t maxhits = atoi(argv[4]);

	//char* unique_mapped = argv[5];

	SamHandler mps_samhandler;

	cout << "read sam " << mps_sam << endl;

	mps_samhandler.ReadAppendSamFile(mps_sam);

	mps_samhandler.CollectStat();

	cout << "sort sam " << mps_sam << endl;

	mps_samhandler.SortSamByTagnameOffset();

	map<string, vector<JuncForSort> > junc_map;

	cout << "ReadSortedJunction"<<endl;

	mps_samhandler.ReadSortedJunction(junc_file, junc_map);

	cout << "mark canon noncanon reads mapping"<<endl;

	mps_samhandler.MarkCanonNoncanon(junc_map);

	cout << "filter canon noncanon " << mps_sam << endl;

	mps_samhandler.FilterCanonNoncaonTag(mps_sam, filtered_canon, filtered_noncanon, fitlered_noncanon_canon);
	
}