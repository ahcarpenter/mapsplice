/*    
 *    filterbyjuncminmis.cpp		
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

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
using namespace std;

int main(int argc, char* argv[])
{
	ifstream ifs(argv[1]);

	cout <<argv[1]<<endl;

	unsigned short min_mismatch = (unsigned short) atoi(argv[2]);

	string canon_junc = argv[3]; ofstream ofs_canon_junc(canon_junc.c_str());

	string noncanon_junc = argv[4]; ofstream ofs_noncanon_junc(noncanon_junc.c_str());

	size_t canon_count = 0, non_canon_count = 0;

	if (ifs.is_open())
	{
		string firstline;
		getline(ifs, firstline);

		ofs_canon_junc <<firstline<<endl;
		ofs_noncanon_junc <<firstline<<endl;

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], lpq[100], il[100], flankchr[100];
			int juncst, juncend, prefixend, suffixst, kinds, hits, flankcase;

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, &flankcase, flankchr, il, lpq, &min_mis, &max_mis, &ave_mis);

			if (min_mis <= min_mismatch)
			{
				++canon_count;
				ofs_canon_junc << line << endl;
			}
			else
			{
				++non_canon_count;
				ofs_noncanon_junc << line << endl;
			}
		}
	}

	cout <<"in min mismatch junction: "<<canon_count<<endl;

	cout <<"not in min mismatch junction: "<<non_canon_count<<endl;
}