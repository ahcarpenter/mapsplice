/*    
 *    merge_sort_sam.cpp		
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

///merge sort sam
int
main(int argc, char** argv)
{
	SamHandler allsam;

	for (int i = 1; i < argc - 3; ++i)
	{
		cout << argv[i] << endl;
		allsam.ReadAppendSamFile(argv[i]);
	}

	allsam.CollectStat();

	cout << "write WriteSamByTagnameOffset "<<argv[argc - 3]<<endl;

	allsam.WriteSamByTagnameOffset(argv[argc - 3]);

	cout << "write WriteSamByChromoOffset "<<argv[argc - 2]<<endl;

	allsam.WriteSamByChromoOffset(argv[argc - 2]);

	allsam.WriteStat(argv[argc - 1]);
}