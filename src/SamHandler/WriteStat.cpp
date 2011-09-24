/*    
 *    WriteStat.cpp		
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

	int unmapped = 0;

	for (int i = 1; i < argc - 1; ++i)
	{
		cout << argv[i] << endl;
		unmapped += allsam.ReadAppendSamFileSimp(argv[i]);
	}

	allsam.CollectStat();

	allsam.WriteStat(argv[argc - 1], unmapped);
}