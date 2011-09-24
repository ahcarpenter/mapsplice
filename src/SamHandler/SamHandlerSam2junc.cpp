/*    
 *    SamHandlerSam2junc.cpp	
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

///sam 2 junc
int
main(int argc, char** argv)
{
	char* junc_file = argv[1];

	size_t read_width = atoi(argv[2]);

	char* chrom_dir = argv[3];

	size_t min_intron = atoi(argv[4]);

	size_t max_intron = atoi(argv[5]);

	SamHandler allsam;

	for (int i = 6; i < argc - 1; ++i)
	{
		allsam.ReadAppendSamFile(argv[i]);
	}

	allsam.Sam2Junc(junc_file, read_width, chrom_dir, min_intron, max_intron);
}