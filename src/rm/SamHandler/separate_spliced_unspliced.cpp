/*    
 *    separate_spliced_unspliced.cpp	
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
	char* sam_file = argv[1];

	char* unspliced = argv[2];

	char* spliced = argv[3];

	SamHandler sam_handler;

	cout << sam_file << endl;

	cout << "number of records before filtered: "<< sam_handler.ReadAppendSamFile(sam_file)<<endl;


	
	pair<size_t, size_t> count_unspliced_spliced = sam_handler.SepSplicedUnspliced(spliced, unspliced);

	cout << "number of unspliced records: "<<count_unspliced_spliced.first<<endl;

	cout << "number of spliced records: "<<count_unspliced_spliced.second<<endl;
}