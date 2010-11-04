/*    
 *    countline.cpp		
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

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int
main(int argc, const char** argv)
{
	if (argc < 2)
	{
		cout << "need input files"<<endl;
		return 0;
	}

	ifstream ifs(argv[1]);

	cout <<argv[1]<<endl;

	size_t count = 0;
	if (ifs.is_open())
	{
		
		while (true)
		{
			if (ifs.eof())
				break;
			string line;
			getline(ifs,line);
			++count;
		}
	}

	cout << count << endl;
	return 0;
}