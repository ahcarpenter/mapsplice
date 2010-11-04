/*    
 *    generate_bash_file.cpp		
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
#include <set>
#include <algorithm>
#include <string.h>
using namespace std;

int main(int argc, char* argv[])
{
	string director_file = argv[1];

	string bash_file = argv[2];

	string abs_path = argv[3];

	ifstream ifs(director_file.c_str());

	ofstream ofs(bash_file.c_str());

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line.empty())
				continue;

			char relative_path[10000], folder_chr[1000];

			sscanf(line.c_str(), "%s\t%s", relative_path, folder_chr);

			//python /homes/kai/code/bin/MapSplice_1.5.0/MapSplice_1.13/mapsplice_segments.py --max-hits 10 -o /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/output 
			//-X 4 -n 10 -w 50 -S fa -Q fq -c /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/ 
			//-u /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/chr2~99735158~99820411~chr2~99989586~99989682.fq 
			//-B /sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/index -L 25 -I 1 -X 4 
			//--full-running --not-rem-temp --fusion --not-rerun-all 2>/sdb1/mapsplice/result/human/2x100bp/GSM_s3_s6_s8_1.13_II/tmp/temp/1~10000/chr2~99735158~99820411~chr2~99989586~99989682/log

			ofs << "python /homes/kai/code/bin/MapSplice_1.5.0/MapSplice_1.13/mapsplice_segments.py --max-hits 10 -o "<<abs_path<<relative_path<<"/output "
				<<"-X 4 -n 10 -w 50 -S fa -Q fq -c "<<abs_path<<relative_path<<" -u "<<abs_path<<relative_path<<'/'<<folder_chr<<".fq "
				<<" -B "<<abs_path<<relative_path<<"/index "<<" -L 25 --not-rem-temp 2>"
				<<abs_path<<relative_path<<".log"<<endl;
		}
	}
	else
		cout <<"can't open file: " << director_file << endl;

}