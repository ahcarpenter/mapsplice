/*    
 *    remove_blankspace_perline.cpp		
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
using namespace std;

string
basename2(string filename) {

	//cout << "bef: "<<filename<<endl;
	const string s(filename.substr(0, filename.find_last_of(".")));
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
	{
		//cout << "aft 1: "<<s.substr(final_slash + 1)<<endl;
		return s.substr(final_slash + 1);
	}
	else
	{
		//cout << "aft 2: "<<s<<endl;
		return s;
	}
}


void convert_ex_fasta(char* extendFile, char* reads_file)
{
	string chromfile = reads_file;

	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	ofstream output_fs(chromfile.c_str());

	string line;
	int count=0;
	while(!input_fs.eof())
	{
		getline(input_fs, line);
		if (line.empty())
			continue;

		if ((line[0] == '@' || line[0] == '>' || line[0] == '+') && (line.find(" ") != string::npos || line.find("\t") != string::npos))
		{
			for (size_t i = 0; i < line.length(); ++i)
			{
				if (line[i] == ' ' || line[i] == '\t')
					line[i] = '_';
			}
		}
		count++;

		output_fs << line << endl;

	}
	input_fs.close();

}

int main(int argc, char* argv[])
{
	if (argc != 3 )
	{
		fprintf(stderr,"./remove_blankspace_perline reads_file formated_reads_file\n");exit(1);
	}

	convert_ex_fasta(argv[1], argv[2]);

	return 0;
}