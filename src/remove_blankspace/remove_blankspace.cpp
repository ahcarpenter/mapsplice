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

void convert_ex_fasta(char* extendFile, string output_path)
{
	output_path.append("/");

	string chromfile = output_path;

	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string firstline;

	getline(input_fs, firstline);

	if (firstline.find(" ") != string::npos || firstline.find("\t") != string::npos)
	{
		for (size_t i = 0; i < firstline.length(); ++i)
		{
			if (firstline[i] == ' ' || firstline[i] == '\t')
				firstline[i] = '_';
		}
	}

	chromfile.append(firstline.c_str() + 1);

	chromfile.append(".fa");

	ofstream output_fs(chromfile.c_str());

	output_fs << firstline<<endl;

	string line;
	int count=0;
	while(!input_fs.eof())
	{
		getline(input_fs, line);
		if (line.empty())
			continue;
		count++;

		output_fs << line << endl;

	}
	input_fs.close();

}

int main(int argc, char* argv[])
{
	if (argc < 3 )
	{
		fprintf(stderr,"./remove_blankspace chromosome_file1 chromosome_file2 ... output_dir\n");exit(1);
	}

	for (int i = 1; i < argc - 1; ++i)
		convert_ex_fasta(argv[i], argv[argc-1]);

	return 0;
}