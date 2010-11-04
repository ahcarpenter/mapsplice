/*    
 *    PERfusion_final.h		
 *    MapPER
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
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




#include <cstdio>
#include <cstring>
#include <cstdlib> 
#include <iostream>
#include <unistd.h>
#include <fstream>

using namespace std;

char chrName[100][100];
int chrNum;

void inputChrName(char *namefile)
{
	ifstream chrFile;
	char filename[500];
	sprintf(filename, "%schrName.txt", namefile);
	chrFile.open(filename);

	char tmp[100];
	for (int i = 0; i < 100; i++)
	{
		tmp[i] = '\0';
	}

	chrFile >> tmp;
	while (tmp[0] != '\0')
	{
		strcpy(chrName[chrNum], tmp);
		chrNum++;
		for (int i = 0; i < 100; i++)
		{
			tmp[i] = '\0';
		}
		chrFile >> tmp;		
	}
	chrFile.close();
}

int main(int argc, char* argv[])
{
	int i, j, category, start, end; //category: 1 all 2 same 3 cross
	char comd[1000];
	char filepath[500], exepath[500];
	
	chrNum = 0;

	if (argc != 4)
	{
		cout << argv[0] << "\t<all/same/cross>\t<exe_path>\t<target_path>" << endl;
		return 1;
	}
	if (strcmp("all", argv[1]) == 0)
	{
		category = 1;
	}
	else if (strcmp("same", argv[1]) == 0)
	{
		category = 2;
	}
	else if (strcmp("cross", argv[1]) == 0)
	{
		category = 3;
	}
	else
	{
		cout << "unrecognized option" << endl;
		return 1;
	}

	strcpy(exepath, argv[2]);
	strcpy(filepath, argv[3]);
	
	inputChrName(filepath);

	for (i = 0; i < chrNum; i++)
	{
		if (category == 1)
		{
			start = i;
			end = chrNum-1;
		}
		else if (category == 2)
		{
			start = i;
			end = i;
		}
		else if (category == 3)
		{
			start = i+1;
			end = chrNum-1;
		}

		for (j = start; j<= end; j++)
		{
			sprintf(comd, "%sMapPER %s %s %s %s cross", exepath, chrName[i], chrName[j], exepath, filepath);
			system(comd);
			cout << "chromosome " << chrName[i] << "   " << chrName[j] << endl;
		}
	}

	return 0;
}