/*    
 *    junc2bed.cpp	
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
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
using namespace std;

int
main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "need two arguments "<<endl; exit(1);
	}

	string outfile_nofirstline = argv[argc - 1]; outfile_nofirstline.append(".nfs");

	ofstream ofs(argv[argc - 1]);

	ofstream ofs_nfs(outfile_nofirstline.c_str());

	string skipline;

	for (size_t i = 1; i < argc - 1; ++i)
	{
		ifstream ifs(argv[i]);

		if (ifs.is_open())
		{
			getline(ifs,skipline);
			
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "")
					continue;

				char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], lpq[100], flankchar[10], il[100];
				int juncst, juncend, prefixend, suffixst, kinds, hits, flankcase;

				unsigned short min_mis, max_mis;

				double ave_mis;

				/*ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<'+'<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< flankcase <<'\t'<<flankstr <<'\t'<<il<<'\t'<<vit->lpq
				<<'\t' << vit->min_mismatch<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<< endl;*/

				//chr20	16408	24646	JUNC_1	125	+	16408	24646	255,0,0	2	58,93,	0,8332,	3.71319	5	GTAG	0.999976	1	1	2	1.056

				sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst,
					juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, &flankcase, flankchar, il, lpq, &min_mis, &max_mis, &ave_mis);

				int prelen, suffixlen;

				sscanf(blocks, "%d,%d", &prelen, &suffixlen);

				string rgbstr = "";
				if (flankcase)
					rgbstr = "255,0,0";
				else
					rgbstr = "0,255,0";

				ofs_nfs<< chromname<<'\t'<<prefixend - prelen<<'\t'<<suffixst + suffixlen - 1<<'\t'<< juncname<<'\t'<< hits<<'\t'<< strand<<'\t'
					<<juncst<<'\t'<<juncend<<'\t'<<rgbstr<<'\t'<<kinds<<'\t'<< prelen<<','<< suffixlen<<'\t' <<0<<','<<suffixst - 1 - (prefixend - prelen)<<endl;
				//string chromstr = chromname;
				//m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset,	rank,  lpq, min_mis, max_mis, ave_mis));
			}
		}
		ifs.close();
	}

	ofs << skipline << endl;

	ofs_nfs.close();

	system("mkdir -p temp");

	string sort_cmd;

	//sort -k3,4 -o test.txt -T ./temp unspliced_map.bwtout.sam
	sort_cmd = "sort -k1,1 -k2,3n -S 2000000 -o ";

	sort_cmd.append(outfile_nofirstline);

	sort_cmd.append(" -T temp ");

	sort_cmd.append(outfile_nofirstline);

	system(sort_cmd.c_str());

	ifstream ifs_nfs(outfile_nofirstline.c_str());

	while (!ifs_nfs.eof() )
	{
		string line;
		getline(ifs_nfs,line);
		if (line == "")
			continue;

		ofs << line << endl;
	}

	ofs.close();

	string rm_command = "rm ";

	rm_command.append(outfile_nofirstline);

	system(rm_command.c_str());
}