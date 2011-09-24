/*    
 *    Compare2Sam.cpp		
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
	char* sam_file1 = argv[1];

	char* sam_file2 = argv[2];

	char* output_path = argv[3];

	SamHandler s1, s2;

	bool open1 = s1.OpenSam(sam_file1);

	bool open2 = s2.OpenSam(sam_file2);

	s1.Compare2Sam(s2, output_path);

	string comp_stat_file = output_path; comp_stat_file.append(s1.base_file_name); comp_stat_file.append("_vs_"); comp_stat_file.append(s2.base_file_name); comp_stat_file.append(".stat");

	s1.WriteCompStat(s2, comp_stat_file.c_str());

	string s1_stat_file = output_path; s1_stat_file.append(s1.base_file_name); s1_stat_file.append(".stat");

	s1.WriteStat(s1_stat_file.c_str());

	string s2_stat_file = output_path; s2_stat_file.append(s2.base_file_name); s2_stat_file.append(".stat");

	s2.WriteStat(s2_stat_file.c_str());

	//char* compared_file = argv[3];

	//size_t range = atoi(argv[4]);

	//char* junc_file = argv[5];

	//char* filtered_multiple = argv[6];

	//size_t maxhits = atoi(argv[7]);

	//char* unique_mapped = argv[8];

	//char* compared_offset = argv[9];

	//char* stat_file = argv[10];

	//int anchor_width = atoi(argv[11]);

	//SamHandler mps_samhandler, compare_samhandler;

	//cout << "read sam " << mps_sam << endl;

	//mps_samhandler.ReadAppendSamFile(mps_sam);

	//mps_samhandler.CollectStat();

	//cout << "sort sam " << mps_sam << endl;

	//mps_samhandler.SortSamByTagnameOffset();

	//cout << "filter multiple mapped " << mps_sam << endl;

	//mps_samhandler.FilterMultipleMappedTag(filtered_multiple, junc_file, maxhits);

	//cout << "write unique mapped " << mps_sam << endl;

	//mps_samhandler.WriteUniqueMappedTagSam(unique_mapped);

	//cout << "read sam " << compare_sam << endl;

	//compare_samhandler.ReadAppendSamFile(compare_sam);

	//compare_samhandler.CollectStat();

	//cout << "sort sam " << compare_sam << endl;

	//compare_samhandler.SortSamByTagnameOffset();

	//cout << "compare sam " << mps_sam <<'\t' <<compare_sam << endl;

	//mps_samhandler.Compare2Sam(compare_samhandler, compared_file, compared_offset, anchor_width, range);

	//mps_samhandler.WriteCompStat(compare_samhandler, stat_file);
}