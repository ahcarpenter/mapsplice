/*    
 *    convert_to_abs_offset.cpp		
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
#include <vector>

#include <string>

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>
using namespace std;

inline char
complement(int i) {
	static const int b2c_size = 20;
	static const char b2c[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
	};
	static const char b2cl[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
	};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

inline string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

size_t parse_spliceway(const string& splicewaystr, vector<pair<size_t, size_t> >& spliceway_vec, size_t prefixst)
{
	size_t index = 0;

	bool exceed = false;
	while (true)
	{
		int maplen, intron;

		if (index == 0)
		{
			sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
			spliceway_vec.push_back(make_pair(prefixst, maplen));

			if (maplen < 0/* || maplen < min_anchor*/)
				exceed =true;

			index = splicewaystr.find("M", index) + 1;
		}
		else
		{
			if (index == splicewaystr.length())
				break;
			sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);

			if (intron >= 0 && maplen >= 0)
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));

			if (maplen < 0/* || maplen < min_anchor*/)
				exceed =true;

			if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
				break;

			index = splicewaystr.find("M", index) + 1;
		}
	}

	return spliceway_vec.size();
}

size_t convert_to_abs_offset(vector<pair<size_t, size_t> >& spliceway_vec, size_t st, size_t end, size_t strand)
{
	if (strand)
	{
		size_t length = end - st;

		for (size_t i = 0; i < spliceway_vec.size(); ++i)
		{
			size_t temp = length - (spliceway_vec[i].first + spliceway_vec[i].second - 1) + 1;

			spliceway_vec[i].first = temp + st;
		}

		reverse(spliceway_vec.begin(), spliceway_vec.end());
	}
	else
	{
		for (size_t i = 0; i < spliceway_vec.size(); ++i)
			spliceway_vec[i].first += st;
	}

	return spliceway_vec.size();
}

int main(int argc, char* argv[])
{
	ifstream ifs_all(argv[1]);

	ofstream ofs(argv[2]);

	string not_converted = argv[2]; not_converted.append(".not_converted");

	ofstream ofs_not_converted(not_converted.c_str());

	ofstream ofs_fusion(argv[3]);

	string abs_path(argv[4]);

	if (ifs_all.is_open())
	{
		while (!ifs_all.eof() )
		{
			string line2;
			getline(ifs_all,line2);

			char path[1000];

			sscanf(line2.c_str(), "%s",  path);

			//cout << line2 << endl;

			//cout << path<<endl;

			string cur_file = abs_path;cur_file.append(path);cur_file.append("/output/tmp/original/ori_paired.sam");

			//cout << cur_file << endl;

			//getchar();

			ifstream ifs(cur_file.c_str());
			
			if (ifs.is_open())
			{
				while (!ifs.eof() )
				{
					string line;
					getline(ifs,line);
					if (line == "")
						continue;

					ofs_not_converted << line << endl;

					char chromname[1000], readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
					char strand = '+';
					size_t prefixst, strand_t, incorrect, mate_offest, readidx;

					int mate_diff;

					unsigned short mis_match;

					char mate_match;

					//12~11400913#0/2	16	chr2~1980062~1980130~chr2~2036858~2155720	109558	255	50M	*	0	0	AAAAAAAAAAAAAATTAGCGGGGTGCGGTGGCGGGCGCCTGTAGTCCCAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	
					//34~15374907#0/2	0	chr2~1980062~1980130~chr2~2036858~2155720	100998	255	50M	*	0	0	CCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCTCTTGAACCTGGGAGGC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0
					sscanf(line.c_str(), "%llu~%s\t%llu\t%s\t%llu\t%llu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu", 
						&readidx, readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, &mate_match, &mate_offest, &mate_diff, chromseq, qualseq, &mis_match);

					prefixst = prefixst - 1;

					char chrom1[1000], chrom2[1000];

					size_t strand1, start1, end1, strand2, start2, end2;

					sscanf(chromname, "%[^~]~%llu~%llu~%llu~%[^~]~%llu~%llu~%llu", chrom1, &strand1, &start1, &end1, chrom2, &strand2, &start2, &end2);

					string chrom1str = chrom1, chrom2str = chrom2;

					size_t length1 = end1 - start1;

					size_t length2 = end2 - start2;

					vector<pair<size_t, size_t> > spliceway_vec;

					parse_spliceway(spliceway, spliceway_vec, prefixst);

					//string splicewaystr = spliceway;
					//size_t index = 0;

					//bool exceed = false;
					//while (true)
					//{
					//	int maplen, intron;

					//	if (index == 0)
					//	{
					//		sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
					//		spliceway_vec.push_back(make_pair(prefixst, maplen));

					//		if (maplen < 0/* || maplen < min_anchor*/)
					//			exceed =true;

					//		index = splicewaystr.find("M", index) + 1;
					//	}
					//	else
					//	{
					//		if (index == splicewaystr.length())
					//			break;
					//		sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);

					//		if (intron >= 0 && maplen >= 0)
					//			spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));

					//		if (maplen < 0/* || maplen < min_anchor*/)
					//			exceed =true;

					//		if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
					//			break;

					//		index = splicewaystr.find("M", index) + 1;
					//	}
					//}

					size_t end_pos = spliceway_vec.back().first + spliceway_vec.back().second - 1;

					vector<pair<size_t, size_t> > spliceway_vec_p1, spliceway_vec_p2;

					if (end_pos <= length1)
					{
						spliceway_vec_p1 = spliceway_vec;
					}
					else if (prefixst > length1)
					{
						spliceway_vec_p2 = spliceway_vec;

						for (size_t i = 0; i < spliceway_vec_p2.size(); ++i)
						{
							spliceway_vec_p2[i].first = spliceway_vec_p2[i].first - length1;
						}
					}
					else if (prefixst <= length1 && end_pos > length1)
					{
						size_t i;

						for (i = 0; i < spliceway_vec.size(); ++i)
						{
							if (spliceway_vec[i].first + spliceway_vec[i].second - 1 <= length1)
								spliceway_vec_p1.push_back(make_pair(spliceway_vec[i].first, spliceway_vec[i].second));
							else
								break;
						}

						for (; i < spliceway_vec.size(); ++i)
						{
							if (spliceway_vec[i].first <= length1 && spliceway_vec[i].first + spliceway_vec[i].second - 1 > length1)
							{
								spliceway_vec_p1.push_back(make_pair(spliceway_vec[i].first, length1 - spliceway_vec[i].first + 1));

								spliceway_vec_p2.push_back(make_pair(1, spliceway_vec[i].second - (length1 - spliceway_vec[i].first + 1)));
							}
							else
								break;
						}

						for (; i < spliceway_vec.size(); ++i)
						{
							if (spliceway_vec[i].first > length1)
								spliceway_vec_p2.push_back(make_pair(spliceway_vec[i].first - length1, spliceway_vec[i].second));
							else
								break;
						}
					}

					if (spliceway_vec_p1.size())
						convert_to_abs_offset(spliceway_vec_p1, start1, end1, strand1);

					if (spliceway_vec_p2.size())
						convert_to_abs_offset(spliceway_vec_p2, start2, end2, strand2);

					if (spliceway_vec_p1.size() && !spliceway_vec_p2.size())
					{
						if (strand1 == 16)
						{
							if (strand_t == 0)
								strand_t = 16;
							else
								strand_t = 0;
						}

						ofs <<readname<<'\t'<<strand_t<<'\t'<<chrom1<<'\t'<<spliceway_vec_p1.front().first<<'\t'<<incorrect<<'\t';

						vector<pair<size_t, size_t> >::iterator spliceway_iter;

						for (spliceway_iter = spliceway_vec_p1.begin(); spliceway_iter != spliceway_vec_p1.end() - 1; ++spliceway_iter)
							ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

						ofs << spliceway_iter->second << 'M';

						string chromstr = chromseq, qualstr = qualseq;

						if (strand1 == 16)
						{
							chromstr = revcomp(chromstr);

							reverse(qualstr.begin(), qualstr.end());
						}

						ofs<<'\t'<<mate_match<<'\t'<<mate_offest<<'\t'<<mate_diff<<'\t'<<chromstr<<'\t'<< qualstr<<"\tNM:i:"<<mis_match<<endl;
					}
					else if (!spliceway_vec_p1.size() && spliceway_vec_p2.size())
					{
						if (strand2 == 16)
						{
							if (strand_t == 0)
								strand_t = 16;
							else
								strand_t = 0;
						}

						ofs <<readname<<'\t'<<strand_t<<'\t'<<chrom2<<'\t'<<spliceway_vec_p2.front().first<<'\t'<<incorrect<<'\t';

						vector<pair<size_t, size_t> >::iterator spliceway_iter;

						for (spliceway_iter = spliceway_vec_p2.begin(); spliceway_iter != spliceway_vec_p2.end() - 1; ++spliceway_iter)
							ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

						ofs << spliceway_iter->second << 'M';	

						string chromstr = chromseq, qualstr = qualseq;

						if (strand2 == 16)
						{
							chromstr = revcomp(chromstr);

							reverse(qualstr.begin(), qualstr.end());
						}

						ofs<<'\t'<<mate_match<<'\t'<<mate_offest<<'\t'<<mate_diff<<'\t'<<chromstr<<'\t'<< qualstr<<"\tNM:i:"<<mis_match<<endl;
					}
					else if (spliceway_vec_p1.size() && spliceway_vec_p2.size())
					{
						//TRAN00000058995:297:67    chr1    +    13794936    13795002    chr15    +    80328563 80328595 TGGACAGCGGCAGGGACATCATAAATACAGAAGATAATTTTAAAGAAGTGTTCCggccgggtgtggtggctcacaactgtaataccaaaactttaggaaataggtagaaggcttgcatcc GTAG    NM:i:4
						if (false && chrom1str == chrom2str && strand1 == 0 && strand2 == 0)
						{
							if (strand1 == 16)
							{
								if (strand_t == 0)
									strand_t = 16;
								else
									strand_t = 0;
							}

							ofs <<readname<<'\t'<<strand_t<<'\t'<<chrom1<<'\t'<<spliceway_vec_p1.front().first<<'\t'<<incorrect<<'\t';//<< length1 - prefixst + 1<<'M'<<start2 - end1 - 1<< 'N'<<end_pos - length1<<'M'

							vector<pair<size_t, size_t> >::iterator spliceway_iter;

							for (spliceway_iter = spliceway_vec_p1.begin(); spliceway_iter != spliceway_vec_p1.end() - 1; ++spliceway_iter)
								ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

							ofs << spliceway_iter->second << 'M';

							ofs << start2 - end1<< 'N';

							for (spliceway_iter = spliceway_vec_p2.begin(); spliceway_iter != spliceway_vec_p2.end() - 1; ++spliceway_iter)
								ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

							ofs << spliceway_iter->second << 'M';

							string chromstr = chromseq, qualstr = qualseq;

							if (strand1 == 16)
							{
								chromstr = revcomp(chromstr);

								reverse(qualstr.begin(), qualstr.end());
							}

							ofs <<'\t'<<mate_match<<'\t'<<mate_offest<<'\t'<<mate_diff<<'\t'<<chromstr<<'\t'<< qualstr<<"\tNM:i:"<<mis_match<<endl;
						}
						else
						{
							if (spliceway_vec_p1.size() > 1 && spliceway_vec_p2.size() > 1)
							{
								cout << line << endl;
								continue;
							}

							char strand_c1 = strand_t ? '-' : '+';

							char strand_c2 = strand_t ? '-' : '+';

							if (strand1 == 16)
							{
								strand_c1 = strand_t ? '+' : '-';
							}

							if (strand2 == 16)
							{
								strand_c2 = strand_t ? '+' : '-';
							}

							size_t prefix_st, prefix_end, suffix_st, suffix_end;
							if (strand_c1 == '+')
							{
								prefix_st = spliceway_vec_p1.front().first;
								prefix_end = spliceway_vec_p1.front().first + spliceway_vec_p1.front().second - 1;
							}
							else
							{
								prefix_end  = spliceway_vec_p1.front().first + spliceway_vec_p1.front().second - 1;
								prefix_st = spliceway_vec_p1.front().first; 
							}

							if (strand_c1 == '+')
							{
								suffix_st = spliceway_vec_p2.front().first;
								suffix_end = spliceway_vec_p2.front().first + spliceway_vec_p2.front().second - 1;
							}
							else
							{
								suffix_end = spliceway_vec_p2.front().first + spliceway_vec_p2.front().second - 1;
								suffix_st = spliceway_vec_p2.front().first;
							}							

							//size_t strand_t1 = strand_t, strand_t2 = strand_t;

							//if (strand1 == 16)
							//{
							//	//if (strand_t1 == 0)
							//	//	strand_t1 = 16;
							//	//else
							//	//	strand_t1 = 0;
							//}

							//if (strand2 == 16)
							//{
							//	//if (strand_t2 == 0)
							//	//	strand_t2 = 16;
							//	//else
							//	//	strand_t2 = 0;
							//}

							//if (strand1 == 16)
							//{
							//	chromstr = revcomp(chromstr);

							//}
							//ofs_fusion <<readname<<'\t'<<strand_t1<<'\t'<<chrom1<<'\t'<<spliceway_vec_p1.front().first<<'\t'<<incorrect<<'\t';

							//vector<pair<size_t, size_t> >::iterator spliceway_iter;

							//size_t prefixlen = 0;

							//for (spliceway_iter = spliceway_vec_p1.begin(); spliceway_iter != spliceway_vec_p1.end() - 1; ++spliceway_iter)
							//{
							//	ofs_fusion << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

							//	prefixlen += spliceway_iter->second;
							//}

							//ofs_fusion << spliceway_iter->second << 'M';

							//prefixlen += spliceway_iter->second;

							//string chromstr = chromseq, qualstr = qualseq;

							//string chromstr1 = chromstr.substr(0, prefixlen);

							//string chromstr2 = chromstr.substr(prefixlen, chromstr.length() - prefixlen);

							//string qualstr1 = qualstr.substr(0, prefixlen);

							//string qualstr2 = qualstr.substr(prefixlen, qualstr.length() - prefixlen);

							//if (strand1 == 16)
							//{
							//	chromstr1 = revcomp(chromstr1);

							//	reverse(qualstr1.begin(), qualstr1.end());
							//}

							//ofs_fusion<<'\t'<<chrom2<<'\t'<<spliceway_vec_p2.front().first<<'\t'<<spliceway_vec_p2.front().first - spliceway_vec_p1.front().first<<'\t'<<chromstr1<<'\t'<< qualstr1<<"\tNM:i:"<<mis_match<<endl;


							////
							//ofs_fusion <<readname<<'\t'<<strand_t2<<'\t'<<chrom2<<'\t'<<spliceway_vec_p2.front().first<<'\t'<<incorrect<<'\t';

							////vector<pair<size_t, size_t> >::iterator spliceway_iter2;

							//size_t suffixlen = 0;

							//for (spliceway_iter = spliceway_vec_p2.begin(); spliceway_iter != spliceway_vec_p2.end() - 1; ++spliceway_iter)
							//{
							//	ofs_fusion << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';

							//	suffixlen += spliceway_iter->second;
							//}

							//ofs_fusion << spliceway_iter->second << 'M';

							//suffixlen += spliceway_iter->second;

							//if (strand2 == 16)
							//{
							//	chromstr2 = revcomp(chromstr2);

							//	reverse(qualstr2.begin(), qualstr2.end());
							//}

							//ofs_fusion<<'\t'<<chrom1<<'\t'<<spliceway_vec_p1.front().first<<'\t'<<spliceway_vec_p1.front().first - spliceway_vec_p2.front().first<<'\t'<<chromstr2<<'\t'<< qualstr2<<"\tNM:i:"<<mis_match<<endl;

							ofs_fusion <<readname<<'\t'<<chrom1<<'\t'<<strand_c1<<'\t'<<prefix_st<<'\t'<<prefix_end <<'\t'
								<<chrom2<<'\t'<<strand_c2<<'\t'<<suffix_st<<'\t'<<suffix_end<<'\t'<< chromseq<<"\tNNNN\tNM:i:"<<mis_match<<endl;

							//ofs_fusion <<readname<<'\t'<<chrom1<<'\t'<<strand<<'\t'<<prefixst + start1<<'\t'<<end1 <<'\t'<<chrom2<<'\t'<<strand<<'\t'<<start2<<'\t'<<start2 + end_pos - length1 - 1
							//	<<'\t'<< chromseq<<"\tNNNN\tNM:i:"<<mis_match<<endl; 
						}

						/*vector<pair<size_t, size_t> >::iterator spliceway_iter;

						if (chrom1str == chrom2str)
						{
							ofs <<readname<<'\t'<<strand_t<<'\t'<<chrom1<<'\t'<<prefixst + start1<<'\t'<<incorrect<<'\t';

							for (spliceway_iter = spliceway_vec.begin(); spliceway_iter != spliceway_vec.end() - 1; ++spliceway_iter)
							{
								if (spliceway_iter->first <= length1 && spliceway_iter->first + spliceway_iter->second > length1)
								{
									ofs << length1 - spliceway_iter->first + 1<<'M'<<start2 - end1 - 1<< 'N'<<spliceway_iter->first + spliceway_iter->second - length1 - 1<<'M'
										<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';
								}
								else if (spliceway_iter->first + spliceway_iter->second <= length1 && (spliceway_iter + 1)->first > length1)
								{
									ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second + start2 - end1 - 1<<'N';
								}
								else
								{
									ofs << spliceway_iter->second << 'M'<< (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second <<'N';
								}
							}

							if (spliceway_iter->first <= length1 && spliceway_iter->first + spliceway_iter->second > length1)
								ofs << length1 - spliceway_iter->first + 1<<'M'<<start2 - end1 - 1<< 'N'<<spliceway_iter->first + spliceway_iter->second - length1 - 1<<'M';
							else
								ofs << spliceway_iter->second << 'M';

							ofs<<'\t'<<mate_match<<'\t'<<mate_offest<<'\t'<<mate_diff<<'\t'<<chromseq<<'\t'<< qualseq<<"\tNM:i:"<<mis_match<<endl;
						}
						else
						{
							if (spliceway_vec.size() > 2)
							{
								cout << line << endl;
							}
							if ((spliceway_vec.front().first <= length1 && spliceway_vec.front().first + spliceway_vec.front().second > length1) ||
								(spliceway_vec.back().first <= length1 && spliceway_vec.back().first + spliceway_vec.back().second > length1))
							{
								cout << line << endl;
							}
							else
							{
								char strand = strand_t ? '-' : '+';

								ofs_fusion <<readname<<'\t'<<chrom1<<'\t'<<strand<<'\t'<<spliceway_vec.front().first + start1<<'\t'<<spliceway_vec.front().first + spliceway_vec.front().second + start1 - 2 <<'\t'
									<<chrom2<<'\t'<<strand<<'\t'<<spliceway_vec.back().first - length1 + start2 - 1<<'\t'<<start2 + end_pos - length1 - 1<<'\t'<< chromseq<<"\tNNNN\tNM:i:"<<mis_match<<endl;
							}
						}*/
						//}
					}
				}

				ifs.close();
			}
		}
	}
}