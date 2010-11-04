/*    
 *    juncdb_bwt2sam.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng, Kai Wang
 *
 *    Authors: Zeng Zheng, Kai Wang
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

#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;
#pragma warning(disable:4996)

const string Is(200, 'I');
const char * jumpcode_spliter = "MN"; 
const char jumptype[] = {'M','N'};

string int_to_str(int numerical)
{
	char c[100];
	sprintf(c,"%d",numerical);
	string str(c);
	return str;
}

struct Differ{
	size_t index;
	char src_char;
	char dst_char;

	Differ(size_t idx, char src, char dst) : index(idx), src_char(src), dst_char(dst) {}
};

struct bowtieRec {
	string tag_name;
	size_t strand;
	size_t offset;
	string tag_seq;
	string alters;
	size_t othermapped;

	//vector<Differ> differs_vec;

	bowtieRec(const string& tname, const size_t& strd, const size_t& oft, const string& tseq, const size_t& omd, const string& alt) : tag_name(tname), strand(strd), offset(oft), tag_seq(tseq), othermapped(omd), alters(alt)
	{
		//if (alters.length() > 1)
		//{
		//	//cout << differs<<endl;
		//	size_t index = 0;
		//	while (true)
		//	{
		//		//0:G>A,1:T>A
		//		int idx;
		//		char src, dst;
		//		sscanf(alters.c_str() + index, "%d:%c>%c", &idx, &src, &dst);

		//		differs_vec.push_back(Differ(idx, src, dst));
		//		index = alters.find(",", index);

		//		if (index == string::npos)
		//			break;

		//		++index;
		//	}
		//}
	}
};

double CalcPSanger(char inchar)
{
	int intchar =  inchar;

	if (intchar >= 64 && intchar <= 126)
		intchar = intchar - 64;

	return pow(double(10), double(-intchar)/double(10));
}

double CalcPSolexa(char inchar)
{
	int intchar =  inchar;

	if (intchar >= 64 && intchar <= 126)
		intchar = intchar - 64;

	double ppow = pow(double(10), double(-intchar)/double(10));

	return ppow / (double(1) + ppow);
}

bool
comp(const bowtieRec& lhs, const bowtieRec& rhs)
{
	return lhs.offset < rhs.offset;
}

void
convert2sam(const map<string, vector<bowtieRec> >& bowtieVECMap, ofstream& ofs, int misinfo)
{
	map<string, vector<bowtieRec> >::const_iterator bVECM_iter;

	size_t chromid = 1;
	for (bVECM_iter = bowtieVECMap.begin(); bVECM_iter != bowtieVECMap.end();  ++bVECM_iter)
	{
		vector<bowtieRec>::const_iterator VEC_iter;

		for (VEC_iter = bVECM_iter->second.begin(); VEC_iter != bVECM_iter->second.end(); ++VEC_iter)
		{
			ofs <<VEC_iter->tag_name<<'\t'<<VEC_iter->strand<<'\t'<<bVECM_iter->first<<'\t'<<VEC_iter->offset + 1<<"\t0\t"<<VEC_iter->tag_seq.length()
				<<"M\t*\t0\t0\t"<<VEC_iter->tag_seq<<'\t'<<Is.substr(0, VEC_iter->tag_seq.length());

			//if (VEC_iter->tag_seq.length() > 80)

			if (misinfo)
				ofs << "\tNM:i:"<<VEC_iter->alters.length() / 5<<'\t' << VEC_iter->alters<<endl;
			else
				ofs << "\tNM:i:"<<VEC_iter->alters.length() / 5<</*'\t' << VEC_iter->alters<<*/endl;
			//else
			//	ofs << endl;
		}

		++chromid;
	}
}

void compute_offset_jumpcode(size_t real_offset, char* real_jumpcode, int offset, int read_len, vector<string>& result)
{
	string final_jumpcode;
	size_t final_offset=real_offset;
	char* p;
	p = strtok(real_jumpcode ,jumpcode_spliter); 
	int step=0;
	int match_len=0;
	while(p != NULL)          //skip offset
	{ 
		if(step%2 == 0)    //match
		{
			match_len=atoi(p);
			if(match_len>offset)    //skipped enough
			{
				final_offset+= offset;
				match_len -= offset;
				break;
			}
			else             // still need to skip
			{
				final_offset+= match_len;
				offset -= match_len;
			}
		} 
		else if(step%2 == 1)  //skip
		{
			final_offset+=atoi(p);
		}
		step++;
		p = strtok(NULL , jumpcode_spliter); 
	} 
	p = strtok(NULL , jumpcode_spliter);     //get next N
	if(match_len>read_len)          //generate jump code
	{
		final_jumpcode+= int_to_str(read_len)+"M";
	}
	else
	{
		final_jumpcode+= int_to_str(match_len)+"M";
		read_len-=match_len;
		step=1;
		while(p != NULL)       
		{ 
			if(step%2 == 0)    //match
			{
				int match_len=atoi(p);
				if(match_len>=read_len)
				{
					final_jumpcode+= int_to_str(read_len)+"M";
					break;
				}
				else
				{
					final_jumpcode+= int_to_str(match_len)+"M";
					read_len-=match_len;
				}
			}                
			else                  //skip
			{
				string new_jump=p;
				final_jumpcode+=new_jump + "N";
			}
			step++;
			p = strtok(NULL , jumpcode_spliter); 
		} 
	}
	result.clear();
	result.push_back(int_to_str((int)final_offset));
	result.push_back(final_jumpcode);
}



void
sortbowtie(const char* bowtiefile, const char* outfile, int misinfo)
{
	ifstream ifs(bowtiefile);
	ofstream ofs(outfile);
	map<string, vector<bowtieRec> > bowtieVECMap;
	if (ifs.is_open())
	{
		cout<<"Converting....."<<endl;
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;
			size_t offset, othermapped;
			char tagname[1000], strand, chromname[1000], tagseq[1000], qualstr[1000], alters[1000];
			int read_count = sscanf(line.c_str(), "%s\t%c\t%s\t%lu\t%s\t%s\t%lu\t%s", tagname,  &strand, chromname, &offset, tagseq, qualstr, &othermapped, alters);
                        

			string alterstr = "";
			if (read_count == 8)
			{
				alterstr = alters;
			}
			else if (read_count != 7)
			{
				cout << "read failed";
				exit(0);
			}
			size_t strandflag = strand == '+' ? 0 : 16;
			string tagseqstr = tagseq;
			string altersstr = alterstr;
			string qualstring = qualstr;
			if (strandflag)
				reverse(qualstring.begin(), qualstring.end());
			if (qualstring.length() != tagseqstr.length())
			{
				cout << "qual str lenght not match "<< qualstring<<endl;
				continue;
			}
			double qual_score = 255;
			if (qualstring.find_first_not_of("I") != string::npos)
			{
				qual_score = 0;
				vector<Differ> differs_vec;
				if (altersstr.length() > 1)
				{
					size_t index = 0;
					while (true)
					{
						int idx;
						char src, dst;
						sscanf(altersstr.c_str() + index, "%d:%c>%c", &idx, &src, &dst);
						differs_vec.push_back(Differ(idx, src, dst));
						index = altersstr.find(",", index);
						if (index == string::npos)
							break;
						++index;
					}
				}
				vector<bool> mismach_vec(tagseqstr.length(), false);
				vector<Differ>::const_iterator dfv_iter;
				for (dfv_iter = differs_vec.begin(); dfv_iter != differs_vec.end(); ++dfv_iter)
				{
					mismach_vec[dfv_iter->index] = true;
					qual_score += log(CalcPSanger(qualstring[dfv_iter->index]) / (double(1) - double(0.25)));
				}
				for (size_t i = 0; i < mismach_vec.size(); ++i)
				{
					if (!mismach_vec[i])
					{
						qual_score += log((double(1) - CalcPSanger(qualstring[i])) / double(0.25));
					}
				}
			}
			char real_chrom[1000], real_jumpcode[1000];
			size_t real_offset;
			sscanf(chromname, "%[^_]_%lu:%s",real_chrom,&real_offset, real_jumpcode);
			vector<string> final_offset_jumpcode;
			compute_offset_jumpcode(real_offset, real_jumpcode, (int)offset, (int)tagseqstr.length(), final_offset_jumpcode);

			size_t count = 0;
			for (size_t i = 0; i < altersstr.length(); ++i)
			{
				if (altersstr[i] == ':')
					++count;
			}

			if (qual_score > 255)
				qual_score = 255;

			if (qual_score < 0)
				qual_score = 0;

			ofs << tagname << '\t' << strandflag <<  '\t' << real_chrom <<  '\t' <<final_offset_jumpcode[0]<<  "\t"<<static_cast<int>(qual_score)<<"\t" <<final_offset_jumpcode[1]<<"\t*\t0\t0\t" <<tagseq <<  '\t' <<qualstr << "\tNM:i:"<<count;
			if (misinfo)
				ofs <<'\t' << altersstr<<endl;
			else
				ofs<<endl;
		}
		ifs.close();
		ofs.close();
		/*cout<<"Sorting....."<<endl;
		system("mkdir -p temp");
		string sort_cmd;
		sort_cmd = "sort -k3,3 -k4,4n -S 2000000 -o ";
		sort_cmd.append(outfile);
		sort_cmd.append(" -T temp ");
		sort_cmd.append(outfile);
		system(sort_cmd.c_str());

		cout <<"sort finished"<<endl;*/
	}
	else
	{
		cout << "can't open file"<<endl;
		exit(0);
	}
}

int 
main(int argc, const char** argv)
{
	if (argc < 4)
	{
	cout << "need two files and misinfo value";
	exit(0);
	}
	int misinfo = atoi(argv[3]);

	sortbowtie(argv[1], argv[2], misinfo);

	//sortbowtie("G://new_mapsplice//1.bwtout", "G://new_mapsplice//2.sam", 1);
	cout <<"sortbowtie finished"<<endl;
	return 0;
}
