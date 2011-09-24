#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
using namespace std;

const string Is(200, 'I');

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
			{
				ofs << "\tNM:i:"<<VEC_iter->alters.length() / 5<<'\t';// << VEC_iter->alters<<endl;

				string differs = VEC_iter->alters;

				size_t index = 0;

				vector<Differ> differs_vec;

				if (differs.length() > 1)
				{
					while (true)
					{
						//0:G>A,1:T>A
						int idx;
						char src, dst;
						sscanf(differs.c_str() + index, "%d:%c>%c", &idx, &src, &dst);
						differs_vec.push_back(Differ(idx, src, dst));
						index = differs.find(",", index);

						if (index == string::npos)
							break;

						++index;
					}

					if (VEC_iter->strand)
					{
						reverse(differs_vec.begin(), differs_vec.end());

						vector<Differ>::iterator v_Differ_iter;
						for (v_Differ_iter = differs_vec.begin(); v_Differ_iter != differs_vec.end(); ++v_Differ_iter)
						{
							v_Differ_iter->index = VEC_iter->tag_seq.length() - v_Differ_iter->index - 1;
						}
					}

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							ofs << ',';
					}
				}

				ofs << endl;

			}
			else
				ofs << "\tNM:i:"<<VEC_iter->alters.length() / 5<</*'\t' << VEC_iter->alters<<*/endl;
			//else
			//	ofs << endl;
		}

		++chromid;
	}
}

void
sortbowtie(const char* bowtiefile, const char* outfile, int misinfo)
{
	ifstream ifs(bowtiefile);

	ofstream ofs(outfile);

	map<string, vector<bowtieRec> > bowtieVECMap;

	size_t count2 = 0;
	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			size_t offset, othermapped;
			char tagname[100], strand, chromname[100], tagseq[1000], qualstr[1000], alters[1000];
			int read_count = sscanf(line.c_str(), "%s\t%c\t%s\t%llu\t%s\t%s\t%d\t%s", tagname,  &strand, chromname, &offset, tagseq, qualstr, &othermapped, alters);

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

			if (qualstring.length() != tagseqstr.length())
			{
				cout << "qual str lenght not match "<< qualstring<<endl;
				continue;
			}

			if (strandflag)
				reverse(qualstring.begin(), qualstring.end());

			double qual_score = 255;

			if (qualstring.find_first_not_of("I") != string::npos)
			{
				//cout << "calculating qual score"<<endl;

				//cout << altersstr << endl;
				qual_score = 0;

				vector<Differ> differs_vec;

				if (altersstr.length() > 1)
				{
					//cout << differs<<endl;
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

				//cout << "mismach_vec"<<endl;

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

			size_t count = 0;
			for (size_t i = 0; i < altersstr.length(); ++i)
			{
				if (altersstr[i] == ':')
					++count;
			}

			string tagnamestr2 = tagname;

			//if (tagnamestr2.find("~") == string::npos)
			//	ofs << count2++<<'~';

			if (qual_score > 255)
				qual_score = 255;

			if (qual_score < 0)
				qual_score = 0;

			ofs << tagname << '\t' << strandflag <<  '\t' << chromname <<  '\t' <<offset + 1 <<  "\t"<<static_cast<size_t>(qual_score)<<"\t" <<tagseqstr.length() <<  "M\t*\t0\t0\t" <<tagseq <<  '\t' <<qualstr << "\tNM:i:"<<count;

			if (misinfo)
				ofs <<'\t' << altersstr<<endl;
			else
				ofs<<endl;

			//bowtieVECMap[chromname].push_back(bowtieRec(tagname, strandflag, offset, tagseq, othermapped, alterstr));

			//ofs <<VEC_iter->tag_name<<'\t'<<VEC_iter->strand<<'\t'<<bVECM_iter->first<<'\t'<<VEC_iter->offset + 1<<"\t0\t"<<VEC_iter->tag_seq.length()
			//	<<"M\t*\t0\t0\t"<<VEC_iter->tag_seq<<'\t'<<Is.substr(0, VEC_iter->tag_seq.length());
		}

		ifs.close();
		ofs.close();

		//system("mkdir temp");

		//string sort_cmd;

		////sort -k3,4 -o test.txt -T ./temp unspliced_map.bwtout.sam
		//sort_cmd = "sort -k3,3 -k4,4n -S 2000000 -o ";

		//sort_cmd.append(outfile);

		//sort_cmd.append(" -T temp ");

		//sort_cmd.append(outfile);

		//system(sort_cmd.c_str());
	}
	else
	{
		cout << "can't open file";
		exit(0);
	}
}

int 
main(int argc, const char** argv)
{
	if (argc < 3)
	{
		cout << "need two files";
		exit(0);
	}

	int misinfo = atoi(argv[3]);

	sortbowtie(argv[1], argv[2], misinfo);

	return 0;
}