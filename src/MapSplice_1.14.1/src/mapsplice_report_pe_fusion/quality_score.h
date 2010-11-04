#ifndef QUALITY_SCORE_H
#define QUALITY_SCORE_H

#include <string>
#include <math.h>
using namespace std;
#pragma warning(disable:4996)


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

int GetQualityScore(const vector<Mismatch>& mismatch_vec, const string quality_str)
{
	double qual_score = 0;
	vector<bool> mismatch_score(quality_str.length(), false);
	vector<Mismatch>::const_iterator mis_iter;
	const   char   *quality_str_c=quality_str.c_str();   
	for (mis_iter = mismatch_vec.begin(); mis_iter != mismatch_vec.end(); ++mis_iter)
	{
		mismatch_score[mis_iter->pos] = true;
		//qual_score += (CalcPSanger(quality_str_c[mis_iter->pos]) / (double(1) - double(0.25)));
		qual_score += log(CalcPSanger(quality_str_c[mis_iter->pos]) /(double(1) - double(0.25)));
	}
	for (size_t i = 0; i < mismatch_score.size(); ++i)
	{
		if (!mismatch_score[i])
		{
			//qual_score += ((double(1) - CalcPSanger(quality_str_c[i])) / double(0.25));
			qual_score += log((double(1) - CalcPSanger(quality_str_c[i])) / double(0.25));
		}
	}
	if(qual_score<0)
		qual_score=0;
	if(qual_score>255)
		qual_score=255;
	return (int)qual_score;
}


#endif