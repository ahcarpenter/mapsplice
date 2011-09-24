#ifndef MISMATCH_INFO_H
#define MISMATCH_INFO_H

using namespace std;

struct Mismatch
{
public:
	int seg_no;
	int pos;
	char src_char;
	char dst_char;

	Mismatch()
	{
	}

	Mismatch(int s, int p, char src, char dst)
	{
		seg_no=s;
		pos=p;
		src_char=src;
		dst_char=dst;
	}

	~Mismatch()
	{
	}
};

#endif