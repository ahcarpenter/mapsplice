/*
 *  mapsplice.cpp
 *  Mapsplice
 *
 *  Created by Kai Wang on 09/09/09.
 *  Author: Kai Wang, Jinze Liu, Jan F. Prins
 *  Copyright 2009 Kai Wang. All rights reserved.
 *
 */

#define _CRT_SECURE_NO_WARNINGS

//#define VS

#include <limits>//vc only
#include <iostream>
#include <vector>
/////////////////////////ReadInfo//////////////////////////

struct ReadInfo {
	ReadInfo(size_t scr,
		size_t chr = std::numeric_limits<size_t>::max( ),
		size_t ls = std::numeric_limits<size_t>::max( ),
		size_t le = std::numeric_limits<size_t>::max( ),
		size_t rs = std::numeric_limits<size_t>::max( ),
		size_t re = std::numeric_limits<size_t>::max( ),
		size_t str = true, size_t unq = true) :
	chrom(chr), lsite_start(ls), lsite_end(le), rsite_start(rs), rsite_end(re), score(scr), strand(str), unique(unq) {}
	void set(size_t scr, size_t chr, size_t ls, size_t le, size_t rs, size_t re, size_t str);
	unsigned chrom  : 15;
	unsigned score  : 15;
	unsigned strand : 1;
	unsigned unique : 1;

	// a read might have multiple matches, so this might be a vector rather than a single item	
	unsigned lsite_start : 32;	
	unsigned lsite_end : 32;	
	unsigned rsite_start : 32;	
	unsigned rsite_end : 32;	

};

inline void 
ReadInfo::set(size_t scr, size_t chr, size_t ls, size_t le, size_t rs, size_t re, size_t str) {
	unique = (scr < score || (lsite_start == ls && unique && chrom == chr));
	chrom = chr;
	lsite_start = ls;
	lsite_end = le;
	rsite_start = rs;
	rsite_end = re;
	score = scr;
	strand = str;
}

static const size_t max_size_t = std::numeric_limits<size_t>::max( );

///////////////////////////////ConjInfo//////////////////////

struct ConjInfo {
	ConjInfo(size_t scr, size_t flk, bool bd, bool rv, size_t pm,
		size_t chr = std::numeric_limits<size_t>::max( ),
		size_t hts = 0, size_t loc = std::numeric_limits<size_t>::max( ),
		size_t read_width = 36) :
	chrom(chr), flank(flk), bad(bd), prim(pm), rev(rv), score(scr), hits(hts), kinds(1), marked(0), loc_count(std::vector<unsigned>(read_width - 1, 0))
	{
		marked = ((size_t)1 << loc);

		if (loc < read_width)
			++loc_count[loc - 1];
		else
			std::cout << "missed loc ?";
	}

	void set_score (size_t scr)
	{
		score = scr;
	}

	void inc_hits ()
	{
		++hits;
	}

	void inc_kinds ()
	{
		++kinds;
	}

	bool mark(const size_t& loc)
	{
		if (loc > loc_count.size())
		{
			std::cout << "too large loc, error?";
		}

		++loc_count[loc - 1];

		if ((marked & ((size_t)1 << loc)) == 0)
		{
			marked |= ((size_t)1 << loc);
			++kinds;
			return true;
		}
		return false;
	}


	unsigned chrom  : 15;
	unsigned score  : 15; // minimum mismatches
	unsigned hits : 15;  //number of hits
	unsigned flank : 8;
	unsigned prim : 8;
	bool rev;
	bool bad;

	size_t marked : 64;
	unsigned kinds : 15;

	std::vector<unsigned> loc_count;
};

#include <string>

#ifdef VS
#include <hash_map> //vc only
#include <hash_set>
//#define _CRT_SECURE_NO_WARNINGS 
//#pragma warning(disable:_CRT_SECURE_NO_WARNINGS)
 
#else
#include <ext/hash_map> //g++ only
#include <ext/hash_set>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;

namespace __gnu_cxx
{
	template<typename Traits, typename Allocator>
	struct hash<std::basic_string<char, Traits, Allocator> >
	{
		size_t operator()(const std::basic_string<char, Traits, Allocator>& __s) const
		{
			return __stl_hash_string(__s.c_str());
		}
	};
}

#endif

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
#ifdef VS
using namespace stdext;
#endif

const size_t THIRTY_TWO = 32;
const size_t ALL_BITS_ON = static_cast<size_t>(-1);
const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

struct JuncInfo{
	JuncInfo(/*int pm, *//*const string fs, */size_t loc, size_t suffix_len, size_t rw, size_t tagidx, unsigned short mis) : /*prim(pm), *//*flankstr(fs),*/ p(rw - 1, 0) 
	{
		++p[loc-1];
		max_prefix_len = loc;
		max_suffix_len = suffix_len;

		max_mismatch = mis;

		min_mismatch = mis;

		sum_mismatch = mis;

		m[tagidx] = 1;
	}
	bool inc_hits(size_t idx, size_t suffix_len, size_t tagidx, unsigned short mis)
	{
		if (m.find(tagidx) != m.end())
			return false;

		++p[idx-1];

		m[tagidx] = 1;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;

		if (mis > max_mismatch)
			max_mismatch = mis;
		
		if (mis < min_mismatch)
			min_mismatch = mis;

		sum_mismatch += mis;

		return true;
	}
	//int prim;
	//string flankstr;
	size_t max_prefix_len;
	size_t max_suffix_len;
	vector<unsigned short> p;

	unsigned short max_mismatch;

	unsigned short min_mismatch;

	unsigned short sum_mismatch;

	map<size_t, int> m;
};

struct JuncInfoFusion{
	JuncInfoFusion( size_t loc, size_t suffix_len, size_t rw, const string& ss, char sd1, char sd2, unsigned short mis, bool bd, size_t fbits, unsigned short pm) : p(rw - 1, 0), max_prefix_len(loc), 
		max_suffix_len(suffix_len),	synstr(ss), strand1(sd1), strand2(sd2), min_mismatch(mis), max_mismatch(mis), sum_mismatch(mis), bads(bd), flank_bits(fbits), prim(pm)
	{
		++p[loc-1];
	}

	bool inc_hits(size_t idx, size_t suffix_len, unsigned short mis)
	{
		++p[idx-1];

		sum_mismatch += mis;

		if (mis > max_mismatch)
			max_mismatch = mis;
		
		if (mis < min_mismatch)
			min_mismatch = mis;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;

		return true;
	}
	//int prim;
	//string flankstr;

	size_t max_prefix_len, max_suffix_len;

	bool bads;
	size_t flank_bits;
	string synstr;
	char strand1, strand2;
	vector<unsigned short> p;
	unsigned short prim;
	unsigned short min_mismatch;
	unsigned short max_mismatch;
	unsigned short sum_mismatch;	
};

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#endif
//g++ only

/////////////////////Globals/////////////////////////////

// This is a command line parameter: globally scoped because needed
// just about everywhere.
static const int VERBOSE = 0;
static const int DEBUG = 0;
static const int CONJ = 1;

//#define TIME

#define TIME2

// This can be changed, if you have the memory and want to use it:
static const size_t MAX_SEED_WIDTH = 31;

static const size_t LEAST_SIG_BIT = static_cast<size_t>(1);
static const size_t SECOND_LSB = (static_cast<size_t>(1) << 1);
//static const size_t ALL_BITS_ON = static_cast<size_t>(-1);

// Set up for 64-bit ONLY!!
static const size_t SIXTY_FOUR = 64;
static const size_t MOST_SIG_BIT = static_cast<size_t>(0x8000000000000000);

//static const size_t THIRTY_TWO = 32;
////static const size_t ALL_BITS_ON = static_cast<size_t>(-1);
//static const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
//static const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

static const size_t bit_GT_upper = 3;
static const size_t bit_GT_lower = 1;

static const size_t bit_TG_upper = 3;
static const size_t bit_TG_lower = 2;

static const size_t bit_AG_upper = 1;
static const size_t bit_AG_lower = 0;

static const size_t bit_GA_upper = 2;
static const size_t bit_GA_lower = 0;

static const size_t bit_GC_upper = 2;
static const size_t bit_GC_lower = 1;

static const size_t bit_CG_upper = 1;
static const size_t bit_CG_lower = 2;

static const size_t bit_AT_upper = 1;
static const size_t bit_AT_lower = 1;

static const size_t bit_TA_upper = 2;
static const size_t bit_TA_lower = 2;

static const size_t bit_AC_upper = 0;
static const size_t bit_AC_lower = 1;

static const size_t bit_CA_upper = 0;
static const size_t bit_CA_lower = 2;

static const size_t bit_CT_upper = 1;
static const size_t bit_CT_lower = 3;

static const size_t bit_TC_upper = 2;
static const size_t bit_TC_lower = 3;

//CTAC

//CTGC

//GTAT

static const size_t bit_GTAG_upper = bit_GT_upper << 2 | bit_AG_upper;
static const size_t bit_GTAG_lower = bit_GT_lower << 2 | bit_AG_lower;
static const size_t bit_GTAG = bit_GTAG_upper << 4 | bit_GTAG_lower;

static const size_t bit_GATG_upper = bit_GA_upper << 2 | bit_TG_upper;
static const size_t bit_GATG_lower = bit_GA_lower << 2 | bit_TG_lower;
static const size_t bit_GATG = bit_GATG_upper << 4 | bit_GATG_lower;

static const size_t bit_GCAG_upper = bit_GC_upper << 2 | bit_AG_upper;
static const size_t bit_GCAG_lower = bit_GC_lower << 2 | bit_AG_lower;
static const size_t bit_GCAG = bit_GCAG_upper << 4 | bit_GCAG_lower;

static const size_t bit_GACG_upper = bit_GA_upper << 2 | bit_CG_upper;
static const size_t bit_GACG_lower = bit_GA_lower << 2 | bit_CG_lower;
static const size_t bit_GACG = bit_GACG_upper << 4 | bit_GACG_lower;

static const size_t bit_ATAC_upper = bit_AT_upper << 2 | bit_AC_upper;
static const size_t bit_ATAC_lower = bit_AT_lower << 2 | bit_AC_lower;
static const size_t bit_ATAC = bit_ATAC_upper << 4 | bit_ATAC_lower;

static const size_t bit_CATA_upper = bit_CA_upper << 2 | bit_TA_upper;
static const size_t bit_CATA_lower = bit_CA_lower << 2 | bit_TA_lower;
static const size_t bit_CATA = bit_CATA_upper << 4 | bit_CATA_lower;

static const size_t bit_CTAC_upper = bit_CT_upper << 2 | bit_AC_upper;
static const size_t bit_CTAC_lower = bit_CT_lower << 2 | bit_AC_lower;
static const size_t bit_CTAC = bit_CTAC_upper << 4 | bit_CTAC_lower;

static const size_t bit_CATC_upper = bit_CA_upper << 2 | bit_TC_upper;
static const size_t bit_CATC_lower = bit_CA_lower << 2 | bit_TC_lower;
static const size_t bit_CATC = bit_CATC_upper << 4 | bit_CATC_lower;

static const size_t bit_CTGC_upper = bit_CT_upper << 2 | bit_GC_upper;
static const size_t bit_CTGC_lower = bit_CT_lower << 2 | bit_GC_lower;
static const size_t bit_CTGC = bit_CTGC_upper << 4 | bit_CTGC_lower;

static const size_t bit_CGTC_upper = bit_CG_upper << 2 | bit_TC_upper;
static const size_t bit_CGTC_lower = bit_CG_lower << 2 | bit_TC_lower;
static const size_t bit_CGTC = bit_CGTC_upper << 4 | bit_CGTC_lower;

static const size_t bit_GTAT_upper = bit_GT_upper << 2 | bit_AT_upper;
static const size_t bit_GTAT_lower = bit_GT_lower << 2 | bit_AT_lower;
static const size_t bit_GTAT = bit_GTAT_upper << 4 | bit_GTAT_lower;

static const size_t bit_TATG_upper = bit_TA_upper << 2 | bit_TG_upper;
static const size_t bit_TATG_lower = bit_TA_lower << 2 | bit_TG_lower;
static const size_t bit_TATG = bit_TATG_upper << 4 | bit_TATG_lower;

static const size_t LAST_FOUR_BIT = 0xf;
static const size_t LAST_TWO_BIT = 3;
static const size_t LAST_THIRD_FOUTH = LAST_TWO_BIT << 2;

// Assumes 4 nucleotide DNA alphabet
static const size_t alphabet_size = 4;


static const string m_canon_in_intron_dir = "canonical/";

static const string m_canon_exceed_intron_dir = "canonical_exceed/";

static const string m_noncanon_in_intron_dir = "noncanonical/";

static const string m_noncanon_exceed_intron_dir = "noncanonical_exceed/";

static const string m_fusion_dir = "fusion/";

static const string m_temp_dir = "tmp/";

static const string m_original_dir = "original/";

static const string m_filteredbest_dir = "filtered_best/";

static const string m_comparison_dir = "comparison/";

string global_qualstr;

inline size_t
base2int(char c) {
	switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 2;
  case 't' : return 3; 
	}
	return 4;
}

inline bool
isvalid(char c) {
	return (base2int(c) != 4);
}

//vc only
#ifndef S_ISDIR
#define S_ISDIR(m) (((m)& 0170000) == 0040000)
#endif

inline bool
isdir(const char *filename) {
	struct stat buffer;
	stat(filename, &buffer);
	return S_ISDIR(buffer.st_mode);
}

struct RMAPException {
	RMAPException(string m) : message(m) {}
	string what() const {return message;}
	string message;
};


template <class T> string toa(T t) {
	ostringstream s;
	s << t;
	return s.str();
}

////////////////////////////////////////////////////////////////////////
// Code for dealing with the DNA alphabet

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

/////////////////////////SeedInfo//////////////////////
struct SeedInfo {
	SeedInfo(const size_t r, /*const size_t sh,*/ const bool st) : 
readID(r), /*shift(sh), */strand(st) {}
inline bool operator<(const SeedInfo &rhs) const {
	return readID < rhs.readID;
	/*shift < rhs.shift;*/
	}

unsigned readID   : 32;
//unsigned shift  : 7;
unsigned strand : 1;

};


///////////////////////////Cmdline/////////////////////////
struct CmdlineArgus
{
	CmdlineArgus(int argc, const char **argv);

	void showhelp();

	const char* chrom_file;
	const char* tissue_file;
	const char* tissuenames_file;
	const char* filenames_file;
	const char* outfile;
	const char* reads_file;
	const char* ambiguous_file;
	const char* fasta_suffix;
	const char* divide_file;
	const char* result_file;
	const char* gene_group;
	const char* TACs_pos;
	const char* path;
	const char* fixhole_file;
	const char* bowtie_segments_file;
	const char* fixtail_file;
	const char* bwtout;
	const char* fastq_file;
	const char* tmp_dir;
	const char* fusion_file;
	const char* fusion_junc;

	size_t input_buffer_size;
	size_t seed_width;
	size_t read_width;
	size_t max_mismatches;
	size_t gap;
	size_t seed_increment;

	size_t max_splice_dist;
	size_t map_size;
	size_t interal;
	size_t max_readID;
	size_t prefix_match;
	size_t anchro_size;
	size_t num_anchor;
	size_t flankstr;
	double rank;



	size_t num_seg;
	size_t seg_len;

	size_t extend_bits;

	size_t total_mismatch;

	size_t total_fusion_mismatch;

	size_t append_score;

	size_t fix_type;

	size_t threads;

	size_t max_insertion;
};

void
CmdlineArgus::showhelp()
{
	cout <<"Usage: rmapspan "<<endl;
	cout << "\t-o, --output=STRING        Name of output file (default: stdout)"<< endl
		 << "\t-c, --chrom=STRING         FASTA file or dir containing chromosome(s)"<< endl
		 << "\t-u, --tissue=STRING        FASTA file or dir containing tissue(s)"<< endl
		 << "\t-H, --fixhole=STRING       file need to fix hole"<< endl
		 << "\t-G, --segment=STRING       number of segments used"<< endl
		 << "\t-r, --junction=STRING      Name of junction output file"<< endl
		 << "\t-b, --buffsize=INT         size of chromosome input buffer(default : 100000000)"<< endl
		 << "\t-s, --suffix=STRING        suffix of FASTA files (assumes -c indicates dir)(default : fa)"<< endl
		 << "\t-F, --filenames=STRING     file listing names of chromosome files"<< endl
		 << "\t-T, --tissuenames=STRING   file listing names of tissue files"<< endl
		 << "\t-t, --exon=STRING          file listing exon positions(if not provide, use whole region of chromosome)"<< endl
		 << "\t-p, --gap=STRING           extend region of exon(default : 10)"<< endl
		 << "\t-h, --hit=INT              width of hit"<< endl
		 << "\t-w, --width=INT            width of reads"<< endl
		 << "\t-m, --mismatch=INT         maximum allowed mismatches(default : 0)"<< endl
		 << "\t-a, --ambiguous=STRING     file to write names of ambiguously mapped reads"<< endl
		 << "\t-z, --size=INT             initial bucket size of hash_map(default : 1024)"<< endl
		 << "\t-i, --min=INT              minimal splice distance(defalut : 0)"<< endl
		 << "\t-x, --max=INT              maximal splice distance(default : 10000)"<< endl
		 << "\t-y, --prefix=INT           option to turn on prefix match(default : 0)"<< endl
		 << "\t-v, --path=INT             Name of output folder (if not provide, use current path)"<< endl
		 << "\t-f, --flankstr=INT         option to decide if output 0 case of flank string"<< endl
		 << "\t-R, --rank=DOUBLE          option to decide if set threshold of rank when output"<< endl
		 << "\t-A, --anchor=INT             option to decide the number of anchor to be used"<< endl

		 <<endl
		 << "Help options:"<<endl
		 << "-?, --help					  Show this help message"<< endl;
	exit(0);
}

CmdlineArgus::CmdlineArgus(int argc, const char **argv) : chrom_file(0), filenames_file(0), tissuenames_file(0), outfile(0), reads_file(0), tissue_file(0), ambiguous_file(0), fasta_suffix("fa"),
						divide_file(0), result_file(0), gene_group(0), TACs_pos(0), input_buffer_size(100000000), seed_width(0), read_width(0), path(0), flankstr(1), rank(0.0), num_anchor(1),
						max_mismatches(0), gap(10), seed_increment(0), max_splice_dist(10000), map_size(1024), interal(0), max_readID(max_size_t), prefix_match(0),anchro_size(6), fusion_file(0), 
						bowtie_segments_file(0), fixtail_file(0), num_seg(0), seg_len(0), fastq_file(0), extend_bits(0), total_mismatch(2), total_fusion_mismatch(2), append_score(0), threads(1), max_insertion(0)
{
	/****************** COMMAND LINE OPTIONS ********************/

	if (argc == 1 || ((argc == 2) && (strcmp(argv[1], "-?") == 0 || strcmp(argv[1], "--help") == 0)))
	{
		showhelp();
	}

	char* stopstring;

	for (int argi = 1; argi < argc - 1; argi = argi + 2)
	{
		if (strcmp(argv[argi], "-o") == 0)
			outfile = argv[argi + 1];
		else if (strcmp(argv[argi], "-c") == 0)
			chrom_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-u") == 0)
			tissue_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-H") == 0)
			fixhole_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-bwtout") == 0)
			bwtout = argv[argi + 1];
		else if (strcmp(argv[argi], "-d") == 0)//d
			divide_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-r") == 0)
			result_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-g") == 0)//d
			gene_group = argv[argi + 1];
		else if (strcmp(argv[argi], "-s") == 0)
			fasta_suffix = argv[argi + 1];
		else if (strcmp(argv[argi], "-S") == 0)
			bowtie_segments_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-FT") == 0)
			fixtail_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-F") == 0)
			filenames_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-T") == 0)
			tissuenames_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-t") == 0)
			TACs_pos = argv[argi + 1];
		else if (strcmp(argv[argi], "-a") == 0)
			ambiguous_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-b") == 0)
			input_buffer_size = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-h") == 0)
			seed_width = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-p") == 0)
			gap = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-w") == 0)
			read_width = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-m") == 0)
			max_mismatches = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-x") == 0)
			max_splice_dist = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-z") == 0)
			map_size = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-i") == 0)
			interal = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-e") == 0)
			max_readID = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-y") == 0)
			prefix_match = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-v") == 0)
			path = argv[argi + 1];
		else if (strcmp(argv[argi], "-n") == 0)
			anchro_size = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-f") == 0)
			flankstr = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-R") == 0)
			rank = strtod(argv[argi + 1], &stopstring);
		else if (strcmp(argv[argi], "-A") == 0)
			num_anchor = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-G") == 0)
			num_seg = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-L") == 0)
			seg_len = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-FQ") == 0)
			fastq_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-E") == 0)
			extend_bits = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-M") == 0)
			total_mismatch = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-FM") == 0)
			total_fusion_mismatch = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-P") == 0)
			append_score = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-FTP") == 0)
			fix_type = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-tmp") == 0)
			tmp_dir = argv[argi + 1];
		else if (strcmp(argv[argi], "-threads") == 0)
			threads = atoi(argv[argi + 1]);
		else if (strcmp(argv[argi], "-fusion") == 0)
			fusion_file = argv[argi + 1];
		else if (strcmp(argv[argi], "-fusionjunc") == 0)
			fusion_junc = argv[argi + 1];
		else if (strcmp(argv[argi], "-max_insertion") == 0)
			max_insertion = atoi(argv[argi + 1]);
		
	}

	/**********************************************************************/

	if (read_width == 0)
		cerr << "ERROR: must specify read width" << endl;

	seed_increment = 
		static_cast<size_t>(floor(double (read_width/(max_mismatches + 1))));

	if (seed_width == 0)
		seed_width = seed_increment;

	if (num_anchor * anchro_size > (read_width / 2) )
		num_anchor = (read_width / 2) / anchro_size;

	if (max_mismatches == 0)
		seed_width = (read_width + 1) / 2;

	if (VERBOSE)
		cerr << "READ WIDTH:     " << read_width << endl
		<< "SEED WIDTH:     " << seed_width << endl
		<< "SEEDS PER READ: " << read_width/seed_increment << endl
		<< "MAX MISMATCH:   " << max_mismatches << endl;
}

//////////////////////////////////GenomicRegion///////////////////////
class GenomicRegion {
public:
	GenomicRegion(string c, const size_t ls, const size_t le, const size_t rs, const size_t re, 
		string n, size_t sc, char str) :
	chrom(c), lstart(ls), lend(le), rstart(rs), rend(re), name(n), score(sc), strand(str) {}
	string tostring() const;
private:
	string chrom, name;
	size_t lstart, lend, rstart, rend;
	unsigned score;
	char strand;
};

ostream& 
operator<<(ostream& s, const GenomicRegion& region) {
  return s << region.tostring();
}

string 
GenomicRegion::tostring() const {
  static const char tab = '\t'; 
  ostringstream s;
  s << chrom << tab << lstart << tab << lend << tab << rstart <<tab << rend << tab;
  if (!name.empty())
    s << tab << name << tab << score << tab << strand;
  return s.str();
}


////////////////////////////////////////////////////////////////////////
// Stuff dealing with FASTA format sequence files

bool
is_valid_filename(const string name, const string& filename_suffix) {
	const string suffix(name.substr(name.find_last_of(".") + 1));
	return (suffix == filename_suffix);
}

string 
path_join(const string& a, const string& b) {
	return a + "/" + b;
}

void 
read_dir(const string& dirname, string filename_suffix,
		 vector<string> &filenames) {
			 DIR *dir;
			 if (!(dir = opendir(dirname.c_str())))
				 throw RMAPException("could not open directory: " + dirname);

			 errno = 0;
			 struct dirent *ent;
			 while ((ent = readdir(dir))) {
				 if (is_valid_filename(ent->d_name, filename_suffix))
					 filenames.push_back(path_join(dirname, string(ent->d_name)));
				 errno = 0;
			 }
			 // check for some errors
			 if (errno)
				 throw RMAPException("error reading directory: " + dirname);
			 if (filenames.empty())
				 throw RMAPException("no valid files found in: " + dirname);
}

//void
//read_fasta_file(const char *filename, vector<string> &names, vector<string> &sequences)
//{
//
//	size_t input_buffer_size = 1000000;
//
//	ifstream in(filename);
//	if (!in) {
//		cout << 1<<endl;
//		getchar();
//		throw RMAPException("cannot open input file " + string(filename));
//	}
//	string s, name;
//	bool first_line = true;
//	while (!in.eof()) {
//		char buffer[1000000 + 1];
//		in.getline(buffer, input_buffer_size);
//		if (in.gcount() == static_cast<int>(input_buffer_size))
//			throw RMAPException("Line in " + name + "\nexceeds max length: " +
//			toa(input_buffer_size));
//		// correct for dos carriage returns before newlines
//		if (buffer[strlen(buffer) - 1] == '\r')
//			buffer[strlen(buffer) - 1] = '\0';
//		if (buffer[0] == '>') {
//			if (first_line == false && s.length() > 0) {
//				names.push_back(name);
//				sequences.push_back(s);
//			}
//			else first_line = false;
//			name = buffer;
//			name = name.substr(name.find_first_not_of("> "));
//			s = "";
//		}
//		else s += buffer;
//	}
//	if (!first_line && s.length() > 0) {
//		names.push_back(name);
//		sequences.push_back(s);
//	}
//}

void
read_filename_file(const char *filename, vector<string> &filenames) {

	static const size_t input_buffer_size = 1000000;

	ifstream in(filename);
	if (!in)
		throw RMAPException("cannot open input file " + string(filename));
	while (!in.eof()) {
		char buffer[input_buffer_size + 1];
		in.getline(buffer, input_buffer_size);
		if (in.gcount() == static_cast<int>(input_buffer_size))
			throw RMAPException("Line in " + string(filename) +
			"\nexceeds max length: " +
			toa(input_buffer_size));
		filenames.push_back(buffer);
		in.peek();
	}
}

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

string
i2mer(size_t n, size_t index) {
	char int2base[] = {'A', 'C', 'G', 'T', 'N'};
	string s(n, ' ');
	do {
		--n;
		s[n] = int2base[index % alphabet_size];
		index /= alphabet_size;
	} while (n > 0);
	return s;
}
//
//bool
//check_seed(const string &read, const bool posstrand,
//		   const size_t seed_width, const size_t mask,
//		   const size_t read_id, const size_t offset,
//		   hash_map<size_t, vector<SeedInfo> > &seed_hash) {
//
//			   size_t seed_key = 0;
//			   bool good_hit = true;
//			   for (size_t k = 0; k < seed_width && good_hit; ++k)
//				   if (isvalid(read[offset + k]))
//					   seed_key = (seed_key << 2) + base2int(read[offset + k]);
//				   else good_hit = false;
//				   if (good_hit) {
//					   seed_key &= mask;
//					   seed_hash[seed_key].push_back(SeedInfo(read_id, offset, posstrand));
//					   return true;
//				   }
//				   return false;
//}

//jinze
//void
//get_read_matches(const vector<string> &reads, const size_t seed_width, 
//				 const size_t seed_increment,
//				 hash_map<size_t, vector<SeedInfo> > &seed_hash) {
//
//					 const size_t mask = (LEAST_SIG_BIT << 2*seed_width) - 1;
//
//					 for (size_t i = 0; i < reads.size(); ++i) {
//						 string read = reads[i];
//
//						 //insert prefix and suffix seeds to the hashtable
//						 //offset will be able to tell whether it is a prefix or a suffix. 
//						 check_seed(read, true, seed_width, mask, i, 0, seed_hash); //prefix
//						 check_seed(read, true, seed_width, mask, i, read.length()-seed_width, seed_hash); //suffix
//
//						 //complementary, reversed, the other strand. Just consider it as a new set of tags.
//						 read = revcomp(read);
//						 check_seed(read, false, seed_width, mask, i, 0, seed_hash); //prefix
//						 check_seed(read, false, seed_width, mask, i, read.length()-seed_width, seed_hash); //suffix
//					 }
//
//					 //for each key entry, sort the reads based on the increasing order of their IDs.  	
//					 for (hash_map<size_t, vector<SeedInfo> >::iterator i = seed_hash.begin(); 
//						 i != seed_hash.end(); ++i) 
//						 if (i->second.size() > 0) {
//							 sort(i->second.begin(), i->second.end());
//						 }
//}

void show_reads_and_seeds(const vector<string> &reads, const size_t seed_width, 
						  hash_map<size_t, vector<SeedInfo> > &seed_hash) {

							  cerr << endl << "Total reads = " << reads.size() << endl;
							  for (size_t ri = 0; ri < reads.size(); ri++) {
								  cerr << ri << ":  " << reads[ri] << endl;
							  }

							  cerr << endl << "Seeds and mappings" << endl;

							  string twobit2base = "ACGT";
							  //char seedstr[seed_width+1];
							  char seedstr[35+1];
							  seedstr[seed_width] = '\0';

							  for (hash_map<size_t, vector<SeedInfo> >::iterator i = seed_hash.begin(); 
								  i != seed_hash.end(); ++i) {

									  for (size_t j = 0; j < i->second.size(); j++) {
										  size_t seed = i->first;
										  SeedInfo seedinfo = i->second[j];

										  for (size_t k = seed_width - 1; k >= 0; k--) {
											  seedstr[k] = twobit2base[seed & 3];
											  seed = seed >> 2;
										  }

										  size_t readid = seedinfo.readID;
										  char strand = seedinfo.strand ? '+' : '-';
										  //size_t offset = seedinfo.shift;
										  string readstr = reads[readid];

										  cerr << "Seed " << setw(seed_width) << seedstr << " -> Read " << setw(3) << readid
											  << " (" << setw(2) /*<< offset */<< "," << strand << ") = " << readstr <<endl;
									  }
							  }
}

size_t get_filesize(string filename)
{
	ifstream f(filename.c_str());
	if (!f.good()) {return 0;}
	size_t begin_pos = f.tellg();
	f.seekg(0, ios_base::end);
	size_t end_pos = f.tellg();
	f.close();
	return end_pos - begin_pos;
}

void
readchrom(const char* filename, string& longseq)
{
	size_t size;  

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);

		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);

			if (line.length() == 0)
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";
}

size_t
reverse_chrom(const char* infile, const char* ourfile)
{
	string longseq;
	size_t size;  

	ifstream longfile(infile);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	string skipline;
	if (longfile.is_open())
	{
		getline(longfile,skipline);
		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);

			if (line.length() == 0)
				continue;

			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";

	size = longseq.length();

	ifstream ifs_exist(ourfile);

	if (!ifs_exist.is_open())
	{
		ofstream ofs(ourfile);

		ofs<<skipline<<endl;
		//for (size_t i = 0; i < longseq.length() / 2; ++i)
		//{
		//	char temp;
		//	temp = longseq[i];
		//	longseq[i] = longseq[longseq.length() - 1 - i];
		//	longseq[longseq.length() - 1 - i] = temp;
		//}
		reverse(longseq.begin(), longseq.end());
		ofs.write(longseq.c_str(), longseq.size());
		ofs.close();
	}

	return size;
}

//void read_chrom_filenames(vector<string>& chrom_files, size_t& total_chromosomes_size, vector<size_t>& chrom_sizes,
//						  const char* filenames_file, const char* chrom_file, const char* fasta_suffix)
//{
//	if (filenames_file)
//		read_filename_file(filenames_file, chrom_files);
//	else if (isdir(chrom_file))
//		read_dir(chrom_file, fasta_suffix, chrom_files);
//	else 
//		chrom_files.push_back(chrom_file);
//
//	//size_t total_chromosomes_size = 0;
//	//vector<size_t> chrom_sizes;
//	if (VERBOSE) 
//	{
//		cerr << endl << "chromosome files found:" << endl;
//		for (vector<string>::const_iterator i = chrom_files.begin();i != chrom_files.end(); ++i) 
//		{
//			chrom_sizes.push_back(get_filesize(*i));
//			total_chromosomes_size += chrom_sizes.back();
//			cerr << *i << endl;
//		}
//	}
//}

//void
//clean_reads(const size_t max_diffs, 
//			const size_t read_width,
//			vector<string> &reads, vector<string> &read_names) {
//				size_t good = 0;
//				for (size_t i = 0; i < reads.size(); ++i) {
//
//					if (reads[i].length() != read_width) {
//						if (reads[i].length() < read_width)
//							throw RMAPException("Incorrect read width");
//						else
//							reads[i] = reads[i].substr(0, read_width);
//					}
//
//					if (count(reads[i].begin(), reads[i].end(), 'N') <= 
//						static_cast<int>(max_diffs)) {
//							reads[good].swap(reads[i]);
//							read_names[good].swap(read_names[i]);
//							++good;
//					}
//				}
//				reads.resize(good);
//				read_names.resize(good);
//}

//void
//sites_to_regions(const vector<string> &chrom, const vector<string> &reads,
//				 const vector<string> &read_names, const vector<ReadInfo> &bests, 
//				 const size_t max_mismatches, vector<GenomicRegion> &hits) {
//					 for (size_t i = 0; i < bests.size(); ++i)
//						 if (bests[i].unique && bests[i].score <= max_mismatches)
//							 hits.push_back(GenomicRegion(chrom[bests[i].chrom], bests[i].lsite_start, bests[i].lsite_end,
//							 bests[i].rsite_start, bests[i].rsite_end,read_names[i], bests[i].score, (bests[i].strand) ? '+' : '-'));
//}

//void
//write_non_uniques(string filename, const vector<string> &read_names,
//				  const vector<ReadInfo> &bests) {
//					  ofstream out(filename.c_str());
//					  for (size_t i = 0; i < read_names.size(); ++i)
//						  if (!bests[i].unique)
//							  out << read_names[i] << endl;
//					  out.close();
//}

///////////////////////////////ChromReader//////////////////////
typedef vector<char> VCHAR;
typedef VCHAR::const_iterator VCHAR_CITER;

class ChromReader {
public:
	ChromReader() {}
	ChromReader(string fn, size_t bc) :
	buffer(VCHAR(bc + SIXTY_FOUR)), filename(fn), 
		buffer_capacity(bc), buffer_size(0), current_offset(0),
		finished_reading(false)
	{
		skip_first_line();
	}
	bool read();
	void get_pointers(VCHAR_CITER &a, VCHAR_CITER &b) 
	{
			a = buffer.begin();
			b = a + buffer_size;
	}
	void skip_first_line();
	size_t get_current_offset(){return current_offset;}

	//return one base in size_t each time
	size_t ScanBase();

	void Clear()
	{
		buffer.clear();
		filename.clear();
		buffer_capacity = 0;
		buffer_size = 0;
		current_offset = 0;
		finished_reading = false;
	}

private:
	VCHAR buffer;
	string filename;
	size_t buffer_capacity;
	size_t buffer_size;
	size_t current_offset;
	bool finished_reading;

	//start and end point for ScanBase
	VCHAR_CITER m_current_base;
	VCHAR_CITER m_final_base_in_buffer;
};

bool 
ChromReader::read() {

	if (finished_reading) 
		return false;

	std::ifstream in(filename.c_str());
	if (!in)
		return false;

	in.seekg(current_offset);
	in.read(&buffer.front(), buffer_capacity);
	buffer_size = in.gcount();
	in.close();

	// This code below is to take care of the case where a read is at
	// the very end of the chromosome, and therefore "dummy" bases need
	// to be pushed through the frame.
	if (buffer_size < buffer_capacity) {
		buffer.insert(buffer.begin() + buffer_size, SIXTY_FOUR, 'N');
		buffer_size += SIXTY_FOUR;
		finished_reading = true;
	}

	current_offset += buffer_size;

	//set start and end point
	get_pointers(m_current_base, m_final_base_in_buffer);

	return true;
}

void
ChromReader::skip_first_line() {
	ifstream fin(filename.c_str());
	char c = '\0';
	do {
		fin.get(c);
		++current_offset;
	}
	while (c != '\n' && c != '\0');
	fin.close();
}

inline size_t
ChromReader::ScanBase()
{
	while (m_current_base != m_final_base_in_buffer)
	{
		if (*m_current_base != '\n' && *m_current_base != '\t' && *m_current_base!= '\r')
			return base2int(*(m_current_base++));
		++m_current_base;
	}
	return -1;
}

struct Kmer{
	bool bad;
	unsigned kmer;
	Kmer(bool bd, unsigned km) : bad(bd), kmer(km) {}
};


////////////////////////WordPair////////////////////////////
// WordPair is developed from Rmap package of Cold Spring Harbor Laboratory 
/*
 * Modified to discover tags covering exon conjuctions. 
 *
 * Copyright (C) 2007 Cold Spring Harbor Laboratory
 *                    Andrew D Smith, Zhenyu Xuan and Michael Q Zhang
 * Author: Andrew D. Smith
 */
struct WordPair {
	WordPair(const string &s);
	WordPair() : upper(0), lower(0), bads(0) {}

	string tostring(size_t mask) const;
	string tostring2(size_t mask) const;
	string tostring3(size_t mask, size_t bitsnum) const;
	char get_char(size_t mask, size_t pos) const;

	size_t score(const WordPair &other, size_t mask) const;

	size_t score(const WordPair &other, size_t mask, size_t& rbit) const;

	size_t score_hmer(const WordPair &other, size_t mask) const;
	size_t ps_score(const WordPair &other,  const WordPair & wp, const size_t mid_width,  const size_t seed_length, const size_t mask, size_t & loc) ;
	void update_key(const size_t update_bit, const size_t mask, size_t &key) const;

	void shift(const size_t i);
	void shift_reserve(const WordPair &other);

	void shift_reserve(const WordPair &other, const size_t reserve_bit);
	void combine(const WordPair &other, size_t shift, WordPair &wp) const;
	void right_shift (const size_t i);
	void left_shift(const size_t i);
	void clear()
	{
		upper = 0;
		lower = 0;
		bads = 0; 
	}
	void ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp);

	void shift_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const size_t leftshift, const size_t rightshift, WordPair suffix_wp, WordPair &wp);
	//the following function is used to Duplicate middle part in a word
	void duplicate(const size_t prefix_mask, const size_t suffix_mask, const size_t mid_buff_mask, const size_t leftshift, const size_t rightshift, WordPair &wp);

	void duplicate_self(const size_t leftshift, WordPair &wp)
	{
		wp.upper = ((upper << leftshift)  + upper) ;
		wp.lower = ((lower << leftshift)  + lower); 

		//why?
		wp.bads =  ((bads << leftshift) + bads);  
	}

		//	vector<Kmer> kmers;

		//for (size_t i = 0; i < m_num_anchors; ++i)
		//{
		//	unsigned kmer = 0;

		//	bool good_hit = true;

		//	size_t offset = readstr.length() - m_anchor_width - (i * m_kmer_step);

		//	for (size_t k = 0; k < m_anchor_width && good_hit; ++k)
		//	{
		//		if (isvalid(readstr[offset + k]))
		//			kmer = (kmer << 2) + (unsigned)base2int(readstr[offset + k]);
		//		else good_hit = false;
		//	}

		//	kmers.push_back(Kmer(!good_hit, kmer));
		//}
		//return kmers;
	Kmer get_kmer(size_t st, size_t end)
	{
		unsigned kmer = 0;

		bool good_hit = true;

		int st_1 = (int)st - 1;

		for (int i = (int)end - 1; i >= st_1; --i)
		{
			size_t cur_bit = LEAST_SIG_BIT << i;
			if ((cur_bit) && bads)
			{
				good_hit = false;
				break;
			}
			else
			{
				kmer = (kmer << 2) + (((cur_bit & upper) != 0) << 1) + ((cur_bit & lower) != 0);
			}				
		}

		return Kmer(good_hit, kmer);
	}

	void get_prefix(const size_t shift, WordPair &prefix_wp);

	void get_suffix(const size_t mask, WordPair &suffix_wp);

	size_t upper;
	size_t lower;
	size_t bads;

	static string bits2string(size_t mask, size_t bits);
	static inline size_t get_upper(const size_t i) {return i > 1;}
	static inline size_t get_lower(const size_t i) {return (i % 2);}
	static inline size_t get_bads(char c) {return !isvalid(c);}
};

inline string 
bits2string2(size_t mask, size_t bits) {
	string s;
	size_t selector = MOST_SIG_BIT;
	for (size_t i = 0; i < SIXTY_FOUR; ++i) {
		s += (selector & bits & mask) ? '1' : '0';
		selector >>= 1;
	}
	return s;
}

inline void
WordPair::combine(const WordPair &other, 
				  const size_t shift, WordPair &wp) const {
					  wp.upper = upper >> shift;
					  wp.lower = lower >> shift;
					  wp.bads  = bads >> shift;
					  if (shift != 0) {
						  const size_t other_shift = (SIXTY_FOUR - shift);
						  wp.upper |= (other.upper << other_shift);
						  wp.lower |= (other.lower << other_shift);
						  wp.bads  |= (other.bads  << other_shift);
					  }  
}

inline void
WordPair::shift(const size_t i) {
	//remove one base on the left and add a new base on the right. 
	upper = ((upper << 1) + (i > 1));
	lower = ((lower << 1) + (i % 2));
}

inline void
WordPair::shift_reserve(const WordPair &other) {
	//the word left shift one bit
	//add another bit which is the most_sig_bit of the other
	//itself and other creates and 64 + 64 buffer on the genome, which is used to align reads within it 

	upper = (upper << 1) + ((other.upper & MOST_SIG_BIT) != 0);
	lower = (lower << 1) + ((other.lower & MOST_SIG_BIT) != 0);
}

inline void
WordPair::shift_reserve(const WordPair &other, const size_t reserve_bit) {
	//the word left shift one bit
	//add another bit which is the most_sig_bit of the other
	//itself and other creates and 64 + 64 buffer on the genome, which is used to align reads within it 

	upper = (upper << 1) + ((other.upper & reserve_bit) != 0);
	lower = (lower << 1) + ((other.lower & reserve_bit) != 0);
}

//return mismatches between two word pairs. 
inline size_t
WordPair::score(const WordPair &other, size_t mask) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
WordPair::score(const WordPair &other, size_t mask, size_t& rbit) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	rbit = bits;

	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}


inline pair<size_t, size_t>
score_string(const string& s1, const string& s2)
{
	pair<size_t, size_t> mis_comp_bits;
	if (s1.length() != s2.length())
	{
		cout << "different length of two strings"<<endl <<s1 <<endl << s2 <<endl ;
		getchar();

		mis_comp_bits = make_pair(-1, -1);
		return mis_comp_bits;
	}

	size_t mask = ALL_BITS_ON >> (SIXTY_FOUR - s2.length());

	WordPair w1(s1), w2(s2);

	register size_t bits = ((w1.upper ^ w2.upper) | 
		(w1.lower ^ w2.lower) | w1.bads | w2.bads) & mask;

	mis_comp_bits.second = bits;
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word

	mis_comp_bits.first = ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF); 
	return mis_comp_bits; 
}

inline size_t
WordPair::score_hmer(const WordPair &other, size_t mask) const {
	register size_t bits = ((other.upper ^ upper) | 
		(other.lower ^ lower) | other.bads | bads) & mask;

	//assume anchor width is less than 8
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	return bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	////cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	//// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	//// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	//return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
socreBits(register size_t bits, size_t anchor_width) 
{
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
	//  cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
	//cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;  
	// do you want to watch it making the sums?  This is the hypercube summation alg, wow, the reult from the score is not right, debugging
	// at this point right here you would have the sum of the top 32 bits and the bottom 32 bits, each sitting in their half of the 64 bit word
	return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);  
}

inline size_t
WordPair::ps_score(const WordPair &other, const WordPair & wp,  const size_t mid_width, const size_t seed_width, const size_t mask_all, size_t& loc) {

	size_t mask_midright = (ALL_BITS_ON >> (SIXTY_FOUR - mid_width))  << seed_width; 
	size_t mask_ps = mask_all - mask_midright; 
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = wp.score (other, mask_ps);

	bool debug = true;
	if (debug) {
		cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, mask_midright) << endl;
		cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, mask_ps) << endl;
		cerr << "suffix_wp" << endl << other.tostring(ALL_BITS_ON) << endl;
		cerr << "prefix_wp" << endl << wp.tostring(ALL_BITS_ON) << endl;
	}
	//   	cerr >> "mask_midright" >> endl >> bits2string2(ALL_BITS_ON, mask_midright);

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((other.upper ^ wp.upper) | 
		(other.lower ^ wp.lower) | other.bads | wp.bads) & mask_all;
	if (debug) {
		cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
	}
	//the following are two pointers each point to the bit to be turned on and turned off
	size_t selector1 = LEAST_SIG_BIT << seed_width ; //to be turned on 
	size_t selector2 = selector1 << mid_width;//to be turned off

	loc = mid_width; // loc is the index to the left
	size_t mins = s;

	if (debug) cerr << "score :" << mins << endl;
	//check out the accumulation score of mismatches when the selectors are moving...
	for (size_t i = 1; i <= mid_width; ++i) {

		if (debug )
		{
			cerr << "loop " << i << "score " << s << endl;
			cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
			cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
		}



		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;     	
		if (mins > s){
			mins = s;
			loc = mid_width - i;
		}
		selector1 <<= 1;
		selector2 <<= 1;

	}
	if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;

	return mins;
}

// this basically merges upper and lower together into one value. 
// the original Key was shifted left 2 bits. 

inline void
WordPair::update_key(const size_t bit, const size_t mask, size_t &key) const {
	key = (((((key << 1) + ((upper & bit) != 0)) << 1) + 
		((lower & bit) != 0)) & mask);
}

ostream& 
operator<<(ostream& s, const WordPair& wp) {
	return s << wp.tostring2(static_cast<size_t>(-1));
}

char
WordPair::get_char(size_t mask, size_t pos) const {
	// 00 -> A, 01 -> C, 10 -> G, 11 -> T
	const size_t selector = (LEAST_SIG_BIT << (pos - 1));

	if ((mask & bads) & selector)
		return 'N';

	const bool upper_bit = ((mask & upper) & selector);
	const bool lower_bit = ((mask & lower) & selector);
	if (upper_bit) return (lower_bit) ? 'T' : 'G';
	else return (lower_bit) ? 'C' : 'A';
}

WordPair::WordPair(const string &s) : upper(0), lower(0), bads(0) {
	string::const_iterator i = s.begin();
	const string::const_iterator limit = s.end();
	while (i != limit) {
		const char c = base2int(*i) & static_cast<size_t>(3);
		upper = ((upper << 1) + get_upper(c));
		lower = ((lower << 1) + get_lower(c));
		bads  = ((bads << 1) + get_bads(*i));
		++i;
	}
}

void
WordPair::get_prefix(const size_t shift, WordPair &prefix_wp)
{
	prefix_wp.upper = upper>>shift;
	prefix_wp.lower = lower>>shift;
	prefix_wp.bads = bads>>shift;
}

void
WordPair::get_suffix(const size_t mask, WordPair &suffix_wp)
{
	suffix_wp.upper = upper & mask;
	suffix_wp.lower = lower & mask;
	suffix_wp.bads = bads & mask;
}

//right most bits corresponds to the combined word

void WordPair::left_shift(const size_t i)
{
	upper = upper << i; //only the left seed_width + mid_buff_size is useful though
	lower = lower << i;
	bads =  bads << i; 
}

void WordPair::right_shift (const size_t i){
	upper = upper >> i; //only the left seed_width + mid_buff_size is useful though
	lower = lower >> i;
	bads =  bads >> i; 
}

inline
void WordPair::ps_combine(const size_t prefix_mask, const size_t suffix_mask, const size_t big_buff_mask, const WordPair& suffix_wp, WordPair &wp)
{
	wp.upper = ((upper &  prefix_mask) + (suffix_wp.upper & suffix_mask) ) & big_buff_mask;
	wp.lower = ((lower &  prefix_mask) + (suffix_wp.lower & suffix_mask) ) & big_buff_mask;
	wp.bads = ((bads &  prefix_mask) + (suffix_wp.bads & suffix_mask) ) & big_buff_mask;
}

inline
void WordPair::duplicate(const size_t prefix_mask, const size_t suffix_mask, const size_t mid_buff_mask, const size_t leftshift, const size_t rightshift, WordPair &wp)
{
	wp.upper = (((upper << leftshift) & prefix_mask)  + ((upper >> rightshift) & suffix_mask)) & mid_buff_mask;
	wp.lower = (((lower << leftshift) & prefix_mask)  + ((lower >> rightshift) & suffix_mask)) & mid_buff_mask; 

	//why?
	wp.bads = (((bads << leftshift) & prefix_mask)  + ((bads >> rightshift) & suffix_mask)) & mid_buff_mask;  
}

string
WordPair::bits2string(size_t mask, size_t bits) {
	string s;
	size_t selector = MOST_SIG_BIT;
	for (size_t i = 0; i < SIXTY_FOUR; ++i) {
		s += (selector & bits & mask) ? '1' : '0';
		selector >>= 1;
	}
	return s;
}

string
WordPair::tostring(size_t mask) const {
	const string s(bits2string(mask, upper) + "\n" +
		bits2string(mask, lower) + "\n" + 
		bits2string(mask, bads) + "\n");
	string seq;
	for (size_t i = SIXTY_FOUR; i > 0; --i)
		seq += get_char(mask, i);
	return s + seq;
}

string
WordPair::tostring2(size_t mask) const {
	string seq;
	for (size_t i = SIXTY_FOUR; i > 0; --i)
		seq += get_char(mask, i);
	return seq;
}

string
WordPair::tostring3(size_t mask, size_t bitsnum) const {
	string seq;
	for (size_t i = bitsnum; i > 0; --i)
		seq += get_char(mask, i);
	return seq;
}

//////////////////////////SeedHash////////////////////

typedef vector<SeedInfo> SEED_VEC;
typedef SEED_VEC::const_iterator SEED_VEC_CONST_ITER;
typedef SEED_VEC::iterator SEED_VEC_ITER;

typedef hash_map<size_t, SEED_VEC > SEEDS;
typedef SEEDS::const_iterator SEEDS_CONST_ITER;
typedef SEEDS::iterator SEEDS_ITER;

class SeedHash{
public:
	SeedHash(size_t seed_length, size_t read_width, size_t max_mismatches);

	void Get_read_matches(const vector<string> &reads);

	void Get_read_matches_rev(const vector<string> &reads);

	bool Check_seed(const string &read, const bool posstrand, const size_t mask, const size_t read_id, const size_t offset);

	bool Check_seedII(const string &read, const bool posstrand, const size_t mask, const size_t read_id, const size_t offset);

	void InsertTagbyPrefix(const string& tag, const size_t mask, const size_t read_id);

	void InsertTagbyPrefixSuffix(const string& tag, const size_t mask, const size_t read_id);

	void InsertTagbykMer(const string& tag, const size_t mask, const size_t read_id);

	SEEDS_CONST_ITER GetTags(size_t seed_key);

	SEEDS_CONST_ITER GetTagsII(size_t seed_key);

	SEEDS_ITER End()
	{
		return m_seed_hash.end();
	}

	SEEDS_ITER EndII()
	{
		return m_seed_hashII.end();
	}

	void clear()
	{
		m_seed_hash.clear();
		m_seed_hashII.clear();
	}

private:
	SEEDS m_seed_hash;

	SEEDS m_seed_hashII;

	size_t m_seed_length;

	size_t m_seed_increment;

	size_t m_read_width;

	size_t m_max_mismatches;
};

SeedHash::SeedHash(size_t seed_length, size_t read_width, size_t max_mismatches) : m_seed_length(12), m_read_width(read_width), m_max_mismatches(max_mismatches)
{ 
	//if (seed_length < 12)
	m_seed_length = seed_length;

	m_seed_increment = static_cast<size_t>(floor(double (m_read_width/(m_max_mismatches + 1))));
}

inline bool
SeedHash::Check_seed(const string &read, const bool posstrand,
					 const size_t mask, const size_t read_id, const size_t offset)
{
	size_t seed_key = 0;
	bool good_hit = true;

	for (size_t k = 0; k < m_seed_length && good_hit; ++k)
	{
		if (isvalid(read[offset + k]))
			seed_key = (seed_key << 2) + base2int(read[offset + k]);
		else good_hit = false;
	}

	if (good_hit) 
	{
		seed_key &= mask;
		m_seed_hash[seed_key].push_back(SeedInfo(read_id, /*offset, */posstrand));
		return true;
	}
	return false;
}

inline bool
SeedHash::Check_seedII(const string &read, const bool posstrand,
					 const size_t mask, const size_t read_id, const size_t offset)
{
	size_t seed_key = 0;
	bool good_hit = true;

	for (size_t k = 0; k < m_seed_length && good_hit; ++k)
	{
		if (isvalid(read[offset + k]))
			seed_key = (seed_key << 2) + base2int(read[offset + k]);
		else good_hit = false;
	}

	if (good_hit) 
	{
		seed_key &= mask;
		m_seed_hashII[seed_key].push_back(SeedInfo(read_id, /*offset, */posstrand));
		return true;
	}
	return false;
}

void
SeedHash::Get_read_matches(const vector<string> &reads) 
{
	const size_t mask = (LEAST_SIG_BIT << 2*m_seed_length) - 1;

	for (size_t i = 0; i < reads.size(); ++i) 
	{
		string read = reads[i];

		//insert prefix and suffix seeds to the hashtable
		//offset will be able to tell whether it is a prefix or a suffix. 
		Check_seed(read, true, mask, i, 0); //prefix
		if (m_max_mismatches)
			Check_seedII(read, true, mask, i, ((read.length() + 1) / 2) - m_seed_length); //second seed

		//complementary, reversed, the other strand. Just consider it as a new set of tags.
		read = revcomp(read);
		Check_seed(read, false, mask, i, 0); //prefix
		if (m_max_mismatches)
			Check_seedII(read, false, mask, i, ((read.length() + 1) / 2) - m_seed_length); //second seed
	}

	//for each key entry, sort the reads based on the increasing order of their IDs.
	//need sort here?
	for (hash_map<size_t, vector<SeedInfo> >::iterator i = m_seed_hash.begin(); i != m_seed_hash.end(); ++i) 
	{
		if (i->second.size() > 0) 
		{
			sort(i->second.begin(), i->second.end());
		}
	}

	for (hash_map<size_t, vector<SeedInfo> >::iterator i = m_seed_hashII.begin(); i != m_seed_hashII.end(); ++i) 
	{
		if (i->second.size() > 0) 
		{
			sort(i->second.begin(), i->second.end());
		}
	}
}

void
SeedHash::Get_read_matches_rev(const vector<string> &reads) 
{
	const size_t mask = (LEAST_SIG_BIT << 2*m_seed_length) - 1;

	for (size_t i = 0; i < reads.size(); ++i) 
	{
		string read = reads[i];
		reverse(read.begin(), read.end());

		//insert prefix and suffix seeds to the hashtable
		//offset will be able to tell whether it is a prefix or a suffix. 
		Check_seed(read, true, mask, i, 0); //prefix
		if (m_max_mismatches)
			Check_seedII(read, true, mask, i, ((read.length() + 1) / 2) - m_seed_length); //second seed

		//complementary, reversed, the other strand. Just consider it as a new set of tags.
		read = revcomp(read);
		Check_seed(read, false, mask, i, 0); //prefix
		if (m_max_mismatches)
			Check_seedII(read, false, mask, i, ((read.length() + 1) / 2) - m_seed_length); //second seed
	}

	//for each key entry, sort the reads based on the increasing order of their IDs.
	//need sort here?
	for (hash_map<size_t, vector<SeedInfo> >::iterator i = m_seed_hash.begin(); i != m_seed_hash.end(); ++i) 
	{
		if (i->second.size() > 0) 
		{
			sort(i->second.begin(), i->second.end());
		}
	}

	for (hash_map<size_t, vector<SeedInfo> >::iterator i = m_seed_hashII.begin(); i != m_seed_hashII.end(); ++i) 
	{
		if (i->second.size() > 0) 
		{
			sort(i->second.begin(), i->second.end());
		}
	}
}

void
SeedHash::InsertTagbyPrefix(const string& tag, const size_t mask, const size_t read_id)
{
	string read = tag;
	Check_seed(read, true, mask, read_id, 0); //prefix

	read = revcomp(read);
	Check_seed(read, false, mask, read_id, 0); //prefix
}

void
SeedHash::InsertTagbyPrefixSuffix(const string& tag, const size_t mask, const size_t read_id)
{
	string read = tag;
	Check_seed(read, true, mask, read_id, 0); //prefix
	Check_seed(read, true, mask, read_id, read.length() - m_seed_length); //suffix

	read = revcomp(read);
	Check_seed(read, false, mask, read_id, 0); //prefix
	Check_seed(read, false, mask, read_id, read.length() - m_seed_length); //suffix
}

void
SeedHash::InsertTagbykMer(const string& tag, const size_t mask, const size_t read_id)
{
	string read = tag;
	const size_t limit = read.length() - m_seed_length;

	size_t j = 0;
	for (j = 0; j <= limit; )
	{
		if (Check_seed(read, true, mask, read_id, j))
			j += m_seed_increment;
		else ++j;
		if (j != read.length())
			Check_seed(read, true, mask, read_id, limit);

		read = revcomp(read);
		j = 0;
		for (j = 0; j <= limit;)
		{
			if (Check_seed(read, false, mask, read_id, j))
				j += m_seed_increment;
			else ++j;
			if (j != read.length())
				Check_seed(read, false, mask, read_id, limit);
		}
	}
}

SEEDS_CONST_ITER 
SeedHash::GetTags(size_t seed_key)
{
	return m_seed_hash.find(seed_key);
}

SEEDS_CONST_ITER 
SeedHash::GetTagsII(size_t seed_key)
{
	return m_seed_hashII.find(seed_key);
}

struct OpenSeed
{
	OpenSeed() : offset(-1), strand(-1),index(-1), bef_offest(-1) {}
	OpenSeed(const size_t& oft, const size_t& idx, const char& strd, size_t boft = -1) : offset(oft), index(idx), strand(strd), bef_offest(boft) {}

	string Tostring()
	{
		char chr[1000];
		sprintf(chr, "%d\t%d\t%c", (int)offset, (int)index, strand);

		string str(chr);
		return str;
	}
	size_t offset;
	//index in the vector
	size_t index;

	size_t bef_offest;
	char strand;
};

//typedef hash_map<unsigned, vector<OpenSeed> > OPEN_READ;
//typedef OPEN_READ::const_iterator OPEN_READ_CITER;
//typedef OPEN_READ::iterator OPEN_READ_ITER;
//
//class OpenRead{
//public:
//	void insert(size_t read_id, size_t chrom_offset);
//
//	void AddRead(const unsigned& readID, const size_t& offset, const unsigned& strand)
//	{
//		m_open_reads[readID].push_back(OpenSeed(offset, strand));
//	}
//
//	OPEN_READ_ITER GetRead(const unsigned& readID)
//	{
//		return m_open_reads.find(readID);
//	}
//
//	OPEN_READ_ITER End()
//	{
//		return m_open_reads.end();
//	}
//
//	void Clear()
//	{
//		m_open_reads.clear();
//	}
//
//private:
//	//to be changed to circular array
//	OPEN_READ m_open_reads;
//};

/////////////////////////////////////////OpenReads/////////////////////////////

//using __gnu_cxx::hash_multimap; //g++ only

//typedef hash_multimap<unsigned int, OpenSeed> ormap_t;
//
//typedef pair<ormap_t::iterator,ormap_t::iterator> ormap_t_pair;
//
//class OpenReadList {
//
//private:
//	// position in genome scan 
//	size_t segment;    // chromosome or separate section thereof
//	size_t offset;     // greatest offset for entry or query in current segment
//	size_t genome_len; // total genome length
//
//	// greatest distance that may be spliced out.  This limits the search.
//	size_t max_splice_dist;
//
//	// map read_id to list of prefix occurrences in the preceding
//	// interval of size max_splice_dist.  
//	size_t map_size;
//	ormap_t map;
//
//	// locn of last map compression to removed expired reads
//	size_t purge_locn; 
//
//	// performance stats
//	size_t num_lookup;
//	size_t num_enter;
//	size_t num_purge;
//
//
//public:
//
//	// create an empty openread map
//	OpenReadList(size_t max_spl_dist = 10000, size_t sz = 1024) : 
//	  max_splice_dist(max_spl_dist), map_size(sz), genome_len(0)/*, map(sz) /*g++ only*/{
//		  segment = offset = purge_locn = 0;
//		  num_lookup = num_enter = num_purge = 0;
//		  map.bucket_size(sz);//for vc
//	  }
//
//	  // start a new segment; deleted all open reads
//	  void segmentStart() {
//		  segment++;
//		  map.clear();
//		  genome_len += offset + 1;
//		  offset = purge_locn = 0;
//	  }
//
//	  /*
//	  * Enter an open read rid starting at offset os.
//	  * Offsets are segment-relative and start at zero.  
//	  * The assumption is that reads are entered in approximately increasing
//	  *   order of offset; they are expired when the largest entered offset 
//	  *   or queried offset follows by distance greater than max_splice_dist
//	  */
//	  void enter(unsigned rid, size_t os, unsigned st) {
//		  num_enter++;
//		  offset = max(offset, os);
//		  if (DEBUG) 
//			  cout << "ENTER " << rid <<" at offset " << os << endl;
//		  map.insert(ormap_t::value_type(rid,OpenSeed(os,st)));
//		  if (offset - purge_locn > 3 * max_splice_dist && map.size() > map_size) {
//			  purge(offset);
//		  }
//	  }
//
//	  // find all open reads for read id rid within max_splice_dist
//	  pair<ormap_t::iterator, ormap_t::iterator> find_prefix_matches(unsigned rid, size_t cur_locn) {
//		  offset = max(offset, cur_locn);
//		  num_lookup++;
//		  if (DEBUG) 
//			  cout << "FIND occurrences of " << rid << " within distance " << max_splice_dist 
//			  << " preceding " << cur_locn << " (offset used = " << offset << "):" << endl;
//
//		  // iterate over occurrences of rid in map deleting expired occurrences
//		  // and determine result range as pair of iterators
//		  ormap_t::iterator itr = map.find(rid);
//		  ormap_t::iterator end = map.end();
//		  ormap_t::iterator res_begin = end;  // result range, initially empty
//		  ormap_t::iterator res_end = end;
//		  while (itr != end && itr->first == rid) {
//			  ormap_t::iterator cur_elt = itr;
//			  itr++;
//			  if (cur_elt->second.offset + max_splice_dist < cur_locn) {
//				  // expired occurrence
//				  if (DEBUG) 
//					  cout << "DELETING expired occurrence of " << cur_elt->first 
//					  << " at posn " << cur_elt->second.offset << endl;
//				  map.erase(cur_elt);
//			  }
//			  else { 
//				  // active occurrence
//				  if (res_begin == res_end) { 
//					  // first occurrence in result range
//					  res_begin = cur_elt;
//				  }
//				  if (DEBUG)
//					  cout << "RETAINED occurrence of " << cur_elt->first
//					  << " at posn " << cur_elt->second.offset << endl;
//			  }
//		  }
//		  // itr holds end of result range
//		  if (res_begin != res_end) {
//			  // nonempty range, update endpoint
//			  res_end = itr;
//		  }
//		  return (make_pair(res_begin,res_end));
//	  }
//
//	  // reduce map size by deleting expired entries
//	  void purge(size_t cur_locn){
//		  num_purge++;
//		  offset = max(offset,cur_locn);
//		  if (DEBUG)
//			  cout << "PURGE at segment " << segment << " offset " << offset 
//			  << " (map size = " << map.size() << ", purge interval = " 
//			  << offset - purge_locn <<")" << endl;
//		  ormap_t::iterator itr = map.begin();
//		  ormap_t::iterator end = map.end();
//		  while (itr != end) {
//			  ormap_t::iterator cur = itr;
//			  itr++;
//			  if (cur->second.offset + max_splice_dist < offset) {
//				  if (DEBUG)
//					  cout << "DELETING expired occurrence of " << cur->first 
//					  << " at posn " << cur->second.offset << endl;
//				  map.erase(cur);
//			  }
//		  }
//		  purge_locn = offset;
//		  if (DEBUG) 
//			  cout << "END PURGE (map size = " << map.size() << ")" << endl;
//	  }
//
//	  void stats() {
//		  if (VERBOSE) {
//			  //size_t entry_size = sizeof(ormap_t::key_type) + sizeof(ormap_t::data_type) + sizeof(ormap_t::pointer);//g++ only
//			  size_t entry_size = sizeof(ormap_t::key_type) + sizeof(ormap_t::mapped_type) + sizeof(ormap_t::pointer);//for vc
//			  size_t ht_size = map.bucket_count() * sizeof(ormap_t::pointer);
//			  size_t map_size = map.size() * entry_size + ht_size;
//			  cout << "STATS Openread:" << endl;
//			  cout << "Genome: " << 1 + segment << " segments, total length = " << genome_len + offset + 1 << endl;
//			  cout << "Reads entered = " << num_enter << ", reads queried = " << num_lookup << endl;
//			  cout << "Table purges = " << num_purge << endl;
//			  cout << "Table size " << map.size() << " entries and " << map.bucket_count() 
//				  << " buckets (about " << map_size << " bytes total)" << endl;
//		  }
//	  }
//
//
//	  // display all open reads currently in map
//	  void show_all(){
//		  cout << "Open Reads at segment " << segment << " offset " << offset << endl;
//		  ormap_t::iterator itr = map.begin();
//		  ormap_t::iterator end = map.end();
//		  while (itr != end) {
//			  cout << "open read of " << itr->first << " at locn " << itr->second.offset << endl;
//			  itr++;
//		  }
//		  cout << "END Open Reads" << endl;
//	  }
//};
//
//void show_matches(OpenReadList &list, unsigned id, size_t locn) {
//	pair<ormap_t::iterator,ormap_t::iterator> p = list.find_prefix_matches(id, locn);
//	ormap_t::iterator begin = p.first;
//	ormap_t::iterator end = p.second;
//	while (begin != end) {
//		cout << "RETURN occurrence of " << begin->first << " at posn " << begin->second.offset << endl;
//		begin++;
//	}
//}
//
//int count_matches(OpenReadList &list, unsigned id, size_t locn) {
//	pair<ormap_t::iterator,ormap_t::iterator> p = list.find_prefix_matches(id, locn);
//	ormap_t::iterator begin = p.first;
//	ormap_t::iterator end = p.second;
//	int result = 0;
//	while (begin != end) {
//		result++;
//		begin++;
//	}
//	return result;
//}

//class CircularList{
//
//public:
//	CircularList(size_t sz) : m_list(sz + 1), m_size(sz), m_used(0) {
//		lst = m_list.end();
//		--lst;
//		lend = lst;
//
//	}
//	CircularList() : m_list(1), m_size(0), m_used(0) {
//		lst = m_list.end();
//		--lst;
//		lend = lst;
//
//	}
//
//	void add_front(const OpenSeed& val){
//
//		if (m_size == m_used)
//		{
//			OpenSeed temp;
//			m_list.insert(lst, m_size / 2 + 1, temp);
//			m_size = m_list.size() - 1;
//		}
//
//		prev(lst);
//		*lst = val;
//		++m_used;
//	}
//
//	void remove_back(list <OpenSeed>::iterator lrm)
//	{
//		while(lend != lrm && lend != lst)
//		{
//			prev(lend);
//			--m_used;
//		}
//	}
//
//	list <OpenSeed>::iterator begin()
//	{
//		return lst;
//	}
//
//	list <OpenSeed>::iterator end()
//	{
//		return lend;
//	}
//
//	void next(list <OpenSeed>::iterator& cur)
//	{
//		if (cur == m_list.end())
//			cur = m_list.begin();
//		else
//			++cur;
//	}
//
//	void prev(list <OpenSeed>::iterator& cur)
//	{
//		if (cur == m_list.begin())
//			cur = m_list.end();
//		--cur;
//	}
//	
//	size_t usedsize()
//	{
//		return m_used;
//	}
//
//	size_t size()
//	{
//		return m_size;
//	}
//
//private:
//	list <OpenSeed> m_list;
//	list <OpenSeed>::iterator lst;
//	list <OpenSeed>::iterator lend;
//	size_t m_size;
//	size_t m_used;
//};

#define MinQueueSize 10


template<typename T>
class CyclicQueue{
public:
	typedef T value_type;
	typedef size_t size_type;
	typedef T& reference;
	typedef const T& const_reference;

	//CyclicQueue() : m_vec(MinQueueSize), m_capacity(MinQueueSize), m_size(0), m_head(0), m_tail(0) {}

	CyclicQueue(size_t sz = MinQueueSize) : m_vec(sz < MinQueueSize ? MinQueueSize : sz), 
		m_capacity(sz < MinQueueSize ? MinQueueSize : sz), m_size(0), m_front(1), m_rear(0), m_expend(0)
	{
		m_vec.reserve(m_capacity * 3);
		m_elem_size = sizeof(value_type);
	}

	~CyclicQueue()
	{
		m_vec.clear();
		m_front = 0;
		m_rear = 0;
		m_size = 0;
		m_capacity = 0;
	}

	bool IsEmpty() 
	{
		return m_size == 0;
	}

	bool IsFull()
	{
		return m_size == m_capacity;
	}

	void MakeEmpty()
	{
		m_size = 0;
		m_front = 1;
		m_rear = 0;
	}

	void Succ(size_t& index)
	{
		if (++index == m_capacity)
			index = 0;
	}

	void Pred(size_t& index)
	{
		if (--index == -1)
			index = m_capacity - 1;
	}

	void Enqueue(const_reference elem)
	{
		if(IsFull())
		{
			m_vec.resize(m_capacity * 4);
			
			if (m_front)
			{
				memcpy(&m_vec[m_front + (m_capacity * 3)], &m_vec[m_front], (m_capacity - m_front) * m_elem_size);
				m_front = m_front + m_capacity * 3;
			}
			else
			{
				if (m_rear != m_capacity - 1)
				{
					cout << "impossible"<<endl;
					exit(0);
				}
				memcpy(&m_vec[m_capacity * 4 - 1], &m_vec[0], m_elem_size);
				memcpy(&m_vec[m_capacity], &m_vec[1], (m_capacity - 1) * m_elem_size);
				memcpy(&m_vec[0], &m_vec[m_capacity], (m_capacity - 1) * m_elem_size);

				m_front = m_capacity * 4 - 1;
				m_rear = m_capacity - 2;
			}
			m_capacity *= 4;
			++m_expend;
		}

		++m_size;
		Succ(m_rear);
		m_vec[m_rear] = elem;
	}

	bool Dequeue()
	{
		if(IsEmpty())
			return false;
		else
		{
			--m_size;
			Succ(m_front);
			return true;
		}
	}

	bool RemoveFront(size_t index)
	{
		if (!Valid(index))
			throw "not valid index";
			//return false;
			
		Succ(index);
		if (index > m_front)
			m_size -= index - m_front;
		else
			m_size -= index - m_front + m_capacity;

		if (m_size > m_capacity)
		{
			throw "not valid index";
			//return false;
			//possible ?
			/*m_size = 0;
			m_front = m_rear;
			Succ(m_front);
			return false;*/
		}
		else
		{
			m_front = index;
			return true;
		}
	}

	reference Front()
	{
		//must exam is empty before use
		//if (!IsEmpty())
		return m_vec[m_front];
	}

	reference At(size_type index)
	{
		if (!Valid(index))
			throw "not valid index";
		return m_vec[index];
	}

	bool Valid(size_type index)
	{
		if (index > m_capacity - 1)
			return false;
		if (m_size == 0)
			return false;
		if (m_size == m_capacity)
			return true;
		if (m_rear >= m_front)
		{
			if (index >= m_front && index <= m_rear)
				return true;
			else
				return false;
		}
		else
		{
			if (index > m_rear && index < m_front)
				return false;
			else
				return true;
		}
	}

	size_type Rear()
	{
		if(IsEmpty())
			return -1;
		else
			return m_rear;
	}

	size_type End()
	{
		if (m_front)
			return m_front - 1;
		else
			return m_capacity - 1;
	}

	size_type Size()
	{
		return m_size;
	}

	size_type Capacity()
	{
		return m_capacity;
	}

	void Disp()
	{
		cout <<"size\t"<<m_size<<endl;
		cout <<"capacity\t"<<m_capacity<<endl;
		cout <<"rear\t"<<m_rear<<endl;
		cout <<"front\t"<<m_front<<endl;
		cout <<"expend\t"<<m_expend<<endl;
		cout <<"element size\t"<<m_elem_size<<endl;
	}

private:
	vector<value_type> m_vec;
	size_type m_front;
	size_type m_rear;
	size_type m_size;
	size_type m_capacity;
	size_type m_expend;
	size_type m_elem_size;
};

void Show(CyclicQueue<OpenSeed>& cqos)
{
	size_t st = cqos.Rear();
	size_t sz = cqos.Size();

	for (size_t i = 1; i <= sz; ++i)
	{
		CyclicQueue<OpenSeed>::reference now =  cqos.At(st);
		cout << "index " << st <<"\t" << now.Tostring() <<endl;
		cqos.Pred(st);
	}
	
	cout << "unused"<<endl;

	size_t cap = cqos.Capacity();
	
	for (size_t i = 1; i <= cap - sz; ++i)
	{
		CyclicQueue<OpenSeed>::reference now =  cqos.At(st);
		cout << "index " << st <<"\t"<< now.Tostring() <<endl;
		cqos.Pred(st);
	}
}

typedef CyclicQueue<OpenSeed> SEEDQUEUE;
typedef SEEDQUEUE::reference SEED_REF;
typedef SEEDQUEUE::const_reference SEED_CONST_REF;
typedef vector<SEEDQUEUE>ormap_t;
typedef ormap_t::iterator ormap_t_iter;

///////////////////////OpenReadList//////////////////////////
/*
 *    OpenReadList is created by Jan F. Prins, improved by Kai Wang with CyclicQueue
 *    Author: Jan F. Prins, Kai Wang
 */
class OpenReadList {

private:
	// position in genome scan 
	size_t segment;    // chromosome or separate section thereof
	size_t offset;     // greatest offset for entry or query in current segment
	size_t genome_len; // total genome length

	// greatest distance that may be spliced out.  This limits the search.
	size_t max_splice_dist;

	// map read_id to list of prefix occurrences in the preceding
	// interval of size max_splice_dist.  
	size_t map_size;
	ormap_t map;

	// locn of last map compression to removed expired reads
	size_t purge_locn; 

	// performance stats
	size_t num_lookup;
	size_t num_enter;
	size_t num_purge;

public:

	// create an empty openread map
	OpenReadList(size_t max_spl_dist = 10000, size_t kmer_sz = 4096) : 
	  max_splice_dist(max_spl_dist), map_size(kmer_sz), genome_len(0), map(kmer_sz, (max_spl_dist * 5 / kmer_sz) + 1 ) /*g++ only*/{
		  segment = offset = purge_locn = 0;
		  num_lookup = num_enter = num_purge = 0;
	  }

	  // start a new segment; deleted all open reads
	  void segmentStart() {
		  segment++;
		  map.clear();
		  map.resize(map_size);
		  genome_len += offset + 1;
		  offset = purge_locn = 0;
	  }

	  /*
	  * Enter an open read rid starting at offset os.
	  * Offsets are segment-relative and start at zero.  
	  * The assumption is that reads are entered in approximately increasing
	  *   order of offset; they are expired when the largest entered offset 
	  *   or queried offset follows by distance greater than max_splice_dist
	  */
	  void enter(size_t kmer, size_t os, size_t pt, unsigned strd, size_t bef_offset = -1) {
		  num_enter++;
		  offset = max(offset, os);
		  if (DEBUG) 
			  cout << "ENTER " << kmer <<" at offset " << os << endl;
		  map[kmer].Enqueue(OpenSeed(os,pt,strd, bef_offset));

		  //if most kmer hashed to one position, check it
		  //if (map[kmer].size() > max_splice_dist / 10 /*&& map.size() > map_size*/) {
			 // purge(offset, kmer);
		  //}
	  }

	  // find all open reads for read id rid within max_splice_dist
	  SEEDQUEUE& find_prefix_matches(size_t kmer, size_t cur_locn) {
		  offset = max(offset, cur_locn);
		  num_lookup++;
		  if (DEBUG) 
			  cout << "FIND occurrences of " << kmer << " within distance " << max_splice_dist 
			  << " preceding " << cur_locn << " (offset used = " << offset << "):" << endl;

		  return map[kmer];
	  }

	  void clear()
	  {
		  map.clear();
	  }

	  //void purge(size_t cur_locn, size_t kmer){
		 // num_purge++;

		 // ormap_t_iter cur_entry = map.begin() + kmer;

		 // SEEDLIST_ITER liter_end = cur_entry->end();

		 // SEEDLIST_ITER liter_back = --liter_end;

		 // while (liter_back != cur_entry->begin() && liter_back->offset + max_splice_dist - 2 < cur_locn)
			//  --liter_back;

		 // ++liter_back;

		 // if (liter_back != liter_end)
			//  cur_entry->erase(liter_back, cur_entry->end());
		 // //offset = max(offset,cur_locn);
		 // /*if (DEBUG)
			//  cout << "PURGE at segment " << segment << " offset " << offset 
			//  << " (map size = " << map.size() << ", purge interval = " 
			//  << offset - purge_locn <<")" << endl;
		 // ormap_t::iterator itr = map.begin();
		 // ormap_t::iterator end = map.end();
		 // while (itr != end) {
			//  ormap_t::iterator cur = itr;
			//  itr++;
			//  if (cur->second.offset + max_splice_dist < offset) {
			//	  if (DEBUG)
			//		  cout << "DELETING expired occurrence of " << cur->first 
			//		  << " at posn " << cur->second.offset << endl;
			//	  map.erase(cur);
			//  }
		 // }
		 // purge_locn = offset;
		 // if (DEBUG) 
			//  cout << "END PURGE (map size = " << map.size() << ")" << endl;*/
	  //}

	  void stats() {
		  if (VERBOSE) {
			  //size_t entry_size = sizeof(ormap_t::key_type) + sizeof(ormap_t::data_type) + sizeof(ormap_t::pointer);//g++ only
			  size_t entry_size = /*sizeof(ormap_t::key_type) + */sizeof(ormap_t::value_type) + sizeof(ormap_t::pointer);//for vc
			  size_t ht_size = map.size() * sizeof(ormap_t::pointer);
			  size_t map_size = map.size() * entry_size + ht_size;
			  cout << "STATS Openread:" << endl;
			  cout << "Genome: " << 1 + segment << " segments, total length = " << genome_len + offset + 1 << endl;
			  cout << "Reads entered = " << num_enter << ", reads queried = " << num_lookup << endl;
			  //cout << "Table purges = " << num_purge << endl;
			  //cout << "Table size " << map.size() << " entries and " << map.bucket_count() 
				//  << " buckets (about " << map_size << " bytes total)" << endl;
		  }
	  }

	  // display all open reads currently in map
	  void show_all(){
		  cout << "Open Reads at segment " << segment << " offset " << offset << endl;
		  for (size_t i = 0; i < map.size(); ++i)
		  {
			  size_t idx = map[i].Rear();
			  size_t sz = map[i].Size();
			   
			  for (size_t k = 1; k <= sz; ++k)
			  {
				  SEED_CONST_REF cur_seed = map[i].At(idx);
				  cout << "open read of " << i << " at locn " << cur_seed.offset << " index "<< cur_seed.index << " strand "<< cur_seed.strand<<endl;
				  map[i].Pred(idx);
			  }
		  }
		  cout << "END Open Reads" << endl;
	  }
};

void show_matches(OpenReadList& orlist, size_t kmer, size_t locn) {
	SEEDQUEUE p = orlist.find_prefix_matches(kmer, locn);
	size_t begin = p.Rear();
	size_t sz = p.Size();

	for (size_t i = 0; i < sz; ++i)
	{
		cout << "RETURN occurrence of " << kmer << " at posn " << p.At(begin).offset << endl;
		p.Pred(begin);
	}
}

size_t count_matches(OpenReadList &orlist, size_t kmer, size_t locn) {
  SEEDQUEUE p = orlist.find_prefix_matches(kmer, locn);
  return p.Size();
}

struct GroupedBwts;

struct TagSeed {
	unsigned score : 6;
	unsigned strand : 1;
	size_t readID;
	size_t pre_offset;
	string spliceway;
	vector<GroupedBwts*>::iterator bwts_ptr;

	TagSeed(unsigned scr, unsigned strd, size_t rid, size_t poft = 0, const string& spway = "") 
		: score(scr), strand(strd), readID(rid), pre_offset(poft), spliceway(spway) {}

	TagSeed(unsigned scr, unsigned strd, size_t rid, vector<GroupedBwts*>::iterator bptr, size_t poft = 0, const string& spway = "") 
		: score(scr), strand(strd), readID(rid), pre_offset(poft), bwts_ptr(bptr), spliceway(spway) {}
};

struct TagSeedReady {
	//unsigned short score : 6;
	char strand;
	vector<pair<WordPair, unsigned short> > tobe_fixed_hmer_segments_wp;
	size_t bef_offset;
	//string spliceway;
	//vector<GroupedBwts*>::iterator bwts_ptr;

	TagSeedReady(char strd, size_t poft, const vector<pair<WordPair, unsigned short> >& tbfhsw) 
		: strand(strd), tobe_fixed_hmer_segments_wp(tbfhsw){}

	//TagSeed(unsigned scr, unsigned strd, size_t rid, vector<GroupedBwts*>::iterator bptr, size_t poft = 0, const string& spway = "") 
	//	: score(scr), strand(strd), readID(rid), pre_offset(poft), bwts_ptr(bptr), spliceway(spway) {}
};

typedef vector<pair<size_t, size_t> > PAIR_SIZE_T_VEC;
typedef PAIR_SIZE_T_VEC::iterator PAIR_SIZE_T_VEC_ITER;
typedef PAIR_SIZE_T_VEC::reverse_iterator PAIR_SIZE_T_VEC_REV_ITER;
typedef map<string, PAIR_SIZE_T_VEC > TAC_POS_MAP;
typedef TAC_POS_MAP::iterator TAC_POS_MAP_ITER;

struct PrefixMatchSeed {
	WordPair prefix_wp;
	WordPair five_prim_suffix;
	
	PAIR_SIZE_T_VEC_ITER cur_tac_pos;
	vector<TagSeed> matchedTagsHalf;
	hash_map<unsigned int, unsigned int> entered_reads;
};

struct TobeFixedHmerReady;

struct PrefixMatchSeedReady {
	WordPair prefix_wp;
	WordPair five_prim_suffix;
	//WordPair prefix_wp2;
	//WordPair five_prim_suffix2;
	PAIR_SIZE_T_VEC_ITER cur_tac_pos;
	vector<TobeFixedHmerReady> matchedTagsHalf;
	//hash_set<size_t> entered_reads;
};

typedef hash_map<size_t, WordPair> PrefixWordsMap;
typedef PrefixWordsMap::iterator PrefixWordsMapIter;
typedef pair<PrefixWordsMapIter, bool> PrefixWordsMapInsertIter;

class PrefixWords{
private:
	PrefixWordsMap prefix_words;
	queue <PrefixWordsMapIter> offsetqueue;
	size_t max_splice_dist;
	size_t max_size;
	PrefixWordsMapInsertIter pr;
public:
	PrefixWords(size_t msd, size_t ms) : max_splice_dist(msd), max_size(ms)
	{
	}

	void set(const size_t chrom_offset, const WordPair& wp)
	{
		pr =  prefix_words.insert(make_pair(chrom_offset, wp));
		offsetqueue.push(pr.first);

		if (max_size < offsetqueue.size())
		{
			while (offsetqueue.back()->first - offsetqueue.front()->first > max_splice_dist)
			{
				prefix_words.erase(offsetqueue.front());
				offsetqueue.pop();
			}
		}
	}

	WordPair get(const size_t chrom_offset)
	{
		return prefix_words[chrom_offset];
	}
};



//////////////////////////Masks////////////////////////////////
struct Masks{
	Masks(const size_t read_width, const size_t seed_width, const size_t max_mismatch, const size_t anchor_width, const size_t num_anchor, const size_t num_seg, const size_t seg_width, const size_t extend_bits);

	void Set(const size_t read_width, const size_t seed_width);

	Masks (size_t seg_len)
	{
		//duplicate reads
		comp_buff_width = 2 * seg_len;

		comb_seg_bits_on = (ALL_BITS_ON >> (SIXTY_FOUR - comp_buff_width)); 

		suffix_seg_bits_on = (ALL_BITS_ON >> (SIXTY_FOUR - seg_len)); ;

		prefix_seg_bits_on = suffix_seg_bits_on << seg_len;

		comb_seg_first_selector_rt = LEAST_SIG_BIT;

		comb_seg_first_selector_lt = LEAST_SIG_BIT << seg_len;

		score_seg_buf_width = seg_len;
	}

	size_t mid_buff_width;

	// marks bits corresponding to the seed_key
	size_t small_mask;

	size_t kmer_mask;

	// marks  rightmost bits numbering the width of a read
	size_t big_mask;

	size_t suffix_mask;

	size_t prefix_mask;

	size_t big_buff_mask;

	size_t mask_ps;

	size_t mask_midright;

	// bit-vector indicating the location of the hit key in the frame
	size_t bad_base_mask;

	size_t bad_base_maskII;

	size_t bad_kmer_mask;

	vector<size_t> bad_kmer_masks;

	vector<size_t> bad_kmer_masks_upper;

	vector<size_t> bad_kmer_masks_lower;

	// marks bit where the key starts in the chromosome frame, the seed_width from the left
	size_t key_update_bit;

	size_t key_update_bitII;

	vector<size_t> kmer_update_bits;

	vector<size_t> kmer_update_bits_upper;

	vector<size_t> kmer_update_bits_lower;

	//vector<size_t> imer_update_bits;

	size_t first_selector_rt;

	size_t first_selector_lt;

	size_t mid_buff_width_ext;

	size_t first_half_mask;

	size_t second_half_mask;

	size_t first_half;

	size_t second_half;

	//////

	size_t comp_buff_width;

	size_t comp_left_shift_width;

	size_t comp_right_shift_width;

	size_t comp_first_half_mask;

	size_t comp_second_half_mask;

	size_t comp_big_buff_mask;

	size_t score_buff_width;

	size_t score_first_half_mask;

	size_t score_second_half_mask;

	size_t score_big_buff_mask;

	size_t score_first_selector_rt;

	size_t score_first_selector_lt;

	size_t mis_first_selector_lt;

	size_t comp_first_selector_rt;

	size_t comp_first_selector_lt;

	size_t comp_flankstr_left_shift_width;

	size_t comp_flankstr_right_shift_width;

	///fix hole

	//duplicate reads
	size_t comb_seg_bits_on;

	size_t suffix_seg_bits_on;

	size_t prefix_seg_bits_on;

	size_t right_shift_seg_width;

	size_t left_shift_seg_width;

	size_t comb_seg_first_selector_rt;

	size_t comb_seg_first_selector_lt;

	size_t score_seg_buf_width;

	//fix tail
	size_t ft_comb_prefix_half_mask;

	size_t ft_comb_suffix_half_mask;

	size_t ft_combined_mask;

	//duplicated reads extend
	size_t comb_seg_bits_on_ext;

	size_t suffix_seg_bits_on_ext;

	size_t prefix_seg_bits_on_ext;

	size_t right_shift_seg_width_ext;

	size_t left_shift_seg_width_ext;

	size_t comb_seg_first_selector_rt_ext;

	size_t comb_seg_first_selector_lt_ext;

	size_t score_seg_buf_width_ext;

	size_t mis_selector_lt;

	size_t mis_selector_rt;

	size_t prefix_ext_mask;

	size_t suffix_ext_mask;

	size_t prefix_append_mask;

	size_t prefix_append_shifted_mask;

	size_t prefix_suffix_append_mask;

	size_t reserve_bit;

};

Masks::Masks(const size_t read_width, const size_t seed_width, const size_t max_mismatch, const size_t anchor_width, const size_t num_anchor, const size_t num_seg, const size_t seg_width, const size_t extend_bits) : 
	bad_kmer_masks(num_anchor), bad_kmer_masks_upper(num_anchor), bad_kmer_masks_lower(num_anchor), kmer_update_bits(num_anchor), kmer_update_bits_upper(num_anchor), kmer_update_bits_lower(num_anchor)/*, imer_update_bits(num_anchor)*/
{
	mid_buff_width = read_width - 2*seed_width + max_mismatch;//+ max_mismatch

	/////

	first_half = read_width;//(read_width + 1) / 2;

	second_half = read_width;// / 2;

	comp_buff_width = read_width;// / 2;// - anchor_width;

	comp_left_shift_width = comp_buff_width;// - anchor_width;

	comp_flankstr_left_shift_width = comp_buff_width - (anchor_width * 2);

	comp_flankstr_right_shift_width = anchor_width;

	comp_right_shift_width = 0;//anchor_width;

	comp_second_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - comp_buff_width)); 

	comp_first_half_mask = comp_second_half_mask << comp_buff_width;

	comp_big_buff_mask = comp_second_half_mask | comp_first_half_mask;

	score_buff_width = seg_width - anchor_width;// read_width - anchor_width;//read_width / 2 - anchor_width;

	score_second_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - score_buff_width)) << anchor_width;

	score_first_half_mask = score_second_half_mask << score_buff_width << anchor_width;

	score_big_buff_mask = score_second_half_mask | score_first_half_mask;

	comp_first_selector_rt = LEAST_SIG_BIT;

	comp_first_selector_lt = comp_first_selector_rt << comp_buff_width;

	score_first_selector_rt = LEAST_SIG_BIT << anchor_width;

	score_first_selector_lt = score_first_selector_rt << seg_width;//score_buff_width << anchor_width;

	mis_first_selector_lt = LEAST_SIG_BIT << (2 * seg_width - 1);

	small_mask = (LEAST_SIG_BIT << 2*seed_width) - 1;

	kmer_mask = (LEAST_SIG_BIT << 2*anchor_width) - 1;

	big_mask = (ALL_BITS_ON >> (SIXTY_FOUR - read_width)); 

	suffix_mask = (ALL_BITS_ON >> (SIXTY_FOUR -  read_width + seed_width - max_mismatch));// - max_mismatch

	prefix_mask = suffix_mask << (read_width - seed_width + max_mismatch);// + max_mismatch

	big_buff_mask = (ALL_BITS_ON >> (SIXTY_FOUR - 2 * read_width + 2 * seed_width - 2 * max_mismatch));//- 2 * max_mismatch

	bad_base_mask = ((LEAST_SIG_BIT << seed_width) - 1) << (read_width - seed_width);

	bad_base_maskII = bad_base_mask >> ((read_width + 1) / 2 - seed_width);

	bad_kmer_mask = ((LEAST_SIG_BIT << anchor_width) - 1);

	for (size_t i = 0; i < num_anchor; ++i)
	{
		bad_kmer_masks[i] = bad_kmer_mask << (i * anchor_width);

		bad_kmer_masks_upper[i] = bad_kmer_mask << (i * anchor_width) << 2;

		bad_kmer_masks_lower[i] = bad_kmer_mask << (seg_width - ((i + 1) * anchor_width)) << 2;
	}

	key_update_bit = LEAST_SIG_BIT << (read_width - seed_width);

	key_update_bitII = key_update_bit >> ((read_width + 1) / 2 - seed_width);

	for (size_t i = 0; i < num_anchor; ++i)
	{
		kmer_update_bits[i] = LEAST_SIG_BIT << (i * anchor_width);

		kmer_update_bits_upper[i] = LEAST_SIG_BIT << (i * anchor_width) << 2;

		kmer_update_bits_lower[i] = LEAST_SIG_BIT << (seg_width - ((i + 1) * anchor_width)) << 2;

		//imer_update_bits[i] = LEAST_SIG_BIT << (seg_width - ((i + 1) * anchor_width)) << 2;
	}

	mask_midright = (ALL_BITS_ON >> (SIXTY_FOUR - mid_buff_width - max_mismatch))  << (seed_width - max_mismatch);//- max_mismatch

	mask_ps = big_buff_mask - mask_midright;

	mid_buff_width_ext = mid_buff_width + max_mismatch;

	first_selector_rt = LEAST_SIG_BIT << (seed_width - max_mismatch);

	first_selector_lt = first_selector_rt << mid_buff_width_ext;//to be turned off

	second_half_mask = ALL_BITS_ON >> (SIXTY_FOUR - ((read_width + 1) / 2));

	first_half_mask = second_half_mask << (read_width / 2);

	second_half_mask = ALL_BITS_ON >> (SIXTY_FOUR - (read_width / 2));

	//

	//size_t seg_width = read_width;// / num_seg;

	///fix hole
	comb_seg_bits_on = ALL_BITS_ON >> (SIXTY_FOUR - (2 * seg_width));

	suffix_seg_bits_on = ALL_BITS_ON >> (SIXTY_FOUR - seg_width);

	prefix_seg_bits_on = suffix_seg_bits_on << seg_width;

	right_shift_seg_width = 0;

	left_shift_seg_width = seg_width;

	comb_seg_first_selector_rt = LEAST_SIG_BIT;

	comb_seg_first_selector_lt = LEAST_SIG_BIT << seg_width;

	score_seg_buf_width = seg_width;

	//fix hole extend
	comb_seg_bits_on_ext = ALL_BITS_ON >> (SIXTY_FOUR - (2 * (seg_width + (2 * extend_bits))));

	suffix_seg_bits_on_ext = ALL_BITS_ON >> (SIXTY_FOUR - seg_width - (2 * extend_bits));

	prefix_seg_bits_on_ext = suffix_seg_bits_on_ext << (seg_width + (2 * extend_bits));

	right_shift_seg_width_ext = 0;

	left_shift_seg_width_ext = seg_width + (2 * extend_bits);

	comb_seg_first_selector_rt_ext = LEAST_SIG_BIT;

	comb_seg_first_selector_lt_ext = LEAST_SIG_BIT << (seg_width + (2 * extend_bits));

	score_seg_buf_width_ext = seg_width + (2 * extend_bits);

    mis_selector_lt = LEAST_SIG_BIT << ((score_seg_buf_width_ext * 2) - 1);

	mis_selector_rt = LEAST_SIG_BIT;

	suffix_ext_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits);

	prefix_ext_mask = suffix_ext_mask << seg_width;

	prefix_append_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits - seg_width);

	prefix_append_shifted_mask = prefix_append_mask << extend_bits;

	prefix_suffix_append_mask = ALL_BITS_ON >> (SIXTY_FOUR - extend_bits - extend_bits - seg_width);

	//fix tail
	ft_comb_suffix_half_mask = (ALL_BITS_ON >> (SIXTY_FOUR - seg_width));

	ft_comb_prefix_half_mask = ft_comb_suffix_half_mask << seg_width;

	ft_combined_mask = ft_comb_suffix_half_mask | ft_comb_prefix_half_mask;

	//fix hmer
	
	reserve_bit = LEAST_SIG_BIT << (seg_width - 1);
}

//////////////////////////////GenomeScan////////////////////////
//typedef hash_map<size_t, ConjInfo> JUNC_SUFFIX;
//typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
//typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;
//
//typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
//typedef JUNC_HASH::iterator JUNC_HASH_ITER;
//typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef hash_map<size_t, JuncInfo> JUNC_SUFFIX;
typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;

typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
typedef JUNC_HASH::iterator JUNC_HASH_ITER;
typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef map<string, JUNC_HASH > CONJ_HASH_MAP;

typedef hash_map<size_t, JuncInfoFusion> JUNCFusion_SUFFIX;
typedef JUNCFusion_SUFFIX::iterator JUNCFusion_SUFFIX_ITER;
typedef JUNCFusion_SUFFIX::const_iterator JUNCFusion_SUFFIX_CITER;

typedef hash_map<size_t, JUNCFusion_SUFFIX > JUNCFusion_HASH;
typedef JUNCFusion_HASH::iterator JUNCFusion_HASH_ITER;
typedef JUNCFusion_HASH::const_iterator JUNCFusion_HASH_CITER;

typedef hash_map<size_t, JuncInfo> JUNC_HASH_COMB;

typedef hash_map<string, JUNC_HASH_COMB> CHROM_JUNC_HASH_COMB;

typedef map<string, JUNCFusion_HASH > CONJFusion_HASH_MAP;

typedef hash_map<size_t, JuncInfoFusion> FUSION_JUNC_HASH_COMB;

typedef FUSION_JUNC_HASH_COMB::iterator FUSION_JUNC_HASH_COMB_ITER ;

typedef hash_map<string, FUSION_JUNC_HASH_COMB> CHROM_FUSION_JUNC_HASH_COMB;

typedef CHROM_FUSION_JUNC_HASH_COMB::iterator CHROM_FUSION_JUNC_HASH_COMB_ITER;

bool compare_pair(const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs)
{
	if (lhs.first == rhs.first)
		return lhs.second < rhs.second;

	return lhs.first < rhs.first;
}

bool compare_pair_rev(const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs)
{
	if (lhs.second == rhs.second)
		return lhs.first > rhs.first;

	return lhs.second > rhs.second;
}

struct CoverageBlock{
	size_t blockst;
	size_t blockend;
	size_t hits;
	CoverageBlock(const size_t& blst, const size_t& blend, const size_t& ht) : blockst(blst), blockend(blend), hits(ht) {}
};

bool
compblk(const CoverageBlock& lhs, const CoverageBlock& rhs)
{
	if (lhs.blockst == rhs.blockst)
		return lhs.blockend < rhs.blockend;
	else
		return lhs.blockst < rhs.blockst;
}

struct BowtieSegment{
	//HWI-EAS217:4:8:196:121#0/1_34212	1	-	chr7	72247608	TTGAGGCCTGGAATGCCAAAATCAC	1	23:A>T,24:C>T	291054
	string tagname;
	size_t segid;
	char strand;
	string chromname;
	size_t offset;
	string segseq;
	size_t othermapped;
	string alters;
	size_t alignid;

	BowtieSegment(const string& tname, const size_t& sid, const char& sd, const string& cname, const size_t& oft, const string& sseq, const size_t& omapped, const string& alt, const size_t& aid)
		: tagname(tname), segid(sid), strand(sd), chromname(cname), offset(oft), segseq(sseq), othermapped(omapped), alters(alt), alignid(aid) {}
};

struct SamFixtail{
	//HWI-EAS217:4:9:774:491#0/1_115277	16	chr1	521896	0	37M30N38M	*	0	0	attcatgcacacccagatcatatat	IIIIIIIIIIIIIIIIIIIIIIIII	2	1	ATAC	1
	string tagname;
	unsigned int strand;
	size_t offset;
	size_t fixtailoffset;
	string splitway;
	vector<pair<size_t, pair<size_t, string> > > flanks;
	size_t readID;

	SamFixtail(const string& tname, const unsigned int& sd, const size_t& oft, const size_t& ftoft, const string& sway, const size_t& rid) :  tagname(tname),
		strand(sd), offset(oft), fixtailoffset(ftoft), splitway(sway), readID(rid) {}
};

bool
compSamFixtail(const SamFixtail& lhs, const SamFixtail& rhs)
{
	return lhs.fixtailoffset < rhs.fixtailoffset;
}

enum FixTypes
{
	mistailandhead_type,
	mistailorhead_type,
	fixhmer_and_hole_type,
	fixhole_and_hmer_type,
	fixfusion_and_hmer_type
};

struct Differ{
	size_t index;
	char src_char;
	char dst_char;

	Differ(size_t idx, char src, char dst) : index(idx), src_char(src), dst_char(dst) {}

};

struct DifferVec{
	vector<Differ> differs_vec;

	DifferVec(const string& differs)
	{
		if (differs.length() > 1)
		{
			//cout << differs<<endl;
			size_t index = 0;
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
		}
	}
};

struct GroupedBwtRec{
	size_t offset;
	size_t seg_id;
	string mapped_chrom_seq;
	size_t other_mapped;
	string differs;
	string chromname;
	char strand;

	vector<Differ> differs_vec;

	GroupedBwtRec(size_t oft, size_t sid, const string& mcs, size_t omd, const string& dfs) : offset(oft),
		seg_id(sid), mapped_chrom_seq(mcs), other_mapped(omd), differs(dfs)
	{
		if (differs.length() > 1)
		{
			//cout << differs<<endl;
			size_t index = 0;
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

			//if (strand == '-')
			//{
			//	reverse(differs_vec.begin(), differs_vec.end());

			//	vector<Differ>::iterator v_Differ_iter;
			//	for (v_Differ_iter = differs_vec.begin(); v_Differ_iter != differs_vec.end(); ++v_Differ_iter)
			//	{
			//		v_Differ_iter->index = mapped_chrom_seq.length() - v_Differ_iter->index - 1;
			//	}
			//}
		}
	}

	GroupedBwtRec(size_t oft, size_t sid, const string& mcs, size_t omd, const string& dfs, const string& cn, char sd) : offset(oft),
		seg_id(sid), mapped_chrom_seq(mcs), other_mapped(omd), differs(dfs), chromname(cn), strand(sd)
	{
		size_t index = 0;
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

			//if (sd == '-')
			//{
			//	reverse(differs_vec.begin(), differs_vec.end());

			//	vector<Differ>::iterator v_Differ_iter;
			//	for (v_Differ_iter = differs_vec.begin(); v_Differ_iter != differs_vec.end(); ++v_Differ_iter)
			//	{
			//		v_Differ_iter->index = mapped_chrom_seq.length() - v_Differ_iter->index - 1;
			//	}
			//}
		}
	}
};

struct WordPairHmer{
	WordPair seg_wp;
	unsigned short seg_len;
	//vector<Kmer> prefix_hmers;
	//vector<Kmer> suffix_hmers;

	WordPairHmer (const WordPair& swp, unsigned short slen/*, size_t hmer_len, size_t num_hmer*/) : seg_wp(swp), seg_len(slen)
	{
		//for (size_t i = 0; i < num_hmer; ++i)
		//{
		//	suffix_hmers.push_back(seg_wp.get_kmer(1 + (i * hmer_len), hmer_len + (i * hmer_len)));

		//	prefix_hmers.push_back(seg_wp.get_kmer(seg_len - hmer_len + 1 - (i * hmer_len), seg_len - (i * hmer_len)));
		//}
	}
};

struct TobeFixedHmerReady{
	string tagname;
	size_t st_offset;
	size_t end_offset;
	size_t st_seg_id;
	size_t end_seg_id;
	vector<WordPairHmer> to_be_fixed_hmer_segments_wp;
	char strand;
	int id1, id2;
	/*size_t other_mapped;*/
	//string differs;
	//string chromname;
	//vector<Differ> differs_vec;

	TobeFixedHmerReady(const string& tn, size_t st_oft, size_t end_oft, size_t st_sid, size_t end_sid, const string& tobe_fixed_hmer_strs, char sd, int i1, int i2, bool is_reverse = false/*, size_t hmer_len, size_t num_hmer*/) : 
						tagname(tn), st_offset(st_oft),	end_offset(end_oft), st_seg_id(st_sid), end_seg_id(end_sid), strand(sd), id1(i1), id2(i2)
	{
		size_t seg_num;

		if (is_reverse)
		{
			if (strand == '-')
				seg_num = end_sid - st_sid + 1;
			else
				seg_num = st_sid - end_sid + 1;
		}
		else
		{
			if (strand == '+')
				seg_num = end_sid - st_sid + 1;
			else
				seg_num = st_sid - end_sid + 1;
		}

		size_t pre_index = 0, index = 0;

		for (size_t i = 0; i < seg_num; ++i)
		{
			index = tobe_fixed_hmer_strs.find(",", pre_index);

			if (index == string::npos)
			{
				cout << "missed a segment?"<< endl;
				cout << tobe_fixed_hmer_strs << endl;
				cout << st_sid <<'\t'<<end_sid<<endl;
				cout <<st_oft<<'\t'<<end_oft<<endl;
				cout << sd << endl;
				getchar();
			}

			string tobe_fixed_hmer_str = tobe_fixed_hmer_strs.substr(pre_index, index - pre_index);

			if (is_reverse)
				reverse(tobe_fixed_hmer_str.begin(), tobe_fixed_hmer_str.end());

			to_be_fixed_hmer_segments_wp.push_back(WordPairHmer(WordPair(tobe_fixed_hmer_str), static_cast<unsigned short>(tobe_fixed_hmer_str.length())) );

			pre_index = index + 1;
		}
	}

	//TobeFixedHmerReady(size_t oft, size_t sid, const string& mcs, size_t omd, const string& dfs, const string& cn, char sd) : offset(oft),
	//	seg_id(sid), mapped_chrom_seq(mcs), other_mapped(omd), differs(dfs), chromname(cn), strand(sd)
	//{
	//	//size_t index = 0;
	//	//if (differs.length() > 1)
	//	//{
	//	//	while (true)
	//	//	{
	//	//		//0:G>A,1:T>A
	//	//		int idx;
	//	//		char src, dst;
	//	//		sscanf(differs.c_str() + index, "%d:%c>%c", &idx, &src, &dst);
	//	//		differs_vec.push_back(Differ(idx, src, dst));
	//	//		index = differs.find(",", index);

	//	//		if (index == string::npos)
	//	//			break;

	//	//		++index;
	//	//	}
	//	//}
	//}
};

struct SideMatch{
	string prefix_str;
	string suffix_str;
	size_t ls_st;
	size_t ls_end;
	size_t rs_st;
	size_t score;

	vector<Differ> differs_vec;

	SideMatch(const string& pstr, const string& sstr, size_t ls, size_t le, size_t rs, size_t sc, const vector<Differ>& d_v) : prefix_str(pstr), 
		suffix_str(sstr),  ls_st(ls), ls_end(le), rs_st(rs), score(sc), differs_vec(d_v) {}
};

struct GroupedBwts
{
	string tagname;
	string chromname;
	char strand;
	FixTypes fixtype;
	size_t hmer_offset;
	size_t hmer_offset_rev;
	size_t tagidx;
	size_t readID;
	size_t readID_rev;
	size_t bef_offset;
	size_t bef_offset_rev;
	vector<GroupedBwtRec> grouped_bwt_recs;

	vector<SideMatch>* prefix_side;
	vector<SideMatch>* suffix_side;

	GroupedBwts(const string& tn, const string& cn, char sd, const size_t tidx, const FixTypes& ft) : tagname(tn), chromname(cn), strand(sd), tagidx(tidx), fixtype(ft), 
		hmer_offset(0), hmer_offset_rev(0), readID(-1), readID_rev(-1), bef_offset(-1), bef_offset_rev(-1) {}

	GroupedBwts(const string& tn, const string& cn, char sd, const size_t tidx, const FixTypes& ft, const vector<GroupedBwtRec>& gbr) : tagname(tn), chromname(cn), strand(sd), tagidx(tidx), 
		fixtype(ft), hmer_offset(0), hmer_offset_rev(0), readID(-1), readID_rev(-1), bef_offset(-1), bef_offset_rev(-1), grouped_bwt_recs(gbr) {} 
};

bool compFixHmer(const GroupedBwts* lhs, const GroupedBwts* rhs)
{
	return lhs->hmer_offset < rhs->hmer_offset;
}

bool compFixHmer_rev(const GroupedBwts* lhs, const GroupedBwts* rhs)
{
	return lhs->hmer_offset_rev < rhs->hmer_offset_rev;
}

typedef vector<CoverageBlock> VEC_COVER_BLOCK;

class GenomeScan{
public:
	GenomeScan(size_t input_buffer_size, size_t seed_width, size_t read_width, size_t max_mismatches);

	GenomeScan(const CmdlineArgus& argulist);

	void GetReads(const char* reads_file);

	void GetReadsFastq(const char* reads_file);

	void GetReadsFastqQualitySeq(const char* reads_file);

	const string& GetQualitySeq(size_t tag_idx, size_t strand_t, const vector<Differ>& differs_vec, double& qual_score);

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

	void CalcReadsQualityScore();

	void CleanReads();

	void GetSeed(const size_t base);

	SEEDS_CONST_ITER GetTags(size_t seed_key);

	size_t PrefixMatch(const WordPair& read_word);

	void AddToOpenReads();

	void TriggerMatch(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed);

	void TriggerMatchFixTail(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed);

	void TriggerMatchFixHmer(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase);

	void TriggerMatchFixHmerReady(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase);

	void TriggerMatchFixHmerReady_Rev(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase);

	void TriggerMatchFixHmerInsertReady(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase);

	void TriggerMatchFixHmerInsertReady_Rev(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase);

	void TriggerMatchFixHole();

	void SetReadWords();

	void SetReadWords_Rev();

	void SetChromFiles(const char* filenames_file, const char* chrom_file, const char* fasta_suffix, string tempchrom = "tempchrom/");

	void LoadChromFiles(const char* filenames_file, const char* chrom_file, const char* fasta_suffix);

	void SetTissueFiles(const char* filenames_file, const char* tissue_file, const char* fasta_suffix);

	void MapReads(size_t chromID);

	void CallMapReads();

	void FixHole(const char* fixhole_file, const char* fixed_hole_file, size_t flankcase, bool avoid = false, bool mark = false, bool exceed = false);

	void FixHoleExt(const char* fixhole_file, const char* fixed_hole_file, size_t flankcase, bool avoid = false, bool mark = false, bool exceed = false);

	void FixHoleReady(const char* fixhole_file, const char* fixed_hole_file, hash_set<size_t>& junc_set, size_t flankcase, bool avoid = false, bool mark = false, bool exceed = false);

	void FixHoleReadyNonCanon(const char* fixhole_file, const char* fixed_hole_file, hash_set<size_t>& junc_set, size_t flankcase, bool avoid = false, bool mark = false, bool exceed = false);

	void FixHole_fusion(const char* fixhole_file, const char* fixed_hole_file, CONJFusion_HASH_MAP& conj_hash_map, bool avoid = false, bool mark = false);

	void FixFusionReady(const char* fusion_file, const char* fixed_fusion_file, CHROM_FUSION_JUNC_HASH_COMB& conj_hash_map, bool avoid, bool mark, size_t flankcase, string& chromseq1, string& chromseq2);

	void UpdateReadCountAndClean()
	{
		for (size_t i = 0; i < m_reads_count.size(); ++i)
		{
			m_reads_count[i] += m_cur_reads_count[i];
			m_cur_reads_count[i] = 0;
		}
	}

	void StartFixHmer(int phase, size_t tissueID, size_t flankcase, bool avoid);

	void FixTail(size_t chromID, const vector<SamFixtail>& cur_fixtail_vec);

	void FixHmer(size_t chromID, vector<GroupedBwts*>& grouped_bwts, size_t flankcase, bool avoid);

	void FixHmerReady(const char* fixhmer_file, const char* fixed_hmer_file, size_t flankcase, bool avoid, bool island);

	void FixHmerReady_Rev(const char* fixhmer_file, const char* fixed_hmer_file, size_t flankcase, bool avoid, bool island);

	void PrepareFixhole(const char* fixhole_file, const char* fixed_hole_file);

	void StartMapReads(int phase, size_t tissueID);

	void OutputJunction();

	void Covert2Junc();

	void SeprateMapreads(vector<string>& comb_mapreads_files);

	void RemoveDupMapreads(vector<string>& comb_mapreads_files);

	void WriteJunc(const CONJ_HASH_MAP& conj_hash_map, ofstream& ofs);

	void WriteJuncFusion(const CONJFusion_HASH_MAP& conj_hash_map, ofstream& ofs);

	void WriteJuncFusion(const CHROM_FUSION_JUNC_HASH_COMB& conj_hash_map, ofstream& ofs);

	void SortJunc(const char* juncfile);

	void SortJuncFusion(const char* juncfile);

	void WriteCoverage(const VEC_COVER_BLOCK& v_coverage_block, const string& chromname, ofstream& ofs, const size_t& chromID);

	void MergeCoverageBlock(const VEC_COVER_BLOCK& v_coverage_block, vector<unsigned short>& v_merged_coverage_block);

	bool Clear();

	bool Release();

    size_t Ps_score(const WordPair & read_word_dup, size_t& loc);

	size_t Ps_score_selective(const WordPair & read_word_dup, size_t& loc, size_t& prim, size_t& rbits);

	size_t Ps_score_selective_rev(const WordPair & read_word_dup, size_t& loc, size_t& prim, size_t& rbits);

	size_t Fixhole_score_selective(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch);

	size_t Fixhole_score_selective_ext(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits);

	size_t Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixhole_score_selective_var_mask_rev(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixhole_score_selective_insert_var_mask_rev(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixfusion_score(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch);

	size_t Fixhole_score_selective_rev(const WordPair & read_word_dup, size_t& loc, size_t& prim);

	void Sites_To_Regions();

	void Write_Non_Uniques();

	void ReadAndMergeTacPos(const char* tac_file);

	void ReadBowtieSegments(const char* bowtie_segments_file);

	void ReadSamFixtail(const char* samfixtailfile);

	void ReadFixHmer();

	void ReadGroupedBwt(const char* infile, const enum FixTypes& fixhmertype);

	void TryAppend(const char* infile, const FixTypes& fixhmertype);

	void Classify_Append(map<string, vector<GroupedBwts> >& fixhmer, const char* infile, const FixTypes& fixhmertype);

	void SortFixtail();

	void WriteMapReads(size_t prefix_width, WordPair prefix_match, const size_t& prim, const size_t& lste_st, const size_t& lste_end, const size_t& rste_st, const size_t& chrom_ID, const size_t& readID, unsigned strand, const string& splitway = "");

	void WriteFixHmer(size_t prefix_width, WordPair& read_word, const size_t& prim, const size_t& lste_st, const size_t& lste_end, const size_t& rste_st, const size_t& chrom_ID, 
		const size_t& readID, unsigned strand, vector<GroupedBwts*>::iterator cur_bwts, const size_t score, vector<Differ>& differs_vec);

	void WriteFixedHeadAndTail(const char* file_fixedHeadAndTail, size_t flankcase);

	bool CheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckFirstTimeRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckBeforeMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckAfterMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	bool CheckAfterGTAGMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s);

	//
	bool FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckFirstTimeRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckBeforeMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterGTAGMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	void FlankString(string& flankstr);

	string FlankString(size_t matched_flank, bool bad);

	void CreateMapReadsFile(size_t tissueID, size_t chromID);

	void CloseCurMapReadsFile()
	{
		m_ofs_map_reads.close();
	}

	void InitReads(size_t tissueID);

	void InitReadsFastq(size_t tissueID);

	void InitReadWords();

	void InitReadWords_Rev();

	void ClearReads();

	void ClearReadWords();

	void CloseMapReadFile()
	{
		m_ofs_map_reads.close();
	}

	vector<Kmer> Str2Kmer(const string& readstr)
	{
		vector<Kmer> kmers;

		for (size_t i = 0; i < m_num_anchors; ++i)
		{
			unsigned kmer = 0;

			bool good_hit = true;

			size_t offset = readstr.length() - m_anchor_width - (i * m_kmer_step);

			for (size_t k = 0; k < m_anchor_width && good_hit; ++k)
			{
				if (isvalid(readstr[offset + k]))
					kmer = (kmer << 2) + (unsigned)base2int(readstr[offset + k]);
				else good_hit = false;
			}

			kmers.push_back(Kmer(!good_hit, kmer));
		}
		return kmers;
	}

	size_t EnterSlideWindow(const SEED_VEC& seed_vec);

	size_t EnterSlideWindowII(const SEED_VEC& seed_vec, size_t prev_hits);

	void EnterSlideWindowFixTail(vector<SamFixtail>::const_iterator ft_seed_begin, vector<SamFixtail>::const_iterator ft_seed_end);

	void EnterSlideWindowFixHmer(vector<GroupedBwts*>::iterator fh_seed_begin, vector<GroupedBwts*>::iterator fh_seed_end, bool avoid);

	void EnterSlideWindowFixHmerReady(vector<TobeFixedHmerReady>& cur_tobe_fixed_hmer_segments, const char* chromname, bool avoid);

	void EnterSlideWindowFixHmerReady_Rev(vector<TobeFixedHmerReady>& cur_tobe_fixed_hmer_segments, const char* chromname, bool avoid);

private:
	
	//global variables, only initialize for one time
	SeedHash m_seed_hash;

	vector<string> m_read_names, m_reads, m_fastq_seq;

	vector<size_t> m_reads_count;

	vector<double> m_quality_scores;

	vector<size_t> m_cur_reads_count;

	vector<ReadInfo> m_bests;

	vector<WordPair> m_read_words, m_read_words_rc;

	vector<BowtieSegment> m_bowtie_segments;

	map<string, vector<GroupedBwts> > m_fixhmer;

	map<string, vector<GroupedBwts*> > m_fixhmer_ptr;

	map<string, vector<GroupedBwts*> > m_fixhmer_ptr_rev;

	map<string, vector<GroupedBwts*> > m_fixheadandtail_ptr;

	map<size_t, map<string, vector<SamFixtail> > > m_sam_fix_tail_map;

	vector<vector<Kmer> > m_read_kmers, m_read_kmers_rc;

	vector<string> m_chrom_names;

	map<string, size_t> m_chrom_names_map;

	vector<string> m_chrom_files;

	vector<string> m_rev_chrom_files;

	vector<size_t> m_chrom_sizes;

	vector<string> m_tissue_files;

	vector<string> m_tissue_names;

	vector<size_t> m_tissue_sizes;

	const Masks m_masks;

	hash_map<size_t, Masks> m_hash_map_masks;

	const char* m_ambiguous_file;

	const char* m_outfile;

	string m_mapreads_file;

	string m_bowtie_segments_file;

	string m_fixtail_file;

	string m_fixhole_file;

	string m_fixhole_exceed_file;

	string m_bwtout;

    string m_fastq_file;

	string m_chrom_seq;

	vector<string> m_mapreads_files;

	size_t m_input_buffer_size;

	size_t m_seed_width;

	size_t m_anchor_width;

	size_t m_read_width;

	size_t m_max_mismatches;

	size_t m_max_mismatches_plus1;

	size_t m_extend_bits;

	size_t m_max_totalmismatch;

	size_t m_append_score;

	size_t m_max_fusion_mismatch;

	size_t m_flank_case;

	double m_max_rank;

	size_t m_left_mismatches;

	size_t m_total_chromosomes_size;

	size_t m_total_tissue_size;
	
	size_t m_max_splice_dist;

	size_t m_max_insert_dist;

	vector<GenomicRegion> m_hits;

	size_t m_interal;

	size_t m_computed_interal;

	size_t m_max_readID;

	size_t m_cur_chrom_size;

	TAC_POS_MAP m_tac_pos_map;

	TAC_POS_MAP m_tac_pos_map_rev;

	size_t m_boundary;

	size_t m_min_boundary;

	double m_total_time;

	double m_map_time;

	double m_checkseedtable_time;

	double m_matainslidewindow_time;

	double m_matainopenreads_time;

	double m_triggermatch_time;

	double m_checkcandclearopenread_time;

	double m_halftag_time;

	double m_output_time;

	double m_init_time;

	double m_any_time;

	bool m_islands_provide;

	string m_islands;

	string m_path;

	size_t m_prefix_match;

	vector<pair<size_t, size_t> >* m_cur_junc_vec_ptr; 
	
	vector<pair<size_t, size_t> >::iterator m_cur_junc_head_iter, m_cur_junc_tail_iter, m_cur_junc_end_iter;

	vector<PrefixMatchSeed> m_slide_window;

	vector<PrefixMatchSeedReady> m_slide_window_ready;

	hash_map<string, string> m_loaded_chromos;

	const size_t m_window_size;

	bool m_reverse_flag;

	size_t m_score;

	size_t m_num_seg;

	size_t m_seg_len;

	//local variables used in MapReads, need to set free after MapReads

	//use new OpenReadList
	vector<OpenReadList> m_open_reads;

	vector<OpenReadList> m_open_reads_prefix;

	vector<OpenReadList> m_open_reads_suffix;

	//hash_map<size_t, WordPair> m_prefix_words;

	//vector <JUNC_HASH> m_conj_hashs;

	ofstream m_ofs_map_reads;

	ofstream m_ofs_allmapreads;

	ofstream m_ofs_fix_hmer;

	ofstream m_ofs_fix_hmer_hole;

	ofstream m_ofs_fix_hole;

	//current local variables
	size_t m_seed_key;

	size_t m_seed_keyII;

	vector<size_t> m_kmer_keys;

	vector<size_t> m_kmer_keys_upper;

	vector<size_t> m_kmer_keys_lower;

	//vector<size_t> m_imer_keys;

	WordPair m_wp, m_wp_reserve, m_shifted_wp;

	WordPair m_five_prim_suffix, m_three_prim_prefix;

	size_t m_chrom_offset, m_file_offset, m_max_doner_offset, m_max_acceptor_offset;

	size_t m_matched_flank, m_matched_bads;

	// bit-vector indicating the location of currently bad bases
	size_t m_bad_bases;

	size_t m_bad_bases_reserve;

	size_t m_candidate_count;

	size_t m_basecount;

	size_t m_enterhalf_count;

	size_t m_passhalf_count;

	size_t m_cur_island_start;

	size_t m_cur_island_end;

	PAIR_SIZE_T_VEC_ITER m_cur_tac_pos, m_tac_pos_end, m_tac_pos_begin;

	string m_qual_str;

	size_t m_num_anchors;

	size_t m_kmer_step;

	size_t m_threads;

	size_t m_segment_mask;

	size_t m_is_fix_tail_head;

};

//GenomeScan::GenomeScan(size_t input_buffer_size, size_t seed_width, size_t read_width, size_t max_mismatches) : m_input_buffer_size(input_buffer_size),
//			m_seed_width(seed_width), m_read_width(read_width), m_max_mismatches(max_mismatches), m_seed_hash(SeedHash(m_seed_width, m_read_width, m_max_mismatches)), 
//			m_masks(Masks(read_width, seed_width, max_mismatches, 6, 1)), m_seed_key(0), m_chrom_offset(0), m_file_offset(0), m_total_chromosomes_size(0), m_ambiguous_file(0), m_bad_bases(ALL_BITS_ON), 
//			m_bad_bases_reserve(ALL_BITS_ON), m_window_size(0)
//{
//}

GenomeScan::GenomeScan(const CmdlineArgus& argulist) : m_input_buffer_size(argulist.input_buffer_size), m_seed_width(argulist.seed_width), m_read_width(argulist.read_width), 
			m_max_mismatches(argulist.max_mismatches), m_seed_hash(SeedHash(argulist.seed_width, argulist.read_width, argulist.max_mismatches)), m_file_offset(0),
			m_masks(Masks(argulist.read_width, argulist.seed_width, argulist.max_mismatches, argulist.anchro_size, argulist.num_anchor, argulist.num_seg, argulist.seg_len, argulist.extend_bits)), m_seed_key(0), m_chrom_offset(0), 
			m_ambiguous_file(argulist.ambiguous_file), m_outfile(argulist.outfile),	m_bad_bases(ALL_BITS_ON), m_bad_bases_reserve(ALL_BITS_ON), m_total_chromosomes_size(0),
			m_max_splice_dist(argulist.max_splice_dist), m_interal(argulist.interal), m_num_anchors(argulist.num_anchor), m_kmer_step(argulist.anchro_size),
			m_max_readID(argulist.max_readID), m_boundary(argulist.seg_len + argulist.seg_len + 2), m_matched_flank(0), m_matched_bads(0), /*m_slide_window(argulist.max_splice_dist), */m_slide_window_ready(argulist.max_splice_dist), 
			m_window_size(argulist.max_splice_dist), m_min_boundary(argulist.seg_len + 2), m_flank_case(argulist.flankstr), m_max_rank(argulist.rank), m_qual_str(argulist.read_width, 'I'),
			m_map_time(0), m_output_time(0), m_total_time(0), m_prefix_match(argulist.prefix_match), m_total_tissue_size(0), m_anchor_width(argulist.anchro_size), m_kmer_keys(argulist.num_anchor, 0),
			m_left_mismatches(argulist.max_mismatches), m_candidate_count(0), m_basecount(0), m_halftag_time(0),m_enterhalf_count(0), m_passhalf_count(0), m_triggermatch_time(0),
			m_checkseedtable_time(0), m_matainslidewindow_time(0), m_matainopenreads_time(0), m_checkcandclearopenread_time(0),m_cur_island_start(0), m_cur_island_end(max_size_t)
			, m_num_seg(argulist.num_seg), m_seg_len(argulist.seg_len), m_bwtout(argulist.bwtout), m_extend_bits(argulist.extend_bits), m_max_totalmismatch(argulist.total_mismatch),
			m_kmer_keys_upper(argulist.num_anchor, 0), m_kmer_keys_lower(argulist.num_anchor, 0), m_max_fusion_mismatch(argulist.total_fusion_mismatch), m_append_score(argulist.append_score),
			m_max_doner_offset(0), m_max_acceptor_offset(0), m_threads(argulist.threads), m_max_mismatches_plus1(argulist.max_mismatches + 1), /*m_imer_keys(argulist.num_anchor, 0),*/
			m_max_insert_dist(argulist.max_insertion)
{
	// Get the reads

	cout <<"genome arguments initialized"<<endl;
	if (argulist.path != 0)
	{
		string mkdir_cmdstr = "mkdir -p ";
		mkdir_cmdstr.append(argulist.path);
		system(mkdir_cmdstr.c_str());
		m_path = argulist.path;
		m_path.append("/");
	}
	else
		m_path = "";

	//ofstream time_ofs;

	clock_t t1, t2;

//#ifdef TIME2
//	string fcounttime = m_path;
//	fcounttime.append("mps_time_log.txt");
//
//	time_ofs.open(fcounttime.c_str(), ios_base::app);
//#endif

#ifdef TIME2
	t1=clock();
#endif

	//m_max_totalmismatch = argulist.num_seg;

	string tempchromdir = argulist.tmp_dir; tempchromdir.append("/tempchrom/");
	string mktempchromdir = "mkdir -p ";mktempchromdir.append(tempchromdir);
	system(mktempchromdir.c_str());

	OpenReadList open_reads(argulist.max_splice_dist, (size_t) pow((double)4,(double)argulist.anchro_size));

	for (size_t i = 0; i < argulist.num_anchor; ++i)
	{
		//m_open_reads.push_back(open_reads);

		m_open_reads_prefix.push_back(open_reads);

		m_open_reads_suffix.push_back(open_reads);
	}

	open_reads.clear();

	cout << "set chrom file\n";

	SetChromFiles(argulist.filenames_file, argulist.chrom_file, argulist.fasta_suffix, tempchromdir);

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "initial time: " << m_any_time << endl;
#endif

#ifdef TIME2
	t1=clock();
#endif

	if (argulist.fusion_file != 0)
	{
		cout << "fix fusion"<<endl;

		vector<string> fusion_files;

		read_dir(argulist.fusion_file, "txt", fusion_files);

		sort(fusion_files.begin(), fusion_files.end());

		string fixed_junc = argulist.fusion_junc;

		CHROM_FUSION_JUNC_HASH_COMB conj_hash_map;

		string chromseq1 = "";
		string chromseq2 = "";
		string prevchrom_file1 = "";
		string prevchrom_file2 = "";

		for (size_t i = 0; i < fusion_files.size(); ++i)
		{
			cout << fusion_files[i] << endl;
			string chromnamebase = basename2(fusion_files[i]);

			char chromname1[1000], chromname2[1000];

			sscanf(chromnamebase.c_str(), "%[^_]_%[^_]", chromname1, chromname2);

			if (m_chrom_names_map.find(chromname1) == m_chrom_names_map.end())
			{
				cout <<"can't find chrom " << chromname1 << endl;
				continue;
			}

			if (m_chrom_names_map.find(chromname2) == m_chrom_names_map.end())
			{
				cout <<"can't find chrom " << chromname1 << endl;
				continue;
			}

			string chromfile1 = m_chrom_files[m_chrom_names_map[chromname1]];
			string chromfile2 = m_chrom_files[m_chrom_names_map[chromname2]];

			//if (m_loaded_chromos.find(chromname1) == m_loaded_chromos.end())
			//	cout << "can't find chromo\t"<<chromname1<<endl;

			//if (m_loaded_chromos.find(chromname2) == m_loaded_chromos.end())
			//	cout << "can't find chromo\t"<<chromname2<<endl;

			//string& chromseq1 = m_loaded_chromos[chromname1];
			//string& chromseq2 = m_loaded_chromos[chromname2];

			cout << "read chrom"<<endl;

			if (prevchrom_file1 == "" || chromfile1 != prevchrom_file1)
			{
				chromseq1.clear();
				cout << "read chromo\t"<<chromfile1<<endl;
				readchrom(chromfile1.c_str(), chromseq1);
				if (chromseq1 == "")
				{
					cout <<"empty chrom: "<<chromfile1<<endl;
					exit(1);
				}
				prevchrom_file1 = chromfile1;
			}

			if (prevchrom_file2 == "" || chromfile2 != prevchrom_file2)
			{
				chromseq2.clear();
				cout << "read chromo\t"<<chromfile2<<endl;
				readchrom(chromfile2.c_str(), chromseq2);
				if (chromseq2 == "")
				{
					cout <<"empty chrom: "<<chromfile2<<endl;
					exit(1);
				}
				prevchrom_file2 = chromfile2;
			}

			string fixed_fusion_file = fusion_files[i]; fixed_fusion_file.append(".fixed");

			//LoadChromFiles(argulist.filenames_file, argulist.chrom_file, argulist.fasta_suffix);

			FixFusionReady(fusion_files[i].c_str(), fixed_fusion_file.c_str(), conj_hash_map, false, false, m_flank_case, chromseq1, chromseq2);			
		}

		ofstream ofs_fixed_junc(fixed_junc.c_str());

		WriteJuncFusion(conj_hash_map, ofs_fixed_junc);

		exit(0);
	}

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "fusion time: " << m_any_time << endl;
#endif

	//SetTissueFiles(argulist.tissuenames_file, argulist.tissue_file, argulist.fasta_suffix);

	//if (argulist.fastq_file != 0)
	//{
	//	m_fastq_file = argulist.fastq_file;

	//	GetReadsFastqQualitySeq(argulist.fastq_file);
	//}

#ifdef TIME2
	t1=clock();
#endif

	if (argulist.TACs_pos != 0)
		ReadAndMergeTacPos(argulist.TACs_pos);

	m_reverse_flag = false;

	string tmp_dir = argulist.tmp_dir;

	string hmer_dir = tmp_dir; hmer_dir.append("single_anchored_middle/");
	string hole_dir = tmp_dir; hole_dir.append("double_anchored/");
	string head_dir = tmp_dir; head_dir.append("single_anchored_head/");
	string tail_dir = tmp_dir; tail_dir.append("single_anchored_tail/");

	vector<string> hmer_tobefixed, hole_tobefixed, head_tobefixed, tail_tobefixed;

	read_dir(hmer_dir, "txt", hmer_tobefixed);
	read_dir(hole_dir, "txt", hole_tobefixed);
	read_dir(head_dir, "txt", head_tobefixed);
	read_dir(tail_dir, "txt", tail_tobefixed);

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "read island and directory time: " << m_any_time << endl;
#endif

	bool use_island = false;

	if (argulist.prefix_match)
		use_island = true;

#ifdef TIME2
	t1=clock();
#endif
	
	hash_map<string, vector<pair<size_t, size_t> > > chrom_junc_map;

	vector<pair<size_t, size_t> > sentinel;

	cout << "fix hole"<<endl;

//typedef hash_map<size_t, JuncInfo> JUNC_HASH_COMB;
//
//typedef hash_map<string, JUNC_HASH_COMB> CHROM_JUNC_HASH_COMB;

	for (size_t i = 0; i < hole_tobefixed.size(); ++i)
	{
		string fixed_fixhole_file = hole_tobefixed[i]; fixed_fixhole_file.append(".fixed");

		string chromnamebase = basename2(hole_tobefixed[i]);

		if (m_chrom_names_map.find(chromnamebase) == m_chrom_names_map.end())
			continue;

		vector<pair<size_t, size_t> >& junc_vec = chrom_junc_map[chromnamebase];

		hash_set<size_t> junc_set;

		if (m_flank_case >= 5)
			FixHoleReady(hole_tobefixed[i].c_str(), fixed_fixhole_file.c_str(), junc_set, m_flank_case, false, true, true);
		else
			FixHoleReadyNonCanon(hole_tobefixed[i].c_str(), fixed_fixhole_file.c_str(), junc_set, m_flank_case, false, true, true);

		hash_set<size_t>::iterator junc_set_iter;

		for (junc_set_iter = junc_set.begin(); junc_set_iter != junc_set.end(); ++junc_set_iter)
		{
			size_t comb_offset = *junc_set_iter;

			size_t prefix_end = comb_offset >> THIRTY_TWO;

			size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

			junc_vec.push_back(make_pair(prefix_end, suffix_st));
		}

		sort(junc_vec.begin(), junc_vec.end(), compare_pair);

		//if (m_threads > 1)
		//{
		//string sort_cmd = "sort -t_ -k1,1n -S 3500000 -o ";

		//sort_cmd.append(fixed_fixhole_file);

		//sort_cmd.append(" -T temp ");

		//sort_cmd.append(fixed_fixhole_file);

		//system(sort_cmd.c_str());
		//}
	}

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "double anchor time: " << m_any_time << endl;
#endif


#ifdef TIME2
	t1=clock();
#endif

	cout << "fix tail"<<endl;

	m_is_fix_tail_head = true;

	for (size_t i = 0; i < tail_tobefixed.size(); ++i)
	{
		string chromnamebase = basename2(hmer_tobefixed[i]);

		if (m_chrom_names_map.find(chromnamebase) == m_chrom_names_map.end())
			continue;

		m_cur_junc_vec_ptr = &sentinel; 

		hash_map<string, vector<pair<size_t, size_t> > >::iterator cur_chrom_iter = chrom_junc_map.find(chromnamebase);

		if (cur_chrom_iter != chrom_junc_map.end())
			m_cur_junc_vec_ptr = &(cur_chrom_iter->second);

		string fixed_fixhmer_file = tail_tobefixed[i]; fixed_fixhmer_file.append(".fixed");

		FixHmerReady(tail_tobefixed[i].c_str(),fixed_fixhmer_file.c_str(), 5, false, true);

		Clear();

		//if (m_threads > 1)
		//{
		string sort_cmd = "sort -t~ -k1,1n -S 3500000 -o ";

		sort_cmd.append(fixed_fixhmer_file);

		sort_cmd.append(" -T temp ");

		sort_cmd.append(fixed_fixhmer_file);

		system(sort_cmd.c_str());
		//}
	}

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "single anchor tail time: " << m_any_time << endl;
#endif

#ifdef TIME2
	t1=clock();
#endif

	cout << "fix hmer"<<endl;

	m_is_fix_tail_head = false;

	for (size_t i = 0; i < hmer_tobefixed.size(); ++i)
	{
		string chromnamebase = basename2(hmer_tobefixed[i]);

		if (m_chrom_names_map.find(chromnamebase) == m_chrom_names_map.end())
			continue;

		m_cur_junc_vec_ptr = &sentinel; 

		hash_map<string, vector<pair<size_t, size_t> > >::iterator cur_chrom_iter = chrom_junc_map.find(chromnamebase);

		if (cur_chrom_iter != chrom_junc_map.end())
			m_cur_junc_vec_ptr = &(cur_chrom_iter->second);

		string fixed_fixhmer_file = hmer_tobefixed[i]; fixed_fixhmer_file.append(".fixed");

		FixHmerReady(hmer_tobefixed[i].c_str(), fixed_fixhmer_file.c_str(), 5, false, false);

		Clear();

		//if (m_threads > 1)
		//{
		string sort_cmd = "sort -t~ -k1,1n -S 3500000 -o ";

		sort_cmd.append(fixed_fixhmer_file);

		sort_cmd.append(" -T temp ");

		sort_cmd.append(fixed_fixhmer_file);

		system(sort_cmd.c_str());
		//}
	}

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "single anchor middle time: " << m_any_time << endl;
#endif

#ifdef TIME2
	t1=clock();
#endif

	cout << "fix head"<<endl;

	m_is_fix_tail_head = true;

	for (size_t i = 0; i < head_tobefixed.size(); ++i)
	{
		string chromnamebase = basename2(hmer_tobefixed[i]);

		if (m_chrom_names_map.find(chromnamebase) == m_chrom_names_map.end())
			continue;

		m_cur_junc_vec_ptr = &sentinel; 

		hash_map<string, vector<pair<size_t, size_t> > >::iterator cur_chrom_iter = chrom_junc_map.find(chromnamebase);

		if (cur_chrom_iter != chrom_junc_map.end())
			m_cur_junc_vec_ptr = &(cur_chrom_iter->second);

		sort(m_cur_junc_vec_ptr->begin(), m_cur_junc_vec_ptr->end(), compare_pair_rev);

		//for (size_t i = 0; i < m_cur_junc_vec_ptr->size(); ++i)
		//	cout << (*m_cur_junc_vec_ptr)[i].first <<'\t'<< (*m_cur_junc_vec_ptr)[i].second<<endl;

		string fixed_fixhmer_file = head_tobefixed[i]; fixed_fixhmer_file.append(".fixed");		

		FixHmerReady_Rev(head_tobefixed[i].c_str(),fixed_fixhmer_file.c_str(), 5, false, true);

		Clear();

		//if (m_threads > 1)
		//{
		string sort_cmd = "sort -t~ -k1,1n -S 3500000 -o ";

		sort_cmd.append(fixed_fixhmer_file);

		sort_cmd.append(" -T temp ");

		sort_cmd.append(fixed_fixhmer_file);

		system(sort_cmd.c_str());
		//}
	}

#ifdef TIME2
	t2=clock();

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cout << "single anchor head time: " << m_any_time << endl;
#endif

	exit(0);

	if (argulist.bwtout != 0)
	{
		m_fixhole_file = argulist.bwtout;
		m_fixhole_file.append(".fixhole");

		m_fixhole_exceed_file = argulist.bwtout;
		m_fixhole_exceed_file.append(".fixholeexceed");

		cout << "TryAppend .mistailorhead"<<endl;
	
		InitReads(0);

		InitReadWords();

		//string fixed_fixhole_file = "test_fix_hole.txt";

		//FixHoleExt(m_fixhole_file.c_str(), fixed_fixhole_file.c_str(), m_flank_case, false, true, false);

		//exit(0);

		//cout << "fix fusion"<<endl;

		//string fixhole_fusion_file = m_bwtout;
		//fixhole_fusion_file.append(".fusion");

		//string fixhole_fusion_fixed_file = m_path;
		//fixhole_fusion_fixed_file.append(m_original_dir);
		//fixhole_fusion_fixed_file.append("fusion.mappedreads");

		//string fusion_junction_file = m_path;
		//fusion_junction_file.append(m_original_dir);
		//fusion_junction_file.append("fusion_junctions.txt");
		//ofstream fusion_junc_ofs(fusion_junction_file.c_str());

		//CONJFusion_HASH_MAP conj_hash_map;

		//FixHole_fusion(fixhole_fusion_file.c_str(), fixhole_fusion_fixed_file.c_str(), conj_hash_map, true, true);

		//exit(0);

		cout << "reads initialed"<<endl;

		string mistailorhead = m_bwtout;
		mistailorhead.append(".mistailorhead");
		TryAppend(mistailorhead.c_str(), mistailorhead_type);

		cout << "TryAppend .mistailandhead"<<endl;
		string mistailandhead = m_bwtout;
		mistailandhead.append(".mistailandhead");
		TryAppend(mistailandhead.c_str(), mistailandhead_type);

		cout << "read fix hmer"<<endl;

		InitReadWords();

		cout << "read words initialed"<<endl;

#ifndef DEBUG_CHECK
		m_reads.clear();
#endif

		ReadFixHmer();
	}

	
}

inline size_t 
GenomeScan::PrefixMatch(const WordPair& read_word)
{
	return read_word.score(m_wp, m_masks.first_half_mask);
}

void 
GenomeScan::SetReadWords()
{
	for (size_t i = 0; i < m_reads.size(); ++i)
	{
		m_read_words.push_back(WordPair(m_reads[i]));
		m_read_kmers.push_back(Str2Kmer(m_reads[i]));
		m_read_words_rc.push_back(WordPair(revcomp(m_reads[i])));
		m_read_kmers_rc.push_back(Str2Kmer(revcomp(m_reads[i])));
	}
}

void 
GenomeScan::SetReadWords_Rev()
{
	for (size_t i = 0; i < m_reads.size(); ++i)
	{
		string revreads = m_reads[i];
		reverse(revreads.begin(), revreads.end());
		m_read_words.push_back(WordPair(revreads));
		m_read_kmers.push_back(Str2Kmer(revreads));
		m_read_words_rc.push_back(WordPair(revcomp(revreads)));
		m_read_kmers_rc.push_back(Str2Kmer(revcomp(revreads)));
	}
}

inline void 
GenomeScan::EnterSlideWindowFixTail(vector<SamFixtail>::const_iterator ft_seed_begin, vector<SamFixtail>::const_iterator ft_seed_end)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;

	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		m_slide_window[modidx].matchedTagsHalf.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	vector<SamFixtail>::const_iterator sv_iter;
	for (sv_iter = ft_seed_begin; sv_iter != ft_seed_end; ++sv_iter)
	{
		unsigned int firstHalfScore = 0;

		if (firstHalfScore <= m_max_mismatches)
		{
			m_slide_window[modidx].matchedTagsHalf.push_back(TagSeed(firstHalfScore, sv_iter->strand, sv_iter->readID, sv_iter->offset, sv_iter->splitway));

			size_t kmer_need = m_max_mismatches - firstHalfScore;

			if (sv_iter->strand)
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers[sv_iter->readID].size()*/m_num_anchors ; ++i)
				{
					if (!m_read_kmers[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}
			else
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
				{
					if (!m_read_kmers_rc[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers_rc[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}
		}
	}

	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window[modidx].five_prim_suffix;

		prefix_wp_ref = m_wp;

		five_prim_suffix_ref = m_wp;
		//- m_max_mismatches
		if (m_anchor_width >=4 && m_anchor_width <= 8)
			five_prim_suffix_ref.right_shift(m_anchor_width - 4);//>>m_anchor_width-2 <<2
		else
		{
			cout << "anchor width must greater or equal to 4 and less or equal to 8";
			exit(0);
		}

		prefix_wp_ref.left_shift(m_masks.left_shift_seg_width);

		if (m_islands_provide)
			m_slide_window[modidx].cur_tac_pos = m_cur_tac_pos;
		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}
}

inline void 
GenomeScan::EnterSlideWindowFixHmerReady(vector<TobeFixedHmerReady>& cur_tobe_fixed_hmer_segments, const char* chromname, bool avoid)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;

	if (m_slide_window_ready[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		m_slide_window_ready[modidx].matchedTagsHalf.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	vector<TobeFixedHmerReady>::iterator sv_iter;

	size_t st_offset;// = cur_tobe_fixed_hmer_segments.front().st_offset;

	if (m_is_fix_tail_head)
	{
		st_offset = cur_tobe_fixed_hmer_segments.front().st_offset;

		while (m_cur_junc_head_iter != m_cur_junc_end_iter && m_cur_junc_head_iter->first < st_offset - 1)
			++m_cur_junc_head_iter;

		while (m_cur_junc_tail_iter != m_cur_junc_end_iter && m_cur_junc_tail_iter->first < st_offset + m_seg_len)
			++m_cur_junc_tail_iter;
	}

	for (sv_iter = cur_tobe_fixed_hmer_segments.begin(); sv_iter != cur_tobe_fixed_hmer_segments.end(); ++sv_iter)
	{
		if (/*sv_iter->end_offset == 0*/m_is_fix_tail_head)
		{
			if (sv_iter->strand == '+' && sv_iter->st_seg_id != m_num_seg)
			{
				continue;
			}
			else if (sv_iter->strand == '-' && sv_iter->st_seg_id != 1)
			{
				continue;
			}
		}

		if (m_is_fix_tail_head)
		{
			bool fixed = false;

			if (m_cur_junc_tail_iter > m_cur_junc_head_iter)
			{
				vector<pair<size_t, size_t> >::iterator junc_iter, min_junc_iter;

				size_t min_mismatch = -1;

				size_t max_loc;

				size_t min_rbits;

				string matched_chromstr;

				for (junc_iter = m_cur_junc_head_iter; junc_iter < m_cur_junc_tail_iter; ++junc_iter)
				{
					size_t prefixlen = junc_iter->first - st_offset + 1;

					size_t suffixlen = m_seg_len - prefixlen;

					string chromseq = m_chrom_seq.substr(st_offset, prefixlen) + m_chrom_seq.substr(junc_iter->second, suffixlen);

					WordPair chromwp(chromseq);

					size_t rbits;

					//size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

					//if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
					//	m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

					//hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

					//cur_masks_ptr = &(cur_hash_map_mask_iter->second);

					size_t score = chromwp.score(sv_iter->to_be_fixed_hmer_segments_wp.front().seg_wp, m_masks.suffix_seg_bits_on, rbits);

					if (score < min_mismatch || (score == min_mismatch && (junc_iter->second + min_junc_iter->first) < (min_junc_iter->second + junc_iter->first)))
					{
						min_mismatch = score;

						min_junc_iter = junc_iter;

						max_loc = prefixlen;

						min_rbits = rbits;

						matched_chromstr = chromseq;
					}
				}

				if (min_mismatch <= m_max_mismatches_plus1)
				{
					size_t strand = sv_iter->strand == '+' ? 0 : 16;

					string segment_str = sv_iter->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, m_seg_len);

					m_ofs_fix_hmer << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
						max_loc<<'M'<<min_junc_iter->second - min_junc_iter->first - 1<< 'N' <<m_seg_len - max_loc << "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					/*cout << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
						max_loc<<'M'<<min_junc_iter->second - min_junc_iter->first - 1<< 'N' <<m_seg_len - max_loc << "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';*/

					//m_ofs_fix_hmer_hole << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
					//	max_loc<<'M'<<min_junc_iter->second - min_junc_iter->first - 1<< 'N' <<m_seg_len - max_loc << "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					size_t selector_bit = LEAST_SIG_BIT << (m_seg_len - 1);

					size_t pre_mapped_len;

					if (sv_iter->strand == '-')
						pre_mapped_len = (m_num_seg - sv_iter->st_seg_id) * m_seg_len;
					else
						pre_mapped_len = (sv_iter->st_seg_id - 1) * m_seg_len;

					for (size_t i = 0; i < m_seg_len; ++i)
					{
						if ((selector_bit >> i) & min_rbits)
						{
							m_ofs_fix_hmer << pre_mapped_len + i<<':'<<matched_chromstr[i]<<'>'<<segment_str[i]<<',';

							//cout << pre_mapped_len + i<<':'<<matched_chromstr[i]<<'>'<<segment_str[i]<<',';

							//m_ofs_fix_hmer_hole << pre_mapped_len + i<<':'<<matched_chromstr[i]<<'>'<<segment_str[i]<<',';
						}
					}

					m_ofs_fix_hmer << endl;

					fixed = true;

					//m_ofs_fix_hmer_hole << endl;

					//m_ofs_fix_hmer.flush();

					//continue;
				}
				
			}

			//else  try append
			{
				string chromseq = m_chrom_seq.substr(st_offset, m_seg_len);

				WordPair chromwp(chromseq);

				size_t rbits;

				size_t score = chromwp.score(sv_iter->to_be_fixed_hmer_segments_wp.front().seg_wp, m_masks.suffix_seg_bits_on, rbits);

				if (score <= m_max_mismatches_plus1)
				{
					size_t strand = sv_iter->strand == '+' ? 0 : 16;

					string segment_str = sv_iter->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, m_seg_len);

					m_ofs_fix_hmer << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
						m_seg_len<<"M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<score<<'\t';

					//cout << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
					//	m_seg_len<<"M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<score<<'\t';

					//m_ofs_fix_hmer_hole << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
					//	max_loc<<'M'<<min_junc_iter->second - min_junc_iter->first - 1<< 'N' <<m_seg_len - max_loc << "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					size_t selector_bit = LEAST_SIG_BIT << (m_seg_len - 1);

					size_t pre_mapped_len;

					if (sv_iter->strand == '-')
						pre_mapped_len = (m_num_seg - sv_iter->st_seg_id) * m_seg_len;
					else
						pre_mapped_len = (sv_iter->st_seg_id - 1) * m_seg_len;

					for (size_t i = 0; i < m_seg_len; ++i)
					{
						if ((selector_bit >> i) & rbits)
						{
							m_ofs_fix_hmer << pre_mapped_len + i<<':'<<chromseq[i]<<'>'<<segment_str[i]<<',';

							//cout << pre_mapped_len + i<<':'<<matched_chromstr[i]<<'>'<<segment_str[i]<<',';

							//m_ofs_fix_hmer_hole << pre_mapped_len + i<<':'<<matched_chromstr[i]<<'>'<<segment_str[i]<<',';
						}
					}

					m_ofs_fix_hmer << endl;

					fixed = true;

					//m_ofs_fix_hmer_hole << endl;

					//m_ofs_fix_hmer.flush();

					//continue;
				}
			}

			//if (fixed)
			//	continue;
		}

		//m_cur_junc_head_iter = m_cur_junc_vec_ptr->begin(); 

		//m_cur_junc_tail_iter = m_cur_junc_vec_ptr->begin();

		//m_cur_junc_end_iter = m_cur_junc_vec_ptr->end();
			
		/*unsigned int firstHalfScore = 0;*/

		//string& tagnamestr = (*sv_iter)->tagname;
		//size_t last_idx = tagnamestr.find_last_of("_");
		//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		//if (avoid && m_reads_count[tagidx])
		//{
		//	//cout << "skip read "<< tagnamestr<<endl;
		//	continue;
		//}

		//if (firstHalfScore <= m_max_mismatches)
		//{
		size_t/* readID,*/ bef_offset;
		//readID = m_reverse_flag ? (*sv_iter)->readID_rev : (*sv_iter)->readID;

		//if (m_reverse_flag)
		//{
		//	bef_offset = sv_iter->st_offset == 0 ? -1 : m_cur_chrom_size - sv_iter->st_offset;
		//}
		//else
		//{
		bef_offset = sv_iter->end_offset == 0 ? -1 : sv_iter->end_offset;

		if (sv_iter->st_offset > m_max_doner_offset)
			m_max_doner_offset = sv_iter->st_offset;

		if (sv_iter->end_offset == 0)
		{
			size_t offset_pass_doner = m_max_doner_offset + m_max_splice_dist + 1000;

			if (m_max_acceptor_offset < offset_pass_doner)
				m_max_acceptor_offset = offset_pass_doner;
		}
		else
		{
			size_t offset_pass_doner = m_max_doner_offset + m_max_splice_dist + 1000;

			size_t offset_pass_end = sv_iter->end_offset + 1000;

			if (offset_pass_end < offset_pass_doner)
			{
				if (m_max_acceptor_offset < offset_pass_end)
					m_max_acceptor_offset = offset_pass_end;
			}
			else
			{
				if (m_max_acceptor_offset < offset_pass_doner)
					m_max_acceptor_offset = offset_pass_doner;
			}
		}

		//}

		m_slide_window_ready[modidx].matchedTagsHalf.push_back(*sv_iter);

		//size_t kmer_need = m_max_mismatches/* - firstHalfScore*/;

		//if (sv_iter->strand == '+')
		//{

		vector<Kmer> prefix_hmers;
		vector<Kmer> suffix_hmers;

		WordPairHmer& seg_wp_hmer = sv_iter->to_be_fixed_hmer_segments_wp.front();

		for (size_t i = 0; i < m_num_anchors; ++i)
		{
			suffix_hmers.push_back(seg_wp_hmer.seg_wp.get_kmer(1 + (i * m_anchor_width), m_anchor_width + (i * m_anchor_width)));
		}

		if (sv_iter->to_be_fixed_hmer_segments_wp.size() > 1)
		{
			WordPairHmer& seg_wp_hmer2 = sv_iter->to_be_fixed_hmer_segments_wp[1];
			for (size_t i = 0; i < m_num_anchors; ++i)
			{
				prefix_hmers.push_back(seg_wp_hmer2.seg_wp.get_kmer(seg_wp_hmer2.seg_len - m_anchor_width + 1 - (i * m_anchor_width), seg_wp_hmer2.seg_len - (i * m_anchor_width)));
			}
		}

		for (size_t i = 0; i <= m_max_mismatches && i < m_num_anchors; ++i)
		{
			if (suffix_hmers[i].bad)
			{
				m_open_reads_suffix[i].enter(suffix_hmers[i].kmer, m_chrom_offset + 1, m_slide_window_ready[modidx].matchedTagsHalf.size() - 1, sv_iter->strand, bef_offset);
			}

			if (sv_iter->to_be_fixed_hmer_segments_wp.size() > 1)
			{
				if (prefix_hmers[i].bad)
				{
					m_open_reads_prefix[i].enter(prefix_hmers[i].kmer, m_chrom_offset + 1, m_slide_window_ready[modidx].matchedTagsHalf.size() - 1, sv_iter->strand, bef_offset);
				}
			}
		}
		//}
		//else
		//{
		//	for (size_t i = 0; i <= m_max_mismatches && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
		//	{
		//		if (!m_read_kmers_rc[readID][i].bad)
		//		{
		//			m_open_reads[i].enter(m_read_kmers_rc[readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter)->strand, bef_offset);
		//		}
		//	}
		//}
		//}
	}

	if (m_slide_window_ready[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window_ready[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window_ready[modidx].five_prim_suffix;

		//WordPair& prefix_wp_ref2 = m_slide_window_ready[modidx].prefix_wp2;

		//WordPair& five_prim_suffix_ref2 = m_slide_window_ready[modidx].five_prim_suffix2;

		//string temp = m_wp_reserve.tostring3(ALL_BITS_ON, SIXTY_FOUR);

		prefix_wp_ref = m_wp_reserve;

		five_prim_suffix_ref = m_wp_reserve;

		//prefix_wp_ref2 = m_wp;

		//five_prim_suffix_ref2 = m_wp;

		//- m_max_mismatches
		//if (m_anchor_width > 8)
		//{
		//	cout << "anchor width must be less or equal to 8";
		//	exit(0);
		//}

		//five_prim_suffix_ref.right_shift(m_anchor_width - 4);//>>m_anchor_width-2 <<2

		prefix_wp_ref.left_shift(m_seg_len - 2);

		//prefix_wp_ref2.left_shift(m_seg_len - 2);

		if (m_islands_provide)
			m_slide_window_ready[modidx].cur_tac_pos = m_cur_tac_pos;
		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}
}


inline void 
GenomeScan::EnterSlideWindowFixHmerReady_Rev(vector<TobeFixedHmerReady>& cur_tobe_fixed_hmer_segments, const char* chromname, bool avoid)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;

	if (m_slide_window_ready[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		m_slide_window_ready[modidx].matchedTagsHalf.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	vector<TobeFixedHmerReady>::iterator sv_iter;

	size_t st_offset;// = cur_tobe_fixed_hmer_segments.front().st_offset;

	if (m_is_fix_tail_head)
	{
		st_offset = cur_tobe_fixed_hmer_segments.front().st_offset;

		while (m_cur_junc_head_iter != m_cur_junc_end_iter && m_cur_junc_head_iter->second > st_offset + 1)
			++m_cur_junc_head_iter;

		while (m_cur_junc_tail_iter != m_cur_junc_end_iter && m_cur_junc_tail_iter->second + m_seg_len > st_offset)
			++m_cur_junc_tail_iter;
	}

	for (sv_iter = cur_tobe_fixed_hmer_segments.begin(); sv_iter != cur_tobe_fixed_hmer_segments.end(); ++sv_iter)
	{
		if (/*sv_iter->end_offset == 0*/m_is_fix_tail_head)
		{
			if (sv_iter->strand == '+' && sv_iter->st_seg_id != 1)
			{
				continue;
			}
			else if (sv_iter->strand == '-' && sv_iter->st_seg_id != m_num_seg)
			{
				continue;
			}
		}

		if (m_is_fix_tail_head)
		{
			bool fixed = false;

			if (m_cur_junc_tail_iter > m_cur_junc_head_iter)
			{
				vector<pair<size_t, size_t> >::iterator junc_iter, min_junc_iter;

				size_t min_mismatch = -1;

				size_t max_loc;

				size_t min_rbits;

				string matched_chromstr;

				for (junc_iter = m_cur_junc_head_iter; junc_iter < m_cur_junc_tail_iter; ++junc_iter)
				{
					size_t suffixlen = st_offset - junc_iter->second + 1;

					size_t prefixlen = m_seg_len - suffixlen;

					size_t rev_end = m_cur_chrom_size - junc_iter->first - 1;

					size_t rev_st = m_cur_chrom_size - junc_iter->second - suffixlen;

					string chromseq = m_chrom_seq.substr(rev_st, suffixlen) + m_chrom_seq.substr(rev_end, prefixlen);

					WordPair chromwp(chromseq);

					size_t rbits;

					//size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

					//if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
					//	m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

					//hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

					//cur_masks_ptr = &(cur_hash_map_mask_iter->second);

					size_t score = chromwp.score(sv_iter->to_be_fixed_hmer_segments_wp.front().seg_wp, m_masks.suffix_seg_bits_on, rbits);

					if (score < min_mismatch || (score == min_mismatch && (junc_iter->second + min_junc_iter->first) < (min_junc_iter->second + junc_iter->first)))
					{
						min_mismatch = score;

						min_junc_iter = junc_iter;

						max_loc = prefixlen;

						min_rbits = rbits;

						matched_chromstr = chromseq;
					}
				}

				if (min_mismatch <= m_max_mismatches_plus1)
				{
					size_t strand = sv_iter->strand == '+' ? 0 : 16;

					string segment_str = sv_iter->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, m_seg_len);

					string bef_rev_segment_str = segment_str;

					reverse(segment_str.begin(), segment_str.end());

					size_t intron_len = min_junc_iter->second - min_junc_iter->first - 1;

					m_ofs_fix_hmer << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset - m_seg_len - intron_len + 1<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
						max_loc<<'M'<<intron_len<< 'N' <<m_seg_len - max_loc<< "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					//m_ofs_fix_hmer_hole << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset - m_seg_len - intron_len + 1<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
					//	max_loc<<'M'<<intron_len<< 'N' <<m_seg_len - max_loc<< "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					size_t selector_bit = LEAST_SIG_BIT << (m_seg_len - 1);

					size_t pre_mapped_len;

					if (sv_iter->strand == '-')
						pre_mapped_len = (m_num_seg - sv_iter->st_seg_id) * m_seg_len;
					else
						pre_mapped_len = (sv_iter->st_seg_id - 1) * m_seg_len;

					for (int i = static_cast<int> (m_seg_len); i >= 0; --i)
					{
						if ((selector_bit >> i) & min_rbits)
						{
							m_ofs_fix_hmer << pre_mapped_len + m_seg_len - i - 1<<':'<<matched_chromstr[i]<<'>'<<bef_rev_segment_str[i]<<',';

							//m_ofs_fix_hmer_hole<< pre_mapped_len + m_seg_len - i - 1<<':'<<matched_chromstr[i]<<'>'<<bef_rev_segment_str[i]<<',';
						}
					}

					m_ofs_fix_hmer<< endl;

					fixed = true;

					//m_ofs_fix_hmer_hole<< endl;

					//continue;
				}
			}

			//else //try append
			{
				size_t rev_st = m_cur_chrom_size - st_offset - 1;

				string chromseq = m_chrom_seq.substr(rev_st, m_seg_len);

				WordPair chromwp(chromseq);

				size_t rbits;

				size_t score = chromwp.score(sv_iter->to_be_fixed_hmer_segments_wp.front().seg_wp, m_masks.suffix_seg_bits_on, rbits);

				if (score <= m_max_mismatches_plus1)
				{
					size_t strand = sv_iter->strand == '+' ? 0 : 16;

					string segment_str = sv_iter->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, m_seg_len);

					string bef_rev_segment_str = segment_str;

					reverse(segment_str.begin(), segment_str.end());

					//size_t intron_len = min_junc_iter->second - min_junc_iter->first - 1;

					m_ofs_fix_hmer << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset - m_seg_len + 1<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
						m_seg_len<< "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<score<<'\t';

					//m_ofs_fix_hmer_hole << sv_iter->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<st_offset - m_seg_len - intron_len + 1<<'\t' << sv_iter->st_seg_id<<'-'<<sv_iter->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
					//	max_loc<<'M'<<intron_len<< 'N' <<m_seg_len - max_loc<< "M\t" <<segment_str<<'\t'<<sv_iter->id1<<'\t'<<sv_iter->id2<<"\tNM:i:"<<min_mismatch<<'\t';

					size_t selector_bit = LEAST_SIG_BIT << (m_seg_len - 1);

					size_t pre_mapped_len;

					if (sv_iter->strand == '-')
						pre_mapped_len = (m_num_seg - sv_iter->st_seg_id) * m_seg_len;
					else
						pre_mapped_len = (sv_iter->st_seg_id - 1) * m_seg_len;

					for (int i = static_cast<int> (m_seg_len); i >= 0; --i)
					{
						if ((selector_bit >> i) & rbits)
						{
							m_ofs_fix_hmer << pre_mapped_len + m_seg_len - i - 1<<':'<<chromseq[i]<<'>'<<bef_rev_segment_str[i]<<',';

							//m_ofs_fix_hmer_hole<< pre_mapped_len + m_seg_len - i - 1<<':'<<matched_chromstr[i]<<'>'<<bef_rev_segment_str[i]<<',';
						}
					}

					m_ofs_fix_hmer<< endl;

					fixed = true;

					//m_ofs_fix_hmer_hole<< endl;

					//continue;
				}
			}

			//if (fixed)
			//	continue;
		}
			
		/*unsigned int firstHalfScore = 0;*/

		//string& tagnamestr = (*sv_iter)->tagname;
		//size_t last_idx = tagnamestr.find_last_of("_");
		//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		//if (avoid && m_reads_count[tagidx])
		//{
		//	//cout << "skip read "<< tagnamestr<<endl;
		//	continue;
		//}

		//if (firstHalfScore <= m_max_mismatches)
		//{
		size_t/* readID,*/ bef_offset;
		//readID = m_reverse_flag ? (*sv_iter)->readID_rev : (*sv_iter)->readID;

		//if (m_reverse_flag)
		//{
		bef_offset = sv_iter->end_offset == 0 ? -1 : m_cur_chrom_size - sv_iter->end_offset - 1;

		if (m_cur_chrom_size - sv_iter->st_offset >  m_max_doner_offset)
			m_max_doner_offset = m_cur_chrom_size - sv_iter->st_offset;

		if (sv_iter->end_offset == 0)
		{
			if (m_max_acceptor_offset < m_max_doner_offset + m_max_splice_dist + 1000)
				m_max_acceptor_offset = m_max_doner_offset + m_max_splice_dist + 1000;
		}
		else
		{
			if (m_max_acceptor_offset < m_cur_chrom_size - sv_iter->end_offset - 1 + 1000)
				m_max_acceptor_offset = m_cur_chrom_size - sv_iter->end_offset - 1 + 1000;
		}

		//}
		//else
		//{
		//	bef_offset = sv_iter->end_offset == 0 ? -1 : sv_iter->end_offset;
		//}

		m_slide_window_ready[modidx].matchedTagsHalf.push_back(*sv_iter);

		//size_t kmer_need = m_max_mismatches/* - firstHalfScore*/;

		//if (sv_iter->strand == '+')
		//{

		vector<Kmer> prefix_hmers;
		vector<Kmer> suffix_hmers;

		WordPairHmer& seg_wp_hmer = sv_iter->to_be_fixed_hmer_segments_wp.front();

		for (size_t i = 0; i < m_num_anchors; ++i)
		{
			suffix_hmers.push_back(seg_wp_hmer.seg_wp.get_kmer(1 + (i * m_anchor_width), m_anchor_width + (i * m_anchor_width)));
		}

		if (sv_iter->to_be_fixed_hmer_segments_wp.size() > 1)
		{
			WordPairHmer& seg_wp_hmer2 = sv_iter->to_be_fixed_hmer_segments_wp[1];
			for (size_t i = 0; i < m_num_anchors; ++i)
			{
				prefix_hmers.push_back(seg_wp_hmer2.seg_wp.get_kmer(seg_wp_hmer2.seg_len - m_anchor_width + 1 - (i * m_anchor_width), seg_wp_hmer2.seg_len - (i * m_anchor_width)));
			}
		}

		for (size_t i = 0; i <= m_max_mismatches && i < m_num_anchors; ++i)
		{
			if (suffix_hmers[i].bad)
			{
				m_open_reads_suffix[i].enter(suffix_hmers[i].kmer, m_chrom_offset + 1, m_slide_window_ready[modidx].matchedTagsHalf.size() - 1, sv_iter->strand, bef_offset);
			}

			if (sv_iter->to_be_fixed_hmer_segments_wp.size() > 1)
			{
				if (prefix_hmers[i].bad)
				{
					m_open_reads_prefix[i].enter(prefix_hmers[i].kmer, m_chrom_offset + 1, m_slide_window_ready[modidx].matchedTagsHalf.size() - 1, sv_iter->strand, bef_offset);
				}
			}
		}
		//}
		//else
		//{
		//	for (size_t i = 0; i <= m_max_mismatches && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
		//	{
		//		if (!m_read_kmers_rc[readID][i].bad)
		//		{
		//			m_open_reads[i].enter(m_read_kmers_rc[readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter)->strand, bef_offset);
		//		}
		//	}
		//}
		//}
	}

	if (m_slide_window_ready[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window_ready[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window_ready[modidx].five_prim_suffix;

		//WordPair& prefix_wp_ref2 = m_slide_window_ready[modidx].prefix_wp2;

		//WordPair& five_prim_suffix_ref2 = m_slide_window_ready[modidx].five_prim_suffix2;

		//string temp = m_wp_reserve.tostring3(ALL_BITS_ON, SIXTY_FOUR);

		prefix_wp_ref = m_wp_reserve;

		five_prim_suffix_ref = m_wp_reserve;

		//prefix_wp_ref2 = m_wp;

		//five_prim_suffix_ref2 = m_wp;

		//- m_max_mismatches
		//if (m_anchor_width > 8)
		//{
		//	cout << "anchor width must be less or equal to 8";
		//	exit(0);
		//}

		//five_prim_suffix_ref.right_shift(m_anchor_width - 4);//>>m_anchor_width-2 <<2

		prefix_wp_ref.left_shift(m_seg_len - 2);

		//prefix_wp_ref2.left_shift(m_seg_len - 2);

		if (m_islands_provide)
			m_slide_window_ready[modidx].cur_tac_pos = m_cur_tac_pos;
		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}
}


inline void 
GenomeScan::EnterSlideWindowFixHmer(vector<GroupedBwts*>::iterator fh_seed_begin, vector<GroupedBwts*>::iterator fh_seed_end, bool avoid)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;

	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		m_slide_window[modidx].matchedTagsHalf.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	vector<GroupedBwts*>::iterator sv_iter;
	for (sv_iter = fh_seed_begin; sv_iter != fh_seed_end; ++sv_iter)
	{
		unsigned int firstHalfScore = 0;

		string& tagnamestr = (*sv_iter)->tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		if (avoid && m_reads_count[tagidx])
		{
			//cout << "skip read "<< tagnamestr<<endl;
			continue;
		}

		if (firstHalfScore <= m_max_mismatches)
		{
			size_t readID, bef_offset;
			readID = m_reverse_flag ? (*sv_iter)->readID_rev : (*sv_iter)->readID;
			bef_offset = m_reverse_flag ? (*sv_iter)->bef_offset_rev : (*sv_iter)->bef_offset;

			m_slide_window[modidx].matchedTagsHalf.push_back(TagSeed(firstHalfScore, (*sv_iter)->strand == '+' ? 1 : 0,  readID, sv_iter));

			size_t kmer_need = m_max_mismatches - firstHalfScore;

			if ((*sv_iter)->strand == '+')
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers[sv_iter->readID].size()*/m_num_anchors ; ++i)
				{
					if (!m_read_kmers[readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers[readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, 1/*(*sv_iter)->strand*/, bef_offset);
					}
				}
			}
			else
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
				{
					if (!m_read_kmers_rc[readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers_rc[readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, 0/*(*sv_iter)->strand*/, bef_offset);
					}
				}
			}
		}
	}

	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window[modidx].five_prim_suffix;

		prefix_wp_ref = m_wp;

		five_prim_suffix_ref = m_wp;
		//- m_max_mismatches
		if (m_anchor_width >=4 && m_anchor_width <= 8)
			five_prim_suffix_ref.right_shift(m_anchor_width - 4);//>>m_anchor_width-2 <<2
		else
		{
			cout << "anchor width must greater or equal to 4 and less or equal to 8";
			exit(0);
		}

		prefix_wp_ref.left_shift(m_masks.left_shift_seg_width);

		if (m_islands_provide)
			m_slide_window[modidx].cur_tac_pos = m_cur_tac_pos;
		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}
}


inline size_t
GenomeScan::EnterSlideWindow(const SEED_VEC& seed_vec)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;
	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		m_slide_window[modidx].matchedTagsHalf.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	if (m_slide_window[modidx].entered_reads.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		m_slide_window[modidx].entered_reads.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	SEED_VEC_CONST_ITER sv_iter;
	for (sv_iter = seed_vec.begin(); sv_iter != seed_vec.end(); ++sv_iter)
	{
		#ifdef TIME
		++m_enterhalf_count;
		#endif
		//first half
		unsigned firstHalfScore;

		if (sv_iter->strand)
		{
			//bool allbad = true;
			//for (size_t i = 0; i < m_read_kmers[sv_iter->readID].size(); ++i)
			//{
			//	if (!m_read_kmers[sv_iter->readID][i].bad)
			//	{
			//		allbad = false;
			//		break;
			//	}
			//}

			//if (allbad == true)
			//	continue;


			#ifdef TIME
			clock_t t1=clock();
			#endif

			if (m_max_mismatches)
				firstHalfScore =(unsigned) m_read_words[sv_iter->readID].score(m_wp, m_masks.first_half_mask);
			else
				firstHalfScore = 0;

			#ifdef TIME
			clock_t t2=clock();

			m_halftag_time += (t2-t1)/(double)CLOCKS_PER_SEC;
			#endif
		}
		else
		{
			//bool allbad = true;
			//for (size_t i = 0; i < m_read_kmers_rc[sv_iter->readID].size(); ++i)
			//{
			//	if (!m_read_kmers[sv_iter->readID][i].bad)
			//	{
			//		allbad = false;
			//		break;
			//	}
			//}

			//if (allbad == true)
			//	continue;

			#ifdef TIME
			clock_t t1=clock();
			#endif

			if (m_max_mismatches)
				firstHalfScore =(unsigned) m_read_words_rc[sv_iter->readID].score(m_wp, m_masks.first_half_mask);
			else
				firstHalfScore = 0;

			#ifdef TIME
			clock_t t2=clock();

			m_halftag_time += (t2-t1)/(double)CLOCKS_PER_SEC;
			#endif
		}

		if (firstHalfScore <= m_max_mismatches)
		{
			#ifdef TIME
			++m_passhalf_count;
			//second half
			clock_t t1=clock();
			#endif

			m_slide_window[modidx].matchedTagsHalf.push_back(TagSeed(firstHalfScore, sv_iter->strand, sv_iter->readID));

			#ifdef TIME
			clock_t t2=clock();
			#endif

			size_t kmer_need = m_max_mismatches - firstHalfScore;

			if (sv_iter->strand)
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers[sv_iter->readID].size()*/ m_num_anchors; ++i)
				{
					if (!m_read_kmers[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}
			else
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
				{
					if (!m_read_kmers_rc[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers_rc[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}

			#ifdef TIME
			clock_t t3=clock();
			#endif

			m_slide_window[modidx].entered_reads.insert(make_pair(sv_iter->readID, sv_iter->strand));

			#ifdef TIME
			clock_t t4=clock();

			m_matainslidewindow_time += (t2-t1 + t4- t3)/(double)CLOCKS_PER_SEC;
			m_matainopenreads_time += (t3-t2)/(double)CLOCKS_PER_SEC;
			#endif
		}
	}

	if (m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window[modidx].five_prim_suffix;

		prefix_wp_ref = m_wp;

		five_prim_suffix_ref = m_wp;
		//- m_max_mismatches
		if (m_anchor_width >=4 && m_anchor_width <= 8)
			five_prim_suffix_ref.right_shift(m_anchor_width - 4);//>>m_anchor_width-2 <<2
		else
		{
			cout << "anchor width must greater or equal to 4 and less or equal to 8";
			exit(0);
		}

		prefix_wp_ref.left_shift(m_masks.comp_left_shift_width);

		if (m_islands_provide)
			m_slide_window[modidx].cur_tac_pos = m_cur_tac_pos;
		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif

		//m_slide_window[modidx].prefix_wp = prefix_wp;

		//m_slide_window[modidx].five_prim_suffix = five_prim_suffix;
	}

	return m_slide_window[modidx].matchedTagsHalf.size();
}

inline size_t
GenomeScan::EnterSlideWindowII(const SEED_VEC& seed_vec, size_t prev_hits)
{
	size_t modidx = (m_chrom_offset + 1) % m_window_size;
	if (prev_hits == 0)
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif

		if (m_slide_window[modidx].matchedTagsHalf.size())
			m_slide_window[modidx].matchedTagsHalf.clear();

		if (m_slide_window[modidx].entered_reads.size())
			m_slide_window[modidx].entered_reads.clear();

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif
	}

	SEED_VEC_CONST_ITER sv_iter;

	for (sv_iter = seed_vec.begin(); sv_iter != seed_vec.end(); ++sv_iter)
	{
		//check if has entered in the first key check

		hash_map<unsigned int, unsigned int>& entry_ref = m_slide_window[modidx].entered_reads;
		if (prev_hits)
		{
			hash_map<unsigned int, unsigned int>::iterator hit = entry_ref.find(sv_iter->readID);
			if (hit != entry_ref.end() && hit->second == sv_iter->strand)
				continue;
		}

		#ifdef TIME
		++m_enterhalf_count;
		#endif

		//first half
		unsigned firstHalfScore;

		if (sv_iter->strand)
		{
			//if (m_read_kmers[sv_iter->readID].bad)
			//	continue;

			//bool allbad = true;
			//for (size_t i = 0; i < m_read_kmers[sv_iter->readID].size(); ++i)
			//{
			//	if (!m_read_kmers[sv_iter->readID][i].bad)
			//	{
			//		allbad = false;
			//		break;
			//	}
			//}

			//if (allbad == true)
			//	continue;

			#ifdef TIME
			clock_t t1=clock();
			#endif

			firstHalfScore =(unsigned) m_read_words[sv_iter->readID].score(m_wp, m_masks.first_half_mask);

			#ifdef TIME
			clock_t t2=clock();

			m_halftag_time += (t2-t1)/(double)CLOCKS_PER_SEC;
			#endif
		}
		else
		{
			//bool allbad = true;
			//for (size_t i = 0; i < m_read_kmers_rc[sv_iter->readID].size(); ++i)
			//{
			//	if (!m_read_kmers[sv_iter->readID][i].bad)
			//	{
			//		allbad = false;
			//		break;
			//	}
			//}

			//if (allbad == true)
			//	continue;
			//if (m_read_kmers_rc[sv_iter->readID].bad)
			//	continue;

			#ifdef TIME
			clock_t t1=clock();
			#endif

			firstHalfScore =(unsigned) m_read_words_rc[sv_iter->readID].score(m_wp, m_masks.first_half_mask);

			#ifdef TIME
			clock_t t2=clock();

			m_halftag_time += (t2-t1)/(double)CLOCKS_PER_SEC;
			#endif
		}

		if (firstHalfScore <= m_max_mismatches)
		{
			#ifdef TIME
			++m_passhalf_count;

			clock_t t1=clock();
			#endif
			//second half
			m_slide_window[modidx].matchedTagsHalf.push_back(TagSeed(firstHalfScore, sv_iter->strand, sv_iter->readID));

			#ifdef TIME
			clock_t t2=clock();
			#endif

			//(sv_iter->strand) ? m_open_reads.enter(m_read_kmers[sv_iter->readID].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand):
			//	m_open_reads.enter(m_read_kmers_rc[sv_iter->readID].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);

			size_t kmer_need = m_max_mismatches - firstHalfScore;

			if (sv_iter->strand)
			{
				for (size_t i = 0; i <= kmer_need && i </* m_read_kmers[sv_iter->readID].size()*/m_num_anchors; ++i)
				{
					if (!m_read_kmers[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}
			else
			{
				for (size_t i = 0; i <= kmer_need && i < /*m_read_kmers_rc[sv_iter->readID].size()*/m_num_anchors; ++i)
				{
					if (!m_read_kmers_rc[sv_iter->readID][i].bad)
					{
						m_open_reads[i].enter(m_read_kmers_rc[sv_iter->readID][i].kmer, m_chrom_offset + 1, m_slide_window[modidx].matchedTagsHalf.size() - 1, sv_iter->strand);
					}
				}
			}

			#ifdef TIME
			clock_t t3=clock();

			m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
			m_matainopenreads_time += (t3-t2)/(double)CLOCKS_PER_SEC;
			#endif
		}
	}

	if (prev_hits == 0 && m_slide_window[modidx].matchedTagsHalf.size())
	{
		#ifdef TIME
		clock_t t1=clock();
		#endif
		//get the sequence at matched prefix_offset of same read ID(same tag)
		WordPair& prefix_wp_ref = m_slide_window[modidx].prefix_wp;

		WordPair& five_prim_suffix_ref = m_slide_window[modidx].five_prim_suffix;

		prefix_wp_ref = m_wp;

		five_prim_suffix_ref = m_wp;
		//- m_max_mismatches

		if (m_anchor_width >=4)
			five_prim_suffix_ref.right_shift(m_anchor_width - 4);
		else
		{
			cout << "anchor width must greater or equal to 4";
			exit(0);
		}

		prefix_wp_ref.left_shift(m_masks.comp_left_shift_width);

		if (m_islands_provide)
			m_slide_window[modidx].cur_tac_pos = m_cur_tac_pos;

		//m_slide_window[modidx].prefix_island_st = m_cur_island_start;

		//m_slide_window[modidx].prefix_island_ed = m_cur_island_end;

		#ifdef TIME
		clock_t t2=clock();

		m_matainslidewindow_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		#endif

		//m_slide_window[modidx].prefix_wp = prefix_wp;

		//m_slide_window[modidx].five_prim_suffix = five_prim_suffix;
	}

	return m_slide_window[modidx].matchedTagsHalf.size();
}

void
GenomeScan::LoadChromFiles(const char* filenames_file, const char* chrom_file, const char* fasta_suffix)
{
	if (filenames_file)
		read_filename_file(filenames_file, m_chrom_files);
	else if (isdir(chrom_file))
		read_dir(chrom_file, fasta_suffix, m_chrom_files);
	else 
		m_chrom_files.push_back(chrom_file);

	//if (VERBOSE) 
	//{
	//	cerr << endl << "chromosome files found:" << endl;
	for (vector<string>::const_iterator i = m_chrom_files.begin();i != m_chrom_files.end(); ++i) 
	{
		string chromseq;

		string& chromseq_ref = m_loaded_chromos.insert(make_pair(basename2(*i), chromseq)).first->second;

		readchrom((*i).c_str(), chromseq_ref);
		
		//cerr << *i << endl;
	}
	//}
}


void
GenomeScan::SetChromFiles(const char* filenames_file, const char* chrom_file, const char* fasta_suffix, string tempchromdir)
{
	if (filenames_file)
		read_filename_file(filenames_file, m_chrom_files);
	else if (isdir(chrom_file))
		read_dir(chrom_file, fasta_suffix, m_chrom_files);
	else 
		m_chrom_files.push_back(chrom_file);

	//if (VERBOSE) 
	//{
	//	cerr << endl << "chromosome files found:" << endl;
	for (vector<string>::const_iterator i = m_chrom_files.begin();i != m_chrom_files.end(); ++i) 
	{
		m_chrom_names.push_back(basename2(*i));
		m_chrom_names_map.insert(make_pair(basename2(*i), i - m_chrom_files.begin()));
		cout <<basename2(*i)<<'\t'<<i - m_chrom_files.begin()<<endl;
		string revchrom(tempchromdir);
		revchrom.append(basename2(*i));
		revchrom.append(".rev.");
		revchrom.append(fasta_suffix);
		m_rev_chrom_files.push_back(revchrom);
		m_chrom_sizes.push_back(reverse_chrom(i->c_str(), revchrom.c_str()));
		cout << basename2(*i) <<'\t'<<m_chrom_sizes.back()<<endl;
		m_total_chromosomes_size += m_chrom_sizes.back();
		//cerr << *i << endl;
	}
	//}
}

void
GenomeScan::SetTissueFiles(const char* filenames_file, const char* tissue_file, const char* fasta_suffix)
{
	if (filenames_file)
		read_filename_file(filenames_file, m_tissue_files);
	else if (isdir(tissue_file))
		read_dir(tissue_file, fasta_suffix, m_tissue_files);
	else 
		m_tissue_files.push_back(tissue_file);

	//if (VERBOSE) 
	//{
	//	cerr << endl << "chromosome files found:" << endl;
	for (vector<string>::const_iterator i = m_tissue_files.begin();i != m_tissue_files.end(); ++i) 
	{
		m_tissue_sizes.push_back(get_filesize(*i));
		m_tissue_names.push_back(basename2(*i));
		m_total_tissue_size += m_tissue_sizes.back();
		//cerr << *i << endl;
	}
	//}
}

void
GenomeScan::GetReads(const char* reads_file)
{
	static const size_t input_buffer_size = 1000000;

	ifstream in(reads_file);
	if (!in) 
	{
		throw RMAPException("cannot open input file " + string(reads_file));
	}

	string s, name;
	bool first_line = true;
	while (!in.eof())
	{
		//if (m_reads.size() >= m_max_readID && m_read_names.size() >= m_max_readID)
		//	return;

		string buffer;
		getline(in, buffer);
		if (buffer == "")
			continue;

		if (in.gcount() == static_cast<int>(input_buffer_size))
			throw RMAPException("Line in " + name + "\nexceeds max length: " + toa(input_buffer_size));

		// correct for dos carriage returns before newlines
		if (buffer[strlen(buffer.c_str()) - 1] == '\r')
			buffer[strlen(buffer.c_str()) - 1] = '\0';

		if (buffer[0] == '>')
		{
			if (first_line == false && s.length() > 0)
			{
				//m_read_names.push_back(name);
				m_reads.push_back(s);
			}
			else first_line = false;
			name = buffer;
			name = name.substr(name.find_first_not_of(">"));
			s = "";
		}
		else s.append(buffer);
	}

	if (!first_line && s.length() > 0)
	{
		//m_read_names.push_back(name);
		m_reads.push_back(s);
	}

	if (m_reads.size() % m_num_seg != 0)
		cout << "total seg number mod num_seg != 0"<<endl;

	m_reads_count.resize(m_reads.size() / m_num_seg, 0);

	m_cur_reads_count.resize(m_reads.size() / m_num_seg, 0);
}

void
GenomeScan::GetReadsFastq(const char* reads_file)
{
	static const size_t input_buffer_size = 1000000;

	ifstream in(reads_file);
	if (!in) 
	{
		throw RMAPException("cannot open input file " + string(reads_file));
	}

	string s, name;
	bool first_line = true;
	while (!in.eof())
	{
		//if (m_reads.size() >= m_max_readID && m_read_names.size() >= m_max_readID)
		//	return;
		string buffer, buffer2, buffer3, buffer4;
		getline(in, buffer);

		getline(in, buffer2);

		getline(in, buffer3);

		getline(in, buffer4);

		if (buffer == "" || buffer2 == "" || buffer3 == "" || buffer4 == "")
			continue;

		if (in.gcount() == static_cast<int>(input_buffer_size))
			throw RMAPException("Line in " + name + "\nexceeds max length: " + toa(input_buffer_size));

		// correct for dos carriage returns before newlines
		if (buffer[strlen(buffer.c_str()) - 1] == '\r')
			buffer[strlen(buffer.c_str()) - 1] = '\0';

		if (buffer2[strlen(buffer2.c_str()) - 1] == '\r')
			buffer2[strlen(buffer2.c_str()) - 1] = '\0';

		m_reads.push_back(buffer2);
	}

	if (m_reads.size() % m_num_seg != 0)
	{
		cout << "total seg number mod num_seg != 0"<<endl;
		cout <<m_reads.size()<<endl;
		cout << m_reads.back() << endl;
		
	}

	m_reads_count.resize(m_reads.size() / m_num_seg, 0);

	m_cur_reads_count.resize(m_reads.size() / m_num_seg, 0);
}

void 
GenomeScan::GetReadsFastqQualitySeq(const char* reads_file)
{
	static const size_t input_buffer_size = 1000000;

	ifstream in(reads_file);
	if (!in) 
	{
		throw RMAPException("cannot open input file " + string(reads_file));
	}

	string s, name;
	bool first_line = true;
	while (!in.eof())
	{
		//if (m_reads.size() >= m_max_readID && m_read_names.size() >= m_max_readID)
		//	return;
		string buffer, buffer2, buffer3, buffer4;
		getline(in, buffer);

		getline(in, buffer2);

		getline(in, buffer3);

		getline(in, buffer4);

		if (buffer == "" || buffer2 == "" || buffer3 == "" || buffer4 == "")
			continue;

		if (in.gcount() == static_cast<int>(input_buffer_size))
			throw RMAPException("Line in " + name + "\nexceeds max length: " + toa(input_buffer_size));

		// correct for dos carriage returns before newlines
		if (buffer[strlen(buffer.c_str()) - 1] == '\r')
			buffer[strlen(buffer.c_str()) - 1] = '\0';

		if (buffer2[strlen(buffer2.c_str()) - 1] == '\r')
			buffer2[strlen(buffer2.c_str()) - 1] = '\0';

		if (buffer4[strlen(buffer4.c_str()) - 1] == '\r')
			buffer4[strlen(buffer4.c_str()) - 1] = '\0';

		m_fastq_seq.push_back(buffer4);
	}
}


const string&
GenomeScan::GetQualitySeq(size_t tag_idx, size_t strand_t, const vector<Differ>& differs_vec, double& qual_score)
{
	qual_score = 0;

	if (m_fastq_file == "")
		return m_qual_str;

	if (tag_idx >= m_fastq_seq.size())
	{
		cerr << "warning exceed m_fastq_seq"<<endl;
		return m_qual_str;
	}

	if (m_fastq_seq[tag_idx].length() != m_read_width)
		cout << "warning qual seq length differ "<<m_fastq_seq[tag_idx];

	vector<bool> mismach_vec(m_read_width, false);

	if ( strand_t > 0)
	{
		global_qualstr = m_fastq_seq[tag_idx];

		reverse(global_qualstr.begin(), global_qualstr.end());

		vector<Differ>::const_iterator dfv_iter;

		for (dfv_iter = differs_vec.begin(); dfv_iter != differs_vec.end(); ++dfv_iter)
		{
			mismach_vec[dfv_iter->index] = true;

			qual_score += log(CalcPSanger(global_qualstr[dfv_iter->index]) / (double(1) - double(0.25)));
		}

		for (size_t i = 0; i < mismach_vec.size(); ++i)
		{
			if (!mismach_vec[i])
			{
				qual_score += log((double(1) - CalcPSanger(global_qualstr[i])) / double(0.25));
			}
		}

		return global_qualstr;
	}
	else
	{
		vector<Differ>::const_iterator dfv_iter;

		for (dfv_iter = differs_vec.begin(); dfv_iter != differs_vec.end(); ++dfv_iter)
		{
			mismach_vec[dfv_iter->index] = true;

			qual_score += log(CalcPSanger(m_fastq_seq[tag_idx][dfv_iter->index]) / (double(1) - double(0.25)));
		}

		for (size_t i = 0; i < mismach_vec.size(); ++i)
		{
			if (!mismach_vec[i])
			{
				qual_score += log((double(1) - CalcPSanger(m_fastq_seq[tag_idx][i])) / double(0.25));
			}
		}
		return m_fastq_seq[tag_idx];
	}
}

void
GenomeScan::CalcReadsQualityScore()
{
	m_quality_scores.resize(m_reads.size() / m_num_seg, 0);
	for (size_t i = 0; i < m_reads.size(); i = i + m_num_seg)
	{
		double quality_score = 0;

		for (size_t j = 0; j < m_num_seg; ++j)
		{
			size_t idx = i + j;
			for (size_t k = 0; k < m_reads[idx].length(); ++k)
			{
				if (m_reads[idx][k] >= 'a' && m_reads[idx][k] <= 'z')
				{
					quality_score += pow(double(10), double(-1)/double(10));
				}
			}
		}

		m_quality_scores[i / m_num_seg] = quality_score;
	}
}

void
GenomeScan::ReadAndMergeTacPos(const char* tac_file)
{
	ifstream in(tac_file);
	if (!in)
	{
		throw RMAPException("cannot open input file " + string(tac_file));
	}

	string skipline;
	getline(in, skipline);
	//string chrom_name_prefix = "Equus_caballus.EquCab2.50.dna_rm.chromosome.";
	string chrom_name_prefix = "";
	while (!in.eof())
	{
		string buffer;
		getline(in, buffer);

		if (buffer == "")
			continue;

		char chrom_name[1000];
		char column2[1000];
		char column3[1000];
		int start, end;
		sscanf(buffer.c_str(), "%s\t%s\t%s\t%d\t%d", chrom_name, column2, column3, &start, &end);
		string chrom_name_str = chrom_name_prefix;
		chrom_name_str.append(chrom_name);
		
		if (chrom_name[0] != 'c')
		{
			if (!m_tac_pos_map[chrom_name_str].empty() && m_tac_pos_map[chrom_name_str].back().second >= start - 2)
				m_tac_pos_map[chrom_name_str].back().second = end - 1;
			else
				m_tac_pos_map[chrom_name_str].push_back(make_pair(start - 1, end - 1));
		}
		else
		{
			if (!m_tac_pos_map[chrom_name].empty() && m_tac_pos_map[chrom_name].back().second >= start - 2)
				m_tac_pos_map[chrom_name].back().second = end - 1;
			else
				m_tac_pos_map[chrom_name].push_back(make_pair(start - 1, end - 1));
		}
	}

	vector<string>::iterator str_iter;// m_chrom_names;

	for (size_t i = 0; i < m_chrom_names.size(); ++i)
	{
		TAC_POS_MAP_ITER TPM_iter = m_tac_pos_map.find(m_chrom_names[i]);

		size_t chrom_size = m_chrom_sizes[i];

		if (TPM_iter == m_tac_pos_map.end())
			continue;

		PAIR_SIZE_T_VEC_REV_ITER PIV_R_iter;
		for (PIV_R_iter = TPM_iter->second.rbegin(); PIV_R_iter != TPM_iter->second.rend(); ++PIV_R_iter)
		{
			m_tac_pos_map_rev[TPM_iter->first].push_back(make_pair(chrom_size - PIV_R_iter->second - 1, chrom_size - PIV_R_iter->first - 1));
		}
	}

	in.close();
}

void
GenomeScan::ReadBowtieSegments(const char* bowtie_segments_file)
{
	ifstream ifs(bowtie_segments_file);

	if( !ifs ) 
	{
		cout <<"error: open extended file error\n";exit(1);
	}

	string line1;
	while(getline(ifs, line1))
	{
		char chromname[1000], tagname[1000], strand, sequence[1000], alters[1000];
		int offset, tagid, othermapped, alignid;
		sscanf(line1.c_str(), "%s\t%d\t%c\t%s\t%d\t%s\t%d\t%s\t%d", tagname, &tagid, &strand, chromname, &offset, sequence, &othermapped, alters, &alignid);

		m_bowtie_segments.push_back(BowtieSegment(tagname, tagid, strand, chromname, offset, sequence, othermapped, alters, alignid));
	}
}

void 
GenomeScan::ReadSamFixtail(const char* samfixtailfile)
{ 
	ifstream ifs(samfixtailfile);

	if( !ifs ) 
	{
		cout <<"error: open extended file error\n";exit(1);
	}

	string line1;
	while(getline(ifs, line1))
	{
		//HWI-EAS217:4:7:1083:1580#0/1_2763	16	chr1	4596	0	75M	*	0	0	GCTGCTGCGGTGGCGGCAGAGGAGGGATGGAGTCTGACACGCGGGCAAAGGCTCCTCCGGGCCCCTCACCAGCCC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	4
		
		char chromname[1000], tagname[1000], sequence[1000], spliceway[1000], Is[1000];
		int offset, strandi, zero, missed4th;
		char star;

		sscanf(line1.c_str(), "%s\t%d\t%s\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\t%d", tagname, &strandi, chromname,  &offset, &zero, spliceway, &star, &zero, &zero, sequence, Is, &missed4th);

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		size_t fixtailoffset = offset;
		if (missed4th != 1)
		{
			int mappedlen;

			sscanf(spliceway, "%dM", &mappedlen);
			fixtailoffset = fixtailoffset + mappedlen - 1;
		}
		else
			fixtailoffset = m_chrom_sizes[m_chrom_names_map[chromname]] - fixtailoffset - 1;

		if (strandi== 0)
			tagidx = tagidx * 4 + missed4th - 1;
		else
			tagidx = tagidx * 4 + 5 - missed4th - 1;

		m_sam_fix_tail_map[missed4th][chromname].push_back(SamFixtail(tagname, strandi?0:1, offset, fixtailoffset, spliceway, tagidx));
	}
}


void 
GenomeScan::ReadFixHmer()
{
	string fixhmer_and_hole = m_bwtout;
	fixhmer_and_hole.append(".fixhmer_and_hole");
	ReadGroupedBwt(fixhmer_and_hole.c_str(), fixhmer_and_hole_type);

	string fixhole_and_hmer = m_bwtout;
	fixhole_and_hmer.append(".fixhole_and_hmer");
	ReadGroupedBwt(fixhole_and_hmer.c_str(), fixhole_and_hmer_type);

	//string mistailandhead = m_bwtout;
	//mistailandhead.append(".fusionhmer");
	//ReadGroupedBwt(mistailandhead.c_str(), fixfusion_and_hmer_type);
	
	//if (m_num_seg >= 4)
	//{
	string mistailandhead = m_bwtout;
	mistailandhead.append(".mistailandhead.notmatched");
	ReadGroupedBwt(mistailandhead.c_str(), mistailandhead_type);
	//}

	string mistailandhead_matched1side = m_bwtout;
	mistailandhead_matched1side.append(".mistailandhead.matched1side");
	ReadGroupedBwt(mistailandhead_matched1side.c_str(), mistailorhead_type);

	string mistailorheadnotmatched = m_bwtout;
	mistailorheadnotmatched.append(".mistailorhead.notmatched");
	ReadGroupedBwt(mistailorheadnotmatched.c_str(), mistailorhead_type);

	map<string, vector<GroupedBwts> >::iterator map_chrom_gbwts_iter;
	for (map_chrom_gbwts_iter = m_fixhmer.begin(); map_chrom_gbwts_iter != m_fixhmer.end(); ++map_chrom_gbwts_iter)
	{
		if (m_chrom_names_map.find(map_chrom_gbwts_iter->first) == m_chrom_names_map.end())
			continue;

		size_t chrom_size = m_chrom_sizes[m_chrom_names_map[map_chrom_gbwts_iter->first]] - 1;

		for (size_t i = 0; i != map_chrom_gbwts_iter->second.size(); ++i)
		{
			//map<string, vector<size_t> > * fixhmer_index_ptr;

			if ((map_chrom_gbwts_iter->second[i].fixtype == mistailorhead_type && map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().seg_id != 1) || 
				(map_chrom_gbwts_iter->second[i].fixtype == fixhole_and_hmer_type && map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().seg_id != 1))
			{
				m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].push_back(&(map_chrom_gbwts_iter->second[i]));
				
				m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->hmer_offset_rev = chrom_size - m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->hmer_offset_rev;

			//	cout <<"chrom size " << chrom_size<<endl;
			//	cout << "hmer_offset_rev " << m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->hmer_offset_rev << endl;
			}
			else if (map_chrom_gbwts_iter->second[i].fixtype == mistailandhead_type || map_chrom_gbwts_iter->second[i].fixtype == fixhmer_and_hole_type)
			{
				m_fixhmer_ptr[map_chrom_gbwts_iter->first].push_back(&(map_chrom_gbwts_iter->second[i]));

				m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].push_back(&(map_chrom_gbwts_iter->second[i]));

				m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->hmer_offset_rev = chrom_size - m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->hmer_offset_rev;

				if (map_chrom_gbwts_iter->second[i].fixtype == fixhmer_and_hole_type)
					m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->bef_offset_rev = chrom_size - m_fixhmer_ptr_rev[map_chrom_gbwts_iter->first].back()->bef_offset_rev;
				else
					m_fixheadandtail_ptr[map_chrom_gbwts_iter->first].push_back(&(map_chrom_gbwts_iter->second[i]));
			}
			else
				m_fixhmer_ptr[map_chrom_gbwts_iter->first].push_back(&(map_chrom_gbwts_iter->second[i]));
		}
	}

	map<string, vector<GroupedBwts*> >::iterator fixhmer_index_iter;
	for (fixhmer_index_iter = m_fixhmer_ptr.begin(); fixhmer_index_iter != m_fixhmer_ptr.end(); ++fixhmer_index_iter)
	{
		sort(fixhmer_index_iter->second.begin(), fixhmer_index_iter->second.end(), compFixHmer);
	}

	for (fixhmer_index_iter = m_fixhmer_ptr_rev.begin(); fixhmer_index_iter != m_fixhmer_ptr_rev.end(); ++fixhmer_index_iter)
	{
		sort(fixhmer_index_iter->second.begin(), fixhmer_index_iter->second.end(), compFixHmer_rev);
	}
}

void
GenomeScan::TryAppend(const char* infile, const FixTypes& fixhmertype)
{
	ifstream input_fs(infile);

	if( !input_fs ) 
	{
		cout <<"error: can't open file error: "<< infile<<endl;return;//exit(1);
	}

	map<string, vector<GroupedBwts> > fixhmer;

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;
		getline(input_fs, line);
		if (line == "")
			continue;

		sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		vector<GroupedBwts> & fixhmer_chrom_ref = fixhmer[chromname];
		fixhmer_chrom_ref.push_back(GroupedBwts(tagname, chromname, strand, tagidx, fixhmertype));
		fixhmer_chrom_ref.back().grouped_bwt_recs.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));

		for (size_t i = 1; i < rec_num; ++i)
		{
			getline(input_fs, line);

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);
			fixhmer_chrom_ref.back().grouped_bwt_recs.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		}

		if (fixhmer_chrom_ref.back().fixtype == mistailandhead_type)
		{
			if (fixhmer_chrom_ref.back().strand == '+')
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg;
			}
			else
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + m_num_seg - 1;
			}
		}
		else if (fixhmer_chrom_ref.back().fixtype == mistailorhead_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id != 1)
		{
			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg;
			else
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + m_num_seg - 1;
		}
		else if (fixhmer_chrom_ref.back().fixtype == mistailorhead_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id == 1)
		{
			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
			else
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg;
		}
		else 
		{
			//missed cases
			cout << "missed case : "<<endl;
			cout << "fixtype : "<<fixhmer_chrom_ref.back().fixtype<<endl;
		}
	}

	input_fs.close();

	Classify_Append(fixhmer, infile, fixhmertype);
}

void 
GenomeScan::Classify_Append(map<string, vector<GroupedBwts> >& fixhmer, const char* infile, const FixTypes& fixhmertype)
{
	string matched = infile;
	matched.append(".matched");
	ofstream matched_ofs(matched.c_str());

	string notmatched = infile;
	notmatched.append(".notmatched");
	ofstream notmatched_ofs(notmatched.c_str());

	string matched1side = infile;
	matched1side.append(".matched1side");
	ofstream matched1side_ofs;

	if (fixhmertype == mistailandhead_type)
	{
		matched1side_ofs.open(matched1side.c_str());
	}

	size_t min_score = m_append_score;

	size_t chrom_id;
	map<string, vector<GroupedBwts> >::iterator map_chrom_gbwts_iter;
	for (map_chrom_gbwts_iter = fixhmer.begin(); map_chrom_gbwts_iter != fixhmer.end(); ++map_chrom_gbwts_iter)
	{
		cout << map_chrom_gbwts_iter->first << endl;
		string chromseq;
		if (m_chrom_names_map.find(map_chrom_gbwts_iter->first) == m_chrom_names_map.end())
		{
			cout <<"can't find chrom " << map_chrom_gbwts_iter->first<<endl;
			continue;
		}
		else
		{
			chrom_id = m_chrom_names_map[map_chrom_gbwts_iter->first];

			string chrom_file = m_chrom_files[chrom_id];

			readchrom(chrom_file.c_str(), chromseq);
		}
		
		for (size_t i = 0; i != map_chrom_gbwts_iter->second.size(); ++i)
		{
			if (map_chrom_gbwts_iter->second[i].fixtype == mistailorhead_type && map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().seg_id != 1)
			{
				string first_read_seg = m_reads[map_chrom_gbwts_iter->second[i].readID_rev];

				if ( map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset < m_seg_len || map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset > chromseq.length())
				{
					cout << "offset exceed 1 "<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len <<endl 
					<< map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset<<endl;
					continue;
				}

				string first_chrom_seq;
				try
				{
					//cout << map_chrom_gbwts_iter->first << " offset 1 " << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len<<endl;
					//cout << chromseq.length() << endl;
					first_chrom_seq = chromseq.substr(map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len, m_seg_len);
				}
				catch(const RMAPException &e)
				{
					cout << "offset exceed 1 catched"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl 
					<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len<<endl;

					continue;
					getchar();
				}
				if (map_chrom_gbwts_iter->second[i].strand == '-')
					first_read_seg = revcomp(first_read_seg);

				pair<size_t, size_t> score_first_bits = score_string(first_read_seg, first_chrom_seq);

				if (score_first_bits.first > min_score)
				{
					notmatched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size()<<endl;
					for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					{
						notmatched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					}
				}
				else if (score_first_bits.first <= min_score)
				{
					size_t score_sum = score_first_bits.first;

					GroupedBwts& cur_bwts_ref = map_chrom_gbwts_iter->second[i];

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						score_sum += cur_bwts_ref.grouped_bwt_recs[j].differs.length() / 5;
					}

					//if (score_sum > m_max_mismatches)
					//	continue;

					size_t lste_st = cur_bwts_ref.grouped_bwt_recs.front().offset - m_seg_len + 1;

					size_t strand_t = 0;
					if (cur_bwts_ref.strand == '-')
						strand_t = 16;

					vector<Differ> differs_vec;

					size_t seg_id = 0;

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] >= 'a' && first_read_seg[k] <= 'z')
							first_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_chrom_seq.length(); ++k)
					{
						if (first_chrom_seq[k] >= 'a' && first_chrom_seq[k] <= 'z')
							first_chrom_seq[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] != first_chrom_seq[k])
							differs_vec.push_back(Differ(seg_id * m_seg_len + k, first_chrom_seq[k], first_read_seg[k]));
					}

					++seg_id;

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						for (size_t k = 0; k < cur_bwts_ref.grouped_bwt_recs[j].differs_vec.size(); ++k)
						{
							differs_vec.push_back(Differ(seg_id * m_seg_len + cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].index, 
								cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].src_char, cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].dst_char));
						}

						++seg_id;
					}

					double qual_score;

					const string& qual_seq = GetQualitySeq(cur_bwts_ref.tagidx, strand_t, differs_vec, qual_score);

					matched_ofs<<cur_bwts_ref.tagname<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_id]<<'\t'<<lste_st <<'\t'<<size_t(qual_score) <<'\t'<<m_read_width << 'M'<<"\t*\t0\t0\t"<<first_read_seg;

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						matched_ofs<<cur_bwts_ref.grouped_bwt_recs[j].mapped_chrom_seq;
					}

					matched_ofs<<'\t'<<qual_seq<<'\t'<< "NM:i:"<<score_sum<<'\t';

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						matched_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							matched_ofs << ',';
					}

					matched_ofs << endl; 

					//matched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size()<<endl;
					//for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					//{
					//	matched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					//}

					//matched_ofs << first_read_seg <<'\t'<< first_chrom_seq << '\t'<<score_first_bits.first<< endl;
				}
			}
			else if (map_chrom_gbwts_iter->second[i].fixtype == mistailorhead_type && map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().seg_id == 1)
			{
				string last_read_seg = m_reads[map_chrom_gbwts_iter->second[i].readID];

				if (map_chrom_gbwts_iter->second[i].strand == '-')
					last_read_seg = revcomp(last_read_seg);

				if ( map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() < 0 || 
					map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len > chromseq.length())
				{
					cout << "offset exceed 2"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl 
					<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len<<endl;

					continue;
				}

				string last_chrom_seq;

				try
				{
					//cout <<map_chrom_gbwts_iter->first << " offset 2 " << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length()<<endl;
					//cout << chromseq.length() << endl;
				
					last_chrom_seq = chromseq.substr(map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length(), m_seg_len);

				}
				catch(const RMAPException &e)
				{
					cout << "offset exceed 2 catched"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl 
					<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len<<endl;

					continue;
					getchar();
				}
				pair<size_t, size_t> score_last_bits = score_string(last_read_seg, last_chrom_seq);

				if (score_last_bits.first > min_score)
				{
					notmatched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size()<<endl;
					for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					{
						notmatched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					}
				}
				else if (score_last_bits.first <= min_score)
				{
					size_t score_sum = score_last_bits.first;

					GroupedBwts& cur_bwts_ref = map_chrom_gbwts_iter->second[i];

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						score_sum += cur_bwts_ref.grouped_bwt_recs[j].differs.length() / 5;
					}

					//if (score_sum > m_max_mismatches)
					//	continue;
				
					size_t lste_st = cur_bwts_ref.grouped_bwt_recs.front().offset + 1;

					size_t strand_t = 0;
					if (cur_bwts_ref.strand == '-')
						strand_t = 16;

					vector<Differ> differs_vec;

					size_t seg_id = 0;

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						for (size_t k = 0; k < cur_bwts_ref.grouped_bwt_recs[j].differs_vec.size(); ++k)
						{
							differs_vec.push_back(Differ(seg_id * m_seg_len + cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].index, 
								cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].src_char, cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].dst_char));
						}

						++seg_id;
					}

					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] >= 'a' && last_read_seg[k] <= 'z')
							last_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < last_chrom_seq.length(); ++k)
					{
						if (last_chrom_seq[k] >= 'a' && last_chrom_seq[k] <= 'z')
							last_chrom_seq[k] += 'A' -'a';
					}

					//last_read_seg, last_chrom_seq);
					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] != last_chrom_seq[k])
							differs_vec.push_back(Differ(seg_id * m_seg_len + k, last_chrom_seq[k] , last_read_seg[k]));
					}

					++seg_id;

					double qual_score;

					const string& qual_seq = GetQualitySeq(cur_bwts_ref.tagidx, strand_t, differs_vec, qual_score);

					matched_ofs<<cur_bwts_ref.tagname<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_id]<<'\t'<<lste_st <<'\t'<<size_t(qual_score) <<'\t'<<m_read_width << 'M'<<"\t*\t0\t0\t";

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						matched_ofs<<cur_bwts_ref.grouped_bwt_recs[j].mapped_chrom_seq;
					}

					matched_ofs<<last_read_seg;
					matched_ofs<<'\t'<< qual_seq<<'\t'<< "NM:i:"<<score_sum<<'\t';
					

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						matched_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							matched_ofs << ',';
					}

					matched_ofs << endl; 

					//matched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size()<<endl;
					//for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					//{
					//	matched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					//}

					//matched_ofs << last_read_seg <<'\t'<< last_chrom_seq << '\t'<<score_last_bits.first<< endl;
				}
			}
			else if (map_chrom_gbwts_iter->second[i].fixtype == mistailandhead_type)
			{
				string first_read_seg = m_reads[map_chrom_gbwts_iter->second[i].readID_rev];
				if (map_chrom_gbwts_iter->second[i].strand == '-')
					first_read_seg = revcomp(first_read_seg);

				if ( map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset < m_seg_len || 
					map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset > chromseq.length())
				{
					cout << "offset exceed 3"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len <<endl 
					<< map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset<<endl;
					continue;
				}

				string first_chrom_seq;

				try
				{
					//cout << map_chrom_gbwts_iter->first << " offset 3 " << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len<<endl;
					//cout << chromseq.length() << endl;

					first_chrom_seq = chromseq.substr(map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len, m_seg_len);
				}
				catch(const RMAPException &e)
				{
					cout << "offset exceed 3 catched"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset - m_seg_len <<endl 
					<< map_chrom_gbwts_iter->second[i].grouped_bwt_recs.front().offset<<endl;
					continue;
					getchar();
				}
				pair<size_t, size_t> score_first_bits = score_string(first_read_seg, first_chrom_seq);

				string last_read_seg = m_reads[map_chrom_gbwts_iter->second[i].readID];
				if (map_chrom_gbwts_iter->second[i].strand == '-')
					last_read_seg = revcomp(last_read_seg);

				if ( map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() < 0 || 
					map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len > chromseq.length())
				{
					cout << "offset exceed 4"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl 
					<< map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len<<endl;
					continue;
				}

				string last_chrom_seq;
				try 
				{
					//cout <<map_chrom_gbwts_iter->first << " offset 4 " << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl;
					//cout << chromseq.length() << endl;
					last_chrom_seq = chromseq.substr(map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length(), m_seg_len);
				}
				catch(const RMAPException &e)
				{
					cout << "offset exceed 4 catched"<<endl<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() <<endl 
					<< map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().mapped_chrom_seq.length() + m_seg_len<<endl;
					continue;
					getchar();
				}

				pair<size_t, size_t> score_last_bits = score_string(last_read_seg, last_chrom_seq);

				if ((score_first_bits.first > min_score && score_last_bits.first > m_max_mismatches) ||
					(score_first_bits.first > m_max_mismatches && score_last_bits.first > min_score))
				{
					notmatched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size()<<endl;
					for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					{
						notmatched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					}
				}
				else if (score_first_bits.first <= min_score && score_last_bits.first <= min_score)
				{
					size_t score_sum = score_first_bits.first + score_last_bits.first;

					GroupedBwts& cur_bwts_ref = map_chrom_gbwts_iter->second[i];

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						score_sum += cur_bwts_ref.grouped_bwt_recs[j].differs.length() / 5;
					}

					//if (score_sum > m_max_mismatches)
					//	continue;

					size_t lste_st = cur_bwts_ref.grouped_bwt_recs.front().offset - m_seg_len + 1;

					size_t strand_t = 0;
					if (cur_bwts_ref.strand == '-')
						strand_t = 16;

					vector<Differ> differs_vec;

					size_t seg_id = 0;

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] >= 'a' && first_read_seg[k] <= 'z')
							first_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_chrom_seq.length(); ++k)
					{
						if (first_chrom_seq[k] >= 'a' && first_chrom_seq[k] <= 'z')
							first_chrom_seq[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] != first_chrom_seq[k])
							differs_vec.push_back(Differ(seg_id * m_seg_len + k, first_chrom_seq[k], first_read_seg[k]));
					}

					++seg_id;

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						for (size_t k = 0; k < cur_bwts_ref.grouped_bwt_recs[j].differs_vec.size(); ++k)
						{
							differs_vec.push_back(Differ(seg_id * m_seg_len + cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].index, 
								cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].src_char, cur_bwts_ref.grouped_bwt_recs[j].differs_vec[k].dst_char));
						}

						++seg_id;
					}

					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] >= 'a' && last_read_seg[k] <= 'z')
							last_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < last_chrom_seq.length(); ++k)
					{
						if (last_chrom_seq[k] >= 'a' && last_chrom_seq[k] <= 'z')
							last_chrom_seq[k] += 'A' -'a';
					}

					//last_read_seg, last_chrom_seq);
					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] != last_chrom_seq[k])
							differs_vec.push_back(Differ(seg_id * m_seg_len + k, last_chrom_seq[k], last_read_seg[k]));
					}

					++seg_id;

					double qual_score;

					const string& qual_seq = GetQualitySeq(cur_bwts_ref.tagidx, strand_t, differs_vec, qual_score);

					matched_ofs<<cur_bwts_ref.tagname<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_id]<<'\t'<<lste_st <<'\t'<<size_t(qual_score) <<'\t'<<m_read_width << 'M'<<"\t*\t0\t0\t"<<first_read_seg;

					for (size_t j = 0; j < cur_bwts_ref.grouped_bwt_recs.size(); ++j)
					{
						matched_ofs<<cur_bwts_ref.grouped_bwt_recs[j].mapped_chrom_seq;
					}

					matched_ofs<<last_read_seg;

					matched_ofs<<'\t'<<qual_seq<<'\t'<< "NM:i:"<<score_sum<<'\t';

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						matched_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							matched_ofs << ',';
					}

					matched_ofs << endl; 

					//matched_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size() <<endl;
					//for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					//{
					//	matched_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
					//		<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					//}

					//matched_ofs << first_read_seg <<'\t'<< first_chrom_seq << '\t'<<score_first_bits.first<< endl;
					//matched_ofs << last_read_seg <<'\t'<< last_chrom_seq << '\t'<<score_last_bits.first<< endl;
				}
				else if (score_first_bits.first <= m_max_mismatches)
				{
					matched1side_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size() + 1<<endl;

					string Ns(score_first_bits.first * 5 + 1, 'N');
					matched1side_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[0].offset - m_seg_len<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[0].seg_id - 1<<'\t'
							<< first_read_seg<<'\t'<<0<<'\t';//<<Ns<<endl;

					vector<Differ> differs_vec;

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] >= 'a' && first_read_seg[k] <= 'z')
							first_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_chrom_seq.length(); ++k)
					{
						if (first_chrom_seq[k] >= 'a' && first_chrom_seq[k] <= 'z')
							first_chrom_seq[k] += 'A' -'a';
					}

					for (size_t k = 0; k < first_read_seg.length(); ++k)
					{
						if (first_read_seg[k] != first_chrom_seq[k])
							differs_vec.push_back(Differ(k, first_chrom_seq[k], first_read_seg[k]));
					}

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						matched1side_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							matched1side_ofs << ',';
					}

					matched1side_ofs << endl; 

					for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					{
						matched1side_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					}
				}
				else if (score_last_bits.first <= m_max_mismatches)
				{
					matched1side_ofs << map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size() + 1<<endl;

					for (size_t j = 0; j < map_chrom_gbwts_iter->second[i].grouped_bwt_recs.size(); ++j)
					{
						matched1side_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].offset<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].seg_id<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].mapped_chrom_seq<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].other_mapped<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs[j].differs<<endl;
					}

					//string Ns(score_last_bits.first * 5 + 1, 'N');
					matched1side_ofs <<  map_chrom_gbwts_iter->second[i].tagname <<'\t'<<map_chrom_gbwts_iter->second[i].chromname<<'\t' << map_chrom_gbwts_iter->second[i].strand<<'\t'
							<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().offset + m_seg_len<<'\t'<<map_chrom_gbwts_iter->second[i].grouped_bwt_recs.back().seg_id + 1<<'\t'
							<< last_read_seg<<'\t'<<0<<'\t';//<<Ns<<endl;

					vector<Differ> differs_vec;

					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] >= 'a' && last_read_seg[k] <= 'z')
							last_read_seg[k] += 'A' -'a';
					}

					for (size_t k = 0; k < last_chrom_seq.length(); ++k)
					{
						if (last_chrom_seq[k] >= 'a' && last_chrom_seq[k] <= 'z')
							last_chrom_seq[k] += 'A' -'a';
					}

					//last_read_seg, last_chrom_seq);
					for (size_t k = 0; k < last_read_seg.length(); ++k)
					{
						if (last_read_seg[k] != last_chrom_seq[k])
							differs_vec.push_back(Differ(k,  last_chrom_seq[k], last_read_seg[k]));
					}

					vector<Differ>::iterator vD_iter;

					for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
					{
						matched1side_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_vec.end() - 1)
							matched1side_ofs << ',';
					}

					matched1side_ofs << endl; 
				}
				else
				{
					cout << "can't arrive here"<<endl;
				}
			}
			else
			{
				//missed cases
				cout << "missed case : "<<endl;
				cout << "fixtype : "<<map_chrom_gbwts_iter->second[i].fixtype<<endl;
			}
		}
	}
	notmatched_ofs.close();
	matched_ofs.close();
	matched1side_ofs.close();
}

void 
GenomeScan::ReadGroupedBwt(const char* infile, const FixTypes& fixhmertype)
{
	ifstream input_fs(infile);

	if( !input_fs ) 
	{
		cout <<"error: can't open file error: "<< infile<<endl;return;//exit(1);
	}

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;
		getline(input_fs, line);
		if (line == "")
			continue;

		sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);

		//map<string, vector<GroupedBwts> > * fixhmer_ptr;
		//if (fixhmertype == mistailorhead_type && seg_id != 1)
		//	fixhmer_ptr = &m_fixhmer_rev;
		//else
		//	fixhmer_ptr = &m_fixhmer;

		//(*fixhmer_ptr)[chromname].push_back(GroupedBwts(tagname, chromname, strand, fixhmertype));
		//(*fixhmer_ptr)[chromname].back().grouped_bwt_recs.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		vector<GroupedBwts> & fixhmer_chrom_ref = m_fixhmer[chromname];
		fixhmer_chrom_ref.push_back(GroupedBwts(tagname, chromname, strand, tagidx, fixhmertype));
		fixhmer_chrom_ref.back().grouped_bwt_recs.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));

		for (size_t i = 1; i < rec_num; ++i)
		{
			getline(input_fs, line);

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);
			fixhmer_chrom_ref.back().grouped_bwt_recs.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		}

		if (fixhmer_chrom_ref.back().fixtype == mistailandhead_type)
		{
			fixhmer_chrom_ref.back().prefix_side = new vector<SideMatch>;
			fixhmer_chrom_ref.back().suffix_side = new vector<SideMatch>;
			fixhmer_chrom_ref.back().hmer_offset = fixhmer_chrom_ref.back().grouped_bwt_recs.back().offset + m_seg_len;
			fixhmer_chrom_ref.back().hmer_offset_rev = fixhmer_chrom_ref.back().grouped_bwt_recs.front().offset - 1;

			if (fixhmer_chrom_ref.back().strand == '+')
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg;
			}
			else
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + m_num_seg - 1;
			}
		}
		else if (fixhmer_chrom_ref.back().fixtype == fixhmer_and_hole_type)
		{
			//possiblly need to be more generalized when num_seg > 4
			fixhmer_chrom_ref.back().hmer_offset = fixhmer_chrom_ref.back().grouped_bwt_recs.front().offset + m_seg_len; 
			fixhmer_chrom_ref.back().hmer_offset_rev = fixhmer_chrom_ref.back().grouped_bwt_recs.back().offset - 1;

			fixhmer_chrom_ref.back().bef_offset = fixhmer_chrom_ref.back().hmer_offset_rev - m_interal;
			fixhmer_chrom_ref.back().bef_offset_rev = fixhmer_chrom_ref.back().hmer_offset + m_interal;
			if (fixhmer_chrom_ref.back().strand == '+')
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + 1;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + m_num_seg - 2;
			}
			else
			{
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 2;
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + 1;
			}
		}
		else if ((fixhmer_chrom_ref.back().fixtype == mistailorhead_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id != 1) || 
				 (fixhmer_chrom_ref.back().fixtype == fixhole_and_hmer_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id != 1))
		{
			fixhmer_chrom_ref.back().hmer_offset_rev = m_fixhmer[chromname].back().grouped_bwt_recs.front().offset - 1;

			//cout << "hmer_offset_rev " << fixhmer_chrom_ref.back().hmer_offset_rev << endl;

			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg;
			else
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg + m_num_seg - 1;
		}
		else if ((fixhmer_chrom_ref.back().fixtype == mistailorhead_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id == 1) || 
				 (fixhmer_chrom_ref.back().fixtype == fixhole_and_hmer_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id == 1))
		{
			fixhmer_chrom_ref.back().hmer_offset = fixhmer_chrom_ref.back().grouped_bwt_recs.back().offset + m_seg_len;

			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
			else
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg;
		}
		else if (fixhmer_chrom_ref.back().fixtype == fixfusion_and_hmer_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id != 1)
		{
			fixhmer_chrom_ref.back().hmer_offset_rev = m_fixhmer[chromname].back().grouped_bwt_recs.front().offset - 1;

			//cout << "hmer_offset_rev " << fixhmer_chrom_ref.back().hmer_offset_rev << endl;

			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID_rev = tagidx * m_num_seg;
			else
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
		}
		else if (fixhmer_chrom_ref.back().fixtype == fixfusion_and_hmer_type && fixhmer_chrom_ref.back().grouped_bwt_recs.front().seg_id == 1)
		{
			fixhmer_chrom_ref.back().hmer_offset = fixhmer_chrom_ref.back().grouped_bwt_recs.back().offset + m_seg_len;

			if (fixhmer_chrom_ref.back().strand == '+')
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg + m_num_seg - 1;
			else
				fixhmer_chrom_ref.back().readID = tagidx * m_num_seg;
		}
		else 
		{
			//missed cases
			cout << "missed case : "<<endl;
			cout << "fixtype : "<<fixhmer_chrom_ref.back().fixtype<<endl;
		}
	}

	input_fs.close();
}

void 
GenomeScan::SortFixtail()
{
	map<size_t, map<string, vector<SamFixtail> > >::iterator missed_id_iter;

	for (missed_id_iter = m_sam_fix_tail_map.begin(); missed_id_iter != m_sam_fix_tail_map.end(); ++missed_id_iter)
	{
		map<string, vector<SamFixtail> >::iterator chrom_iter;

		for (chrom_iter = missed_id_iter->second.begin(); chrom_iter != missed_id_iter->second.end(); ++chrom_iter)
		{
			sort(chrom_iter->second.begin(), chrom_iter->second.end(), compSamFixtail);
		}
	}
}

void
GenomeScan::CleanReads()
{
	size_t good = 0;
	for (size_t i = 0; i < m_reads.size(); ++i)
	{

		if (m_reads[i].length() != m_read_width)
		{
			if (m_reads[i].length() < m_read_width)
			{
				while (m_reads[i].length() < m_read_width)
					m_reads[i].append("N");
				continue;
				cout << i<<'\t'<<m_reads[i].length()<<'\t'<<m_read_width<<endl;
				cout << m_reads[i]<<endl;
				throw RMAPException("Incorrect read width");
			}
			else
				m_reads[i] = m_reads[i].substr(0, m_read_width);
		}

		if (count(m_reads[i].begin(), m_reads[i].end(), 'N') <= static_cast<int>(m_max_mismatches))
		{
			m_reads[good].swap(m_reads[i]);
			m_read_names[good].swap(m_read_names[i]);
			++good;
		}
	}
	m_reads.resize(good);
	m_read_names.resize(good);
}

void
GenomeScan::FixTail(size_t chromID, const vector<SamFixtail>& fixtail_vec)
{
	string chrom_file;

	TAC_POS_MAP* tac_pos_ptr;

	if (m_reverse_flag)
	{
		chrom_file = m_rev_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map_rev;
	}
	else
	{
		chrom_file = m_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map;
	}

	ChromReader chrom_reader(chrom_file, m_input_buffer_size);

	TAC_POS_MAP_ITER chrom_pos_map = tac_pos_ptr->find(m_chrom_names[chromID]);

	m_cur_island_start = m_boundary + 1;
	
	m_cur_island_end = max_size_t - (3 * m_boundary);

	m_islands = "N";

	m_islands_provide = false;

	if (chrom_pos_map != tac_pos_ptr->end())
	{
		m_islands = "Y";
		m_islands_provide = true;
		m_cur_tac_pos = chrom_pos_map->second.begin();
		m_tac_pos_begin = chrom_pos_map->second.begin();
		m_tac_pos_end = chrom_pos_map->second.end();
		m_cur_island_start = m_cur_tac_pos->first;
		m_cur_island_end = m_cur_tac_pos->second;
	}

	cout << "islands: "<<m_islands<<endl;

	m_chrom_offset = -1;

	m_bad_bases = ALL_BITS_ON;
	m_bad_bases_reserve = ALL_BITS_ON;

	m_seed_key = 0;
	m_seed_keyII = 0;

	for (size_t i = 0; i < m_kmer_keys.size(); ++i)
		m_kmer_keys[i] = 0;

	vector<SamFixtail>::const_iterator cur_fixtail_iter = fixtail_vec.begin();

	while (chrom_reader.read())
	{
		size_t base;
		while ((base = chrom_reader.ScanBase()) != -1)
		{
			++m_chrom_offset;

			GetSeed(base);

			m_bad_bases_reserve = ((m_bad_bases_reserve << 1) + static_cast<size_t>(m_bad_bases & MOST_SIG_BIT));

			m_bad_bases = (m_bad_bases << 1) + static_cast<size_t>(base == 4);

			//skip m_chrom_offset not in exons
			while (m_chrom_offset + m_anchor_width > m_cur_island_end + m_boundary + m_seg_len)
			{
				++m_cur_tac_pos;
				if (m_cur_tac_pos == m_tac_pos_end)
					return;

				m_cur_island_start = m_cur_tac_pos->first;
				m_cur_island_end = m_cur_tac_pos->second;

				if (m_chrom_offset + m_anchor_width < m_cur_island_start + m_min_boundary + m_seg_len || 
					(m_chrom_offset + m_anchor_width >= m_cur_island_start + m_min_boundary + m_seg_len && m_chrom_offset + m_anchor_width <= m_cur_island_end + m_boundary + m_seg_len))
					break;
			}

			if (cur_fixtail_iter != fixtail_vec.end() && m_chrom_offset == cur_fixtail_iter->fixtailoffset + m_read_width)
			{
				vector<SamFixtail>::const_iterator first_fixtail_iter = cur_fixtail_iter++;
				while (cur_fixtail_iter != fixtail_vec.end() && m_chrom_offset == cur_fixtail_iter->fixtailoffset + m_read_width)
				{
					//insert first half into slide window
					++cur_fixtail_iter;
				}

				EnterSlideWindowFixTail(first_fixtail_iter, cur_fixtail_iter);
			}


			//check hmers
			if (m_chrom_offset < m_cur_island_start + m_min_boundary)
				continue;

			++m_basecount;

			m_wp_reserve.bads = m_bad_bases_reserve;
			m_wp.bads = m_bad_bases;

			//check suffix kmer
			hash_map<size_t, int> mapped_TagSeed;
			for (size_t i = 0; i < m_num_anchors; ++i)
			{
				if (!(m_masks.bad_kmer_masks[i] & m_bad_bases))
				{
#ifdef TIME
					clock_t t3=clock();
#endif

					SEEDQUEUE& prefix_matched_kmers = m_open_reads[i].find_prefix_matches(m_kmer_keys[i], m_chrom_offset);

#ifdef TIME
					clock_t t4=clock();

					m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

					if (!prefix_matched_kmers.Size())
						continue;

					size_t lit = prefix_matched_kmers.Rear();

					size_t sz = prefix_matched_kmers.Size();

					size_t k;

					for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist - 2 >= m_chrom_offset; ++k)
					{
						//Trigger Match
						//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;
#ifdef TIME
						clock_t t1=clock();
#endif

						TriggerMatchFixTail(prefix_matched_kmers.At(lit), chromID, i, mapped_TagSeed);

#ifdef TIME
						clock_t t2=clock();
#endif

						prefix_matched_kmers.Pred(lit);

#ifdef TIME
						clock_t t5=clock();

						m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

						m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
					}

#ifdef TIME
					clock_t t6=clock();
#endif

					if (k != sz)
						prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
					clock_t t7=clock();

					m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif

					//
				}
			}

		}
	}

}

void
GenomeScan::FixHmerReady(const char* fixhmer_file, const char* fixed_hmer_file, size_t flankcase, bool avoid, bool island)
{
	string sort_cmd = "sort -k4,4n -S 3500000 -o ";

	sort_cmd.append(fixhmer_file);

	sort_cmd.append(" -T temp ");

	sort_cmd.append(fixhmer_file);

	system(sort_cmd.c_str());

	string fixhmer_file_str = fixhmer_file;

	string chromname = basename2(fixhmer_file_str);

	string chromnamestr = chromname;

	cout << chromname << endl;

	if (m_chrom_names_map.find(chromname) == m_chrom_names_map.end())
	{
		cout << "can't find chromosome: " << chromname << endl;
		return;
	}

	string chrom_file;

	m_cur_chrom_size = m_chrom_sizes[m_chrom_names_map[chromname]];

	TAC_POS_MAP* tac_pos_ptr;
	
	//if (m_reverse_flag)
	//{
	//	chrom_file = m_rev_chrom_files[m_chrom_names_map[chromname]];

	//	tac_pos_ptr = &m_tac_pos_map_rev;
	//}
	//else
	//{
	chrom_file = m_chrom_files[m_chrom_names_map[chromname]];

	tac_pos_ptr = &m_tac_pos_map;
	//}

	//load chromo sequence into memory
	m_chrom_seq.clear();
	readchrom(chrom_file.c_str(), m_chrom_seq);

	m_cur_island_start = /*m_boundary + */1;

	m_cur_island_end = max_size_t - (3 * m_boundary);

	m_islands = "N";

	m_islands_provide = false;

	TAC_POS_MAP_ITER chrom_pos_map = tac_pos_ptr->find(chromname);

	if (chrom_pos_map != tac_pos_ptr->end() && island)
	{
		m_islands = "Y";
		m_islands_provide = true;
		m_cur_tac_pos = chrom_pos_map->second.begin();
		m_tac_pos_begin = chrom_pos_map->second.begin();
		m_tac_pos_end = chrom_pos_map->second.end();
		m_cur_island_start = m_cur_tac_pos->first;
		m_cur_island_end = m_cur_tac_pos->second;
	}

	cout << "islands: "<<m_islands<<endl;

	//if (!m_is_fix_tail_head)
	//	cout << "false"<<endl;
	//else
	//	cout << "true"<<endl;

	ChromReader chrom_reader(chrom_file, m_input_buffer_size);

	ifstream input_fs(fixhmer_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhmer_file<<"\n";return;//exit(1);
	}

	m_ofs_fix_hmer.close();
	m_ofs_fix_hmer.open(fixed_hmer_file);

	string fixed_hmer_file_hole = fixed_hmer_file; fixed_hmer_file_hole.append(".hole");
	m_ofs_fix_hmer_hole.close();
	m_ofs_fix_hmer_hole.open(fixed_hmer_file_hole.c_str());

	//if (!m_is_fix_tail_head)
	//{
	//	cout << "file opened "<<endl; 
	//}

	//m_ofs_fix_hmer << "test"<<endl;

	//m_ofs_fix_hmer.close();


	size_t double_seg_len = 2 * m_seg_len;

	vector<TobeFixedHmerReady> cur_tobe_fixed_hmer_segments;

	//m_cur_junc_head_iter, m_cur_junc_tail_iter, m_cur_junc_end_iter;
	if (m_is_fix_tail_head)
	{
		m_cur_junc_head_iter = m_cur_junc_vec_ptr->begin(); 

		m_cur_junc_tail_iter = m_cur_junc_vec_ptr->begin();

		m_cur_junc_end_iter = m_cur_junc_vec_ptr->end();
	}

	//if (!m_is_fix_tail_head)
	//{
	//	cout << "start "<<endl; 
	//}

	//size_t count = 0;
	
	while(!input_fs.eof())
	{
		string line;

		getline(input_fs, line);

		if (line.length() == 0)
			continue;

		//if (!m_is_fix_tail_head)
		//{
		//	cout << line<< endl;
		//	cout << count++<<endl;
		//}

		char chromnamechr[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000], */strand;

		size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id/*, other_mapped*/;
		
		int id1, id2;

		//HWI-EAS217:3:100:1000:1066#0/2	chr15	+	88784825	0	1-1	0-0	CAAGATCGGGGGCATCTTGGCCAAT	-	GAGCTCTCGGTGGCTGAGGCTGNAG,TCCATGCAGCTGTTCTTGCCATCAA,TGAAGCAGTGGAGCGAGGGGTGGTG,

		//40488~TRAN00000035434:439:2818  chr20   +       1168044 1169006 1-1     4-4     CCTGCCTTCTGAGCCTCCCAATAGT       AAAAATGAATGAAACAATATTATTA      TTAGCAGTTGAACCAGATAACTCAG,AGTATATATGTGTAACTGAAAATGA,    0       1

	/*	if (!m_is_fix_tail_head)
			cout << "sscanf "<<endl;*/

		sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%d\t%d", tagname, chromnamechr, &strand, &doner_end, &acceptor_st, 
			&doner_st_id, &doner_end_id, &acceptor_st_id, &acceptor_end_id,	doner_seq, acceptor_seq, tobe_fixed_seq, &id1, &id2);


		//skip 2 segments insertion temp

		//if (!m_is_fix_tail_head)
		//	cout << "sscanf finished"<<endl;

		//if (!m_is_fix_tail_head)
		//	cout << count++<<endl;

		//string doner_str = doner_seq;

		//string acceptor_str = acceptor_seq;
		//
		//cout << doner_str << endl;
		//cout << m_chrom_seq.substr(doner_end - doner_str.length(), doner_str.length())<<endl;
		//cout << acceptor_str << endl;
		//cout << m_chrom_seq.substr(acceptor_st + 1, acceptor_str.length())<<endl;

		//if (acceptor_st && doner_end)
		//	getchar();

		if (strand == '+')
		{
			++doner_end_id;

			if (acceptor_st_id == 0)
				acceptor_st_id = m_num_seg;
			else
				--acceptor_st_id;
		}
		else
		{
			--doner_end_id;

			if (acceptor_st_id == 0)
				acceptor_st_id = 1;
			else
				++acceptor_st_id;
		}

		if (!m_is_fix_tail_head && ((acceptor_st - doner_end + 1) < 2 * m_seg_len))
			continue;

		//if (!m_is_fix_tail_head)
		//	cout << "cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady"<<endl;

		cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2));

		size_t hmer_offset;

		//if (m_reverse_flag)
		//	hmer_offset = m_cur_chrom_size - acceptor_st + double_seg_len;
		//else

		hmer_offset = doner_end + double_seg_len + 1;

		m_chrom_offset = -1;

		m_max_doner_offset = 0, m_max_acceptor_offset = 0;

		m_bad_bases = ALL_BITS_ON;
		m_bad_bases_reserve = ALL_BITS_ON;

		//m_seed_key = 0;
		//m_seed_keyII = 0;

		for (size_t i = 0; i < m_kmer_keys.size(); ++i)
		{
			m_kmer_keys[i] = 0;
			m_kmer_keys_upper[i] = 0;
			m_kmer_keys_lower[i] = 0;
			//m_imer_keys[i] = 0;
		}

		while (chrom_reader.read())
		{
			size_t base;
			while ((base = chrom_reader.ScanBase()) != -1)
			{
				//cout << m_chrom_offset << endl;
				//if (m_chrom_offset == 601)
				//	getchar();
				++m_chrom_offset;

				GetSeed(base);

				m_bad_bases_reserve = ((m_bad_bases_reserve << 1) + static_cast<size_t>(m_bad_bases & m_masks.reserve_bit));

				m_bad_bases = (m_bad_bases << 1) + static_cast<size_t>(base == 4);

				//skip m_chrom_offset not in exons
				while (m_chrom_offset + m_anchor_width > m_cur_island_end + m_boundary + m_seg_len)
				{
					++m_cur_tac_pos;
					if (m_cur_tac_pos == m_tac_pos_end)
						return;

					m_cur_island_start = m_cur_tac_pos->first;
					m_cur_island_end = m_cur_tac_pos->second;

					if (m_chrom_offset + m_anchor_width < m_cur_island_start + m_min_boundary + m_seg_len || 
						(m_chrom_offset + m_anchor_width >= m_cur_island_start + m_min_boundary + m_seg_len && m_chrom_offset + m_anchor_width <= m_cur_island_end + m_boundary + m_seg_len))
						break;
				}

				if (m_chrom_offset == hmer_offset && !input_fs.eof())
				{
					char chromnamechr[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000], */strand;

					size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id/*, other_mapped*/;

					int id1, id2;

					bool initialized = false;

					////HWI-EAS217:3:100:1000:1066#0/2	chr15	+	88784825	0	1-1	0-0	CAAGATCGGGGGCATCTTGGCCAAT	-	GAGCTCTCGGTGGCTGAGGCTGNAG,TCCATGCAGCTGTTCTTGCCATCAA,TGAAGCAGTGGAGCGAGGGGTGGTG,

					while (!input_fs.eof())
					{
						getline(input_fs, line);

						if (line.length() == 0)
							continue;

						sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%d\t%d", tagname, chromnamechr, &strand, &doner_end, &acceptor_st, 
							&doner_st_id, &doner_end_id, &acceptor_st_id, &acceptor_end_id,	doner_seq, acceptor_seq, tobe_fixed_seq, &id1, &id2);

						//if (!m_is_fix_tail_head)
						//	cout << count++<<endl;

						initialized = true;

						if (strand == '+')
						{
							++doner_end_id;

							if (acceptor_st_id == 0)
								acceptor_st_id = m_num_seg;
							else
								--acceptor_st_id;
						}
						else
						{
							--doner_end_id;

							if (acceptor_st_id == 0)
								acceptor_st_id = 1;
							else
								++acceptor_st_id;
						}

						if (!m_is_fix_tail_head && ((acceptor_st - doner_end + 1) < 2 * m_seg_len))
							continue;

						//if (m_reverse_flag)
						//	hmer_offset = m_cur_chrom_size - acceptor_st + double_seg_len;
						//else

						hmer_offset = doner_end + double_seg_len + 1;

						if (m_chrom_offset == hmer_offset)
						{
							cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2));
						}
						else
						{
							break;
						}

						//
						//++cur_fixhmer_iter;	
						//if (cur_fixhmer_iter != grouped_bwts.end())
						//	hmer_offset = m_reverse_flag ? (*cur_fixhmer_iter)->hmer_offset_rev : (*cur_fixhmer_iter)->hmer_offset;
					}

					//cout << "enter slide window " << m_chrom_offset << endl;
					//if (!m_is_fix_tail_head)
					//	cout << "EnterSlideWindowFixHmerReady"<<endl;
					EnterSlideWindowFixHmerReady(cur_tobe_fixed_hmer_segments, chromnamestr.c_str(), avoid);

					//if (!m_is_fix_tail_head)
					//	cout << "EnterSlideWindowFixHmerReady finished"<<endl;

					cur_tobe_fixed_hmer_segments.clear();

					if (initialized)
					{
						if (!m_is_fix_tail_head && ((acceptor_st - doner_end + 1) < 2 * m_seg_len))
							;
						else
							cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2));
					}
				}

				//if (m_chrom_offset == 93)
				//	int tt = 0;
				//check hmers
				if (m_chrom_offset + m_anchor_width < m_cur_island_start + m_min_boundary + m_seg_len || m_chrom_offset + m_anchor_width > m_max_acceptor_offset + m_seg_len)
					continue;				

				++m_basecount;

				m_wp_reserve.bads = m_bad_bases_reserve;
				m_wp.bads = m_bad_bases;

				//check suffix kmer
				hash_map<size_t, int> mapped_TagSeed;

				for (size_t i = 0; i < m_num_anchors; ++i)
				{

					// insertion
//					if (!(m_masks.bad_kmer_masks_lower[i] & m_bad_bases_reserve))
//					{
//#ifdef TIME
//						clock_t t3=clock();
//#endif
//
//						SEEDQUEUE& prefix_matched_kmers = m_open_reads_suffix[i].find_prefix_matches(m_imer_keys[i], m_chrom_offset);
//
//#ifdef TIME
//						clock_t t4=clock();
//
//						m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
//#endif
//
//						if (prefix_matched_kmers.Size())
//						{
//							size_t lit = prefix_matched_kmers.Rear();
//
//							size_t sz = prefix_matched_kmers.Size();
//
//							size_t k;
//
//							for (k = 0; k < sz && m_chrom_offset + m_anchor_width <  prefix_matched_kmers.At(lit).offset + m_seg_len; ++k)
//							{
//								//Trigger Match
//								//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;
//
//								if (m_chrom_offset - 2 > prefix_matched_kmers.At(lit).bef_offest)
//								{
//									//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
//									prefix_matched_kmers.Pred(lit);
//									continue;
//								}
//
//								if (prefix_matched_kmers.At(lit).offset + m_seg_len <= m_chrom_offset + m_anchor_width + m_max_insert_dist)
//								{
//									//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
//									prefix_matched_kmers.Pred(lit);
//									continue;
//								}
//
//#ifdef TIME
//								clock_t t1=clock();
//#endif
//								//cout << "enter triggmer match "<<m_chrom_offset<<endl;
//
//								//skip trigger match
//								//continue;
//								TriggerMatchFixHmerReady(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
//
//#ifdef TIME
//								clock_t t2=clock();
//#endif
//
//								prefix_matched_kmers.Pred(lit);
//
//#ifdef TIME
//								clock_t t5=clock();
//
//								m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;
//
//								m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
//#endif
//							}
//
//#ifdef TIME
//							clock_t t6=clock();
//#endif
//
//							//if (k != sz)
//							//	prefix_matched_kmers.RemoveFront(lit);
//
//#ifdef TIME
//							clock_t t7=clock();
//
//							m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
//#endif
//						}
//
//						//
//					}
					 
					if (!(m_masks.bad_kmer_masks_lower[i] & m_bad_bases_reserve))
					{
#ifdef TIME
						clock_t t3=clock();
#endif

						SEEDQUEUE& prefix_matched_kmers = m_open_reads_suffix[i].find_prefix_matches(m_kmer_keys_upper[i], m_chrom_offset);

#ifdef TIME
						clock_t t4=clock();

						m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

						if (prefix_matched_kmers.Size())
						{
							size_t lit = prefix_matched_kmers.Rear();

							size_t sz = prefix_matched_kmers.Size();

							size_t k;

							for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist + m_seg_len - 2 >= m_chrom_offset + m_anchor_width; ++k)
							{
								//Trigger Match
								//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;

#ifdef TIME
								clock_t t1=clock();
#endif
								
								if (m_chrom_offset + m_anchor_width + 1 < prefix_matched_kmers.At(lit).offset + m_seg_len)//insertion
								{
									//too large insertion
									if (m_chrom_offset + m_anchor_width + m_max_insert_dist < prefix_matched_kmers.At(lit).offset + m_seg_len || !m_is_fix_tail_head)
									{
										prefix_matched_kmers.Pred(lit);
										continue;
									}

									//do insertion here
									TriggerMatchFixHmerInsertReady(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
								}
								//else if (m_chrom_offset + m_anchor_width + 1 == prefix_matched_kmers.At(lit).offset + m_seg_len)//unspliced
								//{
								//}
								else//spliced
								{
									//too large splice
									if (m_chrom_offset + m_anchor_width - 2 - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
									{
										//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;

										if (m_chrom_offset + m_anchor_width - 2 - m_seg_len - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
										{
											prefix_matched_kmers.Pred(lit);
											continue;
										}

										//TriggerMatchFixHmerInsertReady(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
									}

									//cout << "enter triggmer match "<<m_chrom_offset<<endl;

									//skip trigger match
									//continue;
									TriggerMatchFixHmerReady(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
								}

#ifdef TIME
								clock_t t2=clock();
#endif

								prefix_matched_kmers.Pred(lit);

#ifdef TIME
								clock_t t5=clock();

								m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

								m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
							}

#ifdef TIME
							clock_t t6=clock();
#endif

							if (k != sz)
								prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
							clock_t t7=clock();

							m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif
						}

						//
					}

					if (!(m_masks.bad_kmer_masks_lower[i] & m_bad_bases))
					{
#ifdef TIME
						clock_t t3=clock();
#endif

						SEEDQUEUE& prefix_matched_kmers = m_open_reads_prefix[i].find_prefix_matches(m_kmer_keys_lower[i], m_chrom_offset);

#ifdef TIME
						clock_t t4=clock();

						m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

						if (prefix_matched_kmers.Size())
						{
							size_t lit = prefix_matched_kmers.Rear();

							size_t sz = prefix_matched_kmers.Size();

							size_t k;

							for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist + m_seg_len - 2 >= m_chrom_offset + m_anchor_width; ++k)
							{
								//Trigger Match
								//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;

								if (m_chrom_offset + m_anchor_width - 2 - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
								{
									//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
									prefix_matched_kmers.Pred(lit);
									continue;
								}
#ifdef TIME
								clock_t t1=clock();
#endif
								//cout << "enter triggmer match "<<m_chrom_offset<<endl;

								//skip trigger match
								//continue;
								TriggerMatchFixHmerReady(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);

#ifdef TIME
								clock_t t2=clock();
#endif

								prefix_matched_kmers.Pred(lit);

#ifdef TIME
								clock_t t5=clock();

								m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

								m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
							}

#ifdef TIME
							clock_t t6=clock();
#endif

							if (k != sz)
								prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
							clock_t t7=clock();

							m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif
						}
						//
					}


				}
			}
		}
	}
}


void
GenomeScan::FixHmerReady_Rev(const char* fixhmer_file, const char* fixed_hmer_file, size_t flankcase, bool avoid, bool island)
{
	string sort_cmd = "sort -k5,5n -S 3500000 -o ";

	sort_cmd.append(fixhmer_file);

	sort_cmd.append(" -T temp ");

	sort_cmd.append(fixhmer_file);

	system(sort_cmd.c_str());

	string tac_command = "tac ";

	tac_command.append(fixhmer_file);

	tac_command.append(" >");

	tac_command.append(fixhmer_file);

	tac_command.append(".rev");

	system(tac_command.c_str());

	string mv_command("mv ");

	mv_command.append(fixhmer_file);

	mv_command.append(".rev ");

	mv_command.append(fixhmer_file);

	system(mv_command.c_str());

	string fixhmer_file_str = fixhmer_file;

	string chromname = basename2(fixhmer_file_str);

	if (m_chrom_names_map.find(chromname) == m_chrom_names_map.end())
	{
		cout << "can't find chromosome: " << chromname << endl;
		return;
	}

	string chrom_file;

	m_cur_chrom_size = m_chrom_sizes[m_chrom_names_map[chromname]];

	TAC_POS_MAP* tac_pos_ptr;
	
	//if (m_reverse_flag)
	//{
	chrom_file = m_rev_chrom_files[m_chrom_names_map[chromname]];

	tac_pos_ptr = &m_tac_pos_map_rev;
	//}
	//else
	//{
	//	chrom_file = m_chrom_files[m_chrom_names_map[chromname]];

	//	tac_pos_ptr = &m_tac_pos_map;
	//}

	//load chromo sequence into memory
	m_chrom_seq.clear();
	readchrom(chrom_file.c_str(), m_chrom_seq);

	m_cur_island_start = /*m_boundary + */1;

	m_cur_island_end = max_size_t - (3 * m_boundary);

	m_islands = "N";

	m_islands_provide = false;

	TAC_POS_MAP_ITER chrom_pos_map = tac_pos_ptr->find(chromname);

	string chromnamestr = chromname;

	if (chrom_pos_map != tac_pos_ptr->end() && island)
	{
		m_islands = "Y";
		m_islands_provide = true;
		m_cur_tac_pos = chrom_pos_map->second.begin();
		m_tac_pos_begin = chrom_pos_map->second.begin();
		m_tac_pos_end = chrom_pos_map->second.end();
		m_cur_island_start = m_cur_tac_pos->first;
		m_cur_island_end = m_cur_tac_pos->second;
	}

	cout << "islands: "<<m_islands<<endl;

	ChromReader chrom_reader(chrom_file, m_input_buffer_size);

	ifstream input_fs(fixhmer_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhmer_file<<"\n";return;//exit(1);
	}

	m_ofs_fix_hmer.close();

	m_ofs_fix_hmer.open(fixed_hmer_file);

	string fixed_hmer_file_hole = fixed_hmer_file; fixed_hmer_file_hole.append(".hole");
	m_ofs_fix_hmer_hole.close();
	m_ofs_fix_hmer_hole.open(fixed_hmer_file_hole.c_str());

	size_t double_seg_len = 2 * m_seg_len;

	vector<TobeFixedHmerReady> cur_tobe_fixed_hmer_segments;

	if (m_is_fix_tail_head)
	{
		m_cur_junc_head_iter = m_cur_junc_vec_ptr->begin(); 

		m_cur_junc_tail_iter = m_cur_junc_vec_ptr->begin();

		m_cur_junc_end_iter = m_cur_junc_vec_ptr->end();
	}
	
	while(!input_fs.eof())
	{
		string line;

		getline(input_fs, line);

		if (line.length() == 0)
			continue;

		char chromnamechr[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000], */strand;

		size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id/*, other_mapped*/;
		
		int id1, id2;

		//HWI-EAS217:3:100:1000:1066#0/2	chr15	+	88784825	0	1-1	0-0	CAAGATCGGGGGCATCTTGGCCAAT	-	GAGCTCTCGGTGGCTGAGGCTGNAG,TCCATGCAGCTGTTCTTGCCATCAA,TGAAGCAGTGGAGCGAGGGGTGGTG,

		sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%d\t%d", tagname, chromnamechr, &strand, &acceptor_st, &doner_end, 
			&acceptor_end_id, &acceptor_st_id, &doner_end_id, &doner_st_id, acceptor_seq, doner_seq, tobe_fixed_seq, &id1, &id2);

		//string doner_str = doner_seq;

		//string acceptor_str = acceptor_seq;
		//
		//cout << doner_str << endl;
		//cout << m_chrom_seq.substr(doner_end - doner_str.length(), doner_str.length())<<endl;
		//cout << acceptor_str << endl;
		//cout << m_chrom_seq.substr(acceptor_st + 1, acceptor_str.length())<<endl;

		//if (acceptor_st && doner_end)
		//	getchar();

		if (strand == '-')
		{
			++doner_end_id;

			if (acceptor_st_id == 0)
				acceptor_st_id = m_num_seg;
			else
				--acceptor_st_id;
		}
		else
		{
			--doner_end_id;

			if (acceptor_st_id == 0)
				acceptor_st_id = 1;
			else
				++acceptor_st_id;
		}

		//cout << line << endl;
		//cout << doner_end_id << '\t'<<acceptor_st_id << endl;
		//getchar();

		cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2, true));

		size_t hmer_offset;

		//if (m_reverse_flag)
		hmer_offset = m_cur_chrom_size - doner_end + double_seg_len;
		//else
		//	hmer_offset = doner_end + double_seg_len + 1;

		m_chrom_offset = -1;

		m_max_doner_offset = 0, m_max_acceptor_offset = 0;

		m_bad_bases = ALL_BITS_ON;
		m_bad_bases_reserve = ALL_BITS_ON;

		//m_seed_key = 0;
		//m_seed_keyII = 0;

		for (size_t i = 0; i < m_kmer_keys.size(); ++i)
		{
			m_kmer_keys[i] = 0;
			m_kmer_keys_upper[i] = 0;
			m_kmer_keys_lower[i] = 0;
			//m_imer_keys[i] = 0;
		}

		while (chrom_reader.read())
		{
			size_t base;
			while ((base = chrom_reader.ScanBase()) != -1)
			{
				++m_chrom_offset;

				GetSeed(base);

				m_bad_bases_reserve = ((m_bad_bases_reserve << 1) + static_cast<size_t>(m_bad_bases & m_masks.reserve_bit));

				m_bad_bases = (m_bad_bases << 1) + static_cast<size_t>(base == 4);

				//skip m_chrom_offset not in exons
				while (m_chrom_offset + m_anchor_width > m_cur_island_end + m_boundary + m_seg_len)
				{
					++m_cur_tac_pos;
					if (m_cur_tac_pos == m_tac_pos_end)
						return;

					m_cur_island_start = m_cur_tac_pos->first;
					m_cur_island_end = m_cur_tac_pos->second;

					if (m_chrom_offset + m_anchor_width < m_cur_island_start + m_min_boundary + m_seg_len || 
						(m_chrom_offset + m_anchor_width >= m_cur_island_start + m_min_boundary + m_seg_len && m_chrom_offset + m_anchor_width <= m_cur_island_end + m_boundary + m_seg_len))
						break;
				}

				if (m_chrom_offset == hmer_offset && !input_fs.eof())
				{
					char chromnamechr[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000],*/ strand;

					size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id, other_mapped;

					int id1, id2;

					bool initialized = false;

					////HWI-EAS217:3:100:1000:1066#0/2	chr15	+	88784825	0	1-1	0-0	CAAGATCGGGGGCATCTTGGCCAAT	-	GAGCTCTCGGTGGCTGAGGCTGNAG,TCCATGCAGCTGTTCTTGCCATCAA,TGAAGCAGTGGAGCGAGGGGTGGTG,

					while (!input_fs.eof())
					{
						getline(input_fs, line);

						if (line.length() == 0)
							continue;

						sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%d\t%d", tagname, chromnamechr, &strand, &acceptor_st, &doner_end, 
							&acceptor_end_id, &acceptor_st_id, &doner_end_id, &doner_st_id, acceptor_seq, doner_seq, tobe_fixed_seq, &id1, &id2);
						                                  

						initialized = true;

						if (strand == '-')
						{
							++doner_end_id;

							if (acceptor_st_id == 0)
								acceptor_st_id = m_num_seg;
							else
								--acceptor_st_id;
						}
						else
						{
							--doner_end_id;

							if (acceptor_st_id == 0)
								acceptor_st_id = 1;
							else
								++acceptor_st_id;
						}

						//cout << line << endl;
						//cout << doner_end_id << '\t'<<acceptor_st_id << endl;
						//getchar();

						//if (m_reverse_flag)
						hmer_offset = m_cur_chrom_size - doner_end + double_seg_len;
						//else
						//	hmer_offset = doner_end + double_seg_len + 1;

						if (m_chrom_offset == hmer_offset)
						{
							cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, 
								doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2, true));
						}
						else
						{
							break;
						}

						//
						//++cur_fixhmer_iter;	
						//if (cur_fixhmer_iter != grouped_bwts.end())
						//	hmer_offset = m_reverse_flag ? (*cur_fixhmer_iter)->hmer_offset_rev : (*cur_fixhmer_iter)->hmer_offset;
					}

					//cout << "enter slide window " << m_chrom_offset << endl;
					EnterSlideWindowFixHmerReady_Rev(cur_tobe_fixed_hmer_segments, chromnamestr.c_str(), avoid);

					cur_tobe_fixed_hmer_segments.clear();

					if (initialized)
						cur_tobe_fixed_hmer_segments.push_back(TobeFixedHmerReady(tagname, doner_end, acceptor_st, doner_end_id, acceptor_st_id, tobe_fixed_seq, strand, id1, id2, true));
				}

				//check hmers
				if (m_chrom_offset + m_anchor_width < m_cur_island_start + m_min_boundary + m_seg_len || m_chrom_offset + m_anchor_width > m_max_acceptor_offset + m_seg_len)
					continue;

				++m_basecount;

				m_wp_reserve.bads = m_bad_bases_reserve;
				m_wp.bads = m_bad_bases;

				//check suffix kmer
				hash_map<size_t, int> mapped_TagSeed;

				for (size_t i = 0; i < m_num_anchors; ++i)
				{
					if (!(m_masks.bad_kmer_masks_lower[i] & m_bad_bases_reserve))
					{
#ifdef TIME
						clock_t t3=clock();
#endif

						SEEDQUEUE& prefix_matched_kmers = m_open_reads_suffix[i].find_prefix_matches(m_kmer_keys_upper[i], m_chrom_offset);

#ifdef TIME
						clock_t t4=clock();

						m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

						if (prefix_matched_kmers.Size())
						{
							size_t lit = prefix_matched_kmers.Rear();

							size_t sz = prefix_matched_kmers.Size();

							size_t k;

							for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist + m_seg_len - 2 >= m_chrom_offset + m_anchor_width; ++k)
							{
								//Trigger Match
								//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;

#ifdef TIME
								clock_t t1=clock();
#endif

								if (m_chrom_offset + m_anchor_width + 1< prefix_matched_kmers.At(lit).offset + m_seg_len)//insertion
								{
									//too large insertion
									if (m_chrom_offset + m_anchor_width + m_max_insert_dist < prefix_matched_kmers.At(lit).offset + m_seg_len)
									{
										prefix_matched_kmers.Pred(lit);
										continue;
									}

									//do insertion here

									TriggerMatchFixHmerInsertReady_Rev(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
								}
								//else if (m_chrom_offset + m_anchor_width + 1 == prefix_matched_kmers.At(lit).offset + m_seg_len)//unspliced
								//{
								//}
								else//spliced
								{
									//too large splice
									if (m_chrom_offset + m_anchor_width - 2 - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
									{
										//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;

										if (m_chrom_offset + m_anchor_width - 2 - m_seg_len - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
										{
											prefix_matched_kmers.Pred(lit);
											continue;
										}

										//TriggerMatchFixHmerInsertReady_Rev(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
									}

									//cout << "enter triggmer match "<<m_chrom_offset<<endl;

									//skip trigger match
									//continue;
									TriggerMatchFixHmerReady_Rev(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);
								}

								//if (m_chrom_offset - 2 > prefix_matched_kmers.At(lit).bef_offest)
								//{
								//	//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
								//	prefix_matched_kmers.Pred(lit);
								//	continue;
								//}

								//cout << "enter triggmer match "<<m_chrom_offset<<endl;

								//skip trigger match
								//continue;
								//TriggerMatchFixHmerReady_Rev(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);

#ifdef TIME
								clock_t t2=clock();
#endif

								prefix_matched_kmers.Pred(lit);

#ifdef TIME
								clock_t t5=clock();

								m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

								m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
							}

#ifdef TIME
							clock_t t6=clock();
#endif

							if (k != sz)
								prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
							clock_t t7=clock();

							m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif
						}

						//
					}

					if (!(m_masks.bad_kmer_masks_lower[i] & m_bad_bases))
					{
#ifdef TIME
						clock_t t3=clock();
#endif

						SEEDQUEUE& prefix_matched_kmers = m_open_reads_prefix[i].find_prefix_matches(m_kmer_keys_lower[i], m_chrom_offset);

#ifdef TIME
						clock_t t4=clock();

						m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

						if (prefix_matched_kmers.Size())
						{
							size_t lit = prefix_matched_kmers.Rear();

							size_t sz = prefix_matched_kmers.Size();

							size_t k;

							for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist + m_seg_len - 2 >= m_chrom_offset + m_anchor_width; ++k)
							{
								//Trigger Match
								//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;

								if (m_chrom_offset + m_anchor_width - 2 - m_seg_len > prefix_matched_kmers.At(lit).bef_offest)
								{
									//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
									prefix_matched_kmers.Pred(lit);
									continue;
								}
#ifdef TIME
								clock_t t1=clock();
#endif
								//cout << "enter triggmer match "<<m_chrom_offset<<endl;

								//skip trigger match
								//continue;
								TriggerMatchFixHmerReady_Rev(prefix_matched_kmers.At(lit), chromnamestr.c_str(), i, mapped_TagSeed, flankcase);

#ifdef TIME
								clock_t t2=clock();
#endif

								prefix_matched_kmers.Pred(lit);

#ifdef TIME
								clock_t t5=clock();

								m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

								m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
							}

#ifdef TIME
							clock_t t6=clock();
#endif

							if (k != sz)
								prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
							clock_t t7=clock();

							m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif
						}
						//
					}
				}
			}
		}
	}
}



void
GenomeScan::FixHmer(size_t chromID, vector<GroupedBwts*>& grouped_bwts, size_t flankcase, bool avoid)
{
	string chrom_file;

	TAC_POS_MAP* tac_pos_ptr;

	if (m_reverse_flag)
	{
		chrom_file = m_rev_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map_rev;
	}
	else
	{
		chrom_file = m_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map;
	}

	ChromReader chrom_reader(chrom_file, m_input_buffer_size);

	TAC_POS_MAP_ITER chrom_pos_map = tac_pos_ptr->find(m_chrom_names[chromID]);

	m_cur_island_start = m_boundary + 1;
	
	m_cur_island_end = max_size_t - (3 * m_boundary);

	m_islands = "N";

	m_islands_provide = false;

	if (chrom_pos_map != tac_pos_ptr->end())
	{
		m_islands = "Y";
		m_islands_provide = true;
		m_cur_tac_pos = chrom_pos_map->second.begin();
		m_tac_pos_begin = chrom_pos_map->second.begin();
		m_tac_pos_end = chrom_pos_map->second.end();
		m_cur_island_start = m_cur_tac_pos->first;
		m_cur_island_end = m_cur_tac_pos->second;
	}

	cout << "islands: "<<m_islands<<endl;

	m_chrom_offset = -1;

	m_bad_bases = ALL_BITS_ON;
	m_bad_bases_reserve = ALL_BITS_ON;

	m_seed_key = 0;
	m_seed_keyII = 0;

	for (size_t i = 0; i < m_kmer_keys.size(); ++i)
		m_kmer_keys[i] = 0;

	vector<GroupedBwts* >::iterator cur_fixhmer_iter = grouped_bwts.begin();

	while (chrom_reader.read())
	{
		size_t base;
		while ((base = chrom_reader.ScanBase()) != -1)
		{
			++m_chrom_offset;

			GetSeed(base);

			m_bad_bases_reserve = ((m_bad_bases_reserve << 1) + static_cast<size_t>(m_bad_bases & MOST_SIG_BIT));

			m_bad_bases = (m_bad_bases << 1) + static_cast<size_t>(base == 4);

			//skip m_chrom_offset not in exons
			while (m_chrom_offset > m_cur_island_end + m_boundary)
			{
				++m_cur_tac_pos;
				if (m_cur_tac_pos == m_tac_pos_end)
					return;

				m_cur_island_start = m_cur_tac_pos->first;
				m_cur_island_end = m_cur_tac_pos->second;

				if (m_chrom_offset < m_cur_island_start + m_min_boundary || (m_chrom_offset >= m_cur_island_start + m_min_boundary && m_chrom_offset <= m_cur_island_end + m_boundary))
					break;
			}

			size_t hmer_offset;
			
			if (cur_fixhmer_iter != grouped_bwts.end())
			{
				hmer_offset = m_reverse_flag ? (*cur_fixhmer_iter)->hmer_offset_rev : (*cur_fixhmer_iter)->hmer_offset;
				//cout << "hmer_offset "<< hmer_offset<<endl;
			}

			if (cur_fixhmer_iter != grouped_bwts.end() && m_chrom_offset == hmer_offset + m_seg_len - 1)
			{
				vector<GroupedBwts* >::iterator first_fixhmer_iter = cur_fixhmer_iter++;

				if (cur_fixhmer_iter != grouped_bwts.end())
					hmer_offset = m_reverse_flag ? (*cur_fixhmer_iter)->hmer_offset_rev : (*cur_fixhmer_iter)->hmer_offset;
				
				while (cur_fixhmer_iter != grouped_bwts.end() && m_chrom_offset == hmer_offset + m_seg_len - 1)
				{
					
					//insert first half into slide window
					++cur_fixhmer_iter;	
					if (cur_fixhmer_iter != grouped_bwts.end())
						hmer_offset = m_reverse_flag ? (*cur_fixhmer_iter)->hmer_offset_rev : (*cur_fixhmer_iter)->hmer_offset;
				}

				//cout << "enter slide window " << m_chrom_offset << endl;
				EnterSlideWindowFixHmer(first_fixhmer_iter, cur_fixhmer_iter, avoid);
			}


			//check hmers
			if (m_chrom_offset < m_cur_island_start + m_min_boundary)
				continue;

			++m_basecount;

			m_wp_reserve.bads = m_bad_bases_reserve;
			m_wp.bads = m_bad_bases;

			//check suffix kmer
			hash_map<size_t, int> mapped_TagSeed;
			for (size_t i = 0; i < m_num_anchors; ++i)
			{
				if (!(m_masks.bad_kmer_masks[i] & m_bad_bases))
				{
#ifdef TIME
					clock_t t3=clock();
#endif

					SEEDQUEUE& prefix_matched_kmers = m_open_reads[i].find_prefix_matches(m_kmer_keys[i], m_chrom_offset);

#ifdef TIME
					clock_t t4=clock();

					m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

					if (!prefix_matched_kmers.Size())
						continue;

					size_t lit = prefix_matched_kmers.Rear();

					size_t sz = prefix_matched_kmers.Size();

					size_t k;

					for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist - 2 >= m_chrom_offset; ++k)
					{
						//Trigger Match
						//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;

						if (m_chrom_offset > prefix_matched_kmers.At(lit).bef_offest)
						{
							//cout << "skip: " << m_chrom_offset << '\t'<< prefix_matched_kmers.At(lit).bef_offest<<endl;
							prefix_matched_kmers.Pred(lit);
							continue;
						}
#ifdef TIME
						clock_t t1=clock();
#endif
						//cout << "enter triggmer match "<<m_chrom_offset<<endl;

						TriggerMatchFixHmer(prefix_matched_kmers.At(lit), chromID, i, mapped_TagSeed, flankcase);

#ifdef TIME
						clock_t t2=clock();
#endif

						prefix_matched_kmers.Pred(lit);

#ifdef TIME
						clock_t t5=clock();

						m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

						m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
					}

#ifdef TIME
					clock_t t6=clock();
#endif

					if (k != sz)
						prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
					clock_t t7=clock();

					m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif

					//
				}
			}

		}
	}

}

void
GenomeScan::PrepareFixhole(const char* fixhole_file, const char* fixed_hole_file)
{
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;

		for (size_t i = 0; i < rec_num; ++i)
		{
			getline(input_fs, line);

			if (line == "")
				continue;

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);
			fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		}

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		size_t missed_id;

		vector<GroupedBwtRec>::iterator bwt_rec_iter;

		for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		{
			if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			{
				if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
				{
					cout << "missed more than 1 segments"<<endl;
					exit(1);
				}
				else
				{
					missed_id = bwt_rec_iter->seg_id + 1;
					break;
				}
			}
		}

		size_t readID;

		if (strand == '+')
			readID = tagidx * m_num_seg + missed_id - 1;
		else
			readID = tagidx * m_num_seg + m_num_seg - (missed_id - 1);

	}

	input_fs.close();
	ofs.close();
}

void
GenomeScan::MapReads(size_t chromID)
{
	//size_t m_cur_island_start;

	//size_t m_cur_island_end;

	TAC_POS_MAP* tac_pos_ptr;

	string chrom_files;

	if (!m_reverse_flag)
	{
		chrom_files = m_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map;//.find(m_chrom_names[chromID]);
	}
	else
	{
		chrom_files = m_rev_chrom_files[chromID];

		tac_pos_ptr = &m_tac_pos_map_rev;//.find(m_chrom_names[chromID]);
	}

	ChromReader chrom_reader(chrom_files, m_input_buffer_size);

	TAC_POS_MAP_ITER chrom_pos_map = tac_pos_ptr->find(m_chrom_names[chromID]);

	//cout << chrom_pos_map->second.

	m_cur_island_start = m_boundary + 1;
	
	m_cur_island_end = max_size_t - (3 * m_boundary);

	m_islands = "N";

	m_islands_provide = false;

	if (chrom_pos_map != tac_pos_ptr->end())
	{
		m_islands = "Y";
		m_islands_provide = true;
		m_cur_tac_pos = chrom_pos_map->second.begin();
		m_tac_pos_begin = chrom_pos_map->second.begin();
		m_tac_pos_end = chrom_pos_map->second.end();
		m_cur_island_start = m_cur_tac_pos->first;
		m_cur_island_end = m_cur_tac_pos->second;
	}

	cout << "islands: "<<m_islands<<endl;

	//bug fixed 1
	m_chrom_offset = -1;

	m_bad_bases = ALL_BITS_ON;
	m_bad_bases_reserve = ALL_BITS_ON;

	m_seed_key = 0;
	m_seed_keyII = 0;

	for (size_t i = 0; i < m_kmer_keys.size(); ++i)
		m_kmer_keys[i] = 0;

	while (chrom_reader.read())
	{
		size_t base;
		while ((base = chrom_reader.ScanBase()) != -1)
		{
			++m_chrom_offset;

			GetSeed(base);

			//if (( m_chrom_offset == 7944396 && m_seed_key == 1036894) ||  (m_chrom_offset > 7945400 && m_kmer_key == 1228))
			//{
			//	int test = 0;
			//}

			m_bad_bases_reserve = ((m_bad_bases_reserve << 1) + static_cast<size_t>(m_bad_bases & MOST_SIG_BIT));

			m_bad_bases = (m_bad_bases << 1) + static_cast<size_t>(base == 4);

			//skip m_chrom_offset not in exons
			while (m_chrom_offset > m_cur_island_end + m_boundary)
			{
				++m_cur_tac_pos;
				if (m_cur_tac_pos == m_tac_pos_end)
					return;

				m_cur_island_start = m_cur_tac_pos->first;
				m_cur_island_end = m_cur_tac_pos->second;

				if (m_chrom_offset < m_cur_island_start + m_min_boundary || (m_chrom_offset >= m_cur_island_start + m_min_boundary && m_chrom_offset <= m_cur_island_end + m_boundary))
					break;
			}

			if (m_chrom_offset < m_cur_island_start + m_min_boundary)
				continue;

			++m_basecount;

			m_wp_reserve.bads = m_bad_bases_reserve;
			m_wp.bads = m_bad_bases;

			size_t prev_hits = 0;
			//check prefix
			if (!(m_masks.bad_base_mask & m_bad_bases))
			{
				#ifdef TIME
				clock_t t3=clock();
				#endif
				
				SEEDS_CONST_ITER bucket(m_seed_hash.GetTags(m_seed_key));

				#ifdef TIME
				clock_t t4=clock();

				m_checkseedtable_time += (t4-t3)/(double)CLOCKS_PER_SEC;

				#endif

				if (bucket != m_seed_hash.End())
				{
					prev_hits = EnterSlideWindow(bucket->second);
				}
			}

			//check second seedkey
			if (m_max_mismatches && !(m_masks.bad_base_maskII & m_bad_bases))
			{
				#ifdef TIME
				clock_t t3=clock();
				#endif

				SEEDS_CONST_ITER bucket(m_seed_hash.GetTagsII(m_seed_keyII));

				#ifdef TIME
				clock_t t4=clock();

				m_checkseedtable_time += (t4-t3)/(double)CLOCKS_PER_SEC;
				#endif

				if (bucket != m_seed_hash.EndII())
				{
					EnterSlideWindowII(bucket->second, prev_hits);
				}
			}

			//check suffix kmer
			hash_map<size_t, int> mapped_TagSeed;
			for (size_t i = 0; i < m_num_anchors; ++i)
			{
				if (!(m_masks.bad_kmer_masks[i] & m_bad_bases))
				{
#ifdef TIME
					clock_t t3=clock();
#endif

					SEEDQUEUE& prefix_matched_kmers = m_open_reads[i].find_prefix_matches(m_kmer_keys[i], m_chrom_offset);

#ifdef TIME
					clock_t t4=clock();

					m_checkcandclearopenread_time += (t4-t3)/(double)CLOCKS_PER_SEC;
#endif

					if (!prefix_matched_kmers.Size())
						continue;

					size_t lit = prefix_matched_kmers.Rear();

					size_t sz = prefix_matched_kmers.Size();

					size_t k;

					for (k = 0; k < sz && prefix_matched_kmers.At(lit).offset + m_max_splice_dist - 2 >= m_chrom_offset; ++k)
					{
						//Trigger Match
						//WordPair secondHalfTag = m_slide_window[lit->offset % m_window_size].matchedTagsHalf[lit->index].secondHalfT;
#ifdef TIME
						clock_t t1=clock();
#endif

						TriggerMatch(prefix_matched_kmers.At(lit), chromID, i, mapped_TagSeed);

#ifdef TIME
						clock_t t2=clock();
#endif

						prefix_matched_kmers.Pred(lit);

#ifdef TIME
						clock_t t5=clock();

						m_triggermatch_time += (t2-t1)/(double)CLOCKS_PER_SEC;

						m_checkcandclearopenread_time += (t5-t2)/(double)CLOCKS_PER_SEC;
#endif
					}

#ifdef TIME
					clock_t t6=clock();
#endif

					if (k != sz)
						prefix_matched_kmers.RemoveFront(lit);

#ifdef TIME
					clock_t t7=clock();

					m_checkcandclearopenread_time += (t7-t6)/(double)CLOCKS_PER_SEC;
#endif

					//
				}
			}
		}
	}
}

inline void 
GenomeScan::GetSeed(const size_t base)
{
	m_wp_reserve.shift_reserve(m_wp, m_masks.reserve_bit);
	m_wp.shift(base);

	//create the new seed_key 
	//m_wp.update_key(m_masks.key_update_bit, m_masks.small_mask, m_seed_key);

	//m_wp.update_key(m_masks.key_update_bitII, m_masks.small_mask, m_seed_keyII);

	//kmer_update_bits_upper[i] = LEAST_SIG_BIT << (i * anchor_width);

	//kmer_update_bits_lower[i] = LEAST_SIG_BIT << seg_width - ((i + 1) * anchor_width) + 1;

	//vector<size_t> m_kmer_keys_upper;

	//vector<size_t> m_kmer_keys_lower;

	for (size_t i = 0; i < m_num_anchors; ++i)
	{
		//m_wp.update_key(m_masks.kmer_update_bits[i], m_masks.kmer_mask, m_kmer_keys[i]);

		m_wp_reserve.update_key(m_masks.kmer_update_bits_lower[i], m_masks.kmer_mask, m_kmer_keys_upper[i]);

		//m_wp_reserve.update_key(m_masks.kmer_update_bits_lower[i], m_masks.kmer_mask, m_imer_keys[i]);

		//m_wp.update_key(m_masks.kmer_update_bits_lower[i], m_masks.kmer_mask, m_kmer_keys_lower[i]);		
	}
}

inline size_t
GenomeScan::Ps_score(const WordPair & read_word_dup, size_t& loc)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(m_shifted_wp, m_masks.mask_ps);

	bool debug = false;
	if (debug)
	{
		cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_midright) << endl;
		cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_ps) << endl;
		cerr << "suffix_wp" << endl << m_shifted_wp.tostring(ALL_BITS_ON) << endl;
		cerr << "prefix_wp" << endl << read_word_dup.tostring(ALL_BITS_ON) << endl;
	}

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((m_shifted_wp.upper ^ read_word_dup.upper) | 
		(m_shifted_wp.lower ^ read_word_dup.lower) | m_shifted_wp.bads | read_word_dup.bads) & m_masks.big_buff_mask;
	if (debug) cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;

	//the following are two pointers each point to the bit to be turned on and turned off
	size_t selector1 = LEAST_SIG_BIT << m_seed_width ; //to be turned on 
	size_t selector2 = selector1 << m_masks.mid_buff_width;//to be turned off

	loc = m_masks.mid_buff_width; // loc is the index to the left
	size_t mins = s;

	if (debug) cerr << "score :" << mins << endl;

	//check out the accumulation score of mismatches when the selectors are moving...
	for (size_t i = 1; i <= m_masks.mid_buff_width; ++i)
	{
		if (debug )
		{
			cerr << "loop " << i << "score " << s << endl;
			cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
			cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
		}

		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;     	
		if (mins > s)
		{
			mins = s;
			loc = m_masks.mid_buff_width - i;
		}
		selector1 <<= 1;
		selector2 <<= 1;
	}
	if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;

	return mins;
}

inline bool
GenomeScan::CheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = m_masks.score_buff_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);
	
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 6;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::CheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		//combine_upper = ((five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.upper & three_prim_mask)) >> i;
		//combine_lower = ((five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask)) >> i;
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i; 
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 2;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 6;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			if (s < mins)
			{
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = m_masks.score_buff_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::CheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_ATAC:
			{
				if ( s < mins)
				{
					prim = 1;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				if ( s < mins)
				{
					prim = 2;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTGC:
			{
				if ( s < mins)
				{
					prim = 3;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				if ( s < mins)
				{
					prim = 4;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				if ( s < mins)
				{
					prim = 6;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::CheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		if (combine_words == bit_GTAG)
		{
			prim = 5;
			mins = s;
			loc = m_masks.score_buff_width - i;
			m_matched_flank = combine_words;
			m_matched_bads = combine_bads;
			if (mins == 0)
				return true;
		}
	}
	return false;
}

inline bool
GenomeScan::CheckFirstTimeRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = m_masks.score_buff_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);
	
		switch (combine_words)
		{
		case bit_CATA:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CGTC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GACG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_TATG:
			{
				prim = 6;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::CheckBeforeMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		//combine_upper = ((five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.upper & three_prim_mask)) >> i;
		//combine_lower = ((five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask)) >> i;
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i; 
		switch (combine_words)
		{
		case bit_CATA:
			{
				prim = 1;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				prim = 2;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CGTC:
			{
				prim = 3;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GACG:
			{
				prim = 4;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_TATG:
			{
				prim = 6;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		default:
			if (s < mins)
			{
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = m_masks.score_buff_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::CheckAfterMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_CATA:
			{
				if (s < mins)
				{
					prim = 1;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				if (s < mins)
				{
					prim = 2;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CGTC:
			{
				if (s < mins)
				{
					prim = 3;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
					//if (mins == 0)
					//	return true;
				}
			}
			break;
		case bit_GACG:
			{
				if (s < mins)
				{
					prim = 4;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
					//if (mins == 0)
					//	return true;
				}
			}
			break;
		case bit_TATG:
			{
				if (s < mins)
				{
					prim = 6;
					mins = s;
					loc = m_masks.score_buff_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				mins = s;
				loc = m_masks.score_buff_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::CheckAfterGTAGMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		if (combine_words == bit_GATG)
		{
			prim = 5;
			mins = s;
			loc = m_masks.score_buff_width - i;
			m_matched_flank = combine_words;
			m_matched_bads = combine_bads;
			if (mins == 0)
				return true;
		}
	}
	return false;
}

inline size_t
GenomeScan::Ps_score_selective(const WordPair & read_word_dup, size_t& loc, size_t& prim, size_t& rbits)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(m_shifted_wp, m_masks.score_first_half_mask);

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = m_left_mismatches + 1;
	prim = 0;

	register size_t bits = ((m_shifted_wp.upper ^ read_word_dup.upper) | 
		(m_shifted_wp.lower ^ read_word_dup.lower) | m_shifted_wp.bads | read_word_dup.bads) & m_masks.score_big_buff_mask;

	rbits = bits;

	if (s <= m_left_mismatches && CheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s))
		return mins;
	/*{
		//cost 6 bit operations to get combined prims
		//combine_upper = (five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask);
		//combine_lower = (five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask);
		//if (CheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s))
			return mins;
	}*/

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//bool debug = false;
	//if (debug)
	//{
	//	cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_midright) << endl;
	//	cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_ps) << endl;
	//	cerr << "suffix_wp" << endl << m_shifted_wp.tostring(ALL_BITS_ON) << endl;
	//	cerr << "prefix_wp" << endl << read_word_dup.tostring(ALL_BITS_ON) << endl;
	//}

	//generate the bits where 1 indicates the mismatches.

	//if (debug) cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = m_masks.score_first_selector_rt;//score_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = m_masks.score_first_selector_lt;//score_first_selector_lt;//to be turned off
	
	//if (debug) cerr << "score :" << mins << endl;

	//check out the accumulation score of mismatches when the selectors are moving...
	size_t score_buff_width = m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= score_buff_width; ++i)
	{
		//if (debug )
		//{
		//	cerr << "loop " << i << "score " << s << endl;
		//	cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
		//	cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
		//}

		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= m_left_mismatches)
		{
			if (!prim)
			{
				if (CheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}
			else if (prim != 5)
			{
				if (CheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}
			else if (s < mins)
			{
				if (CheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}
	//if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;

	return mins;
}

inline size_t
GenomeScan::Ps_score_selective_rev(const WordPair & read_word_dup, size_t& loc, size_t& prim, size_t& rbits)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(m_shifted_wp, m_masks.score_first_half_mask);

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = m_left_mismatches + 1;
	prim = 0;

	register size_t bits = ((m_shifted_wp.upper ^ read_word_dup.upper) | 
		(m_shifted_wp.lower ^ read_word_dup.lower) | m_shifted_wp.bads | read_word_dup.bads) & m_masks.score_big_buff_mask;

	rbits = bits;

	if (s <= m_left_mismatches && CheckFirstTimeRev(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//bool debug = false;
	//if (debug)
	//{
	//	cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_midright) << endl;
	//	cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_ps) << endl;
	//	cerr << "suffix_wp" << endl << m_shifted_wp.tostring(ALL_BITS_ON) << endl;
	//	cerr << "prefix_wp" << endl << read_word_dup.tostring(ALL_BITS_ON) << endl;
	//}

	//generate the bits where 1 indicates the mismatches.
	
	//if (debug) cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = m_masks.score_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = m_masks.score_first_selector_lt;//to be turned off
	
	//if (debug) cerr << "score :" << mins << endl;

	//check out the accumulation score of mismatches when the selectors are moving...
	size_t score_buff_width = m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= score_buff_width; ++i)
	{
		//if (debug )
		//{
		//	cerr << "loop " << i << "score " << s << endl;
		//	cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
		//	cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
		//}

		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= m_left_mismatches)
		{
			if (!prim)
			{
				if (CheckBeforeMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}
			else if (prim != 5)
			{
				if (CheckAfterMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}
			else if (s < mins)
			{
				if (CheckAfterGTAGMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}
	//if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;

	return mins;
}

inline size_t
GenomeScan::Fixfusion_score(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, m_masks.prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	if (s <= left_mismatch)
		FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, m_masks.score_seg_buf_width);

	if (mins == 0)
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & m_masks.comb_seg_bits_on;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = m_masks.comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = m_masks.comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= m_seg_len/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (/*!prim*/s < mins)
			{
				FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width);

				if (mins == 0)
					return mins;
				//if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
				//	return mins;
			}
			//else if (prim != 5)
			//{
			//	if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
			//		return mins;
			//}
			//else if (s < mins)
			//{
			//	if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
			//		return mins;
			//}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}


inline size_t
GenomeScan::Fixhole_score_selective(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, m_masks.prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	if (s <= left_mismatch && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, m_masks.score_seg_buf_width))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & m_masks.comb_seg_bits_on;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = m_masks.comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = m_masks.comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= m_seg_len/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width))
					return mins;
			}
			else if (prim != 5)
			{
				if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}

inline size_t 
GenomeScan::Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (prim < 5)
			{
				if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}

inline size_t 
GenomeScan::Fixhole_score_selective_var_mask_rev(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch && FixHoleCheckFirstTimeRev(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (prim < 5)
			{
				if (FixHoleCheckAfterMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}



inline size_t 
GenomeScan::Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl; 

	//size_t five_prim_mask = LAST_THIRD_FOUTH;
	//size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;

	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch)// && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
	{
		mins = s;

		loc = mask_ptr->score_seg_buf_width;

		if (!mins) 
			return mins;
	}

	//five_prim_mask <<= 1;
	//three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off

	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on

	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;

	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (mins > s)
		{
			mins = s;

			loc = mask_ptr->score_seg_buf_width - i;
			if (!mins)
				return mins;
		}

		
		//if (s <= left_mismatch)
		//{
		//	if (!prim)
		//	{
		//		if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}
		//	else if (prim < 5)
		//	{
		//		if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}
		//	else if (s < mins)
		//	{
		//		if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		//			return mins;
		//	}

		//}
		selector1 <<= 1;
		selector2 <<= 1;
		//five_prim_mask <<= 1;
		//three_prim_mask <<= 1;
	}

	return mins;
}

//inline size_t 
//GenomeScan::Fixhole_score_selective_insert_var_mask_rev(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
//{
//	//get the score with LPrefix + Rbuf + Rsuffix
//	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);
//
//	//cout << "prefix score: "<<s << endl; 
//
//	size_t five_prim_mask = LAST_THIRD_FOUTH;
//	size_t three_prim_mask = LAST_TWO_BIT;
//
//	size_t mins = left_mismatch + 1;
//	prim = 0;
//
//	//generate the bits where 1 indicates the mismatches.
//	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
//		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;
//
//	rbits = bits;
//
//	if (s <= left_mismatch && FixHoleCheckFirstTimeRev(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
//		return mins;
//
//	five_prim_mask <<= 1;
//	three_prim_mask <<= 1;
//
//	//the following are two pointers each point to the bit to be turned on and turned off
//	                                              // - m_max_mismatch
//	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
//	                                             // + m_max_mismatch
//	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off
//	
//	//check out the accumulation score of mismatches when the selectors are moving...
//	//size_t score_buff_width = m_//m_masks.score_buff_width;
//	                                     // + m_max_mismatch
//	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
//	{
//		s += (selector1 & bits ) ? 1 : 0;
//		s -= (selector2 & bits ) ? 1 : 0;    
//
//		if (s <= left_mismatch)
//		{
//			if (!prim)
//			{
//				if (FixHoleCheckBeforeMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
//					return mins;
//			}
//			else if (prim < 5)
//			{
//				if (FixHoleCheckAfterMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
//					return mins;
//			}
//			else if (s < mins)
//			{
//				if (FixHoleCheckAfterGTAGMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
//					return mins;
//			}
//
//		}
//		selector1 <<= 1;
//		selector2 <<= 1;
//		five_prim_mask <<= 1;
//		three_prim_mask <<= 1;
//	}
//
//	return mins;
//}


inline size_t
GenomeScan::Fixhole_score_selective_ext(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, m_masks.prefix_seg_bits_on_ext);

	//cout << "prefix score: "<<s << endl; 

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) | 
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & m_masks.comb_seg_bits_on_ext;

	rbits = bits;

	if (s <= left_mismatch && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, m_masks.score_seg_buf_width_ext))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = m_masks.comb_seg_first_selector_rt_ext;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
	                                             // + m_max_mismatch
	size_t selector2 = m_masks.comb_seg_first_selector_lt_ext;//to be turned off
	
	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= m_masks.score_seg_buf_width_ext/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;    

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width_ext))
					return mins;
			}
			else if (prim != 5)
			{
				if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width_ext))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, m_masks.score_seg_buf_width_ext))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}

inline bool
GenomeScan::FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = score_buf_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);
	
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		//combine_upper = ((five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.upper & three_prim_mask)) >> i;
		//combine_lower = ((five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask)) >> i;
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i; 
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			if (s < mins)
			{
				mins = s;
				loc = score_buf_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = score_buf_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_ATAC:
			{
				if ( s < mins)
				{
					prim = 1;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				if ( s < mins)
				{
					prim = 6;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				if ( s < mins)
				{
					prim = 3;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				if ( s < mins)
				{
					prim = 4;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				if ( s < mins)
				{
					prim = 2;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;

		switch (combine_words)
		{
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		}
	}
	return false;
}




inline bool
GenomeScan::FixHoleCheckFirstTimeRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = score_buf_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);
	
		switch (combine_words)
		{
		case bit_CATA:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				prim = 6;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CGTC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GACG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_TATG:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::FixHoleCheckBeforeMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		//combine_upper = ((five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.upper & three_prim_mask)) >> i;
		//combine_lower = ((five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask)) >> i;
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i; 
		switch (combine_words)
		{
		case bit_CATA:
			{
				prim = 1;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CGTC:
			{
				prim = 3;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GACG:
			{
				prim = 4;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_TATG:
			{
				prim = 2;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;		
		default:
			if (s < mins)
			{
				mins = s;
				loc = score_buf_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = score_buf_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_CATA:
			{
				if ( s < mins)
				{
					prim = 1;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CATC:
			{
				if ( s < mins)
				{
					prim = 6;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				if (mins == 0)
					return true;
			}
			break;
		case bit_CGTC:
			{
				if ( s < mins)
				{
					prim = 3;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GACG:
			{
				if ( s < mins)
				{
					prim = 4;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_TATG:
			{
				if ( s < mins)
				{
					prim = 2;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GATG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;		
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterGTAGMatchRev(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4 
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;

		switch (combine_words)
		{
		case bit_CATC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		case bit_GATG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		}
	}
	return false;
}


//inline size_t
//GenomeScan::Fixhole_score_selective_rev(const WordPair & read_word_dup, size_t& loc, size_t& prim)
//{
//	//get the score with LPrefix + Rbuf + Rsuffix
//	size_t s = read_word_dup.score(m_shifted_wp, m_masks.score_first_half_mask);
//
//	size_t five_prim_mask = LAST_THIRD_FOUTH;
//	size_t three_prim_mask = LAST_TWO_BIT;
//
//	size_t mins = m_left_mismatches + 1;
//	prim = 0;
//
//	if (s <= m_left_mismatches && CheckFirstTimeRev(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s))
//		return mins;
//	/*{
//		//cost 6 bit operations to get combined prims
//		//combine_upper = (five_prim_suffix.upper & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask);
//		//combine_lower = (five_prim_suffix.lower & five_prim_mask) | (three_prim_prefix.lower & three_prim_mask);
//		if (CheckFirstTimeRev(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s))
//			return mins;
//	}*/
//
//	five_prim_mask <<= 1;
//	three_prim_mask <<= 1;
//
//	//bool debug = false;
//	//if (debug)
//	//{
//	//	cerr << "mask_midright" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_midright) << endl;
//	//	cerr << "mask_ps" << endl << bits2string2(ALL_BITS_ON, m_masks.mask_ps) << endl;
//	//	cerr << "suffix_wp" << endl << m_shifted_wp.tostring(ALL_BITS_ON) << endl;
//	//	cerr << "prefix_wp" << endl << read_word_dup.tostring(ALL_BITS_ON) << endl;
//	//}
//
//	//generate the bits where 1 indicates the mismatches.
//	register size_t bits = ((m_shifted_wp.upper ^ read_word_dup.upper) | 
//		(m_shifted_wp.lower ^ read_word_dup.lower) | m_shifted_wp.bads | read_word_dup.bads) & m_masks.score_big_buff_mask;
//	//if (debug) cerr << "bits " << endl << bits2string2(ALL_BITS_ON, bits) << endl;
//
//	//the following are two pointers each point to the bit to be turned on and turned off
//	                                              // - m_max_mismatch
//	size_t selector1 = m_masks.score_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on 
//	                                             // + m_max_mismatch
//	size_t selector2 = m_masks.score_first_selector_lt;//to be turned off
//	
//	//if (debug) cerr << "score :" << mins << endl;
//
//	//check out the accumulation score of mismatches when the selectors are moving...
//	size_t score_buff_width = m_masks.score_buff_width;
//	                                     // + m_max_mismatch
//	for (size_t i = 1; i <= score_buff_width; ++i)
//	{
//		//if (debug )
//		//{
//		//	cerr << "loop " << i << "score " << s << endl;
//		//	cerr << "selector 1" << endl << bits2string2(ALL_BITS_ON, selector1) << endl;
//		//	cerr << "selector 2" << endl << bits2string2(ALL_BITS_ON, selector2) << endl;
//		//}
//
//		s += (selector1 & bits ) ? 1 : 0;
//		s -= (selector2 & bits ) ? 1 : 0;    
//
//		if (s <= m_left_mismatches)
//		{
//			if (!prim)
//			{
//				if (CheckBeforeMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
//					return mins;
//			}
//			else if (prim != 5)
//			{
//				if (CheckAfterMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
//					return mins;
//			}
//			else if (s < mins)
//			{
//				if (CheckAfterGTAGMatchRev(five_prim_mask, three_prim_mask, i, prim, mins, loc, s))
//					return mins;
//			}
//
//		}
//		selector1 <<= 1;
//		selector2 <<= 1;
//		five_prim_mask <<= 1;
//		three_prim_mask <<= 1;
//	}
//	//if (debug) cerr << "score :" << mins <<  "best at " << loc << endl;
//
//	return mins;
//}


inline void 
GenomeScan::FlankString(string& flankstr)
{
	if (m_matched_bads)
		flankstr = "BADS";
	else
	{
		for (size_t i = 7; i >= 4; --i)
		{
			if (m_matched_flank & (LEAST_SIG_BIT << i))
			{
				if (m_matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "T";
				else
					flankstr += "G";
			}
			else
			{
				if (m_matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "C";
				else
					flankstr += "A";
			}
		}
	}
}

inline string
GenomeScan::FlankString(size_t matched_flank, bool bad)
{
	string flankstr = "";
	if (bad)
		flankstr = "BADS";
	else
	{
		for (size_t i = 7; i >= 4; --i)
		{
			if (matched_flank & (LEAST_SIG_BIT << i))
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "T";
				else
					flankstr += "G";
			}
			else
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "C";
				else
					flankstr += "A";
			}
		}
	}

	return flankstr;
}

//void
//GenomeScan::WriteMapReads(size_t prefix_width, WordPair prefix_match, const size_t& prim, const size_t& lste_end, const size_t& rste_st, const size_t& chrom_ID, const size_t& readID, unsigned strand)
//{
//	size_t suffix_width = m_read_width - prefix_width;
//	size_t suffix_match_mask = ALL_BITS_ON >> (SIXTY_FOUR - suffix_width);
//
//	//WordPair prefix_match, suffix_match;
//	//                           // + m_max_mismatch
//	//m_shifted_wp.get_prefix((m_read_width + m_max_mismatches -  m_seed_width) * 2 - prefix_width, prefix_match);
//	//m_shifted_wp.get_suffix(suffix_match_mask, suffix_match);
//
//	string suffixstr = m_wp.tostring3(ALL_BITS_ON, suffix_width);
//
//	prefix_match.right_shift(suffix_width + m_masks.comp_left_shift_width);
//	string prefixstr = prefix_match.tostring3(ALL_BITS_ON, prefix_width);
//
//	string flankstr = "";
//
//	FlankString(flankstr);
//
//	WordPair matched_read;
//	string strandstr;
//	if (strand)
//	{
//		strandstr = "+";
//		matched_read = m_read_words[readID];
//	}
//	else
//	{
//		strandstr = "-";
//		matched_read = m_read_words_rc[readID];
//	}
//
//	string matched_read_str = matched_read.tostring3(ALL_BITS_ON, m_read_width);
//
//	string dirction = "F";
//
//	if (m_reverse_flag)
//	{
//		reverse(suffixstr.begin(), suffixstr.end());
//		reverse(prefixstr.begin(), prefixstr.end());
//		prefixstr.swap(suffixstr);
//
//		if (!m_matched_bads)
//			reverse(flankstr.begin(), flankstr.end());
//		reverse(matched_read_str.begin(), matched_read_str.end());
//		dirction = "R";
//		prefix_width = m_read_width - prefix_width;
//	}
//
//	m_ofs_map_reads<<m_chrom_names[chrom_ID]<<'\t'<<m_read_names[readID].c_str()<<'\t'<<prim<<'\t'<<flankstr<<'\t'<<strandstr<<'\t'<<matched_read_str<<'\t'
//		<<prefixstr<<'\t'<<suffixstr<<'\t'<<prefix_width<<'\t'<<lste_end<<'\t'<<rste_st<<'\t'<<m_islands<<'\t'<<dirction<<'\t'<<m_score<<endl;
//}

void
GenomeScan::WriteFixHmer(size_t prefix_width, WordPair& read_word, const size_t& prim, const size_t& lste_st, const size_t& lste_end, const size_t& rste_st, 
						 const size_t& chrom_ID, const size_t& readID, unsigned strand, vector<GroupedBwts*>::iterator cur_bwts_iter, const size_t score, vector<Differ>& differs_vec)
{
	//support SAM format
	size_t suffix_width = m_seg_len - prefix_width;
	//size_t suffix_match_mask = ALL_BITS_ON >> (SIXTY_FOUR - suffix_width);

	//string suffixstr = m_wp.tostring3(ALL_BITS_ON, suffix_width);

	//prefix_match.right_shift(m_masks.left_shift_seg_width + suffix_width);
	//string prefixstr = prefix_match.tostring3(ALL_BITS_ON, prefix_width);

	string missed_read_seq = read_word.tostring3(ALL_BITS_ON, m_seg_len);

	string flankstr = "";

	FlankString(flankstr);

	size_t strand_t = 16;
	if (strand)
		strand_t = 0;

	if (m_reverse_flag)
	{
		//reverse(suffixstr.begin(), suffixstr.end());
		//reverse(prefixstr.begin(), prefixstr.end());
		//prefixstr.swap(suffixstr);

		if (!m_matched_bads)
			reverse(flankstr.begin(), flankstr.end());

		suffix_width = prefix_width;
		prefix_width = m_seg_len - prefix_width;

		reverse(missed_read_seq.begin(), missed_read_seq.end());
	}

	string suffixstr = missed_read_seq.substr(prefix_width, suffix_width);

	string prefixstr = missed_read_seq.substr(0, prefix_width);

	size_t score_sum = score;

	string Ns(score * 5 + 1, 'N');

	vector<Differ> differs_all_vec;

	if ((*cur_bwts_iter)->fixtype == mistailorhead_type)
	{
		if (m_reverse_flag)
		{
			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				score_sum += (*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() / 5;
			}

			//if (score_sum > m_max_mismatches)
			//	return;

			differs_all_vec = differs_vec;

			vector<Differ>::iterator vD_iter;

			size_t seg_id = 1;

			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				for (vD_iter = (*cur_bwts_iter)->grouped_bwt_recs[i].differs_vec.begin(); vD_iter != (*cur_bwts_iter)->grouped_bwt_recs[i].differs_vec.end(); ++vD_iter)
					differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

				++seg_id;
			}

			double qual_score;

			const string& qual_seq = GetQualitySeq((*cur_bwts_iter)->tagidx, strand_t, differs_all_vec, qual_score);

			if (score_sum > m_max_totalmismatch)
				return;

			m_ofs_allmapreads<<(*cur_bwts_iter)->tagname<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<lste_st <<'\t'<<size_t(qual_score) <<'\t'<<prefix_width << 'M' <<rste_st - lste_end - 1
				<<'N'<< suffix_width + ((*cur_bwts_iter)->grouped_bwt_recs.size() * m_seg_len)<<'M'<<"\t*\t0\t0\t"<<prefixstr<<suffixstr;

			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				m_ofs_allmapreads<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq;
			}

			m_ofs_allmapreads<<'\t'<< qual_seq <<'\t'<< "NM:i:"<<score_sum<<'\t';// <<'\t'<< <<'\t'<< <<'\t'<< <<

			for (vD_iter = differs_all_vec.begin(); vD_iter != differs_all_vec.end(); ++vD_iter)
			{
				m_ofs_allmapreads << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_all_vec.end() - 1)
					m_ofs_allmapreads << ',';
			}
			m_ofs_allmapreads << endl;
		}
		else
		{
			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				score_sum += int((*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() / 5);
			}

			if (score_sum > m_max_totalmismatch)
				return;

			vector<Differ>::iterator vD_iter;

			size_t seg_id = 0;

			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				for (vD_iter = (*cur_bwts_iter)->grouped_bwt_recs[i].differs_vec.begin(); vD_iter != (*cur_bwts_iter)->grouped_bwt_recs[i].differs_vec.end(); ++vD_iter)
					differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

				++seg_id;
			}

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
				differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

			++seg_id;

			double qual_score;

			const string& qual_seq = GetQualitySeq((*cur_bwts_iter)->tagidx, strand_t, differs_all_vec, qual_score);

			m_ofs_allmapreads<<(*cur_bwts_iter)->tagname<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->offset + 1<<'\t'<<size_t(qual_score) <<'\t'
				<<(*cur_bwts_iter)->grouped_bwt_recs.size() * m_seg_len + prefix_width<< 'M' <<rste_st - lste_end - 1<<'N'<< suffix_width<<'M'<<"\t*\t0\t0\t";

			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				m_ofs_allmapreads<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq;

				//if ((*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() > 1 && (*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() <= 6)
				//	score_sum += 1;
				//else if ((*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() > 7  && (*cur_bwts_iter)->grouped_bwt_recs[i].differs.length() <= 13)
				//	score_sum += 2;
			}

			m_ofs_allmapreads <<prefixstr<<suffixstr;

			m_ofs_allmapreads<<'\t'<< qual_seq <<'\t'<< "NM:i:"<<score_sum<<'\t';

			//differs_all_vec.a = differs_vec;

			

			for (vD_iter = differs_all_vec.begin(); vD_iter != differs_all_vec.end(); ++vD_iter)
			{
				m_ofs_allmapreads << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_all_vec.end() - 1)
					m_ofs_allmapreads << ',';

			}
			m_ofs_allmapreads << endl;
		}

		string& tagnamestr = (*cur_bwts_iter)->tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;
		++m_cur_reads_count[tagidx];
	}
	else if ((*cur_bwts_iter)->fixtype == fixhmer_and_hole_type)
	{
		m_ofs_fix_hole << (*cur_bwts_iter)->grouped_bwt_recs.size() + 1<<endl;
		//need to be more general
		if (m_reverse_flag)
		{
			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.begin()->seg_id 
				<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->mapped_chrom_seq  <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->differs<<endl;

			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<lste_st - 1 <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.back().seg_id - 1 
				<<'\t'<<prefixstr<<'\t'<<score <<'\t';//<<Ns/*flankstr*/<<endl;

			vector<Differ>::iterator vD_iter;

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
			{
				m_ofs_fix_hole << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_vec.end() - 1)
					m_ofs_fix_hole << ',';

			}
			m_ofs_fix_hole << endl;

			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<rste_st - 1 <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.back().seg_id 
				<<'\t'<<suffixstr<<(*cur_bwts_iter)->grouped_bwt_recs.back().mapped_chrom_seq <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().differs<<endl;
		}
		else
		{
			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.begin()->seg_id 
				<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->mapped_chrom_seq <<prefixstr <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.begin()->differs<<endl;

			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<rste_st - 1 <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.begin()->seg_id + 1 
				<<'\t'<<suffixstr<<'\t'<<score <<'\t';//<<Ns/*flankstr*/<<endl;

			vector<Differ>::iterator vD_iter;

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
			{
				m_ofs_fix_hole << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_vec.end() - 1)
					m_ofs_fix_hole << ',';

			}
			m_ofs_fix_hole << endl;

			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.back().seg_id 
				<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().mapped_chrom_seq <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs.back().differs<<endl;
		}
	}
	else if ((*cur_bwts_iter)->fixtype == fixhole_and_hmer_type)
	{
		m_ofs_fix_hole << (*cur_bwts_iter)->grouped_bwt_recs.size() + 1<<endl;
		if (m_reverse_flag)
		{
			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<lste_st - 1 <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.front().seg_id - 1 
				<<'\t'<<prefixstr<<'\t'<<score <<'\t';//<<Ns/*flankstr*/<<endl;

			vector<Differ>::iterator vD_iter;

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
			{
				m_ofs_fix_hole << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_all_vec.end() - 1)
					m_ofs_fix_hole << ',';
			}
			m_ofs_fix_hole << endl;

			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				if (i == 0)
				{
					m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<rste_st - 1<<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs[i].seg_id 
						<<'\t'<<suffixstr<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].differs<<endl;
				}
				else
				{
					m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs[i].seg_id 
						<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].differs<<endl;
				}
			}
		}
		else
		{
			for (size_t i = 0; i < (*cur_bwts_iter)->grouped_bwt_recs.size(); ++i)
			{
				if (i == (*cur_bwts_iter)->grouped_bwt_recs.size() - 1)
				{
					m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs[i].seg_id 
						<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq<<prefixstr  <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].differs<<endl;
				}
				else
				{
					m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].offset <<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs[i].seg_id 
						<<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].mapped_chrom_seq  <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].other_mapped <<'\t'<<(*cur_bwts_iter)->grouped_bwt_recs[i].differs<<endl;
				}
			}

			m_ofs_fix_hole<<(*cur_bwts_iter)->tagname<<'\t'<<m_chrom_names[chrom_ID] <<'\t'<<(*cur_bwts_iter)->strand <<'\t'<<rste_st - 1<<'\t'<< (*cur_bwts_iter)->grouped_bwt_recs.back().seg_id + 1 
				<<'\t'<<suffixstr<<'\t'<<score <<'\t';//<<Ns/*flankstr*/<<endl;

			vector<Differ>::iterator vD_iter;

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
			{
				m_ofs_fix_hole << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_all_vec.end() - 1)
					m_ofs_fix_hole << ',';
			}
			m_ofs_fix_hole << endl;
		}
	}
	else if ((*cur_bwts_iter)->fixtype == mistailandhead_type)
	{
		//cout << prefixstr<<'\t'<< suffixstr<<'\t'<< lste_st<<'\t'<< lste_end<<'\t'<< rste_st<< '\t'<<score<< endl;
		if (m_reverse_flag)
		{
			(*cur_bwts_iter)->prefix_side->push_back(SideMatch(prefixstr, suffixstr, lste_st, lste_end, rste_st, score, differs_vec));
			
		}
		else
		{
			(*cur_bwts_iter)->suffix_side->push_back(SideMatch(prefixstr, suffixstr, lste_st, lste_end, rste_st, score, differs_vec));
		}
	}

	//if (m_reverse_flag)
	//{

		//m_ofs_map_reads<<m_read_names[readID].c_str()<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<lste_st<<'\t'<<0<<'\t'<<prefix_width
		//<<'M'<<int(rste_st-lste_end)<<'N';
		//if (splitway.find("N") != string::npos)
		//{
		//	int firstM, secondM, intronN;
		//	sscanf(splitway.c_str(), "%dM%dN%dM", &firstM, &intronN, &secondM);
		//	m_ofs_map_reads<<suffix_width + firstM<<"M"<<intronN<<"N"<<secondM <<"M"<<"\t*\t0\t0\t"<<prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		//}
		//else
		//{
		//	int firstM;
		//	sscanf(splitway.c_str(), "%dM", &firstM);
		//	m_ofs_map_reads<<suffix_width + firstM<<"M"<<"\t*\t0\t0\t"<<prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		//}
		
	//}
	//else
	//{
		//m_ofs_map_reads<<m_read_names[readID].c_str()<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<lste_st<<'\t'<<0<<'\t';
		//if (splitway.find("N") != string::npos)
		//{
		//	int firstM, secondM, intronN;
		//	sscanf(splitway.c_str(), "%dM%dN%dM", &firstM, &intronN, &secondM);
		//	m_ofs_map_reads<<firstM<<"M"<<intronN<<"N"<<secondM + prefix_width <<"M"<<int(rste_st-lste_end)<<'N'<<suffix_width<<"M\t*\t0\t0\t"<<
		//		prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		//}
		//else
		//{
		//	int firstM;
		//	sscanf(splitway.c_str(), "%dM", &firstM);
		//	m_ofs_map_reads<<firstM + prefix_width<<"M"<<int(rste_st-lste_end)<<'N'<<suffix_width<<"M\t*\t0\t0\t"<<
		//		prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl; 
		//}
		
	//}
}



void
GenomeScan::WriteFixedHeadAndTail(const char* file_fixedHeadAndTail, size_t flankcase)
{

	ofstream fixedHeadAndTail_ofs;
	fixedHeadAndTail_ofs.open(file_fixedHeadAndTail);
	map<string, vector<GroupedBwts*> >::iterator chrom_bwt_iter;

	size_t count = 0;
	size_t count_fixed = 0;

	cout << "write fixed head and tail"<<endl;

	for (chrom_bwt_iter = m_fixheadandtail_ptr.begin(); chrom_bwt_iter != m_fixheadandtail_ptr.end(); ++chrom_bwt_iter)
	{
		vector<GroupedBwts*>::iterator bwt_iter;
		
		for (bwt_iter = chrom_bwt_iter->second.begin(); bwt_iter != chrom_bwt_iter->second.end(); ++bwt_iter)
		{
			vector<SideMatch>::iterator prefix_sm_iter, suffix_sm_iter;

			size_t strand_t = 0;

			if ((*bwt_iter)->strand == '-')
				strand_t = 16;

			count++;

			//if ((*bwt_iter)->prefix_side->size() || (*bwt_iter)->suffix_side->size())
			//{
			//	cout << "prefix size: "<<(*bwt_iter)->prefix_side->size()<<"\tsuffix size: "<<(*bwt_iter)->suffix_side->size()<<endl;

			//	cout << endl;
			//}

			if ((*bwt_iter)->prefix_side->size() && (*bwt_iter)->suffix_side->size())
			{
				//getchar();
			}

			for (prefix_sm_iter = (*bwt_iter)->prefix_side->begin(); prefix_sm_iter != (*bwt_iter)->prefix_side->end(); ++prefix_sm_iter)
			{
				for (suffix_sm_iter = (*bwt_iter)->suffix_side->begin(); suffix_sm_iter != (*bwt_iter)->suffix_side->end(); ++suffix_sm_iter)
				{
					size_t score_sum = prefix_sm_iter->score;

					for (size_t i = 0; i < (*bwt_iter)->grouped_bwt_recs.size(); ++i)
					{
						score_sum += (*bwt_iter)->grouped_bwt_recs[i].differs.length() / 5;
					}

					score_sum += suffix_sm_iter->score;

					if (score_sum > m_max_totalmismatch)
						continue;

					vector<Differ> differs_all_vec;

					size_t seg_id = 0;

					vector<Differ>::iterator vD_iter;

					for (vD_iter = prefix_sm_iter->differs_vec.begin(); vD_iter != prefix_sm_iter->differs_vec.end(); ++vD_iter)
						differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

					++seg_id;

					for (size_t i = 0; i < (*bwt_iter)->grouped_bwt_recs.size(); ++i)
					{
						for (vD_iter = (*bwt_iter)->grouped_bwt_recs[i].differs_vec.begin(); vD_iter != (*bwt_iter)->grouped_bwt_recs[i].differs_vec.end(); ++vD_iter)
							differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

						++seg_id;
					}

					for (vD_iter = suffix_sm_iter->differs_vec.begin(); vD_iter != suffix_sm_iter->differs_vec.end(); ++vD_iter)
						differs_all_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));

					double qual_score;

					const string& qual_seq = GetQualitySeq((*bwt_iter)->tagidx, strand_t, differs_all_vec, qual_score);

					fixedHeadAndTail_ofs<<(*bwt_iter)->tagname<<'\t'<<strand_t<<'\t'<<chrom_bwt_iter->first<<'\t'<<prefix_sm_iter->ls_st <<'\t'<<size_t(qual_score) <<'\t'<<prefix_sm_iter->prefix_str.length() 
						<< 'M' <<prefix_sm_iter->rs_st - prefix_sm_iter->ls_end - 1<<'N'<< prefix_sm_iter->suffix_str.length() + (*bwt_iter)->grouped_bwt_recs.size() * m_seg_len + suffix_sm_iter->prefix_str.length() 
						<<'M'<<suffix_sm_iter->rs_st - suffix_sm_iter->ls_end - 1<<'N'<<suffix_sm_iter->suffix_str.length()<<'M'<<"\t*\t0\t0\t"<<prefix_sm_iter->prefix_str<<prefix_sm_iter->suffix_str;

					for (size_t i = 0; i < (*bwt_iter)->grouped_bwt_recs.size(); ++i)
					{
						fixedHeadAndTail_ofs<<(*bwt_iter)->grouped_bwt_recs[i].mapped_chrom_seq;
					}

					fixedHeadAndTail_ofs <<suffix_sm_iter->prefix_str<<suffix_sm_iter->suffix_str;

					fixedHeadAndTail_ofs<<'\t'<<qual_seq <<'\t' << "NM:i:"<<score_sum<<'\t';

					for (vD_iter = differs_all_vec.begin(); vD_iter != differs_all_vec.end(); ++vD_iter)
					{
						fixedHeadAndTail_ofs << vD_iter->index <<':'<<vD_iter->src_char <<'>'<<vD_iter->dst_char;

						if (vD_iter != differs_all_vec.end() - 1)
							fixedHeadAndTail_ofs << ',';
					}
					fixedHeadAndTail_ofs << endl;

					count_fixed++;

					string& tagnamestr = (*bwt_iter)->tagname;
					size_t last_idx = tagnamestr.find_last_of("_");
					string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
					size_t tagidx = atoi(tagidxstr.c_str()) - 1;
					++m_cur_reads_count[tagidx];
				}
			}

			if (flankcase)
			{
				(*bwt_iter)->prefix_side->clear();
				(*bwt_iter)->suffix_side->clear();
			}
			else
			{
				delete (*bwt_iter)->prefix_side;
				delete (*bwt_iter)->suffix_side;
			}
		}
	}

	cout <<"need fix head tail:"<<count<<endl;
	cout <<"fixed head tail:"<<count_fixed<<endl;
}

void
GenomeScan::WriteMapReads(size_t prefix_width, WordPair prefix_match, const size_t& prim, const size_t& lste_st, const size_t& lste_end, const size_t& rste_st, const size_t& chrom_ID, const size_t& readID, unsigned strand, const string& splitway)
{
	//support SAM format
	size_t suffix_width = m_read_width - prefix_width;
	size_t suffix_match_mask = ALL_BITS_ON >> (SIXTY_FOUR - suffix_width);

	//WordPair prefix_match, suffix_match;
	//                           // + m_max_mismatch
	//m_shifted_wp.get_prefix((m_read_width + m_max_mismatches -  m_seed_width) * 2 - prefix_width, prefix_match);
	//m_shifted_wp.get_suffix(suffix_match_mask, suffix_match);

	string suffixstr = m_wp.tostring3(ALL_BITS_ON, suffix_width);

	prefix_match.right_shift(/*suffix_width + */m_masks.comp_left_shift_width);
	string prefixstr = prefix_match.tostring3(ALL_BITS_ON, prefix_width);

	string flankstr = "";

	FlankString(flankstr);

	//WordPair matched_read;
	//string strandstr;
	size_t strand_t = 16;
	if (strand)
		strand_t = 0;
	//{
		
		//strandstr = "+";
		//matched_read = m_read_words[readID];
	//}
	//else
	//{
	//	strand_t = 16;
	//	//strandstr = "-";
	//	//matched_read = m_read_words_rc[readID];
	//}

	//string matched_read_str = matched_read.tostring3(ALL_BITS_ON, m_read_width);

	//string dirction = "F";

	if (m_reverse_flag)
	{
		reverse(suffixstr.begin(), suffixstr.end());
		reverse(prefixstr.begin(), prefixstr.end());
		prefixstr.swap(suffixstr);

		if (!m_matched_bads)
			reverse(flankstr.begin(), flankstr.end());
		//reverse(matched_read_str.begin(), matched_read_str.end());
		//dirction = "R";
		suffix_width = prefix_width;
		prefix_width = m_read_width - prefix_width;
	}

	if (m_reverse_flag)
	{
		m_ofs_map_reads<<m_read_names[readID].c_str()<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<lste_st<<'\t'<<0<<'\t'<<prefix_width
		<<'M'<<int(rste_st-lste_end)<<'N';
		if (splitway.find("N") != string::npos)
		{
			int firstM, secondM, intronN;
			sscanf(splitway.c_str(), "%dM%dN%dM", &firstM, &intronN, &secondM);
			m_ofs_map_reads<<suffix_width + firstM<<"M"<<intronN<<"N"<<secondM <<"M"<<"\t*\t0\t0\t"<<prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		}
		else
		{
			int firstM;
			sscanf(splitway.c_str(), "%dM", &firstM);
			m_ofs_map_reads<<suffix_width + firstM<<"M"<<"\t*\t0\t0\t"<<prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		}
		
	}
	else
	{
		m_ofs_map_reads<<m_read_names[readID].c_str()<<'\t'<<strand_t<<'\t'<<m_chrom_names[chrom_ID]<<'\t'<<lste_st<<'\t'<<0<<'\t';
		if (splitway.find("N") != string::npos)
		{
			int firstM, secondM, intronN;
			sscanf(splitway.c_str(), "%dM%dN%dM", &firstM, &intronN, &secondM);
			m_ofs_map_reads<<firstM<<"M"<<intronN<<"N"<<secondM + prefix_width <<"M"<<int(rste_st-lste_end)<<'N'<<suffix_width<<"M\t*\t0\t0\t"<<
				prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl;
		}
		else
		{
			int firstM;
			sscanf(splitway.c_str(), "%dM", &firstM);
			m_ofs_map_reads<<firstM + prefix_width<<"M"<<int(rste_st-lste_end)<<'N'<<suffix_width<<"M\t*\t0\t0\t"<<
				prefixstr<<suffixstr<<'\t'<<m_qual_str<<'\t'<<m_score<<'\t'<<prim<<'\t'<<flankstr<<endl; 
		}
		
	}
}

void
GenomeScan::TriggerMatch(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	if ((m_chrom_offset + 2 - cur_prefix_offset) < m_interal)
		return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TagSeed* cur_TagSeed = &(m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	size_t readID = cur_TagSeed->readID;

	unsigned read_strand = cur_TagSeed->strand;

	(read_strand)? 
		m_read_words[readID].duplicate(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, m_masks.comp_left_shift_width, m_masks.comp_right_shift_width, pre_wp) :
	m_read_words_rc[readID].duplicate(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, m_masks.comp_left_shift_width,m_masks.comp_right_shift_width, pre_wp);

	size_t kmer_score = 0;

	if (kmer_ID)
		kmer_score = pre_wp.score_hmer(m_wp, m_masks.bad_kmer_mask);

	m_left_mismatches = m_max_mismatches - m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index].score;

	if (kmer_score > m_left_mismatches)
		return;
	else
		m_left_mismatches -= kmer_score;

	++m_candidate_count;

	m_five_prim_suffix = m_slide_window[slide_window_idx].five_prim_suffix;

	//the following makes the suffix
	WordPair suffix_wp(m_wp);

	////get the right buff_size + right suffix
	//m_wp.combine(m_wp_reserve, m_read_width - m_seed_width, suffix_wp);

	// - m_max_mismatch

	//suffix_wp.right_shift(m_anchor_width);

	m_three_prim_prefix = suffix_wp;

	m_three_prim_prefix.right_shift(m_anchor_width);

	//combine prefix and suffix for scoring  				
	m_slide_window[slide_window_idx].prefix_wp.ps_combine(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, suffix_wp, m_shifted_wp);

	size_t max_loc, prim;
	m_matched_flank = 0;
	m_matched_bads = 0;

	size_t rbits;

	//get match score
	size_t score = (m_reverse_flag) ? Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	if (m_flank_case && !prim)
		return;

	score += m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index].score + kmer_score;

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;

	if (score <= m_max_mismatches)
	{
		size_t lste_st  = cur_prefix_offset - m_read_width + 1;
		size_t lste_end = cur_prefix_offset - m_read_width + (m_read_width + 1) / 2 + max_loc;
		// - m_max_mismatch
		size_t rste_st  = m_chrom_offset + 1 - m_anchor_width - (m_masks.score_buff_width - max_loc) + 1;
		size_t rste_end = m_chrom_offset + 1;

		if (m_islands_provide)
		{
			size_t prev_island_st= m_slide_window[slide_window_idx].cur_tac_pos->first;

			size_t prev_island_ed = m_slide_window[slide_window_idx].cur_tac_pos->second;

			if (lste_end < prev_island_st + 1 ||  rste_st < m_cur_island_start + 1)
				return;

			if (lste_end > prev_island_ed + 1) 
			{
				
				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window[slide_window_idx].cur_tac_pos + 1;
				
				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first - prev_island_ed <= m_boundary && 
					lste_end >= prefix_next_tac_pos->first + 1 && lste_end <= prefix_next_tac_pos->second + 1))
				{
					return;
				}
			}

			if (rste_st > m_cur_island_end + 1)
			{
				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first - m_cur_island_end <= m_boundary && 
					rste_st >= suffix_next_tac_pos->first + 1 && rste_st <= suffix_next_tac_pos->second + 1))
				{
					return;
				}
			}
		}

		size_t prefix_width = max_loc + (m_read_width + 1) / 2;

		m_score = score;

		if (m_reverse_flag)
		{
			size_t temp_lste_st = m_chrom_sizes[chrom_ID] - lste_st + 1;
			size_t temp_lste_end = m_chrom_sizes[chrom_ID] - lste_end + 1;
			size_t temp_rste_st = m_chrom_sizes[chrom_ID] - rste_st + 1;
			size_t temp_rste_end = m_chrom_sizes[chrom_ID] - rste_end + 1;
			lste_end = temp_rste_st;
			lste_st = temp_rste_end;
			rste_st = temp_lste_end;
			rste_end = temp_lste_st;
		}

		//   Also  for conjuction info	
		//JUNC_HASH_ITER iter_conj = m_conj_hashs[chrom_ID].find(rste_st);

		//if (iter_conj == m_conj_hashs[chrom_ID].end())
		//{
		//	JUNC_SUFFIX match_hash_temp;
		//	m_conj_hashs[chrom_ID].insert(JUNC_HASH::value_type(rste_st, match_hash_temp));
		//	iter_conj = m_conj_hashs[chrom_ID].find(rste_st);
		//}

		//JUNC_SUFFIX_ITER iter_conj_pre;

		//size_t real_prefix_width = m_reverse_flag ? (m_read_width- prefix_width) : prefix_width;

		//if ((iter_conj_pre = (iter_conj->second).find(lste_end))!= (iter_conj->second).end())
		//{
		//	if (score <= iter_conj_pre -> second.score)
		//	{
		//		iter_conj_pre->second.set_score(score);
		//		iter_conj_pre->second.inc_hits();
		//		iter_conj_pre->second.mark(real_prefix_width);

		clock_t t1=clock();

		//WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_end, rste_st, chrom_ID, readID, read_strand);

		clock_t t2=clock();

		m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//	}			
		//}	  
		//else 
		//{
			//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
			//(iter_conj->second).insert(JUNC_SUFFIX::value_type(lste_end, ConjInfo(score, m_matched_flank, m_matched_bads > 0, m_reverse_flag, prim, chrom_ID, /*lste_end, rste_st,*/ 1, real_prefix_width, m_read_width)));

			//clock_t t1=clock();

			//WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_end, rste_st, chrom_ID, readID, read_strand);

			//clock_t t2=clock();

			//m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//}

		//const vector<ReadInfo>::iterator current = m_bests.begin() + cur_seed->readID;
		//if (score <= current->score)
		//	current->set(score, chrom_ID, lste_st, lste_end, rste_st, rste_end, cur_seed->strand);
	}
}


void
GenomeScan::TriggerMatchFixHmerReady(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	//if (!m_is_fix_tail_head)
	//	cout << "TriggerMatchFixHmerReady"<<endl;

	WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	//if ((m_chrom_offset - m_seg_len + m_anchor_width + 1 - cur_prefix_offset) < m_interal)
	//	return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TobeFixedHmerReady* cur_TagSeed = &(m_slide_window_ready[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	//if (!m_is_fix_tail_head)
	//	cout << "mapped_TagSeed.find("<<endl;


	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	char read_strand = cur_TagSeed->strand;

	++m_candidate_count;
	//if (!m_is_fix_tail_head)
	//	cout << "Masks* cur_masks_ptr;"<<endl;

	Masks* cur_masks_ptr;

	//if (!m_is_fix_tail_head)
	//	cout << cur_TagSeed->to_be_fixed_hmer_segments_wp.size()<<endl;

	size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

	//if (!m_is_fix_tail_head)
	//	cout << "m_hash_map_masks.find "<<tobe_fixed_len <<endl;

	if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
		m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

	//if (!m_is_fix_tail_head)
	//	cout << "cur_hash_map_mask_iter"<<endl;

	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

	cur_masks_ptr = &(cur_hash_map_mask_iter->second);

	cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp.duplicate_self(tobe_fixed_len, pre_wp);

	//string temp = pre_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	size_t left_mismatches = m_max_mismatches;
	
	m_five_prim_suffix = m_slide_window_ready[slide_window_idx].five_prim_suffix;

	m_five_prim_suffix.left_shift(2);

	//string temp2 = m_five_prim_suffix.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//the following makes the suffix
	WordPair suffix_wp(m_wp_reserve);

	suffix_wp.right_shift(2 + m_seg_len - m_anchor_width);

	m_three_prim_prefix = suffix_wp;

	//string temp3 = suffix_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	WordPair comb_chrom_seq;

	//m_three_prim_prefix.right_shift(m_anchor_width);

	//if (!m_is_fix_tail_head)
	//	cout << "m_slide_window_ready[slide_window_idx].prefix_wp"<<endl;
	//combine prefix and suffix for scoring
	m_slide_window_ready[slide_window_idx].prefix_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, suffix_wp, comb_chrom_seq);

	//string temp4 = comb_chrom_seq.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	size_t max_loc, prim;
	m_matched_flank = 0;
	m_matched_bads = 0;

	size_t rbits;

	//if (!m_is_fix_tail_head)
	//	cout << "score"<<endl;

	//get match score
	size_t score = /*(m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr) :*/
									  Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr);
									  //Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	size_t matched_flank = m_matched_flank, matched_bads = m_matched_bads;

	//test
	//string fl = FlankString(matched_flank, matched_bads > 0); cout << fl << endl;

	WordPair pre_wp2;//(m_wp);

	WordPair comb_chrom_seq2;

	size_t max_loc2, prim2, rbits2, score2;

	size_t tobe_fixed_len2;

	bool second_seg = false;

	Masks* cur_masks_ptr2;

	string matched_flank_str2;

	//if (!m_is_fix_tail_head)
	//	cout << "ur_TagSeed->end_offset != 0"<<endl;

	if (cur_TagSeed->end_offset != 0 && abs(int(cur_TagSeed->end_seg_id - cur_TagSeed->st_seg_id)) == 1 && cur_TagSeed->to_be_fixed_hmer_segments_wp.size() > 1 && score <= m_max_mismatches && (prim >= flankcase || tobe_fixed_len == max_loc))
	{
		//if (!m_is_fix_tail_head)
		//	cout << "matched_flank_str2"<<endl;
		second_seg = true;
		
		tobe_fixed_len2 = cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_len;

		if (m_hash_map_masks.find(tobe_fixed_len2) == m_hash_map_masks.end())
			m_hash_map_masks.insert(make_pair(tobe_fixed_len2, Masks(tobe_fixed_len2)));

		hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter2 = m_hash_map_masks.find(tobe_fixed_len2);

		cur_masks_ptr2 = &(cur_hash_map_mask_iter2->second);

		cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.duplicate_self(tobe_fixed_len2, pre_wp2);

		//string temp5 = pre_wp2.tostring3(ALL_BITS_ON, SIXTY_FOUR);

		size_t left_mismatches2 = m_max_mismatches;

		m_five_prim_suffix = m_wp;

		m_five_prim_suffix.right_shift(m_seg_len - m_anchor_width);

		WordPair prefix_wp2(m_five_prim_suffix);

		prefix_wp2.left_shift(tobe_fixed_len2 - 2);

		m_five_prim_suffix.left_shift(2);

		//the following makes the suffix
		string suffix2_str = m_chrom_seq.substr(cur_TagSeed->end_offset - tobe_fixed_len2 - 1, tobe_fixed_len2 + 2);

		WordPair suffix_wp2(suffix2_str);

		m_three_prim_prefix = suffix_wp2;

		//suffix_wp2.right_shift(2);

		//m_three_prim_prefix.right_shift(m_anchor_width);

		//combine prefix and suffix for scoring
		prefix_wp2.ps_combine(cur_masks_ptr2->prefix_seg_bits_on, cur_masks_ptr2->suffix_seg_bits_on, cur_masks_ptr2->comb_seg_bits_on, suffix_wp2, comb_chrom_seq2);

		//string temp6 = comb_chrom_seq2.tostring3(ALL_BITS_ON, SIXTY_FOUR);
		
		m_matched_flank = 0;
		m_matched_bads = 0;

		//get match score
		score2 = /*(m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2) :*/
									      Fixhole_score_selective_var_mask(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2);
		//Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

		size_t matched_flank2 = m_matched_flank, matched_bads2 = m_matched_bads;

		matched_flank_str2 = FlankString(matched_flank2, matched_bads2 > 0);
	}

	//if ((flankcase == 0 && prim != 0) || (flankcase != 0 && prim == 0))
	//{
	//	/*if (flankcase == 0 && prim != 0)
	//		cout << "not flank 0"<<endl;*/
	//	return;
	//}

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;
	//if (!m_is_fix_tail_head)
	//	cout << "score <= m_max_mismatches"<<endl;

	if (score <= m_max_mismatches && (prim >= flankcase || (second_seg && tobe_fixed_len == max_loc)))
	{
		size_t lste_st  = cur_TagSeed->st_offset;//cur_prefix_offset - m_seg_len + 1;//cur_TagSeed->pre_offset;
		size_t lste_end = cur_TagSeed->st_offset + max_loc - 1;//cur_prefix_offset - m_seg_len/* + (m_read_width + 1) / 2*/ + max_loc;
		size_t rste_st  = m_chrom_offset - tobe_fixed_len + m_anchor_width - tobe_fixed_len - tobe_fixed_len - 1 + max_loc;
		size_t rste_end = m_chrom_offset - tobe_fixed_len + m_anchor_width - tobe_fixed_len - 2;

		if (m_islands_provide)
		{
			size_t prev_island_st = m_slide_window_ready[slide_window_idx].cur_tac_pos->first;

			size_t prev_island_ed = m_slide_window_ready[slide_window_idx].cur_tac_pos->second;

			if (lste_end < prev_island_st || rste_st < m_cur_island_start)
				return;

			if (lste_end > prev_island_ed) 
			{
				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window_ready[slide_window_idx].cur_tac_pos + 1;
				
				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first <= m_seg_len + prev_island_ed && 
					lste_end >= prefix_next_tac_pos->first - 1 && lste_end <= prefix_next_tac_pos->second + 1))
				{
					return;
				}
			}

			if (rste_st > m_cur_island_end)
			{
				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first <= m_seg_len + m_cur_island_end && 
					rste_st >= suffix_next_tac_pos->first - 1 && rste_st <= suffix_next_tac_pos->second + 1))
				{
					return;
				}
			}
		}

		string matched_flank_str1 = FlankString(matched_flank, matched_bads > 0);

		size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

		//if (m_reverse_flag)
		//{
		//	size_t temp_lste_st = m_cur_chrom_size - lste_st + 1;
		//	size_t temp_lste_end = m_cur_chrom_size - lste_end + 1;
		//	size_t temp_rste_st = m_cur_chrom_size - rste_st + 1;
		//	size_t temp_rste_end = m_cur_chrom_size - rste_end + 1;
		//	lste_end = temp_rste_st;
		//	lste_st = temp_rste_end;
		//	rste_st = temp_lste_end;
		//	rste_end = temp_lste_st;
		//}

		//HWI-EAS217:3:100:1000:1310#0/2	0	chr5	167827530	2-4	23	18M986N31M	NGAAGTTCAGTCTTTCAGGGAGAAGATGGCA	NM:i:1
		size_t strand = cur_TagSeed->strand == '+' ? 0 : 16;

		//cout << matched_flank_str1 << '\t' << matched_flank_str2 << endl;

		if (second_seg)
		{ 
			//if (!m_is_fix_tail_head)
			//	cout << "second_seg"<<endl;
			if (score2 <= m_max_mismatches)
			{
				size_t middleM = tobe_fixed_len - max_loc + max_loc2;

				size_t firstN = rste_st - lste_end - 1;

				size_t secondN = cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2);

				if (firstN > m_max_splice_dist || secondN > m_max_splice_dist)
					return;

				if (secondN == 0 && firstN == 0)
				{
					m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
						cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
						tobe_fixed_len + tobe_fixed_len2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
						cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
				}
				else if (middleM == 0)
				{
					string combined_flankstr = matched_flank_str1.substr(0, 2) + matched_flank_str2.substr(2, 2);

					if (combined_flankstr == "GTAG")
						prim = 5;
					else if (combined_flankstr == "CTAC")
						prim = 6;
					else if (combined_flankstr == "CTGC")
						prim = 3;
					else if (combined_flankstr == "GCAG")
						prim = 4;
					else if (combined_flankstr == "ATAC")
						prim = 1;
					else if (combined_flankstr == "GTAT")
						prim = 2;
					else
						prim = 0;

					if (prim < flankcase)
						return;

					if ((rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)) <= m_max_splice_dist)
					{
						m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
							cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
							max_loc<<'M'<<rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
							'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
							cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
					}

					//cout << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
					//	cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
					//	max_loc<<'M'<<rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
					//	'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
					//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
				}
				else if (secondN == 0)
				{
					if (prim < flankcase)
						return;
					m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
						cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
						max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << tobe_fixed_len - max_loc + tobe_fixed_len2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
						cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
				}
				else if (firstN == 0)
				{
					if (prim2 < flankcase)
						return;
					m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
						cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
						tobe_fixed_len + max_loc2 << 'M' << cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
						'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
						cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';

					//cout << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
					//	cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
					//	tobe_fixed_len + max_loc2 << 'M' << cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
					//	'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
					//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
				}
				else
				{
					if (prim < flankcase || prim2 < flankcase)
						return;
					m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
						cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
						max_loc<<'M'<<rste_st - lste_end - 1 <<'N' << tobe_fixed_len - max_loc + max_loc2 <<'M'<<cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
						'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
						cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2 << '\t';

					//cout << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
					//	cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
					//	max_loc<<'M'<<rste_st - lste_end - 1 <<'N' << tobe_fixed_len - max_loc + max_loc2 <<'M'<<cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
					//	'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
					//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2 << '\t';
				}

				if (score + score2)
				{
					string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);
					string comb_chrom_str2 = comb_chrom_seq2.tostring3(ALL_BITS_ON, tobe_fixed_len2 * 2);

					string read_str1 = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);
					string read_str2 = pre_wp2.tostring3(ALL_BITS_ON, tobe_fixed_len2);

					size_t seg1_suffix_len = tobe_fixed_len - max_loc;

					size_t seg2_suffix_len = tobe_fixed_len2 - max_loc2;

					string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len) 
						+ comb_chrom_str2.substr(0, max_loc2) + comb_chrom_str2.substr(comb_chrom_str2.length() - seg2_suffix_len, seg2_suffix_len);

					string combined_read_str = read_str1 + read_str2;

					size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

					size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

					size_t seg2_mask_prefix = cur_masks_ptr2->suffix_seg_bits_on >> seg2_suffix_len << seg2_suffix_len;

					size_t seg2_mask_suffix = cur_masks_ptr2->suffix_seg_bits_on >> max_loc2;

					size_t comb_bits = ((((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix)) << tobe_fixed_len) + 
						(((rbits2 >> tobe_fixed_len2) & seg2_mask_prefix) + (rbits2 & seg2_mask_suffix));

					size_t merged_len = tobe_fixed_len + tobe_fixed_len2;

					size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

					size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					if (cur_TagSeed->strand == '-')
						pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;
					else
						pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					for (size_t i = 0; i < merged_len; ++i)
					{
						if ((selector_bit >> i) & comb_bits)
						{
							m_ofs_fix_hmer << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
						}
					}
				}

				m_ofs_fix_hmer /*<<"\t"<< matched_flank_str1<<'\t' <<matched_flank_str2*/<< endl;
			}
		}
		else
		{
			//if (m_reverse_flag)
			//{
				//m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
				//	cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
				//	tobe_fixed_len - max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << max_loc << "M\t" << 
				//	cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<endl;
			//}
			//else
			//{
			//if (flankcase && !prim )
			//	return;

			if (rste_st - lste_end - 1)
				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
				lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
				max_loc<<'M'<<rste_st - lste_end - 1<< 'N' <<tobe_fixed_len - max_loc << "M\t" << 
				cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';
			else
				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
				lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<tobe_fixed_len<<"M\t"<<
				cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';

			//cout << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
			//	lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
			//	max_loc<<'M'<<rste_st - lste_end - 1<< 'N' <<tobe_fixed_len - max_loc << "M\t" << 
			//	cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';

			if (score)
			{
				string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

				string combined_read_str = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

				size_t seg1_suffix_len = tobe_fixed_len - max_loc;

				string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

				size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

				size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

				size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

				size_t merged_len = tobe_fixed_len;

				size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

				size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

				if (cur_TagSeed->strand == '-')
					pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;
				else
					pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

				for (size_t i = 0; i < merged_len; ++i)
				{
					if ((selector_bit >> i) & comb_bits)
					{
						m_ofs_fix_hmer << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
					}
				}
			}

			m_ofs_fix_hmer <</*"\t"<< matched_flank_str1<<*/ endl;
			//}
		}
	}
}


void
GenomeScan::TriggerMatchFixHmerReady_Rev(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	//if ((m_chrom_offset - m_seg_len + m_anchor_width + 1 - cur_prefix_offset) < m_interal)
	//	return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TobeFixedHmerReady* cur_TagSeed = &(m_slide_window_ready[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	char read_strand = cur_TagSeed->strand;

	++m_candidate_count;

	Masks* cur_masks_ptr;

	size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

	if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
		m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

	cur_masks_ptr = &(cur_hash_map_mask_iter->second);

	cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp.duplicate_self(tobe_fixed_len, pre_wp);

	//string temp = pre_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	size_t left_mismatches = m_max_mismatches;
	
	m_five_prim_suffix = m_slide_window_ready[slide_window_idx].five_prim_suffix;

	m_five_prim_suffix.left_shift(2);

	//string temp2 = m_five_prim_suffix.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//the following makes the suffix
	WordPair suffix_wp(m_wp_reserve);
	
	suffix_wp.right_shift(2 + m_seg_len - m_anchor_width);

	m_three_prim_prefix = suffix_wp;

	//string temp3 = suffix_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	WordPair comb_chrom_seq;

	//m_three_prim_prefix.right_shift(m_anchor_width);

	//combine prefix and suffix for scoring
	m_slide_window_ready[slide_window_idx].prefix_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, suffix_wp, comb_chrom_seq);

	//string temp4 = comb_chrom_seq.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	size_t max_loc, prim;
	m_matched_flank = 0;
	m_matched_bads = 0;

	size_t rbits;

	//get match score
	size_t score = /*(m_reverse_flag) ?*/ Fixhole_score_selective_var_mask_rev(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr)/* :
									  Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr)*/;
									  //Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	size_t matched_flank = m_matched_flank, matched_bads = m_matched_bads;

	//string matched_flank_str1 = FlankString(matched_flank, matched_bads > 0);

	//cout << matched_flank_str1<< endl;

	//WordPair pre_wp2(m_wp);

	//size_t max_loc2, prim2, rbits2, score2;

	//size_t tobe_fixed_len2;

	//bool second_seg = false;

	//if (cur_TagSeed->end_offset != 0 && abs(int(cur_TagSeed->end_seg_id - cur_TagSeed->st_seg_id)) == 1 && cur_TagSeed->to_be_fixed_hmer_segments_wp.size() > 1 && score <= m_max_mismatches)
	//{
	//	second_seg = true;

	//	Masks* cur_masks_ptr2;

	//	tobe_fixed_len2 = cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_len;

	//	if (m_hash_map_masks.find(tobe_fixed_len2) == m_hash_map_masks.end())
	//		m_hash_map_masks.insert(make_pair(tobe_fixed_len2, Masks(tobe_fixed_len2)));

	//	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter2 = m_hash_map_masks.find(tobe_fixed_len2);

	//	cur_masks_ptr2 = &(cur_hash_map_mask_iter2->second);

	//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.duplicate_self(tobe_fixed_len2, pre_wp2);

	//	string temp5 = pre_wp2.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//	size_t left_mismatches2 = m_max_mismatches;

	//	m_five_prim_suffix = m_wp;

	//	WordPair prefix_wp2(m_five_prim_suffix);

	//	prefix_wp2.left_shift(tobe_fixed_len2 - 2);

	//	//the following makes the suffix
	//	string suffix2_str = m_chrom_seq.substr(cur_TagSeed->end_offset - tobe_fixed_len2 + 1, tobe_fixed_len2 + 2);

	//	WordPair suffix_wp2(suffix2_str);

	//	m_three_prim_prefix = suffix_wp2;

	//	suffix_wp2.right_shift(2);

	//	WordPair comb_chrom_seq2;		

	//	//m_three_prim_prefix.right_shift(m_anchor_width);

	//	//combine prefix and suffix for scoring
	//	prefix_wp2.ps_combine(cur_masks_ptr2->prefix_seg_bits_on, cur_masks_ptr2->suffix_seg_bits_on, cur_masks_ptr2->comb_seg_bits_on, suffix_wp2, comb_chrom_seq2);

	//	string temp6 = comb_chrom_seq2.tostring3(ALL_BITS_ON, SIXTY_FOUR);
	//	
	//	m_matched_flank = 0;
	//	m_matched_bads = 0;

	//	//get match score
	//	score2 = (m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2) :
	//								      Fixhole_score_selective_var_mask(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2);
	//	//Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	//	size_t matched_flank2 = m_matched_flank, matched_bads2 = m_matched_bads;
	//}

	//if ((flankcase == 0 && prim != 0) || (flankcase != 0 && prim == 0))
	//{
	//	/*if (flankcase == 0 && prim != 0)
	//		cout << "not flank 0"<<endl;*/
	//	return;
	//}

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;

	if (score <= m_max_mismatches && prim >= flankcase)
	{
		//if (flankcase && !prim)
		//	return;

		string matched_flank_str1 = FlankString(matched_flank, matched_bads > 0);

		reverse(matched_flank_str1.begin(), matched_flank_str1.end());

		size_t lste_st  = m_cur_chrom_size - cur_TagSeed->st_offset - 1;//cur_prefix_offset - m_seg_len + 1;//cur_TagSeed->pre_offset;
		size_t lste_end = lste_st + max_loc - 1;//cur_prefix_offset - m_seg_len/* + (m_read_width + 1) / 2*/ + max_loc;
		size_t rste_st  = m_chrom_offset - tobe_fixed_len + m_anchor_width - tobe_fixed_len - tobe_fixed_len - 1 + max_loc;
		size_t rste_end = m_chrom_offset - tobe_fixed_len + m_anchor_width - tobe_fixed_len - 2;

		if (m_islands_provide)
		{
			size_t prev_island_st = m_slide_window_ready[slide_window_idx].cur_tac_pos->first;

			size_t prev_island_ed = m_slide_window_ready[slide_window_idx].cur_tac_pos->second;

			if (lste_end < prev_island_st ||  rste_st < m_cur_island_start)
				return;

			if (lste_end > prev_island_ed) 
			{
				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window_ready[slide_window_idx].cur_tac_pos + 1;
				
				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first <= m_seg_len + prev_island_ed && 
					lste_end >= prefix_next_tac_pos->first - 1 && lste_end <= prefix_next_tac_pos->second + 1))
				{
					return;
				}
			}

			if (rste_st > m_cur_island_end)
			{
				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first <= m_seg_len + m_cur_island_end && 
					rste_st >= suffix_next_tac_pos->first - 1 && rste_st <= suffix_next_tac_pos->second + 1))
				{
					return;
				}
			}
		}

		size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

		//if (m_reverse_flag)
		//{
		size_t temp_lste_st = m_cur_chrom_size - lste_st - 1;
		size_t temp_lste_end = m_cur_chrom_size - lste_end - 1;
		size_t temp_rste_st = m_cur_chrom_size - rste_st - 1;
		size_t temp_rste_end = m_cur_chrom_size - rste_end - 1;
		lste_end = temp_rste_st;
		lste_st = temp_rste_end;
		rste_st = temp_lste_end;
		rste_end = temp_lste_st;
		//}

		//HWI-EAS217:3:100:1000:1310#0/2	0	chr5	167827530	2-4	23	18M986N31M	NGAAGTTCAGTCTTTCAGGGAGAAGATGGCA	NM:i:1
		size_t strand = cur_TagSeed->strand == '+' ? 0 : 16;

		size_t prefixlen = tobe_fixed_len - max_loc;

		size_t suffixlen = max_loc;

		size_t intronlen = rste_st - lste_end - 1;

		string fixed_str = cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

		reverse(fixed_str.begin(), fixed_str.end());

		//if (second_seg)
		//{ 
		//	if (score2 <= m_max_mismatches)
		//	{
		//		size_t middleM = tobe_fixed_len - max_loc + max_loc2;

		//		size_t firstN = rste_st - lste_end - 1;

		//		size_t secondN = cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2);

		//		if (middleM == 0)
		//		{
		//			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		//				max_loc<<'M'<<rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
		//				'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		//				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		//		}
		//		else if (secondN == 0)
		//		{
		//			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		//				max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << tobe_fixed_len - max_loc + tobe_fixed_len2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		//				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		//		}
		//		else if (firstN == 0)
		//		{
		//			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		//				tobe_fixed_len + max_loc2 << 'M' << cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
		//				'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		//				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		//		}				
		//	}
		//}
		//else
		//{
		//	if (m_reverse_flag)
		//	{

		if (intronlen)
			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
			cur_TagSeed->st_offset - tobe_fixed_len - intronlen + 1 <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
			tobe_fixed_len - max_loc<<'M'<<intronlen << 'N' << max_loc << "M\t" << 
			fixed_str<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';//endl;
		else
			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
			cur_TagSeed->st_offset - tobe_fixed_len - intronlen + 1 <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
			tobe_fixed_len<<"M\t"<<fixed_str<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';//endl;

		//cout << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//	cur_TagSeed->st_offset - tobe_fixed_len - intronlen + 1 <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//	tobe_fixed_len - max_loc<<'M'<<intronlen << 'N' << max_loc << "M\t" << 
		//	fixed_str<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';//endl;

		if (score)
		{
			string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

			string combined_read_str = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

			size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

			int merged_len = (int)tobe_fixed_len;

			size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

			size_t pre_mapped_len;
			
			if (cur_TagSeed->strand == '+')
				pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len; 
			else
				pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;

			for (int i = merged_len; i >= 0; --i)
			{
				if ((selector_bit >> i) & comb_bits)
				{
					m_ofs_fix_hmer << pre_mapped_len + merged_len - i - 1<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
				}
			}
		}

		m_ofs_fix_hmer <</*"\t"<<matched_flank_str1 <<*/ endl;
		//	}
		//	else
		//	{
		//		m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//			lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//			max_loc<<'M'<<rste_st - lste_end - 1<< 'N' <<tobe_fixed_len - max_loc << "M\t" << 
		//			cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<endl;
		//	}
		//}
	}
}



void
GenomeScan::TriggerMatchFixHmerInsertReady(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	//if (!m_is_fix_tail_head)
	//	cout << "TriggerMatchFixHmerReady"<<endl;

	size_t cur_prefix_offset = cur_seed.offset;

	//if ((m_chrom_offset + 1 - cur_prefix_offset) < m_interal)
	//	return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TobeFixedHmerReady* cur_TagSeed = &(m_slide_window_ready[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	//if (!m_is_fix_tail_head)
	//	cout << "mapped_TagSeed.find("<<endl;

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	char read_strand = cur_TagSeed->strand;

	++m_candidate_count;
	//if (!m_is_fix_tail_head)
	//	cout << "Masks* cur_masks_ptr;"<<endl;

	Masks* cur_masks_ptr;

	//if (!m_is_fix_tail_head)
	//	cout << cur_TagSeed->to_be_fixed_hmer_segments_wp.size()<<endl;

	size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

	size_t insert_len = cur_prefix_offset + m_seg_len - m_chrom_offset - m_anchor_width - 1;

	tobe_fixed_len = tobe_fixed_len - insert_len;

	//if (!m_is_fix_tail_head)
	//	cout << "m_hash_map_masks.find "<<tobe_fixed_len <<endl;

	if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
		m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

	//if (!m_is_fix_tail_head)
	//	cout << "cur_hash_map_mask_iter"<<endl;

	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

	cur_masks_ptr = &(cur_hash_map_mask_iter->second);

	//calc insert length
	
//******************

	WordPair read_wp = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp;

	WordPair shift_read_wp(read_wp);

	shift_read_wp.left_shift(tobe_fixed_len - insert_len);

	WordPair comb_read_wp;

	shift_read_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, read_wp, comb_read_wp);

	//cout << "sub 3:"<<endl;

	//string chromstr = chromseq.substr(doner_end, tobe_fixed_len);

	//cout << "sub 3 finished:"<<endl;

	WordPair dup_chrom_wp, chrom_wp(m_wp_reserve);

	string test = chrom_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	chrom_wp.left_shift(SIXTY_FOUR - (m_seg_len - m_anchor_width + 2 + tobe_fixed_len));

	chrom_wp.right_shift(SIXTY_FOUR - tobe_fixed_len);

	chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

	size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;

	size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);

	if (score <= m_max_mismatches)
	{
		int strandi = (cur_TagSeed->strand == '+' ? 0 : 16);

		m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strandi <<  '\t' << chromname << '\t' <<cur_TagSeed->st_offset<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
			max_loc<<'M'<<insert_len<< 'I' <<tobe_fixed_len - max_loc << "M\t" <<cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_len)<<
			'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<<"\tNM:i:"<<score<<'\t';

		//cout << cur_TagSeed->tagname << '\t' << strandi <<  '\t' << chromname << '\t' <<cur_TagSeed->st_offset<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//	max_loc<<'M'<<insert_len<< 'I' <<tobe_fixed_len - max_loc << "M\t" <<cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_len)<<
		//	'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<<"\tNM:i:"<<score<<'\t';
		//m_ofs_fix_hmer.flush();
		//ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/max_loc
		//	<<'M'<< insert_len<<'I'<<tobe_fixed_len - max_loc << 'M'<<'\t'<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

		if (score)
		{
			string combined_read_str = comb_read_wp.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

			string chromstr = chrom_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_read_str = combined_read_str.substr(0, max_loc) + combined_read_str.substr(combined_read_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

			size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

			size_t merged_len = tobe_fixed_len;

			size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

			size_t pre_mapped_len;

			if (cur_TagSeed->strand == '-')
				pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;
			else
				pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

			for (size_t i = 0; i < merged_len; ++i)
			{
				if ((selector_bit >> i) & comb_bits)
				{
					if (i > max_loc)
						m_ofs_fix_hmer << pre_mapped_len + i + insert_len<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
					else
						m_ofs_fix_hmer << pre_mapped_len + i<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
				}
			}
		}

		m_ofs_fix_hmer /*<<"\t"<< flankstr*//*<<'\t' <<matched_flank_str2*/<< endl;
	}

//******************

	//cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp.duplicate_self(tobe_fixed_len, pre_wp);

	//string temp = pre_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//size_t left_mismatches = m_max_mismatches;
	//
	//m_five_prim_suffix = m_slide_window_ready[slide_window_idx].five_prim_suffix;

	//m_five_prim_suffix.left_shift(2);

	//string temp2 = m_five_prim_suffix.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	////the following makes the suffix
	//WordPair suffix_wp(m_wp_reserve);

	//suffix_wp.right_shift(2);

	//m_three_prim_prefix = suffix_wp;

	////string temp3 = suffix_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//WordPair comb_chrom_seq;

	////m_three_prim_prefix.right_shift(m_anchor_width);

	////if (!m_is_fix_tail_head)
	////	cout << "m_slide_window_ready[slide_window_idx].prefix_wp"<<endl;
	////combine prefix and suffix for scoring
	//m_slide_window_ready[slide_window_idx].prefix_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, suffix_wp, comb_chrom_seq);

	////string temp4 = comb_chrom_seq.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//size_t max_loc, prim;
	//m_matched_flank = 0;
	//m_matched_bads = 0;

	//size_t rbits;

	////if (!m_is_fix_tail_head)
	////	cout << "score"<<endl;

	////get match score
	//size_t score = /*(m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr) :*/
	//								  Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr);
	//								  //Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	//size_t matched_flank = m_matched_flank, matched_bads = m_matched_bads;

	//WordPair pre_wp2(m_wp);

	//WordPair comb_chrom_seq2;

	//size_t max_loc2, prim2, rbits2, score2;

	//size_t tobe_fixed_len2;

	//bool second_seg = false;

	//Masks* cur_masks_ptr2;

	//string matched_flank_str2;

	////if (!m_is_fix_tail_head)
	////	cout << "ur_TagSeed->end_offset != 0"<<endl;

	//if (cur_TagSeed->end_offset != 0 && abs(int(cur_TagSeed->end_seg_id - cur_TagSeed->st_seg_id)) == 1 && cur_TagSeed->to_be_fixed_hmer_segments_wp.size() > 1 && score <= m_max_mismatches)
	//{
	//	//if (!m_is_fix_tail_head)
	//	//	cout << "matched_flank_str2"<<endl;
	//	second_seg = true;
	//	
	//	tobe_fixed_len2 = cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_len;

	//	if (m_hash_map_masks.find(tobe_fixed_len2) == m_hash_map_masks.end())
	//		m_hash_map_masks.insert(make_pair(tobe_fixed_len2, Masks(tobe_fixed_len2)));

	//	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter2 = m_hash_map_masks.find(tobe_fixed_len2);

	//	cur_masks_ptr2 = &(cur_hash_map_mask_iter2->second);

	//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.duplicate_self(tobe_fixed_len2, pre_wp2);

	//	//string temp5 = pre_wp2.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//	size_t left_mismatches2 = m_max_mismatches;

	//	m_five_prim_suffix = m_wp;

	//	WordPair prefix_wp2(m_five_prim_suffix);

	//	prefix_wp2.left_shift(tobe_fixed_len2 - 2);

	//	m_five_prim_suffix.left_shift(2);

	//	//the following makes the suffix
	//	string suffix2_str = m_chrom_seq.substr(cur_TagSeed->end_offset - tobe_fixed_len2 - 1, tobe_fixed_len2 + 2);

	//	WordPair suffix_wp2(suffix2_str);

	//	m_three_prim_prefix = suffix_wp2;

	//	//suffix_wp2.right_shift(2);

	//	//m_three_prim_prefix.right_shift(m_anchor_width);

	//	//combine prefix and suffix for scoring
	//	prefix_wp2.ps_combine(cur_masks_ptr2->prefix_seg_bits_on, cur_masks_ptr2->suffix_seg_bits_on, cur_masks_ptr2->comb_seg_bits_on, suffix_wp2, comb_chrom_seq2);

	//	//string temp6 = comb_chrom_seq2.tostring3(ALL_BITS_ON, SIXTY_FOUR);
	//	
	//	m_matched_flank = 0;
	//	m_matched_bads = 0;

	//	//get match score
	//	score2 = /*(m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2) :*/
	//								      Fixhole_score_selective_var_mask(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2);
	//	//Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	//	size_t matched_flank2 = m_matched_flank, matched_bads2 = m_matched_bads;

	//	matched_flank_str2 = FlankString(matched_flank2, matched_bads2 > 0);
	//}

	////if ((flankcase == 0 && prim != 0) || (flankcase != 0 && prim == 0))
	////{
	////	/*if (flankcase == 0 && prim != 0)
	////		cout << "not flank 0"<<endl;*/
	////	return;
	////}

	////if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	////	return;
	////if (!m_is_fix_tail_head)
	////	cout << "score <= m_max_mismatches"<<endl;

	//if (score <= m_max_mismatches && (prim >= flankcase || (second_seg && tobe_fixed_len == max_loc)))
	//{
	//	size_t lste_st  = cur_TagSeed->st_offset;//cur_prefix_offset - m_seg_len + 1;//cur_TagSeed->pre_offset;
	//	size_t lste_end = cur_TagSeed->st_offset + max_loc - 1;//cur_prefix_offset - m_seg_len/* + (m_read_width + 1) / 2*/ + max_loc;
	//	size_t rste_st  = m_chrom_offset - tobe_fixed_len - tobe_fixed_len - 1 + max_loc;
	//	size_t rste_end = m_chrom_offset - tobe_fixed_len - 2;

	//	if (m_islands_provide)
	//	{
	//		size_t prev_island_st = m_slide_window_ready[slide_window_idx].cur_tac_pos->first;

	//		size_t prev_island_ed = m_slide_window_ready[slide_window_idx].cur_tac_pos->second;

	//		if (lste_end < prev_island_st || rste_st < m_cur_island_start)
	//			return;

	//		if (lste_end > prev_island_ed) 
	//		{
	//			PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window_ready[slide_window_idx].cur_tac_pos + 1;
	//			
	//			if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first <= m_seg_len + prev_island_ed && 
	//				lste_end >= prefix_next_tac_pos->first - 1 && lste_end <= prefix_next_tac_pos->second + 1))
	//			{
	//				return;
	//			}
	//		}

	//		if (rste_st > m_cur_island_end)
	//		{
	//			PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

	//			if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first <= m_seg_len + m_cur_island_end && 
	//				rste_st >= suffix_next_tac_pos->first - 1 && rste_st <= suffix_next_tac_pos->second + 1))
	//			{
	//				return;
	//			}
	//		}
	//	}

	//	string matched_flank_str1 = FlankString(matched_flank, matched_bads > 0);

	//	size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

	//	//if (m_reverse_flag)
	//	//{
	//	//	size_t temp_lste_st = m_cur_chrom_size - lste_st + 1;
	//	//	size_t temp_lste_end = m_cur_chrom_size - lste_end + 1;
	//	//	size_t temp_rste_st = m_cur_chrom_size - rste_st + 1;
	//	//	size_t temp_rste_end = m_cur_chrom_size - rste_end + 1;
	//	//	lste_end = temp_rste_st;
	//	//	lste_st = temp_rste_end;
	//	//	rste_st = temp_lste_end;
	//	//	rste_end = temp_lste_st;
	//	//}

	//	//HWI-EAS217:3:100:1000:1310#0/2	0	chr5	167827530	2-4	23	18M986N31M	NGAAGTTCAGTCTTTCAGGGAGAAGATGGCA	NM:i:1
	//	size_t strand = cur_TagSeed->strand == '+' ? 0 : 16;

	//	if (second_seg)
	//	{ 
	//		//if (!m_is_fix_tail_head)
	//		//	cout << "second_seg"<<endl;
	//		if (score2 <= m_max_mismatches)
	//		{
	//			size_t middleM = tobe_fixed_len - max_loc + max_loc2;

	//			size_t firstN = rste_st - lste_end - 1;

	//			size_t secondN = cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2);

	//			if (middleM == 0)
	//			{
	//				string combined_flankstr = matched_flank_str1.substr(0, 2) + matched_flank_str2.substr(2, 2);

	//				if (combined_flankstr == "GTAG")
	//					prim = 5;
	//				else if (combined_flankstr == "CTAC")
	//					prim = 6;
	//				else if (combined_flankstr == "CTGC")
	//					prim = 3;
	//				else if (combined_flankstr == "GCAG")
	//					prim = 4;
	//				else if (combined_flankstr == "ATAC")
	//					prim = 1;
	//				else if (combined_flankstr == "GTAT")
	//					prim = 2;
	//				else
	//					prim = 0;

	//				if (prim < flankcase)
	//					return;

	//				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//					cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
	//					max_loc<<'M'<<rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
	//					'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
	//					cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
	//			}
	//			else if (secondN == 0)
	//			{
	//				if (prim < flankcase)
	//					return;
	//				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//					cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
	//					max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << tobe_fixed_len - max_loc + tobe_fixed_len2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
	//					cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
	//			}
	//			else if (firstN == 0)
	//			{
	//				if (prim2 < flankcase)
	//					return;
	//				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//					cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
	//					tobe_fixed_len + max_loc2 << 'M' << cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
	//					'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
	//					cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<< '\t';
	//			}
	//			else
	//			{
	//				if (prim < flankcase || prim2 < flankcase)
	//					return;
	//				m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//					cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
	//					max_loc<<'M'<<rste_st - lste_end - 1 <<'N' << tobe_fixed_len - max_loc + max_loc2 <<'M'<<cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
	//					'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
	//					cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2 << '\t';
	//			}

	//			if (score + score2)
	//			{
	//				string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);
	//				string comb_chrom_str2 = comb_chrom_seq2.tostring3(ALL_BITS_ON, tobe_fixed_len2 * 2);

	//				string read_str1 = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);
	//				string read_str2 = pre_wp2.tostring3(ALL_BITS_ON, tobe_fixed_len2);

	//				size_t seg1_suffix_len = tobe_fixed_len - max_loc;

	//				size_t seg2_suffix_len = tobe_fixed_len2 - max_loc2;

	//				string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len) 
	//					+ comb_chrom_str2.substr(0, max_loc2) + comb_chrom_str2.substr(comb_chrom_str2.length() - seg2_suffix_len, seg2_suffix_len);

	//				string combined_read_str = read_str1 + read_str2;

	//				size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

	//				size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

	//				size_t seg2_mask_prefix = cur_masks_ptr2->suffix_seg_bits_on >> seg2_suffix_len << seg2_suffix_len;

	//				size_t seg2_mask_suffix = cur_masks_ptr2->suffix_seg_bits_on >> max_loc2;

	//				size_t comb_bits = ((((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix)) << tobe_fixed_len) + 
	//					(((rbits2 >> tobe_fixed_len2) & seg2_mask_prefix) + (rbits2 & seg2_mask_suffix));

	//				size_t merged_len = tobe_fixed_len + tobe_fixed_len2;

	//				size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

	//				size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

	//				if (cur_TagSeed->strand == '-')
	//					pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;
	//				else
	//					pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

	//				for (size_t i = 0; i < merged_len; ++i)
	//				{
	//					if ((selector_bit >> i) & comb_bits)
	//					{
	//						m_ofs_fix_hmer << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
	//					}
	//				}
	//			}

	//			m_ofs_fix_hmer /*<<"\t"<< matched_flank_str1<<'\t' <<matched_flank_str2*/<< endl;
	//		}
	//	}
	//	else
	//	{
	//		//if (m_reverse_flag)
	//		//{
	//			//m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//			//	cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
	//			//	tobe_fixed_len - max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << max_loc << "M\t" << 
	//			//	cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<endl;
	//		//}
	//		//else
	//		//{
	//		//if (flankcase && !prim )
	//		//	return;

	//		m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
	//			lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
	//			max_loc<<'M'<<rste_st - lste_end - 1<< 'N' <<tobe_fixed_len - max_loc << "M\t" << 
	//			cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';

	//		if (score)
	//		{
	//			string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

	//			string combined_read_str = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

	//			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

	//			string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

	//			size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

	//			size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

	//			size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

	//			size_t merged_len = tobe_fixed_len;

	//			size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

	//			size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

	//			if (cur_TagSeed->strand == '-')
	//				pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;
	//			else
	//				pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

	//			for (size_t i = 0; i < merged_len; ++i)
	//			{
	//				if ((selector_bit >> i) & comb_bits)
	//				{
	//					m_ofs_fix_hmer << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
	//				}
	//			}
	//		}

	//		m_ofs_fix_hmer <</*"\t"<< matched_flank_str1<<*/ endl;
	//		//}
	//	}
	//}
}


void
GenomeScan::TriggerMatchFixHmerInsertReady_Rev(SEED_CONST_REF cur_seed, const string& chromname, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	//WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	//if ((m_chrom_offset + 1 - cur_prefix_offset) < m_interal)
	//	return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TobeFixedHmerReady* cur_TagSeed = &(m_slide_window_ready[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	char read_strand = cur_TagSeed->strand;

	++m_candidate_count;

	Masks* cur_masks_ptr;

	size_t tobe_fixed_len = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_len;

	size_t insert_len = cur_prefix_offset + m_seg_len - m_chrom_offset - m_anchor_width - 1;

	tobe_fixed_len = tobe_fixed_len - insert_len;

	if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
		m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

	cur_masks_ptr = &(cur_hash_map_mask_iter->second);

//******************

	WordPair read_wp = cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp;

	WordPair shift_read_wp(read_wp);

	shift_read_wp.left_shift(tobe_fixed_len - insert_len);

	WordPair comb_read_wp;

	shift_read_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, read_wp, comb_read_wp);

	//cur_TagSeed->to_be_fixed_hmer_segments_wp.front().seg_wp.duplicate_self(tobe_fixed_len, pre_wp);

	//string temp = pre_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//size_t 
	
	//m_five_prim_suffix = m_slide_window_ready[slide_window_idx].five_prim_suffix;

	//m_five_prim_suffix.left_shift(2);

	//string temp2 = m_five_prim_suffix.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//the following makes the suffix
	WordPair dup_chrom_wp, chrom_wp(m_wp_reserve);
	
	string test = chrom_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	chrom_wp.left_shift(SIXTY_FOUR - (m_seg_len - m_anchor_width + 2 + tobe_fixed_len));

	chrom_wp.right_shift(SIXTY_FOUR - tobe_fixed_len);

	chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

	//suffix_wp.right_shift(2);

	//m_three_prim_prefix = suffix_wp;

	////string temp3 = suffix_wp.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//WordPair comb_chrom_seq;

	////m_three_prim_prefix.right_shift(m_anchor_width);

	////combine prefix and suffix for scoring
	//m_slide_window_ready[slide_window_idx].prefix_wp.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, suffix_wp, comb_chrom_seq);

	//string temp4 = comb_chrom_seq.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	size_t max_loc, prim, left_mismatches = m_max_mismatches, rbits;

	//get match score
	size_t score = /*(m_reverse_flag) ?*/ Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatches, rbits, cur_masks_ptr)/* :
									  Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, cur_masks_ptr)*/;
									  //Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	//size_t matched_flank = m_matched_flank, matched_bads = m_matched_bads;

	

	//WordPair pre_wp2(m_wp);

	//size_t max_loc2, prim2, rbits2, score2;

	//size_t tobe_fixed_len2;

	//bool second_seg = false;

	//if (cur_TagSeed->end_offset != 0 && abs(int(cur_TagSeed->end_seg_id - cur_TagSeed->st_seg_id)) == 1 && cur_TagSeed->to_be_fixed_hmer_segments_wp.size() > 1 && score <= m_max_mismatches)
	//{
	//	second_seg = true;

	//	Masks* cur_masks_ptr2;

	//	tobe_fixed_len2 = cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_len;

	//	if (m_hash_map_masks.find(tobe_fixed_len2) == m_hash_map_masks.end())
	//		m_hash_map_masks.insert(make_pair(tobe_fixed_len2, Masks(tobe_fixed_len2)));

	//	hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter2 = m_hash_map_masks.find(tobe_fixed_len2);

	//	cur_masks_ptr2 = &(cur_hash_map_mask_iter2->second);

	//	cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.duplicate_self(tobe_fixed_len2, pre_wp2);

	//	string temp5 = pre_wp2.tostring3(ALL_BITS_ON, SIXTY_FOUR);

	//	size_t left_mismatches2 = m_max_mismatches;

	//	m_five_prim_suffix = m_wp;

	//	WordPair prefix_wp2(m_five_prim_suffix);

	//	prefix_wp2.left_shift(tobe_fixed_len2 - 2);

	//	//the following makes the suffix
	//	string suffix2_str = m_chrom_seq.substr(cur_TagSeed->end_offset - tobe_fixed_len2 + 1, tobe_fixed_len2 + 2);

	//	WordPair suffix_wp2(suffix2_str);

	//	m_three_prim_prefix = suffix_wp2;

	//	suffix_wp2.right_shift(2);

	//	WordPair comb_chrom_seq2;		

	//	//m_three_prim_prefix.right_shift(m_anchor_width);

	//	//combine prefix and suffix for scoring
	//	prefix_wp2.ps_combine(cur_masks_ptr2->prefix_seg_bits_on, cur_masks_ptr2->suffix_seg_bits_on, cur_masks_ptr2->comb_seg_bits_on, suffix_wp2, comb_chrom_seq2);

	//	string temp6 = comb_chrom_seq2.tostring3(ALL_BITS_ON, SIXTY_FOUR);
	//	
	//	m_matched_flank = 0;
	//	m_matched_bads = 0;

	//	//get match score
	//	score2 = (m_reverse_flag) ? Fixhole_score_selective_var_mask_rev(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2) :
	//								      Fixhole_score_selective_var_mask(pre_wp2, comb_chrom_seq2, max_loc2, prim2, left_mismatches2, rbits2, cur_masks_ptr2);
	//	//Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	//	size_t matched_flank2 = m_matched_flank, matched_bads2 = m_matched_bads;
	//}

	//if ((flankcase == 0 && prim != 0) || (flankcase != 0 && prim == 0))
	//{
	//	/*if (flankcase == 0 && prim != 0)
	//		cout << "not flank 0"<<endl;*/
	//	return;
	//}

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;

	if (score <= m_max_mismatches/* && prim >= flankcase*/)
	{
		//if (flankcase && !prim)
		//	return;

		//string matched_flank_str1 = FlankString(matched_flank, matched_bads > 0);

		//reverse(matched_flank_str1.begin(), matched_flank_str1.end());

		//size_t lste_st  = m_cur_chrom_size - cur_TagSeed->st_offset - 1;//cur_prefix_offset - m_seg_len + 1;//cur_TagSeed->pre_offset;
		//size_t lste_end = lste_st + max_loc - 1;//cur_prefix_offset - m_seg_len/* + (m_read_width + 1) / 2*/ + max_loc;
		//size_t rste_st  = m_chrom_offset - tobe_fixed_len - tobe_fixed_len - 1 + max_loc;
		//size_t rste_end = m_chrom_offset - tobe_fixed_len - 2;

		//if (m_islands_provide)
		//{
		//	size_t prev_island_st = m_slide_window_ready[slide_window_idx].cur_tac_pos->first;

		//	size_t prev_island_ed = m_slide_window_ready[slide_window_idx].cur_tac_pos->second;

		//	if (lste_end < prev_island_st ||  rste_st < m_cur_island_start)
		//		return;

		//	if (lste_end > prev_island_ed) 
		//	{
		//		PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window_ready[slide_window_idx].cur_tac_pos + 1;
		//		
		//		if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first <= m_seg_len + prev_island_ed && 
		//			lste_end >= prefix_next_tac_pos->first - 1 && lste_end <= prefix_next_tac_pos->second + 1))
		//		{
		//			return;
		//		}
		//	}

		//	if (rste_st > m_cur_island_end)
		//	{
		//		PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

		//		if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first <= m_seg_len + m_cur_island_end && 
		//			rste_st >= suffix_next_tac_pos->first - 1 && rste_st <= suffix_next_tac_pos->second + 1))
		//		{
		//			return;
		//		}
		//	}
		//}

		//size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

		////if (m_reverse_flag)
		////{
		//size_t temp_lste_st = m_cur_chrom_size - lste_st - 1;
		//size_t temp_lste_end = m_cur_chrom_size - lste_end - 1;
		//size_t temp_rste_st = m_cur_chrom_size - rste_st - 1;
		//size_t temp_rste_end = m_cur_chrom_size - rste_end - 1;
		//lste_end = temp_rste_st;
		//lste_st = temp_rste_end;
		//rste_st = temp_lste_end;
		//rste_end = temp_lste_st;
		////}

		////HWI-EAS217:3:100:1000:1310#0/2	0	chr5	167827530	2-4	23	18M986N31M	NGAAGTTCAGTCTTTCAGGGAGAAGATGGCA	NM:i:1
		//size_t strand = cur_TagSeed->strand == '+' ? 0 : 16;

		//size_t prefixlen = tobe_fixed_len - max_loc;

		//size_t suffixlen = max_loc;

		//size_t intronlen = rste_st - lste_end - 1;

		//string fixed_str = cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

		//reverse(fixed_str.begin(), fixed_str.end());

		////if (second_seg)
		////{ 
		////	if (score2 <= m_max_mismatches)
		////	{
		////		size_t middleM = tobe_fixed_len - max_loc + max_loc2;

		////		size_t firstN = rste_st - lste_end - 1;

		////		size_t secondN = cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2);

		////		if (middleM == 0)
		////		{
		////			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		////				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		////				max_loc<<'M'<<rste_st - lste_end - 1 + cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
		////				'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		////				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		////		}
		////		else if (secondN == 0)
		////		{
		////			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		////				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		////				max_loc<<'M'<<rste_st - lste_end - 1 << 'N' << tobe_fixed_len - max_loc + tobe_fixed_len2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		////				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		////		}
		////		else if (firstN == 0)
		////		{
		////			m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		////				cur_TagSeed->st_offset <<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len + 1<<'\t'<<
		////				tobe_fixed_len + max_loc2 << 'M' << cur_TagSeed->end_offset - (rste_end + tobe_fixed_len2)<<
		////				'N'<<tobe_fixed_len2 - max_loc2 << "M\t" << cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<
		////				cur_TagSeed->to_be_fixed_hmer_segments_wp[1].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len2)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score + score2<<endl;
		////		}				
		////	}
		////}
		////else
		////{
		////	if (m_reverse_flag)
		////	{
		//m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//	cur_TagSeed->st_offset - tobe_fixed_len + insert_len + 1 <<'\t' << cur_TagSeed->st_seg_id
		//	<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//	tobe_fixed_len - max_loc<<'M'<<intronlen << 'N' << max_loc << "M\t" << 
		//	fixed_str<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<'\t';//endl;

		int strandi = (cur_TagSeed->strand == '+' ? 0 : 16);

		string rev_read = cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_len);

		reverse(rev_read.begin(), rev_read.end());

		m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strandi <<  '\t' << chromname << '\t' 
			<<cur_TagSeed->st_offset - tobe_fixed_len + 1<<'\t' 
			<< cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
			tobe_fixed_len - max_loc <<'M'<<insert_len<< 'I' <<max_loc<< "M\t" <<rev_read
			<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<<"\tNM:i:"<<score<<'\t';

		//cout << cur_TagSeed->tagname << '\t' << strandi <<  '\t' << chromname << '\t' 
		//	<<cur_TagSeed->st_offset - tobe_fixed_len + 1<<'\t' 
		//	<< cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//	tobe_fixed_len - max_loc <<'M'<<insert_len<< 'I' <<max_loc << "M\t" <<rev_read
		//	<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<<"\tNM:i:"<<score<<'\t';

		if (score)
		{
			//string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

			//string combined_read_str = pre_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

			//size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

			//size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			//size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

			//size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

			string combined_read_str = comb_read_wp.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

			string chromstr = chrom_wp.tostring3(ALL_BITS_ON, tobe_fixed_len);

			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_read_str = combined_read_str.substr(0, max_loc) + combined_read_str.substr(combined_read_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

			size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

			int merged_len = (int)tobe_fixed_len;

			size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

			size_t pre_mapped_len;
			
			if (cur_TagSeed->strand == '+')
				pre_mapped_len = (cur_TagSeed->st_seg_id - 1) * m_seg_len; 
			else
				pre_mapped_len = (m_num_seg - cur_TagSeed->st_seg_id) * m_seg_len;

			for (int i = merged_len; i >= 0; --i)
			{
				if ((selector_bit >> i) & comb_bits)
				{
					if (i < max_loc)
						m_ofs_fix_hmer << pre_mapped_len + merged_len - i - 1 + insert_len<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
					else
						m_ofs_fix_hmer << pre_mapped_len + merged_len - i - 1<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
				}
			}
		}

		m_ofs_fix_hmer <</*"\t"<<matched_flank_str1 <<*/ endl;
		//	}
		//	else
		//	{
		//		m_ofs_fix_hmer << cur_TagSeed->tagname << '\t' << strand <<  '\t' << chromname << '\t' <<
		//			lste_st<<'\t' << cur_TagSeed->st_seg_id<<'-'<<cur_TagSeed->end_seg_id<<'\t'<<m_seg_len+1<<'\t'<<
		//			max_loc<<'M'<<rste_st - lste_end - 1<< 'N' <<tobe_fixed_len - max_loc << "M\t" << 
		//			cur_TagSeed->to_be_fixed_hmer_segments_wp[0].seg_wp.tostring3(ALL_BITS_ON, tobe_fixed_len)<<'\t'<<cur_TagSeed->id1<<'\t'<<cur_TagSeed->id2<</*endl;*/"\tNM:i:"<<score<<endl;
		//	}
		//}
	}
}



void
GenomeScan::TriggerMatchFixHmer(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed, size_t flankcase)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	if ((m_chrom_offset + 1 - cur_prefix_offset) < m_interal)
		return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TagSeed* cur_TagSeed = &(m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	size_t readID = cur_TagSeed->readID;

	vector<GroupedBwts*>::iterator bwts_ptr = cur_TagSeed->bwts_ptr;

	unsigned read_strand = cur_TagSeed->strand;

	(read_strand)? 
		m_read_words[readID].duplicate(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, pre_wp) :
	m_read_words_rc[readID].duplicate(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, pre_wp);

	WordPair read_word = (read_strand) ? m_read_words[readID] : m_read_words_rc[readID];

	size_t kmer_score = 0;

	if (kmer_ID)
		kmer_score = pre_wp.score_hmer(m_wp, m_masks.bad_kmer_mask);

	m_left_mismatches = m_max_mismatches - cur_TagSeed->score;

	if (kmer_score > m_left_mismatches)
		return;
	else
		m_left_mismatches -= kmer_score;

	++m_candidate_count;

	m_five_prim_suffix = m_slide_window[slide_window_idx].five_prim_suffix;

	//the following makes the suffix
	WordPair suffix_wp(m_wp);

	////get the right buff_size + right suffix
	//m_wp.combine(m_wp_reserve, m_read_width - m_seed_width, suffix_wp);

	// - m_max_mismatch

	//suffix_wp.right_shift(m_anchor_width);

	m_three_prim_prefix = suffix_wp;

	m_three_prim_prefix.right_shift(m_anchor_width);

	//combine prefix and suffix for scoring  				
	m_slide_window[slide_window_idx].prefix_wp.ps_combine(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, suffix_wp, m_shifted_wp);

	size_t max_loc, prim;
	m_matched_flank = 0;
	m_matched_bads = 0;

	size_t rbits;

	//get match score
	size_t score = (m_reverse_flag) ? Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	if ((flankcase == 0 && prim != 0) || (flankcase != 0 && prim == 0))
	{
		/*if (flankcase == 0 && prim != 0)
			cout << "not flank 0"<<endl;*/
		return;
	}

	score += cur_TagSeed->score + kmer_score;

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;

	if (score <= m_max_mismatches)
	{
		size_t lste_st  = cur_prefix_offset - m_seg_len + 1;//cur_TagSeed->pre_offset;
		size_t lste_end = cur_prefix_offset - m_seg_len/* + (m_read_width + 1) / 2*/ + max_loc;
		// - m_max_mismatch
		size_t rste_st  = m_chrom_offset + 1 - m_anchor_width - (m_masks.score_buff_width - max_loc) + 1;
		size_t rste_end = m_chrom_offset + 1;

		if (m_islands_provide)
		{
			size_t prev_island_st= m_slide_window[slide_window_idx].cur_tac_pos->first;

			size_t prev_island_ed = m_slide_window[slide_window_idx].cur_tac_pos->second;

			if (lste_end < prev_island_st + 1 ||  rste_st < m_cur_island_start + 1)
				return;

			if (lste_end > prev_island_ed + 1) 
			{
				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window[slide_window_idx].cur_tac_pos + 1;
				
				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first - prev_island_ed <= m_boundary && 
					lste_end >= prefix_next_tac_pos->first + 1 && lste_end <= prefix_next_tac_pos->second + 1))
				{
					return;
				}
			}

			if (rste_st > m_cur_island_end + 1)
			{
				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first - m_cur_island_end <= m_boundary && 
					rste_st >= suffix_next_tac_pos->first + 1 && rste_st <= suffix_next_tac_pos->second + 1))
				{
					return;
				}
			}
		}

		size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

		vector<Differ> differs_vec;

		size_t selector_lt = m_masks.mis_first_selector_lt;

		//cout << (2 * m_seg_len - 1) << endl << selector_lt << endl;

		for (size_t i = 0; i < max_loc; ++i)
		{
			if (selector_lt & rbits)
			{
				char dst = pre_wp.get_char(ALL_BITS_ON, 2 * m_seg_len - i);

				char src = m_shifted_wp.get_char(ALL_BITS_ON, 2 * m_seg_len - i);//prefixflankstr[i];

				if (dst == src)
				{

					cout << max_loc << endl <<i <<endl<< rbits << endl << selector_lt << endl << src << endl << dst<<endl;
					continue;
					getchar();
				}

				differs_vec.push_back(Differ(i, src, dst));
			}

			selector_lt >>= 1;
		}

		//cout << selector_lt << endl;

		selector_lt >>= m_seg_len;

		//cout << selector_lt << endl;

		for (size_t j = max_loc; j < m_seg_len; ++j)
		{
			if (selector_lt & rbits)
			{
				char dst = pre_wp.get_char(ALL_BITS_ON, m_seg_len - j);

				char src = m_shifted_wp.get_char(ALL_BITS_ON, m_seg_len - j);

				if (dst == src)
				{
					cout << max_loc << endl <<j <<endl<<rbits << endl << selector_lt << endl << src << endl << dst<<endl;
					continue;
					getchar();
				}

				differs_vec.push_back(Differ(j, src, dst));
			}

			selector_lt >>= 1;
		}

		m_score = score;

		if (m_reverse_flag)
		{
			size_t temp_lste_st = m_chrom_sizes[chrom_ID] - lste_st + 1;
			size_t temp_lste_end = m_chrom_sizes[chrom_ID] - lste_end + 1;
			size_t temp_rste_st = m_chrom_sizes[chrom_ID] - rste_st + 1;
			size_t temp_rste_end = m_chrom_sizes[chrom_ID] - rste_end + 1;
			lste_end = temp_rste_st;
			lste_st = temp_rste_end;
			rste_st = temp_lste_end;
			rste_end = temp_lste_st;

			reverse(differs_vec.begin(), differs_vec.end());

			vector<Differ>::iterator vD_iter;

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
				vD_iter->index = m_seg_len - vD_iter->index - 1;
		}

		//   Also  for conjuction info	
		//JUNC_HASH_ITER iter_conj = m_conj_hashs[chrom_ID].find(rste_st);

		//if (iter_conj == m_conj_hashs[chrom_ID].end())
		//{
		//	JUNC_SUFFIX match_hash_temp;
		//	m_conj_hashs[chrom_ID].insert(JUNC_HASH::value_type(rste_st, match_hash_temp));
		//	iter_conj = m_conj_hashs[chrom_ID].find(rste_st);
		//}

		//JUNC_SUFFIX_ITER iter_conj_pre;

		//size_t real_prefix_width = m_reverse_flag ? (m_read_width- prefix_width) : prefix_width;

		//if ((iter_conj_pre = (iter_conj->second).find(lste_end))!= (iter_conj->second).end())
		//{
		//	if (score <= iter_conj_pre -> second.score)
		//	{
		//		iter_conj_pre->second.set_score(score);
		//		iter_conj_pre->second.inc_hits();
		//		iter_conj_pre->second.mark(real_prefix_width);

		clock_t t1=clock();

		WordPair prefix_wp = m_slide_window[slide_window_idx].prefix_wp;
		WriteFixHmer(prefix_width, read_word, prim, lste_st, lste_end, rste_st, chrom_ID, readID, read_strand, bwts_ptr, m_score, differs_vec);

		clock_t t2=clock();

		m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//	}			
		//}	  
		//else 
		//{
			//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
			//(iter_conj->second).insert(JUNC_SUFFIX::value_type(lste_end, ConjInfo(score, m_matched_flank, m_matched_bads > 0, m_reverse_flag, prim, chrom_ID, /*lste_end, rste_st,*/ 1, real_prefix_width, m_read_width)));

			//clock_t t1=clock();

			//WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_end, rste_st, chrom_ID, readID, read_strand);

			//clock_t t2=clock();

			//m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//}

		//const vector<ReadInfo>::iterator current = m_bests.begin() + cur_seed->readID;
		//if (score <= current->score)
		//	current->set(score, chrom_ID, lste_st, lste_end, rste_st, rste_end, cur_seed->strand);
	}
}



void
GenomeScan::TriggerMatchFixTail(SEED_CONST_REF cur_seed, size_t chrom_ID, size_t kmer_ID, hash_map<size_t, int>& mapped_TagSeed)
{
    //prepare read for scoring, move down to avoid unnecessary operation
	WordPair pre_wp;

	size_t cur_prefix_offset = cur_seed.offset;

	if ((m_chrom_offset + 2 - cur_prefix_offset) < m_interal)
		return;

	size_t slide_window_idx = cur_prefix_offset % m_window_size;

	TagSeed* cur_TagSeed = &(m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index]);

	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
		return;
	else
		mapped_TagSeed[(size_t)cur_TagSeed] = 1;

	size_t readID = cur_TagSeed->readID;

	unsigned read_strand = cur_TagSeed->strand;

	(read_strand)? 
		m_read_words[readID].duplicate(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, pre_wp) :
	m_read_words_rc[readID].duplicate(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, pre_wp);

	size_t kmer_score = 0;

	if (kmer_ID)
		kmer_score = pre_wp.score_hmer(m_wp, m_masks.bad_kmer_mask);

	m_left_mismatches = m_max_mismatches - cur_TagSeed->score;

	if (kmer_score > m_left_mismatches)
		return;
	else
		m_left_mismatches -= kmer_score;

	++m_candidate_count;

	m_five_prim_suffix = m_slide_window[slide_window_idx].five_prim_suffix;

	//the following makes the suffix
	WordPair suffix_wp(m_wp);

	////get the right buff_size + right suffix
	//m_wp.combine(m_wp_reserve, m_read_width - m_seed_width, suffix_wp);

	// - m_max_mismatch

	//suffix_wp.right_shift(m_anchor_width);

	m_three_prim_prefix = suffix_wp;

	m_three_prim_prefix.right_shift(m_anchor_width);

	//combine prefix and suffix for scoring  				
	m_slide_window[slide_window_idx].prefix_wp.ps_combine(m_masks.ft_comb_prefix_half_mask, m_masks.ft_comb_suffix_half_mask, m_masks.ft_combined_mask, suffix_wp, m_shifted_wp);

	size_t max_loc, prim;
	m_matched_flank = 0;
	m_matched_bads = 0;

	size_t rbits;

	//get match score
	size_t score = (m_reverse_flag) ? Ps_score_selective_rev(pre_wp, max_loc, prim, rbits) : Ps_score_selective(pre_wp, max_loc, prim, rbits);

	if (m_flank_case && !prim)
		return;

	score += cur_TagSeed->score + kmer_score;

	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
	//	return;

	if (score <= m_max_mismatches)
	{
		size_t lste_st  = cur_TagSeed->pre_offset;//cur_prefix_offset - m_read_width + 1 - 
		size_t lste_end = cur_prefix_offset - m_read_width/* + (m_read_width + 1) / 2*/ + max_loc;
		// - m_max_mismatch
		size_t rste_st  = m_chrom_offset + 1 - m_anchor_width - (m_masks.score_buff_width - max_loc) + 1;
		size_t rste_end = m_chrom_offset + 1;

		if (m_islands_provide)
		{
			size_t prev_island_st= m_slide_window[slide_window_idx].cur_tac_pos->first;

			size_t prev_island_ed = m_slide_window[slide_window_idx].cur_tac_pos->second;

			if (lste_end < prev_island_st + 1 ||  rste_st < m_cur_island_start + 1)
				return;

			if (lste_end > prev_island_ed + 1) 
			{
				
				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window[slide_window_idx].cur_tac_pos + 1;
				
				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first - prev_island_ed <= m_boundary && 
					lste_end >= prefix_next_tac_pos->first + 1 && lste_end <= prefix_next_tac_pos->second + 1))
				{
					return;
				}
			}

			if (rste_st > m_cur_island_end + 1)
			{
				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;

				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first - m_cur_island_end <= m_boundary && 
					rste_st >= suffix_next_tac_pos->first + 1 && rste_st <= suffix_next_tac_pos->second + 1))
				{
					return;
				}
			}
		}

		size_t prefix_width = max_loc;//lste_end - lste_st + 1;//max_loc + (m_read_width + 1) / 2;

		m_score = score;

		if (m_reverse_flag)
		{
			size_t temp_lste_st = lste_st;//m_chrom_sizes[chrom_ID] - lste_st + 1;
			size_t temp_lste_end = m_chrom_sizes[chrom_ID] - lste_end + 1;
			size_t temp_rste_st = m_chrom_sizes[chrom_ID] - rste_st + 1;
			size_t temp_rste_end = m_chrom_sizes[chrom_ID] - rste_end + 1;
			lste_end = temp_rste_st;
			lste_st = temp_rste_end;
			rste_st = temp_lste_end;
			rste_end = temp_lste_st;
		}

		//   Also  for conjuction info	
		//JUNC_HASH_ITER iter_conj = m_conj_hashs[chrom_ID].find(rste_st);

		//if (iter_conj == m_conj_hashs[chrom_ID].end())
		//{
		//	JUNC_SUFFIX match_hash_temp;
		//	m_conj_hashs[chrom_ID].insert(JUNC_HASH::value_type(rste_st, match_hash_temp));
		//	iter_conj = m_conj_hashs[chrom_ID].find(rste_st);
		//}

		//JUNC_SUFFIX_ITER iter_conj_pre;

		//size_t real_prefix_width = m_reverse_flag ? (m_read_width- prefix_width) : prefix_width;

		//if ((iter_conj_pre = (iter_conj->second).find(lste_end))!= (iter_conj->second).end())
		//{
		//	if (score <= iter_conj_pre -> second.score)
		//	{
		//		iter_conj_pre->second.set_score(score);
		//		iter_conj_pre->second.inc_hits();
		//		iter_conj_pre->second.mark(real_prefix_width);

		clock_t t1=clock();

		WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_st, lste_end, rste_st, chrom_ID, readID, read_strand, cur_TagSeed->spliceway);

		clock_t t2=clock();

		m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//	}			
		//}	  
		//else 
		//{
			//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
			//(iter_conj->second).insert(JUNC_SUFFIX::value_type(lste_end, ConjInfo(score, m_matched_flank, m_matched_bads > 0, m_reverse_flag, prim, chrom_ID, /*lste_end, rste_st,*/ 1, real_prefix_width, m_read_width)));

			//clock_t t1=clock();

			//WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_end, rste_st, chrom_ID, readID, read_strand);

			//clock_t t2=clock();

			//m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
		//}

		//const vector<ReadInfo>::iterator current = m_bests.begin() + cur_seed->readID;
		//if (score <= current->score)
		//	current->set(score, chrom_ID, lste_st, lste_end, rste_st, rste_end, cur_seed->strand);
	}
}


void
GenomeScan::CreateMapReadsFile(size_t tissueID, size_t chromID)
{
	char allmapreads_file[1000];

	sprintf(allmapreads_file, "%smapreads_%s.txt", m_path.c_str(), m_tissue_names[tissueID].c_str());

	m_ofs_allmapreads.open(allmapreads_file);

	//char map_reads[1000];
	//
	//string prefixmatch;
	//if (m_prefix_match)
	//	prefixmatch = "Y";
	//else
	//	prefixmatch = "N";

	//sprintf(map_reads, "%smapreads_%s_%s.txt", m_path.c_str(), m_tissue_names[tissueID].c_str(), m_chrom_names[chromID].c_str());

	//m_mapreads_file = map_reads;

	//cout << map_reads<<endl;

	//if (m_reverse_flag)
	//	m_ofs_map_reads.open(map_reads, ios_base::app);
	//else
	//	m_ofs_map_reads.open(map_reads);

	cout << allmapreads_file<<endl;

	if (m_reverse_flag)
		m_ofs_allmapreads.open(allmapreads_file, ios_base::app);
	else
		m_ofs_allmapreads.open(allmapreads_file);

	if (!m_ofs_allmapreads.is_open())
	{
		return;
		getchar();
	}

	//if (!m_reverse_flag)
	//{
	//	m_mapreads_files.push_back(m_mapreads_file);

	//	m_ofs_map_reads << "chromname\tReadID\tcase\tG[end+1:end+2]+G[start1-2, start1-1]\tstrand\tReadSequence\tG[start:end]\tG[start1:end1]\tstart-end+1\tend\tstart1\tisland\tmismatch\n";
	//}

	if (!m_reverse_flag)
	{
		m_mapreads_files.push_back(allmapreads_file);

		//m_ofs_allmapreads << "chromname\tReadID\tcase\tG[end+1:end+2]+G[start1-2, start1-1]\tstrand\tReadSequence\tG[start:end]\tG[start1:end1]\tstart-end+1\tend\tstart1\tisland\tmismatch\n";
	}
}

void
GenomeScan::InitReads(size_t tissueID)
{
	cout << "get reads"<<endl;

	if (m_fastq_file.length())
		GetReadsFastq(m_tissue_files[tissueID].c_str());
	else
		GetReads(m_tissue_files[tissueID].c_str());

	//cout << "clean reads\n";
	//// remove any bad reads
	//CleanReads();

	if (VERBOSE)
		cerr << "READS AFTER QC: " << m_reads.size() << endl;
}

void
GenomeScan::InitReadsFastq(size_t tissueID)
{
	cout << "get reads\n";
	GetReads(m_tissue_files[tissueID].c_str());

	//cout << "clean reads\n";
	//// remove any bad reads
	//CleanReads();

	if (VERBOSE)
		cerr << "READS AFTER QC: " << m_reads.size() << endl;
}

void
GenomeScan::InitReadWords()
{
	clock_t t1=clock();

	// jinze: get_read_matches has been changed.
	cout << "hash seeds\n";
	//m_seed_hash.Get_read_matches(m_reads);

	cout << "set read words\n";

	SetReadWords();

	clock_t t2=clock();

	m_init_time += (t2-t1)/(double)CLOCKS_PER_SEC;

	/*cout <<"init m_best\n";

	m_bests.insert(m_bests.begin(), m_reads.size(), ReadInfo(m_max_mismatches + 1));*/

	cout << "init end\n";
}

void
GenomeScan::InitReadWords_Rev()
{
	clock_t t1=clock();

	// jinze: get_read_matches has been changed.
	cout << "hash seeds reverse\n";
	// m_seed_hash.Get_read_matches_rev(m_reads);

	cout << "set read words reverse\n";

	SetReadWords_Rev();

	clock_t t2=clock();

	m_init_time += (t2-t1)/(double)CLOCKS_PER_SEC;

	/*cout <<"init m_best\n";

	m_bests.insert(m_bests.begin(), m_reads.size(), ReadInfo(m_max_mismatches + 1));*/

	cout << "reverse init end\n";
}

bool
GenomeScan::Clear()
{
	for (size_t i = 0; i < m_open_reads_prefix.size(); ++i)
	{
		//m_open_reads[i].segmentStart();
		m_open_reads_prefix[i].segmentStart();
		m_open_reads_suffix[i].segmentStart();
	}
	//m_slide_window.clear();
	//m_slide_window.resize(m_window_size);
	m_slide_window_ready.clear();
	m_slide_window_ready.resize(m_window_size);
	//m_prefix_words.clear();
	m_seed_key = 0;
	for (size_t i = 0; i < m_kmer_keys.size(); ++i)
	{
		m_kmer_keys[i] = 0;
		m_kmer_keys_lower[i] = 0;
		m_kmer_keys_upper[i] = 0;
	}

	m_candidate_count = 0;
	m_basecount = 0;
	m_enterhalf_count = 0;
	m_passhalf_count = 0;
	m_wp.clear(); 
	m_wp_reserve.clear();
	m_shifted_wp.clear();
	m_chrom_offset = 0;
	m_file_offset = 0;
	m_bad_bases = ALL_BITS_ON;
	m_bad_bases_reserve = ALL_BITS_ON;
	return true;
}

bool
GenomeScan::Release()
{
	ClearReads();
	m_slide_window.clear();
	m_open_reads.clear();
	m_tac_pos_map.clear();
	m_tac_pos_map_rev.clear();
	return true;
}

void
GenomeScan::ClearReads()
{
	m_reads.clear();
	m_read_names.clear();
	m_seed_hash.clear();
	m_read_words.clear();
	m_read_words_rc.clear();
	m_read_kmers.clear();
	m_read_kmers_rc.clear();
	//m_bests.clear();
}

void
GenomeScan::ClearReadWords()
{
	m_seed_hash.clear();
	m_read_words.clear();
	m_read_words_rc.clear();
	m_read_kmers.clear();
	m_read_kmers_rc.clear();
	//m_bests.clear();
}

//void
//GenomeScan::FixHole()
//{
//	ifstream ifs_fixhole(m_fixhole_file.c_str());
//
//	string fixed_file = m_path + "mapreads.fixed";
//	string unmapped_file = m_path + "mapreads.unmapped";
//	string fixedtail_file = m_path + "maprads.fixtail";
//	ofstream ofs_fixed(fixed_file.c_str());
//	ofstream ofs_unmapped(unmapped_file.c_str());
//	ofstream ofs_fixtail(fixedtail_file.c_str());
//
//	if (!ifs_fixhole)
//	{
//		cout <<"error: unable to open file: " << m_fixhole_file<<endl;
//	}
//
//	if (!ofs_fixed)
//	{
//		cout <<"error: unable to open file: " << ofs_fixed<<endl;
//	}
//
//	cout << "fixing hole"<<endl;
//	
//	string line;
//	while(getline(ifs_fixhole, line))
//	{
//		//chr1	HWI-EAS217:4:8:196:1268#0/1_34222	-	1	291064	5786	0	3	291066	6495	0	2	TGAAGGAGACGTGCTGCTGCTGCTG	TGGAGCCGGGAGGTGGGGAACAGGGCA	CCTGAAGGAGACGCGCTGCTGCTGCTG
//		char chromname[1000], tagname[1000], strand, alters1[1000], alters2[1000], missedseq[1000], prim5flank[1000], prim3flank[1000];
//		int offset1, tagid1, alignid1, offset2, tagid2, alignid2, missedtagid;
//		sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s", chromname, tagname, &strand, &tagid1, &alignid1, &offset1, alters1, &tagid2, &alignid2, &offset2, alters2, 
//			&missedtagid, missedseq, prim5flank, prim3flank);
//
//		string tagnamestr = tagname;
//		size_t last_idx = tagnamestr.find_last_of("_");
//		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
//		int tagidx = atoi(tagidxstr.c_str()) - 1;
//
//		//WordPair missed_read
//
//		m_five_prim_suffix = WordPair(prim5flank);
//		m_three_prim_prefix = WordPair(prim3flank);
//		//WordPair prim5flank_wp(prim5flank);
//
//		//WordPair prim3flank_wp(prim3flank);
//
//		WordPair dup_missed_read;
//
//		if (strand == '+')
//		{
//			m_read_words[tagidx * 4 + missedtagid - 1].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
//		}
//		else
//		{
//			m_read_words_rc[tagidx * 4 + 4 - missedtagid].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
//		}
//
//		//(strand == '+') ? 
//		//	m_read_words[tagidx * 4 + missedtagid].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read) : 
//		//m_read_words_rc[tagidx * 4 + 5 - missedtagid].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
//
//		//WordPair prim5flank_shifted = prim5flank_wp;
//
//		WordPair prefix_chrom = m_five_prim_suffix;
//
//		m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string
//
//		prefix_chrom.left_shift(m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs
//
//		WordPair comb_chrom_seq;
//
//		prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);
//
//		size_t max_loc, prim, left_mismatch = m_max_mismatches;
//		m_matched_flank = 0;
//		m_matched_bads = 0;
//
//		size_t score = Fixhole_score_selective(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);
//
//		if (score <= m_max_mismatches)
//		{
//			int strandi = 0;
//			if (strand == '-')
//				strandi = 16;
//
//			string prim5flankstr = prim5flank;
//			string prim3flankstr = prim3flank;
//			string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_read_width - max_loc), m_read_width - max_loc);
//
//			string qual_str(comb_mapped_chromseq.length(), 'I');
//
//			string flankstr;
//
//			FlankString(flankstr);
//
//			size_t missed4th = 0;
//			size_t prefixext = 0;
//			size_t suffixext = 0;
//
//			//check if 1 or 4 can append
//			if (missedtagid == 2)
//			{
//				if (m_bowtie_segments[alignid2 + 1].chromname == string(chromname) && 
//					m_bowtie_segments[alignid2 + 1].tagname == string(tagname) && 
//					m_bowtie_segments[alignid2 + 1].strand == strand && 
//					m_bowtie_segments[alignid2 + 1].segid == tagid2 + 1 &&
//					m_bowtie_segments[alignid2 + 1].offset == offset2 + m_read_width)
//					suffixext = m_read_width;
//				else
//					missed4th = 4;
//				
//			}
//			else if (missedtagid == 3)
//			{
//				
//				if (m_bowtie_segments[alignid1 - 1].chromname == string(chromname) && 
//					m_bowtie_segments[alignid1 - 1].tagname == string(tagname) && 
//					m_bowtie_segments[alignid1 - 1].strand == strand && 
//					m_bowtie_segments[alignid1 - 1].segid == tagid1 - 1 &&
//					m_bowtie_segments[alignid1 - 1].offset == offset1 - m_read_width)
//					prefixext = m_read_width;
//				else
//					missed4th = 1;
//			}
//			else
//			{
//				cout <<"error : missedtagid not 2 or 3 : "<<endl<<line<<endl;
//			}
//
//			if (!missed4th)
//				ofs_fixed<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<offset1 - prefixext<<'\t'<<0<<'\t'<<m_read_width + max_loc + prefixext<<'M'<<offset2 - m_read_width - offset1 - m_read_width  <<'N'
//				<<m_read_width + m_read_width - max_loc + suffixext<<'M'<<"\t*\t0\t0"<<'\t'<<comb_mapped_chromseq<<'\t'<< qual_str<<'\t'<< score<<'\t'<< prim<<'\t'<< flankstr<<'\t'<<missed4th<<endl;
//			else
//			{
//				size_t prefix_matched_len = m_read_width + max_loc;
//				size_t intron_size = offset2 - m_read_width - offset1 - m_read_width;
//				size_t suffix_matched_len = m_read_width + m_read_width - max_loc;
//
//				ofs_fixtail<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<offset1 - prefixext<<'\t'<<0<<'\t'<<prefix_matched_len <<'M'<<intron_size <<'N'
//					<<suffix_matched_len <<'M'<<"\t*\t0\t0"<<'\t'<<comb_mapped_chromseq<<'\t'<< qual_str<<'\t'<< score<<'\t'<< prim<<'\t'<< flankstr<<'\t'<<missed4th<<endl;
//
//				string tagnamestr = tagname;
//				size_t last_idx = tagnamestr.find_last_of("_");
//				string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
//				size_t tagidx = atoi(tagidxstr.c_str()) - 1;
//
//				size_t fixtailoffset = offset1;
//				if (missed4th != 1)
//					fixtailoffset = offset2 + m_read_width - 1;
//				else
//					fixtailoffset = m_chrom_sizes[m_chrom_names_map[chromname]] - fixtailoffset - 1;
//
//				char spliceway[1000];
//
//				sprintf(spliceway, "%dM%dN%dM", prefix_matched_len, intron_size, suffix_matched_len);
//
//				if (strandi== 0)
//					tagidx = tagidx * 4 + missed4th - 1;
//				else
//					tagidx = tagidx * 4 + 5 - missed4th - 1;
//
//				//if (missed4th == 1)
//
//				m_sam_fix_tail_map[missed4th][chromname].push_back(SamFixtail(tagname, strandi?0:1, offset1, fixtailoffset, spliceway, tagidx));
//
//				m_sam_fix_tail_map[missed4th][chromname].back().flanks.push_back(make_pair(score, make_pair(prim, flankstr)));
//
//			}
//		}
//		else
//		{
//			//int strandi = 0;
//			//if (strand == '-')
//			//	strandi = 16;
//
//			//string prim5flankstr = prim5flank;
//			//string prim3flankstr = prim3flank;
//			//string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_read_width - max_loc), m_read_width - max_loc);
//
//			//string qual_str(comb_mapped_chromseq.length(), 'I');
//
//			//string flankstr;
//
//			//FlankString(flankstr);
//
//			ofs_unmapped<<line<<endl;
//		}
//	}
//
//	ofs_fixtail.close();
//	ofs_fixed.close();
//	ofs_unmapped.close();
//	ifs_fixhole.close();
//}


void
GenomeScan::FixHole(const char* fixhole_file, const char* fixed_hole_file, size_t flankcase, bool avoid, bool mark, bool exceed)
{
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string chromseq = "";
	string prevchrom_file = "";

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		//cout <<"read number"<<endl;
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;

		//cout <<"read line"<<endl;
		for (size_t i = 0; i < rec_num; ++i)
		{
			getline(input_fs, line);

			if (line == "")
				continue;

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);
			fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		}

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		if (avoid && m_reads_count[tagidx])
			continue;

		size_t missed_id;

		vector<GroupedBwtRec>::iterator bwt_rec_iter;

		for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		{
			if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			{
				if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
				{
					cout << "missed more than 1 segments"<<endl;
					exit(1);
				}
				else
				{
					missed_id = bwt_rec_iter->seg_id + 1;
					break;
				}
			}
		}
		if (exceed)
			;
		else if (bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length()/*m_seg_len*/ + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
			continue;

		size_t readID;

		WordPair dup_missed_read;

		if (strand == '+')
		{
			readID = tagidx * m_num_seg + missed_id - 1;
			m_read_words[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		}
		else
		{
			readID = tagidx * m_num_seg + m_num_seg - missed_id;
			m_read_words_rc[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		}

		string chromfile = m_chrom_files[m_chrom_names_map[chromname]];

		//if (string(chromname) != "chr1")
		//	continue;
		
		if (prevchrom_file == "" || chromfile != prevchrom_file)
		{
			chromseq.clear();
			//cout << "read chromo\t"<<chromfile<<endl;
			readchrom(chromfile.c_str(), chromseq);
			if (chromseq == "")
			{
				cout <<"empty chrom: "<<chromfile<<endl;
				exit(1);
			}
			prevchrom_file = chromfile;
		}

		//cout << "get substring"<<endl;
		string prefixflankstr = chromseq.substr(bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length(), m_seg_len + 2);
		string suffixflankstr = chromseq.substr((bwt_rec_iter + 1)->offset - m_seg_len - 2, m_seg_len + 2);

		m_five_prim_suffix = WordPair(prefixflankstr);
		m_three_prim_prefix = WordPair(suffixflankstr);

		WordPair prefix_chrom = m_five_prim_suffix;

		m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

		prefix_chrom.left_shift(m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

		WordPair comb_chrom_seq;

		//cout <<"combine chrom"<<endl;
		prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

		size_t max_loc, prim, left_mismatch = m_max_mismatches;
		m_matched_flank = 0;
		m_matched_bads = 0;

		//cout << "scoring"<<endl;
		size_t score = Fixhole_score_selective(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);

		//ofstream* cur_ofs_ptr;

		//if (m_flank_case && !prim)
		//	continue;
		//cout <<"score: "<<score<<endl;
		//cout <<"max_loc: "<<max_loc<<endl;
		if (score <= m_max_mismatches)
		{
			if ((prim > 0 && flankcase > 0) || (prim == 0 && flankcase == 0))
			{
				if (mark)
				{
					++m_cur_reads_count[tagidx];
				}
			}
			else
				continue;

			size_t score_sum = score;

			vector<GroupedBwtRec>::iterator sp_bwt_rec_iter;

			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			{
				score_sum += sp_bwt_rec_iter->differs.length() / 5;
			}

			for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				score_sum += sp_bwt_rec_iter->differs.length() / 5;
			}

			if (score_sum > m_max_totalmismatch)
				continue;

			int strandi = 0;
			if (strand == '-')
				strandi = 16;

			string prim5flankstr = prefixflankstr;
			string prim3flankstr = suffixflankstr;
			string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

			string flankstr;

			FlankString(flankstr);

			ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t';
			//cout<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t'<<endl;

			//cout <<"prefix part"<<endl;
			size_t cur_map_len = fix_hole_vec.begin()->mapped_chrom_seq.length();

			size_t outputlen = 0;
			bool outputprefix = false;

			vector <pair<size_t, size_t> > splice_way;

			splice_way.push_back(make_pair(fix_hole_vec.begin()->offset, fix_hole_vec.begin()->mapped_chrom_seq.length()));
			for (sp_bwt_rec_iter = fix_hole_vec.begin() + 1; sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			{
				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				{
					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"before expend : "<<cur_map_len<<endl;
					//cout <<"expend prefix part"<<endl;
					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();

					//cout <<"after expend : "<<cur_map_len<<endl;
				}
				else
				{
					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
					//outputlen += cur_map_len;
					//outputprefix = true;
				}
			}

			splice_way.back().second += max_loc;

			splice_way.push_back(make_pair((bwt_rec_iter + 1)->offset - m_seg_len + max_loc, m_seg_len - max_loc));

			//cout <<"middle part"<<endl;
			//if (outputprefix)
			//{
			//	ofs<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M';
			//	cout<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M'<<endl;
			//	outputlen += sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc;
			//}
			//else
			//{
			//	ofs<<cur_map_len + max_loc<<'M';
			//	cout<<cur_map_len + max_loc<<'M'<<endl;
			//	outputlen += cur_map_len + max_loc;
			//}

			//ofs<<(bwt_rec_iter + 1)->offset - bwt_rec_iter->offset - m_seg_len - m_seg_len<<'N';

			//cur_map_len = m_seg_len - max_loc + (bwt_rec_iter + 1)->mapped_chrom_seq.length();
			
			//cout <<"suffix part"<<endl;

			splice_way.back().second += (bwt_rec_iter + 1)->mapped_chrom_seq.length();

			//bool outputsuffix = false;
			for (sp_bwt_rec_iter = bwt_rec_iter + 2; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				{
					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"before expend : "<<cur_map_len<<endl;
					//cout <<"expend suffix part"<<endl;
					//
					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"after expend : "<<cur_map_len<<endl;
				}
				else
				{
					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
					//outputlen += cur_map_len;
					//outputsuffix = true;
				}
			}

			//cout << "tail part"<<endl;

			vector <pair<size_t, size_t> >::iterator sp_iter;

			size_t lensum = 0;

			string combseq = "";

			for (sp_iter = splice_way.begin(); sp_iter != splice_way.end(); ++sp_iter)
			{
				if (sp_iter != splice_way.begin())
				{
					ofs << sp_iter->first - (sp_iter - 1)->first - (sp_iter - 1)->second<<'N';
				}
				ofs << sp_iter->second<<'M';

				combseq += chromseq.substr(sp_iter->first, sp_iter->second);

				//cout << "combine position: "<<sp_iter->first<<'\t'<< sp_iter->second<<endl;

				lensum += sp_iter->second;
			}

			if (lensum != m_read_width)
				cout << "too short len"<<endl;

			if (strandi)
				combseq = revcomp(combseq);
			for (size_t ridx = 0; ridx < m_num_seg; ++ridx)
			{
				int errsum = 0;

#ifdef DEBUG_CHECK
				for (size_t j = 0; j < m_seg_len; ++j)
				{
					if (m_reads[tagidx * m_num_seg + ridx][j] != combseq[ridx * m_seg_len + j] && 
						m_reads[tagidx * m_num_seg + ridx][j] + 'A' - 'a' != combseq[ridx * m_seg_len + j] &&
						m_reads[tagidx * m_num_seg + ridx][j] - ('A' - 'a') != combseq[ridx * m_seg_len + j])
						errsum++;
				}


				if (errsum > m_max_mismatches)
				{
					for (size_t j = 0; j < splice_way.size(); ++j)
						cout <<splice_way[j].first<<'\t'<<splice_way[j].second<<endl;
					cout <<ridx<<endl;
					cout <<m_reads[tagidx * m_num_seg + ridx]<<endl;
					cout << combseq<<endl;
					cout << "check why"<<endl;
					//getchar();
				}
#endif
			}
			//if (outputsuffix)
			//{
			//	ofs<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M';
			//	cout<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M'<<endl;
			//	outputlen += (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length();
			//	
			//}
			//else
			//{
			//	ofs<<cur_map_len<<'M';
			//	cout <<cur_map_len<<'M'<<endl;
			//	outputlen += cur_map_len;
			//}

			ofs<<"\t*\t0\t0\t";

			//cout << "mapped first sequence"<<endl;

			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			{
				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
			}

			//cout << "mapped middle sequence"<<endl;
			ofs<<comb_mapped_chromseq;

			//cout << "mapped suffix sequence"<<endl;
			for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
			}

			//cout <<"output finished"<<endl;

			ofs<<'\t'<</*GetQualitySeq(tagidx, strandi)*/m_qual_str<<'\t'<<"NM:i:"<<score_sum<<endl;
		}
	}
	input_fs.close();
	ofs.close();
	//ofs_f0.close();

}

void
GenomeScan::FixHoleExt(const char* fixhole_file, const char* fixed_hole_file, size_t flankcase, bool avoid, bool mark, bool exceed)
{
	cout << fixhole_file << endl;
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string chromseq = "";
	string prevchrom_file = "";

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		//cout <<"read number"<<endl;
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		//if (rec_num != m_num_seg)
		//{
		//	cout << "missed segments in all mapped"<<endl;
		//	continue;
		//}

		vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;

		//cout <<"read line"<<endl;
		for (size_t i = 0; i < rec_num; ++i)
		{
			getline(input_fs, line);

			if (line == "")
				continue;

			//cout << line << endl;

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);
			fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs));
		}

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		if (avoid && m_reads_count[tagidx])
			continue;

		size_t missed_id, fixid1, fixid2, offset1, offset2;

		vector<GroupedBwtRec>::iterator bwt_rec_iter;

		for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		{
			//if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			//{
			//	cout << "missed more 1 segments"<<endl;
			//	exit(1);
			//}
			//else if (bwt_rec_iter->offset + m_seg_len != (bwt_rec_iter + 1)->offset)
			//{
			//	
			//}

			if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			{
				if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
				{
					cout << "missed more than 1 segments"<<endl;
					exit(1);
				}
				else
				{
					//missed_id = bwt_rec_iter->seg_id + 1;
					//break;
					fixid1 = bwt_rec_iter->seg_id;
					fixid2 = (bwt_rec_iter + 1)->seg_id;

					offset1 = bwt_rec_iter->offset;
					offset2 = (bwt_rec_iter + 1)->offset;

					missed_id = bwt_rec_iter->seg_id + 1;
					break;
				}
			}
		}
		if (exceed)
			;
		else if (bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length()/*m_seg_len*/ + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
			continue;

		size_t readID1, readID2, readID;

		WordPair dup_missed_read, prefix_suffix_append;

		//cout << "prepare read word"<<endl;

		//cout << strand << endl;

		if (strand == '+')
		{
			readID1 = tagidx * m_num_seg + fixid1 - 1;

			readID2 = tagidx * m_num_seg + fixid2 - 1;

			readID = tagidx * m_num_seg + missed_id - 1;

			WordPair wp1(m_read_words[readID1]);

			WordPair wp2(m_read_words[readID2]);

			//cout << "wp1 wp2"<<endl;

			//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			WordPair prefixappend;

			wp1.left_shift(m_seg_len);

			wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words[readID], prefixappend);

			prefixappend.left_shift(m_extend_bits);

			wp2.right_shift(m_seg_len - m_extend_bits);

			//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//cout << "prefix_suffix_append"<<endl;

			//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			//cout << m_masks.prefix_seg_bits_on_ext << endl << m_masks.suffix_seg_bits_on_ext << endl << m_masks.comb_seg_bits_on_ext<<endl;

			prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
				m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
		}
		else
		{
			readID1 = tagidx * m_num_seg + m_num_seg - fixid1;

			readID2 = tagidx * m_num_seg + m_num_seg - fixid2;
			
			readID = tagidx * m_num_seg + m_num_seg - missed_id;

			WordPair wp1(m_read_words_rc[readID1]);

			WordPair wp2(m_read_words_rc[readID2]);

			//cout << "wp1 wp2"<<endl;

			//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			WordPair prefixappend;

			wp1.left_shift(m_seg_len);

			wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words_rc[readID], prefixappend);

			prefixappend.left_shift(m_extend_bits);

			wp2.right_shift(m_seg_len - m_extend_bits);

			//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//cout << "prefix_suffix_append"<<endl;

			//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
				m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
		}

		//cout << "dup_missed_read"<<endl;

		//cout << dup_missed_read.tostring3(ALL_BITS_ON, 2 * (2*m_max_mismatches + m_seg_len))<<endl;

		if (m_chrom_names_map.find(chromname) == m_chrom_names_map.end())
			continue;

		string chromfile = m_chrom_files[m_chrom_names_map[chromname]];

		//cout << "prepare chrom sequence"<<endl;
		
		if (prevchrom_file == "" || chromfile != prevchrom_file)
		{
			chromseq.clear();
			//cout << "read chromo\t"<<chromfile<<endl;
			readchrom(chromfile.c_str(), chromseq);
			if (chromseq == "")
			{
				cout <<"empty chrom: "<<chromfile<<endl;
				exit(1);
			}
			prevchrom_file = chromfile;
		}

		//cout << "get substring"<<endl;
		string prefixflankstr = chromseq.substr(bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length() - m_extend_bits, m_seg_len + 2 + (2 * m_extend_bits));
		string suffixflankstr = chromseq.substr((bwt_rec_iter + 1)->offset - m_seg_len - 2 - m_extend_bits, m_seg_len + 2 + (2 * m_extend_bits));

		//cout << "prefix chrom " << endl;

		//cout << prefixflankstr << endl;

		//cout << "suffix chrom " << endl;

		//cout << suffixflankstr << endl;

		m_five_prim_suffix = WordPair(prefixflankstr);
		m_three_prim_prefix = WordPair(suffixflankstr);

		WordPair prefix_chrom = m_five_prim_suffix;

		m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

		prefix_chrom.left_shift(m_masks.left_shift_seg_width_ext - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

		WordPair comb_chrom_seq;

		//cout <<"combine chrom"<<endl;
		prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, m_three_prim_prefix, comb_chrom_seq);

		//cout << comb_chrom_seq.tostring3(ALL_BITS_ON, 2 * (2*m_max_mismatches + m_seg_len))<<endl;

		size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;
		m_matched_flank = 0;
		m_matched_bads = 0;

		//cout << "scoring"<<endl;
		size_t score = Fixhole_score_selective_ext(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits);

		//cout << score <<'\t' << max_loc<< '\t' << prim<<endl;

		//getchar();

		//ofstream* cur_ofs_ptr;

		//if (m_flank_case && !prim)
		//	continue;
		//cout <<"score: "<<score<<endl;
		//cout <<"max_loc: "<<max_loc<<endl;
		if (score <= m_max_mismatches)
		{
			//string flankstrtest;

			//FlankString(flankstrtest);

			//cout << flankstrtest << endl;

			//getchar();
			if ((prim > 0 && flankcase > 0) || (prim == 0 && flankcase == 0))
			{
				if (mark)
				{
					++m_cur_reads_count[tagidx];
				}
			}
			else
				continue;

			size_t score_sum = score;

			vector<GroupedBwtRec>::iterator sp_bwt_rec_iter;

			vector<Differ> differs_vec;

			size_t seg_id = 0;

			//cout << "first part 1"<<endl;

			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter; ++sp_bwt_rec_iter)
			{
				vector<Differ>::iterator vD_iter;

				for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				{
					differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
					++score_sum;

					//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				}

				++seg_id;
			}

			vector<Differ>::iterator vD_iter;

			//cout << "first part 2"<<endl;

			for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
			{
				if (vD_iter->index < m_seg_len - m_extend_bits)
				{
					differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
					++score_sum;

					//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				}
			}

			++sp_bwt_rec_iter;

			++seg_id;

			size_t selector_lt = m_masks.mis_selector_lt;

			//cout << "middle part 1"<<endl;

			for (size_t i = 0; i < max_loc; ++i)
			{
				if (selector_lt & rbits)
				{
					char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - i);

					char src = prefixflankstr[i];

					if (src >= 'a' && src <= 'z')
						src += 'A' - 'a';

					if (dst == src)
					{
						cout << prefixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<i<<endl;
						continue;
						getchar();
					}

					differs_vec.push_back(Differ(seg_id * m_seg_len + i - m_extend_bits, src, dst));

					//cout << seg_id * m_seg_len + i - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				}

				selector_lt >>= 1;
			}

			selector_lt >>= m_masks.score_seg_buf_width_ext;

			//cout << "middle part 2"<<endl;

			for (size_t j = max_loc; j < m_masks.score_seg_buf_width_ext; ++j)
			{
				if (selector_lt & rbits)
				{
					char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - j);

					char src = suffixflankstr[j + suffixflankstr.length() - m_masks.score_seg_buf_width_ext];

					if (src >= 'a' && src <= 'z')
						src += 'A' - 'a';

					if (dst == src)
					{
						cout << suffixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<j<<endl;
						continue;
						getchar();
					}

					differs_vec.push_back(Differ(seg_id * m_seg_len + j - m_extend_bits, src, dst));

					//cout << seg_id * m_seg_len + j - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				}

				selector_lt >>= 1;
			}

			++seg_id;

			//cout << "third part 1"<<endl;

			for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
			{
				if (vD_iter->index >= m_extend_bits)
				{
					differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
					++score_sum;

					//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				}
			}

			++sp_bwt_rec_iter;

			++seg_id;

			//cout << "third part 2"<<endl;

			for (; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				vector<Differ>::iterator vD_iter;

				for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				{
					differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
					++score_sum;

					//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				}

				++seg_id;
			}

			//size_t prefix_mask = ALL_BITS_ON >> (SIXTY_FOUR - max_loc)<< (m_masks.score_seg_buf_width_ext - max_loc);

			//size_t suffix_mask = ALL_BITS_ON >> (m_masks.score_seg_buf_width_ext - max_loc);

			//size_t mbits = ((rbits >> m_masks.score_seg_buf_width_ext) & prefix_mask) & (rbits & suffix_mask);

			//for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			//{
			//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
			//}

			//for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			//{
			//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
			//}

			if (score_sum > m_max_totalmismatch)
				continue;

			int strandi = 0;
			if (strand == '-')
				strandi = 16;

			double qual_score;

			const string& qual_seq = GetQualitySeq(tagidx, strandi, differs_vec, qual_score);

			//string prim5flankstr = prefixflankstr;
			//string prim3flankstr = suffixflankstr;
			//string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

			string flankstr;

			FlankString(flankstr);

			ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<size_t(qual_score)<<'\t';
			//cout<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t'<<endl;

			//cout <<"prefix part"<<endl;
			size_t cur_map_len = fix_hole_vec.begin()->mapped_chrom_seq.length();

			size_t outputlen = 0;
			bool outputprefix = false;

			vector <pair<size_t, size_t> > splice_way;

			splice_way.push_back(make_pair(fix_hole_vec.begin()->offset, fix_hole_vec.begin()->mapped_chrom_seq.length()));
			for (sp_bwt_rec_iter = fix_hole_vec.begin() + 1; sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			{
				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				{
					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"before expend : "<<cur_map_len<<endl;
					//cout <<"expend prefix part"<<endl;
					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();

					//cout <<"after expend : "<<cur_map_len<<endl;
				}
				else
				{
					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
					//outputlen += cur_map_len;
					//outputprefix = true;
				}
			}

			splice_way.back().second += max_loc - m_extend_bits;

			splice_way.push_back(make_pair((bwt_rec_iter + 1)->offset - m_seg_len + max_loc - m_extend_bits, m_seg_len - max_loc + m_extend_bits));

			//cout <<"middle part"<<endl;
			//if (outputprefix)
			//{
			//	ofs<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M';
			//	cout<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M'<<endl;
			//	outputlen += sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc;
			//}
			//else
			//{
			//	ofs<<cur_map_len + max_loc<<'M';
			//	cout<<cur_map_len + max_loc<<'M'<<endl;
			//	outputlen += cur_map_len + max_loc;
			//}

			//ofs<<(bwt_rec_iter + 1)->offset - bwt_rec_iter->offset - m_seg_len - m_seg_len<<'N';

			//cur_map_len = m_seg_len - max_loc + (bwt_rec_iter + 1)->mapped_chrom_seq.length();
			
			//cout <<"suffix part"<<endl;

			splice_way.back().second += (bwt_rec_iter + 1)->mapped_chrom_seq.length();

			//bool outputsuffix = false;
			for (sp_bwt_rec_iter = bwt_rec_iter + 2; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				{
					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"before expend : "<<cur_map_len<<endl;
					//cout <<"expend suffix part"<<endl;
					//
					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
					//cout <<"after expend : "<<cur_map_len<<endl;
				}
				else
				{
					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
					//outputlen += cur_map_len;
					//outputsuffix = true;
				}
			}

			//cout << "tail part"<<endl;

			vector <pair<size_t, size_t> >::iterator sp_iter;

			size_t lensum = 0;

			string combseq = "";

			for (sp_iter = splice_way.begin(); sp_iter != splice_way.end(); ++sp_iter)
			{
				if (sp_iter != splice_way.begin())
				{
					ofs << sp_iter->first - (sp_iter - 1)->first - (sp_iter - 1)->second<<'N';
				}
				ofs << sp_iter->second<<'M';

				combseq += chromseq.substr(sp_iter->first, sp_iter->second);

				//cout << "combine position: "<<sp_iter->first<<'\t'<< sp_iter->second<<endl;

				lensum += sp_iter->second;
			}

			if (lensum != m_read_width)
				cout << "too short len"<<endl;

			if (strandi)
				combseq = revcomp(combseq);
			for (size_t ridx = 0; ridx < m_num_seg; ++ridx)
			{
				int errsum = 0;

#ifdef DEBUG_CHECK
				for (size_t j = 0; j < m_seg_len; ++j)
				{
					if (m_reads[tagidx * m_num_seg + ridx][j] != combseq[ridx * m_seg_len + j] && 
						m_reads[tagidx * m_num_seg + ridx][j] + 'A' - 'a' != combseq[ridx * m_seg_len + j] &&
						m_reads[tagidx * m_num_seg + ridx][j] - ('A' - 'a') != combseq[ridx * m_seg_len + j])
						errsum++;
				}

				if (errsum > m_max_mismatches)
				{
					for (size_t j = 0; j < splice_way.size(); ++j)
						cout <<splice_way[j].first<<'\t'<<splice_way[j].second<<endl;
					cout <<ridx<<endl;
					cout <<m_reads[tagidx * m_num_seg + ridx]<<endl;
					cout << combseq<<endl;
					cout << "check why"<<endl;
					//getchar();
				}
#endif
			}
			//if (outputsuffix)
			//{
			//	ofs<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M';
			//	cout<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M'<<endl;
			//	outputlen += (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length();
			//	
			//}
			//else
			//{
			//	ofs<<cur_map_len<<'M';
			//	cout <<cur_map_len<<'M'<<endl;
			//	outputlen += cur_map_len;
			//}

			ofs<<"\t*\t0\t0\t";

			//cout << "mapped first sequence"<<endl;

			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
			{
				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
			}

			//cout << "mapped middle sequence"<<endl;

			string missed_read_seq;

			if (strand == '+')
				missed_read_seq = m_read_words[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);
			else
				missed_read_seq = m_read_words_rc[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);

			ofs<<missed_read_seq;

			//cout << "mapped suffix sequence"<<endl;
			for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
			{
				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
			}

			//cout <<"output finished"<<endl;

			ofs<<'\t'<<qual_seq <<'\t'<<"NM:i:"<<score_sum<<'\t';

			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
			{
				ofs << vD_iter->index << ':' <<vD_iter->src_char <<'>'<<vD_iter->dst_char;

				if (vD_iter != differs_vec.end() - 1)
					ofs << ',';
			}

			ofs << endl;
		}
	}
	input_fs.close();
	ofs.close();
	//ofs_f0.close();

}

void
GenomeScan::FixHoleReady(const char* fixhole_file, const char* fixed_hole_file, hash_set<size_t>& junc_set, size_t flankcase, bool avoid, bool mark, bool exceed)
{
	cout << fixhole_file << endl;
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string chromseq = "";
	string prevchrom_file = "";

	string line;
	while(getline(input_fs, line))
	{
		if (line.empty())
			continue;
		//cout <<"read number"<<endl;

		//char *endptr;

		//if (rec_num != m_num_seg)
		//{
		//	cout << "missed segments in all mapped"<<endl;
		//	continue;
		//}

		//vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000], */strand;

		size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id, /*other_mapped, */id1, id2;

		//char s1[1000], s2[1000];

		//size_t sz = sizeof(doner_end);

		//string line;

		//cout <<"read line"<<endl;
		//HWI-EAS217:3:100:1000:1066#0/2	chr1	-	154799716	154799979	4-3	2-1	TCCATGCAGCTGT	GCTGAGGCTGNAG	TCTTGCCATCAAGAGCTCTCGGTG

		sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%llu\t%llu", tagname, chromname, &strand, &doner_end, &acceptor_st, 
			&doner_st_id, &doner_end_id, &acceptor_st_id, &acceptor_end_id,	doner_seq, acceptor_seq, tobe_fixed_seq, &id1, &id2);

		//cout << line << endl;

		//cout << line << endl;

		//if (strand == '-')
		//	continue;

		//cout << doner_end << '\t'<< acceptor_st << '\t'<< doner_st_id<< '\t'<< doner_end_id<< '\t'<< acceptor_st_id<< '\t'<< acceptor_end_id<<endl; getchar();

		string doner_str = doner_seq;

		string acceptor_str = acceptor_seq;

		string tobe_fixed_str = tobe_fixed_seq;

		size_t tobe_fixed_len = tobe_fixed_str.length();

		bool isinsertion = false;

		if (doner_end >= acceptor_st)
		{
			continue;
		}

		//indel
		if (acceptor_st + 1 < tobe_fixed_len + doner_end)
		{
			isinsertion = true;
			tobe_fixed_len = acceptor_st - doner_end + 1;
		}

		Masks* cur_masks_ptr;

		if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
			m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

		hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

		cur_masks_ptr = &(cur_hash_map_mask_iter->second);

		if (m_chrom_names_map.find(chromname) == m_chrom_names_map.end())
		{
			cout << "can't find chromosome: " << chromname << endl;
			continue;
		}

		string chromfile = m_chrom_files[m_chrom_names_map[chromname]];

		//cout << "prepare chrom sequence"<<endl;

		if (prevchrom_file == "" || chromfile != prevchrom_file)
		{
			chromseq.clear();
			//cout << "read chromo\t"<<chromfile<<endl;
			readchrom(chromfile.c_str(), chromseq);
			if (chromseq == "")
			{
				cout <<"empty chrom: "<<chromfile<<endl;
				exit(1);
			}
			prevchrom_file = chromfile;
		}

		if (isinsertion)
		{
			//cout << "sub 1:"<<endl;

			string prefix_read = tobe_fixed_str.substr(0, tobe_fixed_len);

			//cout << "sub 1 finished:"<<endl;

			//cout << "sub 2:"<<endl;

			string suffix_read = tobe_fixed_str.substr(tobe_fixed_str.length() - tobe_fixed_len, tobe_fixed_len);

			//cout << "sub 2 finished:"<<endl;

			string comb_read = prefix_read + suffix_read;

			WordPair comb_read_wp(comb_read);

			//cout << "sub 3:"<<endl;

			string chromstr = chromseq.substr(doner_end, tobe_fixed_len);

			//cout << "sub 3 finished:"<<endl;

			WordPair dup_chrom_wp, chrom_wp(chromstr);

			chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

			size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;

			size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);

			if (score <= m_max_mismatches)
			{
				int strandi = (strand == '+' ? 0 : 16);

				size_t insert_len = tobe_fixed_str.length() - tobe_fixed_len;

				ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/max_loc
					<<'M'<< insert_len<<'I'<<tobe_fixed_len - max_loc << 'M'<<'\t'<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

				if (score)
				{
					size_t seg1_suffix_len = tobe_fixed_len - max_loc;

					//cout << "sub 4:"<<endl;

					string mapped_read_str = comb_read.substr(0, max_loc) + comb_read.substr(comb_read.length() - seg1_suffix_len, seg1_suffix_len);

					//cout << "sub 4 finished:"<<endl;

					size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

					size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

					size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

					size_t merged_len = tobe_fixed_len;

					size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

					size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					if (strand == '-')
						pre_mapped_len = (m_num_seg - doner_end_id) * m_seg_len + doner_str.length();
					else
						pre_mapped_len = (doner_end_id - 1) * m_seg_len + doner_str.length();

					for (size_t i = 0; i < merged_len; ++i)
					{
						if ((selector_bit >> i) & comb_bits)
						{
							if (i > max_loc)
								ofs << pre_mapped_len + i + insert_len<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
							else
								ofs << pre_mapped_len + i<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
						}
					}
				}

				ofs /*<<"\t"<< flankstr*//*<<'\t' <<matched_flank_str2*/<< endl;
			}
		}
		else
		{
			//string tagnamestr = tagname;
			//size_t last_idx = tagnamestr.find_last_of("_");
			//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
			//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

			//if (avoid && m_reads_count[tagidx])
			//	continue;

			//size_t missed_id, fixid1, fixid2, offset1, offset2;

			////vector<GroupedBwtRec>::iterator bwt_rec_iter;

			//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
			//{
			//	//if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			//	//{
			//	//	cout << "missed more 1 segments"<<endl;
			//	//	exit(1);
			//	//}
			//	//else if (bwt_rec_iter->offset + m_seg_len != (bwt_rec_iter + 1)->offset)
			//	//{
			//	//	
			//	//}

			//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			//	{
			//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
			//		{
			//			cout << "missed more than 1 segments"<<endl;
			//			exit(1);
			//		}
			//		else
			//		{
			//			//missed_id = bwt_rec_iter->seg_id + 1;
			//			//break;
			//			fixid1 = bwt_rec_iter->seg_id;
			//			fixid2 = (bwt_rec_iter + 1)->seg_id;

			//			offset1 = bwt_rec_iter->offset;
			//			offset2 = (bwt_rec_iter + 1)->offset;

			//			missed_id = bwt_rec_iter->seg_id + 1;
			//			break;
			//		}
			//	}
			//}

			//if (exceed)
			//	;
			//else if (doner_end + tobe_fixed_len - 1 + m_interal > acceptor_st || doner_end + tobe_fixed_len + m_max_splice_dist - 1 < acceptor_st)
			//	continue;

			WordPair dup_tobe_fixed_read, tobe_fixed_read(tobe_fixed_str);//, prefix_suffix_append;

			tobe_fixed_read.duplicate_self(tobe_fixed_len, dup_tobe_fixed_read);

			//cout << "dup_missed_read"<<endl;

			//cout << dup_tobe_fixed_read.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;
			/*prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);*/

			//cout << "prepare read word"<<endl;

			//cout << strand << endl;

			//if (strand == '+')
			//{
			//	readID1 = tagidx * m_num_seg + fixid1 - 1;

			//	readID2 = tagidx * m_num_seg + fixid2 - 1;

			//	readID = tagidx * m_num_seg + missed_id - 1;

			//	WordPair wp1(m_read_words[readID1]);

			//	WordPair wp2(m_read_words[readID2]);

			//	//cout << "wp1 wp2"<<endl;

			//	//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	WordPair prefixappend;

			//	wp1.left_shift(m_seg_len);

			//	wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words[readID], prefixappend);

			//	prefixappend.left_shift(m_extend_bits);

			//	wp2.right_shift(m_seg_len - m_extend_bits);

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			//	prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//	//cout << "prefix_suffix_append"<<endl;

			//	//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			//	//cout << m_masks.prefix_seg_bits_on_ext << endl << m_masks.suffix_seg_bits_on_ext << endl << m_masks.comb_seg_bits_on_ext<<endl;

			//	prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			//		m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
			//}
			//else
			//{
			//	readID1 = tagidx * m_num_seg + m_num_seg - fixid1;

			//	readID2 = tagidx * m_num_seg + m_num_seg - fixid2;
			//	
			//	readID = tagidx * m_num_seg + m_num_seg - missed_id;

			//	WordPair wp1(m_read_words_rc[readID1]);

			//	WordPair wp2(m_read_words_rc[readID2]);

			//	//cout << "wp1 wp2"<<endl;

			//	//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	WordPair prefixappend;

			//	wp1.left_shift(m_seg_len);

			//	wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words_rc[readID], prefixappend);

			//	prefixappend.left_shift(m_extend_bits);

			//	wp2.right_shift(m_seg_len - m_extend_bits);

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			//	prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//	//cout << "prefix_suffix_append"<<endl;

			//	//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			//	prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			//		m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
			//}

			//cout << "dup_missed_read"<<endl;

			//cout << dup_missed_read.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;

			

			//cout << "get substring"<<endl;

			string prefixflankstr = chromseq.substr(doner_end, tobe_fixed_len + 2);

			string prefixmatchstr = chromseq.substr(doner_end - doner_str.length(), doner_str.length());

			//string prefixchromstr = chromseq.substr(doner_end, tobe_fixed_len);

			string suffixflankstr = chromseq.substr(acceptor_st - tobe_fixed_len - 1, tobe_fixed_len + 2);

			string suffixmatchstr = chromseq.substr(acceptor_st + 1, acceptor_str.length());

			//string suffixchromstr = chromseq.substr(acceptor_st - tobe_fixed_len, tobe_fixed_len);

			//cout << "prefix chrom " << endl;

			//cout << prefixflankstr << endl;

			//cout << "suffix chrom " << endl;

			//cout << suffixflankstr << endl;

			m_five_prim_suffix = WordPair(prefixflankstr);
			m_three_prim_prefix = WordPair(suffixflankstr);

			WordPair prefix_chrom = m_five_prim_suffix;

			m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

			prefix_chrom.left_shift(tobe_fixed_len - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

			WordPair comb_chrom_seq;

			//cout <<"combine chrom"<<endl;

			prefix_chrom.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

			//cout << comb_chrom_seq.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;

			size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;

			m_matched_flank = 0;

			m_matched_bads = 0;

			//cout << "scoring"<<endl;
			size_t score = Fixhole_score_selective_var_mask(dup_tobe_fixed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);

			//cout << score <<'\t' << max_loc<< '\t' << prim<<endl;

			//getchar();

			//ofstream* cur_ofs_ptr;

			//if (m_flank_case && !prim)
			//	continue;
			//cout <<"score: "<<score<<endl;
			//cout <<"max_loc: "<<max_loc<<endl;

			size_t intron_len = acceptor_st - doner_end + 1 - tobe_fixed_len;

			if (score <= m_max_mismatches && (prim >= flankcase || intron_len == 0) /*(!flankcase || (flankcase && prim))*/)
			{
				//string flankstrtest;

				//FlankString(flankstrtest);

				//cout << flankstrtest << endl;

				//getchar();
				//if ((prim > 0 && flankcase > 0) || (prim == 0 && flankcase == 0))
				//{
				//	if (mark)
				//	{
				//		//++m_cur_reads_count[tagidx];
				//	}
				//}
				//else
				//	continue;

				//size_t score_sum = score;

				//vector<GroupedBwtRec>::iterator sp_bwt_rec_iter;

				//vector<Differ> differs_vec;

				//size_t seg_id = 0;

				//cout << "first part 1"<<endl;

				//for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter; ++sp_bwt_rec_iter)
				//{
				//	vector<Differ>::iterator vD_iter;

				//	for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}

				//	++seg_id;
				//}

				//vector<Differ>::iterator vD_iter;

				//cout << "first part 2"<<endl;

				//for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//{
				//	if (vD_iter->index < m_seg_len - m_extend_bits)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}
				//}

				//++sp_bwt_rec_iter;

				//++seg_id;

				//size_t selector_lt = m_masks.mis_selector_lt;

				//cout << "middle part 1"<<endl;

				//for (size_t i = 0; i < max_loc; ++i)
				//{
				//	if (selector_lt & rbits)
				//	{
				//		char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - i);

				//		char src = prefixflankstr[i];

				//		if (src >= 'a' && src <= 'z')
				//			src += 'A' - 'a';

				//		if (dst == src)
				//		{
				//			cout << prefixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<i<<endl;
				//			continue;
				//			getchar();
				//		}

				//		differs_vec.push_back(Differ(seg_id * m_seg_len + i - m_extend_bits, src, dst));

				//		//cout << seg_id * m_seg_len + i - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				//	}

				//	selector_lt >>= 1;
				//}

				//selector_lt >>= m_masks.score_seg_buf_width_ext;

				//cout << "middle part 2"<<endl;

				//for (size_t j = max_loc; j < m_masks.score_seg_buf_width_ext; ++j)
				//{
				//	if (selector_lt & rbits)
				//	{
				//		char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - j);

				//		char src = suffixflankstr[j + suffixflankstr.length() - m_masks.score_seg_buf_width_ext];

				//		if (src >= 'a' && src <= 'z')
				//			src += 'A' - 'a';

				//		if (dst == src)
				//		{
				//			cout << suffixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<j<<endl;
				//			continue;
				//			getchar();
				//		}

				//		differs_vec.push_back(Differ(seg_id * m_seg_len + j - m_extend_bits, src, dst));

				//		//cout << seg_id * m_seg_len + j - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				//	}

				//	selector_lt >>= 1;
				//}

				//++seg_id;

				//cout << "third part 1"<<endl;

				//for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//{
				//	if (vD_iter->index >= m_extend_bits)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}
				//}

				//++sp_bwt_rec_iter;

				//++seg_id;

				//cout << "third part 2"<<endl;

				//for (; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//{
				//	vector<Differ>::iterator vD_iter;

				//	for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}

				//	++seg_id;
				//}

				//size_t prefix_mask = ALL_BITS_ON >> (SIXTY_FOUR - max_loc)<< (m_masks.score_seg_buf_width_ext - max_loc);

				//size_t suffix_mask = ALL_BITS_ON >> (m_masks.score_seg_buf_width_ext - max_loc);

				//size_t mbits = ((rbits >> m_masks.score_seg_buf_width_ext) & prefix_mask) & (rbits & suffix_mask);

				//for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//{
				//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
				//}

				//for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//{
				//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
				//}

				//if (score_sum > m_max_totalmismatch)
				//	continue;

				int strandi = 0;
				if (strand == '-')
					strandi = 16;

				double qual_score = 0;

				//const string& qual_seq = GetQualitySeq(tagidx, strandi, differs_vec, qual_score);

				//string prim5flankstr = prefixflankstr;
				//string prim3flankstr = suffixflankstr;
				//string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

				string flankstr;

				FlankString(flankstr);

				//size_t
				

				if (intron_len)
				{
					ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/max_loc
						<<'M'<< intron_len<<'N'<<tobe_fixed_len - max_loc << 'M'<<'\t'<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

					size_t prefixend = doner_end + max_loc - 1;

					size_t suffixst = prefixend + intron_len + 1;

					size_t combined_offset = (prefixend << THIRTY_TWO) + suffixst;

					junc_set.insert(combined_offset);
				}
				else				
					ofs<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/tobe_fixed_len<<'M'<<'\t'
						<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

				if (score)
				{
					string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

					string combined_read_str = tobe_fixed_str;

					size_t seg1_suffix_len = tobe_fixed_len - max_loc;

					string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

					size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

					size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

					size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

					//comb_bits += (rbits >> tobe_fixed_len) & seg1_mask_prefix;// + 

					size_t merged_len = tobe_fixed_len;

					size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

					size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					if (strand == '-')
						pre_mapped_len = (m_num_seg - doner_end_id) * m_seg_len + doner_str.length();
					else
						pre_mapped_len = (doner_end_id - 1) * m_seg_len + doner_str.length();

					for (size_t i = 0; i < merged_len; ++i)
					{
						if ((selector_bit >> i) & comb_bits)
						{
							ofs << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
						}
					}
				}

				ofs /*<<"\t"<< flankstr*//*<<'\t' <<matched_flank_str2*/<< endl;


				//<<endl;
				//cout<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t'<<endl;

				//cout <<"prefix part"<<endl;
				//size_t cur_map_len = fix_hole_vec.begin()->mapped_chrom_seq.length();

				//			size_t outputlen = 0;
				//			bool outputprefix = false;
				//
				//			vector <pair<size_t, size_t> > splice_way;
				//
				//			splice_way.push_back(make_pair(fix_hole_vec.begin()->offset, fix_hole_vec.begin()->mapped_chrom_seq.length()));
				//			for (sp_bwt_rec_iter = fix_hole_vec.begin() + 1; sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//			{
				//				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				//				{
				//					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"before expend : "<<cur_map_len<<endl;
				//					//cout <<"expend prefix part"<<endl;
				//					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//
				//					//cout <<"after expend : "<<cur_map_len<<endl;
				//				}
				//				else
				//				{
				//					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
				//					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
				//					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
				//					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//outputlen += cur_map_len;
				//					//outputprefix = true;
				//				}
				//			}
				//
				//			splice_way.back().second += max_loc - m_extend_bits;
				//
				//			splice_way.push_back(make_pair((bwt_rec_iter + 1)->offset - m_seg_len + max_loc - m_extend_bits, m_seg_len - max_loc + m_extend_bits));
				//
				//			//cout <<"middle part"<<endl;
				//			//if (outputprefix)
				//			//{
				//			//	ofs<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M';
				//			//	cout<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M'<<endl;
				//			//	outputlen += sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc;
				//			//}
				//			//else
				//			//{
				//			//	ofs<<cur_map_len + max_loc<<'M';
				//			//	cout<<cur_map_len + max_loc<<'M'<<endl;
				//			//	outputlen += cur_map_len + max_loc;
				//			//}
				//
				//			//ofs<<(bwt_rec_iter + 1)->offset - bwt_rec_iter->offset - m_seg_len - m_seg_len<<'N';
				//
				//			//cur_map_len = m_seg_len - max_loc + (bwt_rec_iter + 1)->mapped_chrom_seq.length();
				//			
				//			//cout <<"suffix part"<<endl;
				//
				//			splice_way.back().second += (bwt_rec_iter + 1)->mapped_chrom_seq.length();
				//
				//			//bool outputsuffix = false;
				//			for (sp_bwt_rec_iter = bwt_rec_iter + 2; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//			{
				//				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				//				{
				//					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"before expend : "<<cur_map_len<<endl;
				//					//cout <<"expend suffix part"<<endl;
				//					//
				//					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"after expend : "<<cur_map_len<<endl;
				//				}
				//				else
				//				{
				//					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
				//					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
				//					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
				//					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//outputlen += cur_map_len;
				//					//outputsuffix = true;
				//				}
				//			}
				//
				//			//cout << "tail part"<<endl;
				//
				//			vector <pair<size_t, size_t> >::iterator sp_iter;
				//
				//			size_t lensum = 0;
				//
				//			string combseq = "";
				//
				//			for (sp_iter = splice_way.begin(); sp_iter != splice_way.end(); ++sp_iter)
				//			{
				//				if (sp_iter != splice_way.begin())
				//				{
				//					ofs << sp_iter->first - (sp_iter - 1)->first - (sp_iter - 1)->second<<'N';
				//				}
				//				ofs << sp_iter->second<<'M';
				//
				//				combseq += chromseq.substr(sp_iter->first, sp_iter->second);
				//
				//				//cout << "combine position: "<<sp_iter->first<<'\t'<< sp_iter->second<<endl;
				//
				//				lensum += sp_iter->second;
				//			}
				//
				//			if (lensum != m_read_width)
				//				cout << "too short len"<<endl;
				//
				//			if (strandi)
				//				combseq = revcomp(combseq);
				//			for (size_t ridx = 0; ridx < m_num_seg; ++ridx)
				//			{
				//				int errsum = 0;
				//
				//#ifdef DEBUG_CHECK
				//				for (size_t j = 0; j < m_seg_len; ++j)
				//				{
				//					if (m_reads[tagidx * m_num_seg + ridx][j] != combseq[ridx * m_seg_len + j] && 
				//						m_reads[tagidx * m_num_seg + ridx][j] + 'A' - 'a' != combseq[ridx * m_seg_len + j] &&
				//						m_reads[tagidx * m_num_seg + ridx][j] - ('A' - 'a') != combseq[ridx * m_seg_len + j])
				//						errsum++;
				//				}
				//
				//				if (errsum > m_max_mismatches)
				//				{
				//					for (size_t j = 0; j < splice_way.size(); ++j)
				//						cout <<splice_way[j].first<<'\t'<<splice_way[j].second<<endl;
				//					cout <<ridx<<endl;
				//					cout <<m_reads[tagidx * m_num_seg + ridx]<<endl;
				//					cout << combseq<<endl;
				//					cout << "check why"<<endl;
				//					//getchar();
				//				}
				//#endif
				//			}
				//			//if (outputsuffix)
				//			//{
				//			//	ofs<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M';
				//			//	cout<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M'<<endl;
				//			//	outputlen += (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length();
				//			//	
				//			//}
				//			//else
				//			//{
				//			//	ofs<<cur_map_len<<'M';
				//			//	cout <<cur_map_len<<'M'<<endl;
				//			//	outputlen += cur_map_len;
				//			//}
				//
				//			ofs<<"\t*\t0\t0\t";
				//
				//			//cout << "mapped first sequence"<<endl;
				//
				//			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//			{
				//				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
				//			}
				//
				//			//cout << "mapped middle sequence"<<endl;
				//
				//			string missed_read_seq;
				//
				//			if (strand == '+')
				//				missed_read_seq = m_read_words[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);
				//			else
				//				missed_read_seq = m_read_words_rc[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);
				//
				//			ofs<<missed_read_seq;
				//
				//			//cout << "mapped suffix sequence"<<endl;
				//			for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//			{
				//				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
				//			}
				//
				//			//cout <<"output finished"<<endl;
				//
				//			ofs<<'\t'<<qual_seq <<'\t'<<"NM:i:"<<score_sum<<'\t';
				//
				//			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
				//			{
				//				ofs << vD_iter->index << ':' <<vD_iter->src_char <<'>'<<vD_iter->dst_char;
				//
				//				if (vD_iter != differs_vec.end() - 1)
				//					ofs << ',';
				//			}

				//ofs << endl;
			}
		}
	}
	input_fs.close();
	ofs.close();
	//ofs_f0.close();

}


void
GenomeScan::FixHoleReadyNonCanon(const char* fixhole_file, const char* fixed_hole_file, hash_set<size_t>& junc_set, size_t flankcase, bool avoid, bool mark, bool exceed)
{
	cout << "noncanon: " << fixhole_file << endl;
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string chromseq = "";
	string prevchrom_file = "";

	//1 unspliced, 2 canon-spliced,  3 noncanon-spliced,  4 insertion
	vector<vector<string> > mapped_segments(4);

	string prev_tagnme = "";

	string line;
	while(getline(input_fs, line))
	{
		if (line.empty())
			continue;
		//cout <<"read number"<<endl;

		//char *endptr;

		//if (rec_num != m_num_seg)
		//{
		//	cout << "missed segments in all mapped"<<endl;
		//	continue;
		//}

		//vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], /*differs[1000], */strand;

		size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id, /*other_mapped, */id1, id2;

		//char s1[1000], s2[1000];

		//size_t sz = sizeof(doner_end);

		//string line;

		//cout <<"read line"<<endl;
		//HWI-EAS217:3:100:1000:1066#0/2	chr1	-	154799716	154799979	4-3	2-1	TCCATGCAGCTGT	GCTGAGGCTGNAG	TCTTGCCATCAAGAGCTCTCGGTG

		sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%llu\t%llu", tagname, chromname, &strand, &doner_end, &acceptor_st, 
			&doner_st_id, &doner_end_id, &acceptor_st_id, &acceptor_end_id,	doner_seq, acceptor_seq, tobe_fixed_seq, &id1, &id2);

		if (prev_tagnme != tagname)
		{
			//1 unspliced, 2 canon-spliced,  3 noncanon-spliced,  4 insertion
 
			if (mapped_segments[2].size() > 50)
				mapped_segments[2].clear();

			if (!mapped_segments[2].empty() && (!mapped_segments[0].empty() || !mapped_segments[1].empty()))
				mapped_segments[2].clear();

			vector<vector<string> >::iterator mapped_segments_iter;
			for (mapped_segments_iter = mapped_segments.begin(); mapped_segments_iter != mapped_segments.end(); ++mapped_segments_iter)
			{
				vector<string>::iterator mapped_segments_cat_iter;
				for (mapped_segments_cat_iter = mapped_segments_iter->begin(); mapped_segments_cat_iter != mapped_segments_iter->end(); ++mapped_segments_cat_iter)
				{
					ofs << *mapped_segments_cat_iter;
				}				
			}

			mapped_segments.clear();
			mapped_segments.resize(4);		

			prev_tagnme = tagname;
		}

		//cout << line << endl;

		//cout << line << endl;

		//if (strand == '-')
		//	continue;

		//cout << doner_end << '\t'<< acceptor_st << '\t'<< doner_st_id<< '\t'<< doner_end_id<< '\t'<< acceptor_st_id<< '\t'<< acceptor_end_id<<endl; getchar();

		string doner_str = doner_seq;

		string acceptor_str = acceptor_seq;

		string tobe_fixed_str = tobe_fixed_seq;

		size_t tobe_fixed_len = tobe_fixed_str.length();

		bool isinsertion = false;

		if (doner_end >= acceptor_st)
		{
			continue;
		}

		//indel
		if (acceptor_st + 1 < tobe_fixed_len + doner_end)
		{
			isinsertion = true;
			tobe_fixed_len = acceptor_st - doner_end + 1;
		}

		Masks* cur_masks_ptr;

		if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
			m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

		hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

		cur_masks_ptr = &(cur_hash_map_mask_iter->second);

		if (m_chrom_names_map.find(chromname) == m_chrom_names_map.end())
		{
			cout << "can't find chromosome: " << chromname << endl;
			continue;
		}

		string chromfile = m_chrom_files[m_chrom_names_map[chromname]];

		//cout << "prepare chrom sequence"<<endl;

		if (prevchrom_file == "" || chromfile != prevchrom_file)
		{
			chromseq.clear();
			//cout << "read chromo\t"<<chromfile<<endl;
			readchrom(chromfile.c_str(), chromseq);
			if (chromseq == "")
			{
				cout <<"empty chrom: "<<chromfile<<endl;
				exit(1);
			}
			prevchrom_file = chromfile;
		}

		if (isinsertion)
		{
			//cout << "sub 1:"<<endl;
			ostringstream s;

			string prefix_read = tobe_fixed_str.substr(0, tobe_fixed_len);

			//cout << "sub 1 finished:"<<endl;

			//cout << "sub 2:"<<endl;

			string suffix_read = tobe_fixed_str.substr(tobe_fixed_str.length() - tobe_fixed_len, tobe_fixed_len);

			//cout << "sub 2 finished:"<<endl;

			string comb_read = prefix_read + suffix_read;

			WordPair comb_read_wp(comb_read);

			//cout << "sub 3:"<<endl;

			string chromstr = chromseq.substr(doner_end, tobe_fixed_len);

			//cout << "sub 3 finished:"<<endl;

			WordPair dup_chrom_wp, chrom_wp(chromstr);

			chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

			size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;

			size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);

			if (score <= m_max_mismatches)
			{
				int strandi = (strand == '+' ? 0 : 16);

				size_t insert_len = tobe_fixed_str.length() - tobe_fixed_len;

				s<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/max_loc
					<<'M'<< insert_len<<'I'<<tobe_fixed_len - max_loc << 'M'<<'\t'<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

				if (score)
				{
					size_t seg1_suffix_len = tobe_fixed_len - max_loc;

					//cout << "sub 4:"<<endl;

					string mapped_read_str = comb_read.substr(0, max_loc) + comb_read.substr(comb_read.length() - seg1_suffix_len, seg1_suffix_len);

					//cout << "sub 4 finished:"<<endl;

					size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

					size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

					size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

					size_t merged_len = tobe_fixed_len;

					size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

					size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					if (strand == '-')
						pre_mapped_len = (m_num_seg - doner_end_id) * m_seg_len + doner_str.length();
					else
						pre_mapped_len = (doner_end_id - 1) * m_seg_len + doner_str.length();

					for (size_t i = 0; i < merged_len; ++i)
					{
						if ((selector_bit >> i) & comb_bits)
						{
							if (i > max_loc)
								s << pre_mapped_len + i + insert_len<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
							else
								s << pre_mapped_len + i<<':'<<chromstr[i]<<'>'<<mapped_read_str[i]<<',';
						}
					}
				}

				s /*<<"\t"<< flankstr*//*<<'\t' <<matched_flank_str2*/<< endl;

				//1 unspliced, 2 canon-spliced,  3 noncanon-spliced,  4 insertion
				mapped_segments[3].push_back(s.str());
			}
		}
		else
		{
			//string tagnamestr = tagname;
			//size_t last_idx = tagnamestr.find_last_of("_");
			//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
			//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

			//if (avoid && m_reads_count[tagidx])
			//	continue;

			//size_t missed_id, fixid1, fixid2, offset1, offset2;

			////vector<GroupedBwtRec>::iterator bwt_rec_iter;

			//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
			//{
			//	//if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			//	//{
			//	//	cout << "missed more 1 segments"<<endl;
			//	//	exit(1);
			//	//}
			//	//else if (bwt_rec_iter->offset + m_seg_len != (bwt_rec_iter + 1)->offset)
			//	//{
			//	//	
			//	//}

			//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			//	{
			//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
			//		{
			//			cout << "missed more than 1 segments"<<endl;
			//			exit(1);
			//		}
			//		else
			//		{
			//			//missed_id = bwt_rec_iter->seg_id + 1;
			//			//break;
			//			fixid1 = bwt_rec_iter->seg_id;
			//			fixid2 = (bwt_rec_iter + 1)->seg_id;

			//			offset1 = bwt_rec_iter->offset;
			//			offset2 = (bwt_rec_iter + 1)->offset;

			//			missed_id = bwt_rec_iter->seg_id + 1;
			//			break;
			//		}
			//	}
			//}

			//if (exceed)
			//	;
			//else if (doner_end + tobe_fixed_len - 1 + m_interal > acceptor_st || doner_end + tobe_fixed_len + m_max_splice_dist - 1 < acceptor_st)
			//	continue;

			ostringstream s;

			WordPair dup_tobe_fixed_read, tobe_fixed_read(tobe_fixed_str);//, prefix_suffix_append;

			tobe_fixed_read.duplicate_self(tobe_fixed_len, dup_tobe_fixed_read);

			//cout << "dup_missed_read"<<endl;

			//cout << dup_tobe_fixed_read.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;
			/*prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);*/

			//cout << "prepare read word"<<endl;

			//cout << strand << endl;

			//if (strand == '+')
			//{
			//	readID1 = tagidx * m_num_seg + fixid1 - 1;

			//	readID2 = tagidx * m_num_seg + fixid2 - 1;

			//	readID = tagidx * m_num_seg + missed_id - 1;

			//	WordPair wp1(m_read_words[readID1]);

			//	WordPair wp2(m_read_words[readID2]);

			//	//cout << "wp1 wp2"<<endl;

			//	//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	WordPair prefixappend;

			//	wp1.left_shift(m_seg_len);

			//	wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words[readID], prefixappend);

			//	prefixappend.left_shift(m_extend_bits);

			//	wp2.right_shift(m_seg_len - m_extend_bits);

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			//	prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//	//cout << "prefix_suffix_append"<<endl;

			//	//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			//	//cout << m_masks.prefix_seg_bits_on_ext << endl << m_masks.suffix_seg_bits_on_ext << endl << m_masks.comb_seg_bits_on_ext<<endl;

			//	prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			//		m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
			//}
			//else
			//{
			//	readID1 = tagidx * m_num_seg + m_num_seg - fixid1;

			//	readID2 = tagidx * m_num_seg + m_num_seg - fixid2;
			//	
			//	readID = tagidx * m_num_seg + m_num_seg - missed_id;

			//	WordPair wp1(m_read_words_rc[readID1]);

			//	WordPair wp2(m_read_words_rc[readID2]);

			//	//cout << "wp1 wp2"<<endl;

			//	//cout << wp1.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_seg_len) <<endl;

			//	WordPair prefixappend;

			//	wp1.left_shift(m_seg_len);

			//	wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words_rc[readID], prefixappend);

			//	prefixappend.left_shift(m_extend_bits);

			//	wp2.right_shift(m_seg_len - m_extend_bits);

			//	//cout << wp2.tostring3(ALL_BITS_ON, m_max_mismatches)<<endl;

			//	prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

			//	//cout << "prefix_suffix_append"<<endl;

			//	//cout << prefix_suffix_append.tostring3(ALL_BITS_ON, m_seg_len + (2 * m_max_mismatches)) << endl;

			//	prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
			//		m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);
			//}

			//cout << "dup_missed_read"<<endl;

			//cout << dup_missed_read.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;

			

			//cout << "get substring"<<endl;

			string prefixflankstr = chromseq.substr(doner_end, tobe_fixed_len + 2);

			string prefixmatchstr = chromseq.substr(doner_end - doner_str.length(), doner_str.length());

			//string prefixchromstr = chromseq.substr(doner_end, tobe_fixed_len);

			string suffixflankstr = chromseq.substr(acceptor_st - tobe_fixed_len - 1, tobe_fixed_len + 2);

			string suffixmatchstr = chromseq.substr(acceptor_st + 1, acceptor_str.length());

			//string suffixchromstr = chromseq.substr(acceptor_st - tobe_fixed_len, tobe_fixed_len);

			//cout << "prefix chrom " << endl;

			//cout << prefixflankstr << endl;

			//cout << "suffix chrom " << endl;

			//cout << suffixflankstr << endl;

			m_five_prim_suffix = WordPair(prefixflankstr);
			m_three_prim_prefix = WordPair(suffixflankstr);

			WordPair prefix_chrom = m_five_prim_suffix;

			m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

			prefix_chrom.left_shift(tobe_fixed_len - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

			WordPair comb_chrom_seq;

			//cout <<"combine chrom"<<endl;

			prefix_chrom.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

			//cout << comb_chrom_seq.tostring3(ALL_BITS_ON, 2 * tobe_fixed_len)<<endl;

			size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;

			m_matched_flank = 0;

			m_matched_bads = 0;

			//cout << "scoring"<<endl;
			size_t score = Fixhole_score_selective_var_mask(dup_tobe_fixed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);

			//cout << score <<'\t' << max_loc<< '\t' << prim<<endl;

			//getchar();

			//ofstream* cur_ofs_ptr;

			//if (m_flank_case && !prim)
			//	continue;
			//cout <<"score: "<<score<<endl;
			//cout <<"max_loc: "<<max_loc<<endl;

			size_t intron_len = acceptor_st - doner_end + 1 - tobe_fixed_len;

			if (score <= m_max_mismatches && (prim >= flankcase || intron_len == 0) /*(!flankcase || (flankcase && prim))*/)
			{
				//string flankstrtest;

				//FlankString(flankstrtest);

				//cout << flankstrtest << endl;

				//getchar();
				//if ((prim > 0 && flankcase > 0) || (prim == 0 && flankcase == 0))
				//{
				//	if (mark)
				//	{
				//		//++m_cur_reads_count[tagidx];
				//	}
				//}
				//else
				//	continue;

				//size_t score_sum = score;

				//vector<GroupedBwtRec>::iterator sp_bwt_rec_iter;

				//vector<Differ> differs_vec;

				//size_t seg_id = 0;

				//cout << "first part 1"<<endl;

				//for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter; ++sp_bwt_rec_iter)
				//{
				//	vector<Differ>::iterator vD_iter;

				//	for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}

				//	++seg_id;
				//}

				//vector<Differ>::iterator vD_iter;

				//cout << "first part 2"<<endl;

				//for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//{
				//	if (vD_iter->index < m_seg_len - m_extend_bits)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}
				//}

				//++sp_bwt_rec_iter;

				//++seg_id;

				//size_t selector_lt = m_masks.mis_selector_lt;

				//cout << "middle part 1"<<endl;

				//for (size_t i = 0; i < max_loc; ++i)
				//{
				//	if (selector_lt & rbits)
				//	{
				//		char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - i);

				//		char src = prefixflankstr[i];

				//		if (src >= 'a' && src <= 'z')
				//			src += 'A' - 'a';

				//		if (dst == src)
				//		{
				//			cout << prefixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<i<<endl;
				//			continue;
				//			getchar();
				//		}

				//		differs_vec.push_back(Differ(seg_id * m_seg_len + i - m_extend_bits, src, dst));

				//		//cout << seg_id * m_seg_len + i - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				//	}

				//	selector_lt >>= 1;
				//}

				//selector_lt >>= m_masks.score_seg_buf_width_ext;

				//cout << "middle part 2"<<endl;

				//for (size_t j = max_loc; j < m_masks.score_seg_buf_width_ext; ++j)
				//{
				//	if (selector_lt & rbits)
				//	{
				//		char dst = prefix_suffix_append.get_char(m_masks.prefix_suffix_append_mask, m_masks.score_seg_buf_width_ext - j);

				//		char src = suffixflankstr[j + suffixflankstr.length() - m_masks.score_seg_buf_width_ext];

				//		if (src >= 'a' && src <= 'z')
				//			src += 'A' - 'a';

				//		if (dst == src)
				//		{
				//			cout << suffixflankstr << endl << prefix_suffix_append.tostring3(ALL_BITS_ON, m_masks.score_seg_buf_width_ext)<<endl<<j<<endl;
				//			continue;
				//			getchar();
				//		}

				//		differs_vec.push_back(Differ(seg_id * m_seg_len + j - m_extend_bits, src, dst));

				//		//cout << seg_id * m_seg_len + j - m_max_mismatches<<'\t'<< src<<'\t'<< dst<<endl;
				//	}

				//	selector_lt >>= 1;
				//}

				//++seg_id;

				//cout << "third part 1"<<endl;

				//for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//{
				//	if (vD_iter->index >= m_extend_bits)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}
				//}

				//++sp_bwt_rec_iter;

				//++seg_id;

				//cout << "third part 2"<<endl;

				//for (; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//{
				//	vector<Differ>::iterator vD_iter;

				//	for (vD_iter = sp_bwt_rec_iter->differs_vec.begin(); vD_iter != sp_bwt_rec_iter->differs_vec.end(); ++vD_iter)
				//	{
				//		differs_vec.push_back(Differ(seg_id * m_seg_len + vD_iter->index, vD_iter->src_char, vD_iter->dst_char));
				//		++score_sum;

				//		//cout << seg_id * m_seg_len + vD_iter->index<<'\t'<< vD_iter->src_char<<'\t'<< vD_iter->dst_char<<endl;
				//	}

				//	++seg_id;
				//}

				//size_t prefix_mask = ALL_BITS_ON >> (SIXTY_FOUR - max_loc)<< (m_masks.score_seg_buf_width_ext - max_loc);

				//size_t suffix_mask = ALL_BITS_ON >> (m_masks.score_seg_buf_width_ext - max_loc);

				//size_t mbits = ((rbits >> m_masks.score_seg_buf_width_ext) & prefix_mask) & (rbits & suffix_mask);

				//for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//{
				//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
				//}

				//for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//{
				//	score_sum += sp_bwt_rec_iter->differs.length() / 5;
				//}

				//if (score_sum > m_max_totalmismatch)
				//	continue;

				int strandi = 0;
				if (strand == '-')
					strandi = 16;

				double qual_score = 0;

				//const string& qual_seq = GetQualitySeq(tagidx, strandi, differs_vec, qual_score);

				//string prim5flankstr = prefixflankstr;
				//string prim3flankstr = suffixflankstr;
				//string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

				string flankstr;

				FlankString(flankstr);

				//size_t
				

				if (intron_len)
				{
					s<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/max_loc
						<<'M'<< intron_len<<'N'<<tobe_fixed_len - max_loc << 'M'<<'\t'<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

					size_t prefixend = doner_end + max_loc - 1;

					size_t suffixst = prefixend + intron_len + 1;

					size_t combined_offset = (prefixend << THIRTY_TWO) + suffixst;

					junc_set.insert(combined_offset);
				}
				else				
					s<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<doner_end<<'\t'<<doner_end_id<<'-'<<acceptor_st_id <<'\t'<<doner_str.length() + 1 <<'\t'<</*size_t(qual_score)<<'\t'<<*/tobe_fixed_len<<'M'<<'\t'
						<<tobe_fixed_str<<'\t'<<id1<<'\t'<<id2<<"\tNM:i:"<<score<<'\t'/*<<flankstr*/;

				if (score)
				{
					string comb_chrom_str1 = comb_chrom_seq.tostring3(ALL_BITS_ON, tobe_fixed_len * 2);

					string combined_read_str = tobe_fixed_str;

					size_t seg1_suffix_len = tobe_fixed_len - max_loc;

					string mapped_chrom_str = comb_chrom_str1.substr(0, max_loc) + comb_chrom_str1.substr(comb_chrom_str1.length() - seg1_suffix_len, seg1_suffix_len);

					size_t seg1_mask_prefix = cur_masks_ptr->suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

					size_t seg1_mask_suffix = cur_masks_ptr->suffix_seg_bits_on >> max_loc;

					size_t comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);

					//comb_bits += (rbits >> tobe_fixed_len) & seg1_mask_prefix;// + 

					size_t merged_len = tobe_fixed_len;

					size_t selector_bit = LEAST_SIG_BIT << (merged_len - 1);

					size_t pre_mapped_len;// = (cur_TagSeed->st_seg_id - 1) * m_seg_len;

					if (strand == '-')
						pre_mapped_len = (m_num_seg - doner_end_id) * m_seg_len + doner_str.length();
					else
						pre_mapped_len = (doner_end_id - 1) * m_seg_len + doner_str.length();

					for (size_t i = 0; i < merged_len; ++i)
					{
						if ((selector_bit >> i) & comb_bits)
						{
							s << pre_mapped_len + i<<':'<<mapped_chrom_str[i]<<'>'<<combined_read_str[i]<<',';
						}
					}
				}

				s /*<<"\t"<< flankstr*//*<<'\t' <<matched_flank_str2*/<< endl;

				if (intron_len)
				{
					if (prim >= 5)
					{
						//1 unspliced, 2 canon-spliced,  3 noncanon-spliced,  4 insertion
						mapped_segments[1].push_back(s.str());
					}
					else
					{
						mapped_segments[2].push_back(s.str());
					}
				}
				else
					mapped_segments[0].push_back(s.str());

				//<<endl;
				//cout<<tagname<<'\t'<<strandi<<'\t'<< chromname<<'\t'<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t'<<endl;

				//cout <<"prefix part"<<endl;
				//size_t cur_map_len = fix_hole_vec.begin()->mapped_chrom_seq.length();

				//			size_t outputlen = 0;
				//			bool outputprefix = false;
				//
				//			vector <pair<size_t, size_t> > splice_way;
				//
				//			splice_way.push_back(make_pair(fix_hole_vec.begin()->offset, fix_hole_vec.begin()->mapped_chrom_seq.length()));
				//			for (sp_bwt_rec_iter = fix_hole_vec.begin() + 1; sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//			{
				//				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				//				{
				//					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"before expend : "<<cur_map_len<<endl;
				//					//cout <<"expend prefix part"<<endl;
				//					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//
				//					//cout <<"after expend : "<<cur_map_len<<endl;
				//				}
				//				else
				//				{
				//					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
				//					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
				//					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
				//					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//outputlen += cur_map_len;
				//					//outputprefix = true;
				//				}
				//			}
				//
				//			splice_way.back().second += max_loc - m_extend_bits;
				//
				//			splice_way.push_back(make_pair((bwt_rec_iter + 1)->offset - m_seg_len + max_loc - m_extend_bits, m_seg_len - max_loc + m_extend_bits));
				//
				//			//cout <<"middle part"<<endl;
				//			//if (outputprefix)
				//			//{
				//			//	ofs<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M';
				//			//	cout<<sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc<<'M'<<endl;
				//			//	outputlen += sp_bwt_rec_iter->mapped_chrom_seq.length() + max_loc;
				//			//}
				//			//else
				//			//{
				//			//	ofs<<cur_map_len + max_loc<<'M';
				//			//	cout<<cur_map_len + max_loc<<'M'<<endl;
				//			//	outputlen += cur_map_len + max_loc;
				//			//}
				//
				//			//ofs<<(bwt_rec_iter + 1)->offset - bwt_rec_iter->offset - m_seg_len - m_seg_len<<'N';
				//
				//			//cur_map_len = m_seg_len - max_loc + (bwt_rec_iter + 1)->mapped_chrom_seq.length();
				//			
				//			//cout <<"suffix part"<<endl;
				//
				//			splice_way.back().second += (bwt_rec_iter + 1)->mapped_chrom_seq.length();
				//
				//			//bool outputsuffix = false;
				//			for (sp_bwt_rec_iter = bwt_rec_iter + 2; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//			{
				//				if ((sp_bwt_rec_iter - 1)->offset + (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length() == (sp_bwt_rec_iter)->offset)
				//				{
				//					splice_way.back().second += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"before expend : "<<cur_map_len<<endl;
				//					//cout <<"expend suffix part"<<endl;
				//					//
				//					//cur_map_len += sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//cout <<"after expend : "<<cur_map_len<<endl;
				//				}
				//				else
				//				{
				//					splice_way.push_back(make_pair(sp_bwt_rec_iter->offset, sp_bwt_rec_iter->mapped_chrom_seq.length()));
				//					//ofs <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N';
				//					//cout <<cur_map_len<<'M'<<sp_bwt_rec_iter->offset - (sp_bwt_rec_iter - 1)->offset - (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'N'<<endl;
				//					//cur_map_len = sp_bwt_rec_iter->mapped_chrom_seq.length();
				//					//outputlen += cur_map_len;
				//					//outputsuffix = true;
				//				}
				//			}
				//
				//			//cout << "tail part"<<endl;
				//
				//			vector <pair<size_t, size_t> >::iterator sp_iter;
				//
				//			size_t lensum = 0;
				//
				//			string combseq = "";
				//
				//			for (sp_iter = splice_way.begin(); sp_iter != splice_way.end(); ++sp_iter)
				//			{
				//				if (sp_iter != splice_way.begin())
				//				{
				//					ofs << sp_iter->first - (sp_iter - 1)->first - (sp_iter - 1)->second<<'N';
				//				}
				//				ofs << sp_iter->second<<'M';
				//
				//				combseq += chromseq.substr(sp_iter->first, sp_iter->second);
				//
				//				//cout << "combine position: "<<sp_iter->first<<'\t'<< sp_iter->second<<endl;
				//
				//				lensum += sp_iter->second;
				//			}
				//
				//			if (lensum != m_read_width)
				//				cout << "too short len"<<endl;
				//
				//			if (strandi)
				//				combseq = revcomp(combseq);
				//			for (size_t ridx = 0; ridx < m_num_seg; ++ridx)
				//			{
				//				int errsum = 0;
				//
				//#ifdef DEBUG_CHECK
				//				for (size_t j = 0; j < m_seg_len; ++j)
				//				{
				//					if (m_reads[tagidx * m_num_seg + ridx][j] != combseq[ridx * m_seg_len + j] && 
				//						m_reads[tagidx * m_num_seg + ridx][j] + 'A' - 'a' != combseq[ridx * m_seg_len + j] &&
				//						m_reads[tagidx * m_num_seg + ridx][j] - ('A' - 'a') != combseq[ridx * m_seg_len + j])
				//						errsum++;
				//				}
				//
				//				if (errsum > m_max_mismatches)
				//				{
				//					for (size_t j = 0; j < splice_way.size(); ++j)
				//						cout <<splice_way[j].first<<'\t'<<splice_way[j].second<<endl;
				//					cout <<ridx<<endl;
				//					cout <<m_reads[tagidx * m_num_seg + ridx]<<endl;
				//					cout << combseq<<endl;
				//					cout << "check why"<<endl;
				//					//getchar();
				//				}
				//#endif
				//			}
				//			//if (outputsuffix)
				//			//{
				//			//	ofs<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M';
				//			//	cout<<(sp_bwt_rec_iter - 1)->mapped_chrom_seq.length()<<'M'<<endl;
				//			//	outputlen += (sp_bwt_rec_iter - 1)->mapped_chrom_seq.length();
				//			//	
				//			//}
				//			//else
				//			//{
				//			//	ofs<<cur_map_len<<'M';
				//			//	cout <<cur_map_len<<'M'<<endl;
				//			//	outputlen += cur_map_len;
				//			//}
				//
				//			ofs<<"\t*\t0\t0\t";
				//
				//			//cout << "mapped first sequence"<<endl;
				//
				//			for (sp_bwt_rec_iter = fix_hole_vec.begin(); sp_bwt_rec_iter != bwt_rec_iter + 1; ++sp_bwt_rec_iter)
				//			{
				//				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
				//			}
				//
				//			//cout << "mapped middle sequence"<<endl;
				//
				//			string missed_read_seq;
				//
				//			if (strand == '+')
				//				missed_read_seq = m_read_words[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);
				//			else
				//				missed_read_seq = m_read_words_rc[readID].tostring3(m_masks.suffix_seg_bits_on, m_seg_len);
				//
				//			ofs<<missed_read_seq;
				//
				//			//cout << "mapped suffix sequence"<<endl;
				//			for (sp_bwt_rec_iter = bwt_rec_iter + 1; sp_bwt_rec_iter != fix_hole_vec.end(); ++sp_bwt_rec_iter)
				//			{
				//				ofs<<sp_bwt_rec_iter->mapped_chrom_seq;
				//			}
				//
				//			//cout <<"output finished"<<endl;
				//
				//			ofs<<'\t'<<qual_seq <<'\t'<<"NM:i:"<<score_sum<<'\t';
				//
				//			for (vD_iter = differs_vec.begin(); vD_iter != differs_vec.end(); ++vD_iter)
				//			{
				//				ofs << vD_iter->index << ':' <<vD_iter->src_char <<'>'<<vD_iter->dst_char;
				//
				//				if (vD_iter != differs_vec.end() - 1)
				//					ofs << ',';
				//			}

				//ofs << endl;
			}
		}
	}

	if (mapped_segments[2].size() > 50)
		mapped_segments[2].clear();

	if (!mapped_segments[2].empty() && (!mapped_segments[0].empty() || !mapped_segments[1].empty()))
		mapped_segments[2].clear();

	vector<vector<string> >::iterator mapped_segments_iter;
	for (mapped_segments_iter = mapped_segments.begin(); mapped_segments_iter != mapped_segments.end(); ++mapped_segments_iter)
	{
		vector<string>::iterator mapped_segments_cat_iter;
		for (mapped_segments_cat_iter = mapped_segments_iter->begin(); mapped_segments_cat_iter != mapped_segments_iter->end(); ++mapped_segments_cat_iter)
		{
			ofs << *mapped_segments_cat_iter;
		}				
	}

	input_fs.close();
	ofs.close();
	//ofs_f0.close();

}


void
GenomeScan::FixFusionReady(const char* fusion_file, const char* fixed_fusion_file, CHROM_FUSION_JUNC_HASH_COMB& conj_hash_map, bool avoid, bool mark, size_t flankcase, string& chromseq1, string& chromseq2)
{
	ifstream input_fs(fusion_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fusion_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_fusion_file);

	//string chromseq1 = "";
	//string chromseq2 = "";
	//string prevchrom_file1 = "";
	//string prevchrom_file2 = "";

	string line;
	while(getline(input_fs, line))
	{
		//cout <<"read number"<<endl;

		if (line.empty())
			continue;

		char chromname1[1000], chromname2[1000], tagname[1000], differs[1000], doner_seq[1000], acceptor_seq[1000], tobe_fixed_seq[1000], strand1, strand2, quality_seq[1000], doner_splice[1000], acceptor_splice[1000];
		size_t doner_end, acceptor_st, doner_st_id, doner_end_id, acceptor_st_id, acceptor_end_id;

		//string line;

		//HWI-EAS217:4:34:972:1977#0/1	chr11	chr9	+-	118979705	101107896 1-1 3-4	seg1+ seg3- seg2+_comb_seg2-
		//TRAN00000105895:153:9 chr20 chr20 -+ 25614644 18234950 1-1 2-2 AGAACTACAGCAA TTCTAAACCAGAA AGTGAGACCAGGGGGAATTTCATT CBCCCDCCCCCCCBBBBBBAAAAAAAAAAAA@@?@@@@?@@@@??>>>=> 38:G>T,25:A>G, 

		//cout <<"read line"<<endl;

		int read_count = sscanf(line.c_str(), "%s\t%s\t%s\t%c%c\t%llu\t%llu\t%llu-%llu\t%llu-%llu\t%s\t%s\t%s\t%s\t%s\t%s\t%s", tagname, chromname1, chromname2, &strand1, &strand2, 
			&doner_end, &acceptor_st, &doner_st_id, &doner_end_id, &acceptor_st_id, &acceptor_end_id, doner_seq, acceptor_seq, tobe_fixed_seq, quality_seq, doner_splice, acceptor_splice, differs);

		string differstr;

		if (read_count == 18)
			differstr = differs;

		//string chr1, chr2;

		string doner_str = doner_seq;

		string acceptor_str = acceptor_seq;

		string tobe_fixed_str = tobe_fixed_seq;

		size_t tobe_fixed_len = tobe_fixed_str.length();

		Masks* cur_masks_ptr;

		if (m_hash_map_masks.find(tobe_fixed_len) == m_hash_map_masks.end())
			m_hash_map_masks.insert(make_pair(tobe_fixed_len, Masks(tobe_fixed_len)));

		hash_map<size_t, Masks>::iterator cur_hash_map_mask_iter = m_hash_map_masks.find(tobe_fixed_len);

		cur_masks_ptr = &(cur_hash_map_mask_iter->second);

		//string tagnamestr = tagname;
		//size_t last_idx = tagnamestr.find_last_of("_");
		//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		//if (avoid && m_reads_count[tagidx])
		//	continue;

		//size_t missed_id, fixid1, fixid2, offset1, offset2;

		////vector<GroupedBwtRec>::iterator bwt_rec_iter;

		//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		//{
		//	//if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
		//	//{
		//	//	cout << "missed more 1 segments"<<endl;
		//	//	exit(1);
		//	//}
		//	//else if (bwt_rec_iter->offset + m_seg_len != (bwt_rec_iter + 1)->offset)
		//	//{
		//	//	
		//	//}

		//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
		//	{
		//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
		//		{
		//			cout << "missed more than 1 segments"<<endl;
		//			exit(1);
		//		}
		//		else
		//		{
		//			//missed_id = bwt_rec_iter->seg_id + 1;
		//			//break;
		//			fixid1 = bwt_rec_iter->seg_id;
		//			fixid2 = (bwt_rec_iter + 1)->seg_id;

		//			offset1 = bwt_rec_iter->offset;
		//			offset2 = (bwt_rec_iter + 1)->offset;

		//			missed_id = bwt_rec_iter->seg_id + 1;
		//			break;
		//		}
		//	}
		//}
		//if (exceed)
		//	;
		//else if (doner_end + tobe_fixed_len - 1 + m_interal > acceptor_st || doner_end + tobe_fixed_len + m_max_splice_dist - 1 < acceptor_st)
		//	continue;

		//cout << line << endl;

		//string tagnamestr = tagname;
		//size_t last_idx = tagnamestr.find_last_of("_");
		//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		//size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		//if (avoid && m_reads_count[tagidx])
		//	continue;

		

		//char strand1, strand2;

		//size_t offset1, offset2;

		//cout << "find missed seg"<<endl;

		//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		//{
		//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
		//	{
		//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
		//		{
		//			cout << "missed more than 1 segments"<<endl;

		//			for (size_t i = 0; i < fix_hole_vec.size(); ++i)
		//			{
		//				cout << tagname <<'\t' << fix_hole_vec[i].offset << '\t' <<fix_hole_vec[i].seg_id <<endl;
		//			}
		//			exit(1);
		//		}
		//		else
		//		{
		//			chr1 = bwt_rec_iter->chromname;

		//			chr2 = (bwt_rec_iter + 1)->chromname;

		//			strand1 =  bwt_rec_iter->strand;

		//			strand2 = (bwt_rec_iter + 1)->strand;

		//			offset1 = bwt_rec_iter->offset;

		//			offset2 = (bwt_rec_iter + 1)->offset;

		//			missed_id = bwt_rec_iter->seg_id + 1;
		//			break;
		//		}
		//	}
		//}

		//vector<GroupedBwtRec>::iterator bwt_rec_cur_iter;

		//size_t other_mismatch = 0;

		//for (bwt_rec_cur_iter = fix_hole_vec.begin(); bwt_rec_cur_iter != fix_hole_vec.end(); ++bwt_rec_cur_iter)
		//{
		//	vector<Differ>::iterator differs_vec_iter;
		//	if (bwt_rec_cur_iter == bwt_rec_iter)
		//	{
		//		for (differs_vec_iter = bwt_rec_cur_iter->differs_vec.begin(); differs_vec_iter != bwt_rec_cur_iter->differs_vec.end(); ++differs_vec_iter)
		//		{
		//			if (differs_vec_iter->index < m_seg_len - m_extend_bits)
		//				++other_mismatch;
		//		}
		//	}
		//	else if (bwt_rec_cur_iter == bwt_rec_iter + 1)
		//	{
		//		for (differs_vec_iter = bwt_rec_cur_iter->differs_vec.begin(); differs_vec_iter != bwt_rec_cur_iter->differs_vec.end(); ++differs_vec_iter)
		//		{
		//			if (differs_vec_iter->index >= m_extend_bits)
		//				++other_mismatch;
		//		}
		//	}
		//	else
		//		other_mismatch += bwt_rec_cur_iter->differs_vec.size();
		//}

		//string seg1readstr = "";

		//string seg3readstr = "";

		//vector<GroupedBwtRec>::iterator bwt_mapped_iter;

		//for (bwt_mapped_iter = fix_hole_vec.begin(); bwt_mapped_iter != bwt_rec_iter + 1; ++bwt_mapped_iter)
		//{
		//	seg1readstr.append(bwt_mapped_iter->mapped_chrom_seq);
		//}

		//for (bwt_mapped_iter = bwt_rec_iter + 1; bwt_mapped_iter != fix_hole_vec.end(); ++bwt_mapped_iter)
		//{
		//	seg3readstr.append(bwt_mapped_iter->mapped_chrom_seq);
		//}

		//if (exceed)
		//	;
		//else if (bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length()/*m_seg_len*/ + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
		//	continue;

		//size_t readID, readID1, readID2;

		//WordPair dup_missed_read;//, prefix_suffix_append;

		WordPair dup_tobe_fixed_read, tobe_fixed_read(tobe_fixed_str);//, prefix_suffix_append;

		tobe_fixed_read.duplicate_self(tobe_fixed_len, dup_tobe_fixed_read);

		//if (strand == '+')
		//{
		//readID = tagidx * m_num_seg + missed_id - 1;

		//readID1 = readID - 1;

		//readID2 = readID + 1;

		//WordPair wp1(m_read_words[readID1]);

		//WordPair wp2(m_read_words[readID2]);

		//WordPair prefixappend;

		//wp1.left_shift(m_seg_len);

		//wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words[readID], prefixappend);

		//prefixappend.left_shift(m_extend_bits);

		//wp2.right_shift(m_seg_len - m_extend_bits);

		//prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

		//prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
		//		m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);

		//cout << "readID:"<<readID<<endl;

		//cout << "readID1:"<<readID1<<endl;

		//cout << "readID2:"<<readID2<<endl;

		//cout << "m_read_words[readID]:"<<m_read_words[readID].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//cout << "m_read_words[readID1]:"<<m_read_words[readID1].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//cout << "m_read_words[readID2]:"<<m_read_words[readID2].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//m_read_words[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		//}
		//else
		//{
		//	readID = tagidx * m_num_seg + m_num_seg - missed_id;
		//	m_read_words_rc[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		//}

		//cout << "rid : "<<readID<<endl;

		

		//cout << "combine mapped seg"<<endl;

//#ifdef DEBUG_CHECK
		//for (size_t i = missed_id - 1; i >= 1 ; --i)
		//{
		//	cout << readID - i << endl;
		//	cout << m_reads[readID - i] << endl;
		//	seg1readstr.append(m_reads[readID - i]);
		//}

		//cout << "seg1readstr.length(): " << seg1readstr.length()<<endl;
//#endif

		

//#ifdef DEBUG_CHECK
		//for (size_t i = missed_id; i < m_num_seg ; ++i)
		//{
		//	//cout << readID + i - missed_id + 1 << endl;
		//	//cout << m_reads[readID + i - missed_id + 1] << endl;
		//	seg3readstr.append(m_reads[readID + i - missed_id + 1]);
		//}

		//cout << "seg3readstr.length(): " << seg3readstr.length()<<endl;
//#endif
		string prefixflankstr;
		string suffixflankstr;

		string doner_chromstr;
		string acceptor_chromstr;

		size_t tobe_fixed_len_ext_flank = tobe_fixed_len + 2;

		//cout << "prepare to fix"<<endl;

		if (strand1 == '+' && strand2 == '+')
		{
			if (doner_end + tobe_fixed_len_ext_flank> chromseq1.length() || acceptor_st < tobe_fixed_len_ext_flank)
				continue;

			prefixflankstr = chromseq1.substr(doner_end, tobe_fixed_len_ext_flank);
			suffixflankstr = chromseq2.substr(acceptor_st - tobe_fixed_len_ext_flank + 1, tobe_fixed_len_ext_flank);

			doner_chromstr = chromseq1.substr(doner_end - doner_str.length(), doner_str.length());
			acceptor_chromstr = chromseq2.substr(acceptor_st + 1, acceptor_str.length());
		}
		else if (strand1 == '+' && strand2 == '-')
		{
			if (doner_end + tobe_fixed_len_ext_flank > chromseq1.length() || acceptor_st + tobe_fixed_len_ext_flank > chromseq2.length())
				continue;

			prefixflankstr = chromseq1.substr(doner_end, tobe_fixed_len_ext_flank );
			suffixflankstr = chromseq2.substr(acceptor_st, tobe_fixed_len_ext_flank);

			suffixflankstr = revcomp(suffixflankstr);

			doner_chromstr = chromseq1.substr(doner_end - doner_str.length(), doner_str.length());
			acceptor_chromstr = chromseq2.substr(acceptor_st - acceptor_str.length(), acceptor_str.length());

			acceptor_chromstr = revcomp(acceptor_chromstr);

		}
		else if (strand1 == '-' && strand2 == '-')
		{
			if (doner_end < tobe_fixed_len_ext_flank || acceptor_st + tobe_fixed_len_ext_flank> chromseq2.length())
				continue;
			prefixflankstr = chromseq1.substr(doner_end - tobe_fixed_len_ext_flank + 1, tobe_fixed_len_ext_flank);
			suffixflankstr = chromseq2.substr(acceptor_st, tobe_fixed_len_ext_flank);

			//prefixflankstr = revcomp(prefixflankstr);

			//suffixflankstr = revcomp(suffixflankstr);

			//seg1chromstr = chromseq1.substr(offset1, m_seg_len);
			//seg3chromstr = chromseq2.substr(offset2, m_seg_len);
			doner_chromstr = chromseq1.substr(doner_end + 1, doner_str.length());
			acceptor_chromstr = chromseq2.substr(acceptor_st - acceptor_str.length(), acceptor_str.length());

			doner_chromstr = revcomp(doner_chromstr);
			acceptor_chromstr = revcomp(acceptor_chromstr);
		}
		else
		{
			if (doner_end < tobe_fixed_len_ext_flank || acceptor_st < tobe_fixed_len_ext_flank)
				continue;
			prefixflankstr = chromseq1.substr(doner_end - tobe_fixed_len_ext_flank + 1, tobe_fixed_len_ext_flank);
			suffixflankstr = chromseq2.substr(acceptor_st - tobe_fixed_len_ext_flank + 1, tobe_fixed_len_ext_flank);

			//prefixflankstr = revcomp(prefixflankstr);

			//seg1chromstr = chromseq1.substr(offset1, m_seg_len);
			//seg3chromstr = chromseq2.substr(offset2, m_seg_len);
			doner_chromstr = chromseq1.substr(doner_end + 1, doner_str.length());
			acceptor_chromstr = chromseq2.substr(acceptor_st + 1, acceptor_str.length());

			doner_chromstr = revcomp(doner_chromstr);
		}

		//cout << "extracted sequence "<<endl;

		//cout << m_reads[readID] << endl;
		//cout << prefixflankstr << endl;
		//cout << suffixflankstr << endl;
		//cout << seg1chromstr << endl;
		//cout << seg1readstr << endl;
		//cout << seg3chromstr << endl;
		//cout << seg3readstr << endl;

#ifdef DEBUG_CHECK
		if (seg1readstr.length() != seg1chromstr.length() || seg3readstr.length() != seg3chromstr.length())
		{
			cout << "read length does not match chrom length "<<endl;
			cout <<seg1readstr<<endl;
			cout <<seg1chromstr<<endl; 
			cout <<seg3readstr<<endl;
			cout <<seg3chromstr<<endl; 
		}

		for (size_t i = 0; i < seg1readstr.length(); ++i)
		{
			if (seg1readstr[i] >= 'a' && seg1readstr[i] <= 'z')
				seg1readstr[i] = seg1readstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg1chromstr.length(); ++i)
		{
			if (seg1chromstr[i] >= 'a' && seg1chromstr[i] <= 'z')
				seg1chromstr[i] = seg1chromstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg3readstr.length(); ++i)
		{
			if (seg3readstr[i] >= 'a' && seg3readstr[i] <= 'z')
				seg3readstr[i] = seg3readstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg3chromstr.length(); ++i)
		{
			if (seg3chromstr[i] >= 'a' && seg3chromstr[i] <= 'z')
				seg3chromstr[i] = seg3chromstr[i] + 'A' - 'a';
		}

		size_t count1 = 0;
		for (size_t i = 0; i < seg3readstr.length(); ++i)
		{
			if (seg3readstr[i] != seg3chromstr[i])
				++count1;
		}

		if (count1 > m_max_mismatches * seg3readstr.length() / m_seg_len)
		{
			cout << "seg 3 mismatch exceed "<<endl;
			cout << seg3readstr<<endl;
			cout << seg3chromstr << endl;

			getchar();
		}

		size_t count2 = 0;

		for (size_t i = 0; i < seg1readstr.length(); ++i)
		{
			if (seg1readstr[i] != seg1chromstr[i])
				++count2;
		}

		if (count2 > m_max_mismatches * seg1readstr.length() / m_seg_len)
		{
			cout << "seg 1 mismatch exceed "<<endl;
			cout << seg1readstr<<endl;
			cout << seg1chromstr << endl;

			getchar();
		}
#endif
		//getchar();

		m_five_prim_suffix = WordPair(prefixflankstr);
		m_three_prim_prefix = WordPair(suffixflankstr);

		WordPair prefix_chrom = m_five_prim_suffix;

		m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

		prefix_chrom.left_shift(tobe_fixed_len - 2);//m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

		WordPair comb_chrom_seq;

		//cout <<"combine chrom"<<endl;
		prefix_chrom.ps_combine(cur_masks_ptr->prefix_seg_bits_on, cur_masks_ptr->suffix_seg_bits_on, cur_masks_ptr->comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);
		//                     (m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

		size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;
		m_matched_flank = 0;
		m_matched_bads = 0;

		//cout << "scoring"<<endl;

		//cout << "chromwp " << comb_chrom_seq.tostring2(m_masks.comb_seg_bits_on_ext)<<endl;
		//cout << "readswp " << dup_missed_read.tostring2(m_masks.comb_seg_bits_on_ext)<<endl;

		size_t score = Fixhole_score_selective_var_mask(dup_tobe_fixed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits, cur_masks_ptr);
					 //Fixhole_score_selective_ext(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits);
					 //Fixfusion_score(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);

		//cout << "score:" << score<<endl;
		//cout << "max_loc:" << max_loc<<endl;
		//cout << "prim:" << prim<<endl;

		//getchar();
		//size_t total_mismatch = other_mismatch + score;

		//cout <<  "total score:" <<total_mismatch << endl;

		//getchar();

		//cout << "flank case: "<<flankcase << '\t' << "prim: " << prim << endl;
		if (score <= m_max_mismatches && prim >= flankcase /*(!flankcase || (flankcase && prim))*//* && max_loc <= m_seg_len + (2 * m_extend_bits) && total_mismatch <= m_max_fusion_mismatch*/)
		{
			//if (mark)
			//{
			//	++m_cur_reads_count[tagidx];
			//}
			//int strandi = 0;
			//if (strand == '-')
			//	strandi = 16;

			//string prim5flankstr = prefixflankstr;
			//string prim3flankstr = suffixflankstr;

			//string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len + (2 * m_extend_bits) - max_loc), m_seg_len + (2 * m_extend_bits) - max_loc);

			//comb_mapped_chromseq = comb_mapped_chromseq.substr(m_extend_bits, comb_mapped_chromseq.length() - m_extend_bits);

			//cout <<"comb_mapped_chromseq:"<<comb_mapped_chromseq<<endl; 

			string flankstr;

			FlankString(flankstr);	

			//cout << "flankstr:" << flankstr<<endl;

			//if (prim == 0)
			//	getchar();
			
			//<<strand1<<strand2<<'\t'<< chr1<< '_'<<chr2<<'\t';

			//revcomp flank string

			//string flankstr1 = flankstr.substr(0, 2);

			//if (strand1 == '-')
			//	flankstr1 = revcomp(flankstr1);

			//string flankstr2 = flankstr.substr(2, 2);

			//if (strand2 == '-')
			//	flankstr2 = revcomp(flankstr2);

			//flankstr = flankstr1 + flankstr2;

			size_t suffixst, prefixend;

			string prefixsyn, suffixsyn;

			size_t synlen = 60;

			DifferVec differs_vec(differstr);

			unsigned short total_mismatch = (unsigned short) (differs_vec.differs_vec.size() + score);

			if (strand1 == '+' && strand2 == '+')
			{
				prefixend = doner_end + max_loc - 1;
				suffixst = acceptor_st - tobe_fixed_len + max_loc + 1;

				if (prefixend < synlen || suffixst + synlen >= chromseq2.length())
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend <<'\t' <<suffixst + synlen<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - synlen + 1, synlen);
				suffixsyn = chromseq2.substr(suffixst, synlen);

				//ofs << offset + m_seg_len + max_loc << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;
				ofs <<tagname<<'\t'<< chromname1 << '\t' << strand1 << '\t' << prefixend - max_loc - doner_str.length() + 1 << '\t' <<prefixend 
					<<'\t'<<chromname2 << '\t' << strand2 << '\t' << suffixst << '\t' << suffixst + tobe_fixed_len - max_loc + acceptor_str.length() - 1
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr <<'\t' <<  "NM:i:"<<total_mismatch <<endl;

			}
			else if (strand1 == '+' && strand2 == '-')
			{
				prefixend = doner_end + max_loc - 1;
				suffixst = acceptor_st + tobe_fixed_len - max_loc - 1;
				//ofs << offset + m_seg_len + max_loc << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

				if (prefixend < synlen || suffixst < synlen)
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend <<'\t' <<suffixst<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - synlen + 1, synlen);
				suffixsyn = chromseq2.substr(suffixst - synlen + 1, synlen);

				suffixsyn = revcomp(suffixsyn);

				ofs <<tagname<<'\t'<< chromname1 << '\t' << strand1 << '\t' << prefixend - max_loc - doner_str.length() + 1 << '\t' <<prefixend 
					<<'\t'<<chromname2 << '\t' << strand2 << '\t' << suffixst - (tobe_fixed_len - max_loc + acceptor_str.length() - 1) << '\t' << suffixst
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/  << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;
				
			}
			else if (strand1 == '-' && strand2 == '-')
			{
				prefixend = doner_end - max_loc + 1;
				suffixst = acceptor_st + tobe_fixed_len - max_loc + 1;

				if (prefixend + synlen >= chromseq1.length() || suffixst < synlen)
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend + synlen <<'\t' <<suffixst<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend, synlen);
				suffixsyn = chromseq2.substr(suffixst - synlen + 1, synlen);

				prefixsyn = revcomp(prefixsyn);
				suffixsyn = revcomp(suffixsyn);

				//ofs << offset - max_loc + 1 << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;
				ofs <<tagname<<'\t'<< chromname1 << '\t' << strand1 << '\t' << prefixend << '\t' <<prefixend + max_loc + doner_str.length() - 1 
					<<'\t'<<chromname2 << '\t' << strand2 << '\t' << suffixst - (tobe_fixed_len - max_loc + acceptor_str.length() - 1) << '\t' << suffixst
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;

			}
			else
			{
				prefixend  = doner_end - max_loc + 1;
				suffixst = acceptor_st - tobe_fixed_len + max_loc + 1;

				if (prefixend + synlen >= chromseq1.length() || suffixst + synlen >= chromseq2.length())
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend + synlen <<'\t' <<suffixst + synlen<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend, synlen);
				suffixsyn = chromseq2.substr(suffixst, synlen);

				prefixsyn = revcomp(prefixsyn);

				//ofs << offset - max_loc + 1 << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;
				ofs <<tagname<<'\t' << chromname1 << '\t' << strand1 << '\t' << prefixend << '\t' <<prefixend + max_loc + doner_str.length() - 1 
					<<'\t'<<chromname2 << '\t' << strand2 << '\t' << suffixst << '\t' << suffixst + tobe_fixed_len - max_loc + acceptor_str.length() - 1
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;

			}

			//ofs << '\t' <<flankstr <<endl;

			//cout <<"insert junction"<<endl;

			string chrom_ID = chromname1;
			chrom_ID.append("_");
			chrom_ID.append(chromname2);

			size_t prefixlen = max_loc + doner_str.length(), suffixlen = tobe_fixed_len - max_loc + acceptor_str.length();

			//cout << "find(suffixst) "<<endl;

			size_t combined_offset = (prefixend << THIRTY_TWO) + suffixst;

			CHROM_FUSION_JUNC_HASH_COMB_ITER chrom_junc_hash_iter = conj_hash_map.find(chrom_ID);

			if (chrom_junc_hash_iter == conj_hash_map.end())
			{
				FUSION_JUNC_HASH_COMB junc_hash_comb;

				chrom_junc_hash_iter = (conj_hash_map.insert(CHROM_FUSION_JUNC_HASH_COMB::value_type(chrom_ID, junc_hash_comb))).first;							
			}

			FUSION_JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

			FUSION_JUNC_HASH_COMB_ITER junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

			if (junc_hash_comb_iter != junc_hash_comb.end())
			{
				junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, total_mismatch);
			}
			else
			{
				string synstr = prefixsyn;
				synstr.append(suffixsyn);
				
				//string strand12 = "++";
				//strand12[0] = strand1;
				//strand12[1] = strand2;

				bool isbads = m_matched_bads ? true : false;

				unsigned short pm = static_cast<unsigned short>(prim);

				chrom_junc_hash_iter->second.insert(FUSION_JUNC_HASH_COMB::value_type(combined_offset, JuncInfoFusion(prefixlen, suffixlen, m_read_width, synstr, strand1, strand2, total_mismatch, isbads, m_matched_flank, pm)));
			}

			//JUNCFusion_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

			//if (iter_conj == conj_hash_map[chrom_ID].end())
			//{
			//	//cout << "insert(JUNC_HASH: "<<endl;
			//	JUNCFusion_SUFFIX match_hash_temp;
			//	conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
			//	iter_conj = conj_hash_map[chrom_ID].find(suffixst);
			//}

			//JUNCFusion_SUFFIX_ITER iter_conj_pre;

			//
			//if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
			//{
			//	//if (iter_conj_pre->second.flankstr != flankstr)
			//	//{
			//	//	cout <<"flank does not match"<<endl;
			//	//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
			//	//}
			//	//cout << "inc_hits(max_loc) "<<endl;
			//	iter_conj_pre->second.inc_hits(max_loc + seg1chromstr.length() - m_extend_bits, (unsigned short)total_mismatch);	
			//}	  
			//else 
			//{
			//	//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
			//	//cout << "insert(JUNC_SUFFIX "<<endl;
			//	//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

			//	string synstr = prefixsyn;
			//	synstr.append(suffixsyn);
			//	
			//	string strand12 = "++";
			//	strand12[0] = strand1;
			//	strand12[1] = strand2;
			//	(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankstr, max_loc + seg1chromstr.length() - m_extend_bits, m_read_width, synstr, strand12, (unsigned short)total_mismatch)));
			//}
		}

	}

	input_fs.close();
	ofs.close();

	//while(!input_fs.eof())
	//{
	//	//cout <<"read number"<<endl;
	//	//if (rec_num_str == "")
	//	//	continue;
	//	//char *endptr;
	//	//size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);
	//	string line1, line2;

	//	getline(input_fs, line1);
	//	getline(input_fs, line2);

	//	if (line1 == "" || line2 == "")
	//		continue;

	//	vector<GroupedBwtRec> fix_hole_vec;

	//	char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
	//	int offset, seg_id, other_mapped;

	//	sscanf(line1.c_str(), "%s\t%d\t%c\t%s\t%d", tagname, &seg_id, &strand, chromname, &offset);
	//	fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, "0", 0, "0"));

	//	char chromname2[1000], tagname2[1000], mapped_seq2[1000], differs2[1000], strand2;
	//	int offset2, seg_id2, other_mapped2;
	//	sscanf(line2.c_str(), "%s\t%d\t%c\t%s\t%d", tagname2, &seg_id2, &strand2, chromname2, &offset2);
	//	fix_hole_vec.push_back(GroupedBwtRec(offset2, seg_id2, "0", 0, "0"));
	//	//}

	//	//string tagnamestr = tagname;
	//	//size_t last_idx = tagnamestr.find_last_of("_");
	//	//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);

	//	if (m_reads_name_map_id.find(tagname) == m_reads_name_map_id.end())
	//	{
	//		cout << "can't find tagname "<<endl<< line1<<endl<<line2<<endl;
	//		continue;
	//	}

	//	string chr1 = chromname;

	//	string chr2 = chromname2;

	//	if (chr1 == chr2)
	//		return;

	//	size_t tagidx = m_reads_name_map_id[tagname] + 1;

	//	size_t missed_id = 2;

	//	//vector<GroupedBwtRec>::iterator bwt_rec_iter;

	//	//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
	//	//{
	//	//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
	//	//	{
	//	//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
	//	//		{
	//	//			cout << "missed more than 1 segments"<<endl;
	//	//			exit(1);
	//	//		}
	//	//		else
	//	//		{
	//	//			missed_id = bwt_rec_iter->seg_id + 1;
	//	//			break;
	//	//		}
	//	//	}
	//	//}

	//	//if (bwt_rec_iter->offset + m_seg_len + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
	//	//	continue;

	//	size_t readID;

	//	WordPair dup_missed_read;

	//	readID = tagidx;
	//	m_read_words[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);

	//	string seg1readstr = m_reads[readID - 1];
	//	string seg3readstr = m_reads[readID + 1];

	//	string chromfile1 = m_chrom_files[m_chrom_names_map[chromname]];
	//	string chromfile2 = m_chrom_files[m_chrom_names_map[chromname2]];
	//	
	//	if (prevchrom_file1 == "" || chromfile1 != prevchrom_file1)
	//	{
	//		chromseq1.clear();
	//		//cout << "read chromo\t"<<chromfile<<endl;
	//		readchrom(chromfile1.c_str(), chromseq1);
	//		if (chromseq1 == "")
	//		{
	//			cout <<"empty chrom: "<<chromfile1<<endl;
	//			exit(1);
	//		}
	//		prevchrom_file1 = chromfile1;
	//	}

	//	if (prevchrom_file2 == "" || chromfile2 != prevchrom_file2)
	//	{
	//		chromseq2.clear();
	//		//cout << "read chromo\t"<<chromfile<<endl;
	//		readchrom(chromfile2.c_str(), chromseq2);
	//		if (chromseq2 == "")
	//		{
	//			cout <<"empty chrom: "<<chromfile2<<endl;
	//			exit(1);
	//		}
	//		prevchrom_file2 = chromfile2;
	//	}
	//	

	//	//cout << "get substring"<<endl;
	//	string prefixflankstr;// = chromseq.substr(bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length(), m_seg_len + 2);
	//	string suffixflankstr;// = chromseq.substr((bwt_rec_iter + 1)->offset - m_seg_len - 2, m_seg_len + 2);

	//	string seg1chromstr;
	//	string seg3chromstr;

	//	if (strand == '+' && strand2 == '+')
	//	{
	//		prefixflankstr = chromseq1.substr(offset + m_seg_len, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2, m_seg_len + 2);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);
	//	}
	//	else if (strand == '+' && strand2 == '-')
	//	{
	//		prefixflankstr = chromseq1.substr(offset + m_seg_len, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 + m_seg_len, m_seg_len + 2);

	//		suffixflankstr = revcomp(suffixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg3chromstr = revcomp(seg3chromstr);

	//	}
	//	else if (strand == '-' && strand2 == '-')
	//	{
	//		prefixflankstr = chromseq1.substr(offset - m_seg_len - 2, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 + m_seg_len, m_seg_len + 2);

	//		prefixflankstr = revcomp(prefixflankstr);

	//		suffixflankstr = revcomp(suffixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg1chromstr = revcomp(seg1chromstr);
	//		seg3chromstr = revcomp(seg3chromstr);
	//	}
	//	else
	//	{
	//		prefixflankstr = chromseq1.substr(offset - m_seg_len - 2, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2, m_seg_len + 2);

	//		prefixflankstr = revcomp(prefixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg1chromstr = revcomp(seg1chromstr);
	//	}

	//	if (seg1readstr.length() != seg1chromstr.length() || seg3readstr.length() != seg3chromstr.length())
	//	{
	//		cout << "read length does not match chrom length "<<endl;
	//		cout <<seg1readstr<<endl;
	//		cout <<seg1chromstr<<endl; 
	//		cout <<seg3readstr<<endl;
	//		cout <<seg3chromstr<<endl; 
	//	}

	//	for (size_t i = 0; i < seg1readstr.length(); ++i)
	//	{
	//		if (seg1readstr[i] >= 'a' && seg1readstr[i] <= 'z')
	//			seg1readstr[i] = seg1readstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg1chromstr.length(); ++i)
	//	{
	//		if (seg1chromstr[i] >= 'a' && seg1chromstr[i] <= 'z')
	//			seg1chromstr[i] = seg1chromstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg3readstr.length(); ++i)
	//	{
	//		if (seg3readstr[i] >= 'a' && seg3readstr[i] <= 'z')
	//			seg3readstr[i] = seg3readstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg3chromstr.length(); ++i)
	//	{
	//		if (seg3chromstr[i] >= 'a' && seg3chromstr[i] <= 'z')
	//			seg3chromstr[i] = seg3chromstr[i] + 'A' - 'a';
	//	}

	//	size_t count = 0;
	//	for (size_t i = 0; i < seg3readstr.length(); ++i)
	//	{
	//		if (seg3readstr[i] != seg3chromstr[i])
	//			++count;
	//	}

	//	if (count > m_max_mismatches)
	//	{
	//		cout << "seg 3 mismatch exceed "<<endl;
	//		cout << seg3readstr<<endl;
	//		cout << seg3chromstr << endl;

	//		getchar();
	//	}

	//	size_t count2 = 0;

	//	for (size_t i = 0; i < seg1readstr.length(); ++i)
	//	{
	//		if (seg1readstr[i] != seg1chromstr[i])
	//			++count2;
	//	}

	//	if (count2 > m_max_mismatches)
	//	{
	//		cout << "seg 1 mismatch exceed "<<endl;
	//		cout << seg1readstr<<endl;
	//		cout << seg1chromstr << endl;

	//		getchar();
	//	}

	//	m_five_prim_suffix = WordPair(prefixflankstr);
	//	m_three_prim_prefix = WordPair(suffixflankstr);

	//	WordPair prefix_chrom = m_five_prim_suffix;

	//	m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

	//	prefix_chrom.left_shift(m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

	//	WordPair comb_chrom_seq;

	//	//cout <<"combine chrom"<<endl;
	//	prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	//	size_t max_loc, prim, left_mismatch = m_max_mismatches;
	//	m_matched_flank = 0;
	//	m_matched_bads = 0;

	//	//cout << "scoring"<<endl;
	//	size_t score = Fixhole_score_selective(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);

	//	//cout << "5 prim and 3 prim"<<endl;
	//	//cout << prefixflankstr << endl;
	//	//cout << suffixflankstr << endl;

	//	//cout <<"score: "<<score<<endl;
	//	//cout <<"max_loc: "<<max_loc<<endl;
	//	if (score <= m_max_mismatches && max_loc <= m_seg_len)
	//	{
	//		int strandi = 0;
	//		if (strand == '-')
	//			strandi = 16;

	//		string prim5flankstr = prefixflankstr;
	//		string prim3flankstr = suffixflankstr;
	//		string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

	//		string flankstr;

	//		FlankString(flankstr);

	//		ofs<<tagname<<'\t'<<strand<<strand2<<'\t'<< chromname<< '_'<<chromname2<<'\t';//<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t';

	//		size_t suffixst, prefixend;

	//		string prefixsyn, suffixsyn;

	//		if (strand == '+' && strand2 == '+')
	//		{
	//			prefixend = offset + m_seg_len + max_loc;
	//			suffixst = offset2 - m_seg_len + max_loc + 1;
	//			ofs << offset + m_seg_len + max_loc << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 30, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 1, 30);
	//		}
	//		else if (strand == '+' && strand2 == '-')
	//		{
	//			prefixend = offset + m_seg_len + max_loc;
	//			suffixst = offset2 + m_seg_len + m_seg_len - max_loc;
	//			ofs << offset + m_seg_len + max_loc << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 30, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 30, 30);

	//			suffixsyn = revcomp(suffixsyn);
	//		}
	//		else if (strand == '-' && strand2 == '-')
	//		{
	//			prefixend = offset - max_loc + 1;
	//			suffixst = offset2 + m_seg_len + m_seg_len - max_loc;
	//			ofs << offset - max_loc + 1 << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 1, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 30, 30);

	//			prefixsyn = revcomp(prefixsyn);
	//			suffixsyn = revcomp(suffixsyn);
	//		}
	//		else
	//		{
	//			prefixend  = offset - max_loc + 1;
	//			suffixst = offset2 - m_seg_len + max_loc + 1;
	//			ofs << offset - max_loc + 1 << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 1, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 1, 30);

	//			prefixsyn = revcomp(prefixsyn);
	//		}

	//		ofs << '\t' <<flankstr <<endl;

	//		string chrom_ID = chromname;
	//		chrom_ID.append("_");
	//		chrom_ID.append(chromname2);

	//		//cout << "find(suffixst) "<<endl;
	//		JUNCFusion_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

	//		if (iter_conj == conj_hash_map[chrom_ID].end())
	//		{
	//			//cout << "insert(JUNC_HASH: "<<endl;
	//			JUNCFusion_SUFFIX match_hash_temp;
	//			conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
	//			iter_conj = conj_hash_map[chrom_ID].find(suffixst);
	//		}

	//		JUNCFusion_SUFFIX_ITER iter_conj_pre;

	//		
	//		if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
	//		{
	//			//if (iter_conj_pre->second.flankstr != flankstr)
	//			//{
	//			//	cout <<"flank does not match"<<endl;
	//			//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
	//			//}
	//			//cout << "inc_hits(max_loc) "<<endl;
	//			iter_conj_pre->second.inc_hits(max_loc);	
	//		}	  
	//		else 
	//		{
	//			//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
	//			//cout << "insert(JUNC_SUFFIX "<<endl;
	//			//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

	//			string synstr = prefixsyn;
	//			synstr.append(suffixsyn);
	//			
	//			string strand12 = "++";
	//			strand12[0] = strand;
	//			strand12[1] = strand2;
	//			(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankstr, max_loc, m_read_width, synstr, strand12)));
	//		}
	//	}
	//}
	
}

void
GenomeScan::FixHole_fusion(const char* fixhole_file, const char* fixed_hole_file, CONJFusion_HASH_MAP& conj_hash_map, bool avoid, bool mark)
{
	ifstream input_fs(fixhole_file);

	if( !input_fs ) 
	{
		cout <<"error: open file failed: "<<fixhole_file<<"\n";return;//exit(1);
	}

	ofstream ofs(fixed_hole_file);

	string chromseq1 = "";
	string chromseq2 = "";
	string prevchrom_file1 = "";
	string prevchrom_file2 = "";

	string rec_num_str;
	while(getline(input_fs, rec_num_str))
	{
		//cout <<"read number"<<endl;
		if (rec_num_str == "")
			continue;
		char *endptr;
		size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);

		vector<GroupedBwtRec> fix_hole_vec;

		char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
		int offset, seg_id, other_mapped;

		string line;

		//cout <<"read line"<<endl;
		for (size_t i = 0; i < rec_num; ++i)
		{
			getline(input_fs, line);

			if (line == "")
				continue;

			sscanf(line.c_str(), "%s\t%s\t%c\t%d\t%d\t%s\t%d\t%s", tagname, chromname, &strand, &offset, &seg_id, mapped_seq, &other_mapped, differs);

			//cout << line << endl;
			fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, mapped_seq, other_mapped, differs, chromname, strand));
		}

		string tagnamestr = tagname;
		size_t last_idx = tagnamestr.find_last_of("_");
		string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
		size_t tagidx = atoi(tagidxstr.c_str()) - 1;

		if (avoid && m_reads_count[tagidx])
			continue;

		size_t missed_id;

		vector<GroupedBwtRec>::iterator bwt_rec_iter;

		string chr1, chr2;

		char strand1, strand2;

		size_t offset1, offset2;

		//cout << "find missed seg"<<endl;

		for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
		{
			if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
			{
				if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
				{
					cout << "missed more than 1 segments"<<endl;

					for (size_t i = 0; i < fix_hole_vec.size(); ++i)
					{
						cout << tagname <<'\t' << fix_hole_vec[i].offset << '\t' <<fix_hole_vec[i].seg_id <<endl;
					}
					exit(1);
				}
				else
				{
					chr1 = bwt_rec_iter->chromname;

					chr2 = (bwt_rec_iter + 1)->chromname;

					strand1 =  bwt_rec_iter->strand;

					strand2 = (bwt_rec_iter + 1)->strand;

					offset1 = bwt_rec_iter->offset;

					offset2 = (bwt_rec_iter + 1)->offset;

					missed_id = bwt_rec_iter->seg_id + 1;
					break;
				}
			}
		}

		vector<GroupedBwtRec>::iterator bwt_rec_cur_iter;

		size_t other_mismatch = 0;

		for (bwt_rec_cur_iter = fix_hole_vec.begin(); bwt_rec_cur_iter != fix_hole_vec.end(); ++bwt_rec_cur_iter)
		{
			vector<Differ>::iterator differs_vec_iter;
			if (bwt_rec_cur_iter == bwt_rec_iter)
			{
				for (differs_vec_iter = bwt_rec_cur_iter->differs_vec.begin(); differs_vec_iter != bwt_rec_cur_iter->differs_vec.end(); ++differs_vec_iter)
				{
					if (differs_vec_iter->index < m_seg_len - m_extend_bits)
						++other_mismatch;
				}
			}
			else if (bwt_rec_cur_iter == bwt_rec_iter + 1)
			{
				for (differs_vec_iter = bwt_rec_cur_iter->differs_vec.begin(); differs_vec_iter != bwt_rec_cur_iter->differs_vec.end(); ++differs_vec_iter)
				{
					if (differs_vec_iter->index >= m_extend_bits)
						++other_mismatch;
				}
			}
			else
				other_mismatch += bwt_rec_cur_iter->differs_vec.size();
		}

		string seg1readstr = "";

		string seg3readstr = "";

		vector<GroupedBwtRec>::iterator bwt_mapped_iter;

		for (bwt_mapped_iter = fix_hole_vec.begin(); bwt_mapped_iter != bwt_rec_iter + 1; ++bwt_mapped_iter)
		{
			seg1readstr.append(bwt_mapped_iter->mapped_chrom_seq);
		}

		for (bwt_mapped_iter = bwt_rec_iter + 1; bwt_mapped_iter != fix_hole_vec.end(); ++bwt_mapped_iter)
		{
			seg3readstr.append(bwt_mapped_iter->mapped_chrom_seq);
		}

		//if (exceed)
		//	;
		//else if (bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length()/*m_seg_len*/ + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
		//	continue;

		size_t readID, readID1, readID2;

		WordPair dup_missed_read, prefix_suffix_append;

		//if (strand == '+')
		//{
		readID = tagidx * m_num_seg + missed_id - 1;

		readID1 = readID - 1;

		readID2 = readID + 1;

		WordPair wp1(m_read_words[readID1]);

		WordPair wp2(m_read_words[readID2]);

		WordPair prefixappend;

		wp1.left_shift(m_seg_len);

		wp1.ps_combine(m_masks.prefix_ext_mask, m_masks.suffix_seg_bits_on, m_masks.prefix_append_mask, m_read_words[readID], prefixappend);

		prefixappend.left_shift(m_extend_bits);

		wp2.right_shift(m_seg_len - m_extend_bits);

		prefixappend.ps_combine(m_masks.prefix_append_shifted_mask, m_masks.suffix_ext_mask, m_masks.prefix_suffix_append_mask, wp2, prefix_suffix_append);

		prefix_suffix_append.duplicate(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, 
				m_masks.left_shift_seg_width_ext, m_masks.right_shift_seg_width_ext, dup_missed_read);

		//cout << "readID:"<<readID<<endl;

		//cout << "readID1:"<<readID1<<endl;

		//cout << "readID2:"<<readID2<<endl;

		//cout << "m_read_words[readID]:"<<m_read_words[readID].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//cout << "m_read_words[readID1]:"<<m_read_words[readID1].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//cout << "m_read_words[readID2]:"<<m_read_words[readID2].tostring3(m_masks.comb_seg_bits_on_ext, m_seg_len) << endl;

		//m_read_words[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		//}
		//else
		//{
		//	readID = tagidx * m_num_seg + m_num_seg - missed_id;
		//	m_read_words_rc[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);
		//}

		//cout << "rid : "<<readID<<endl;

		

		//cout << "combine mapped seg"<<endl;

//#ifdef DEBUG_CHECK
		//for (size_t i = missed_id - 1; i >= 1 ; --i)
		//{
		//	cout << readID - i << endl;
		//	cout << m_reads[readID - i] << endl;
		//	seg1readstr.append(m_reads[readID - i]);
		//}

		//cout << "seg1readstr.length(): " << seg1readstr.length()<<endl;
//#endif

		

//#ifdef DEBUG_CHECK
		//for (size_t i = missed_id; i < m_num_seg ; ++i)
		//{
		//	//cout << readID + i - missed_id + 1 << endl;
		//	//cout << m_reads[readID + i - missed_id + 1] << endl;
		//	seg3readstr.append(m_reads[readID + i - missed_id + 1]);
		//}

		//cout << "seg3readstr.length(): " << seg3readstr.length()<<endl;
//#endif
		if (m_chrom_names_map.find(chr1) == m_chrom_names_map.end())
		{
			cout <<"can't find chrom " << chr1 << endl;
			continue;
		}

		if (m_chrom_names_map.find(chr2) == m_chrom_names_map.end())
		{
			cout <<"can't find chrom " << chr2 << endl;
			continue;
		}

		string chromfile1 = m_chrom_files[m_chrom_names_map[chr1]];
		string chromfile2 = m_chrom_files[m_chrom_names_map[chr2]];
		
		//cout << "read chrom"<<endl;
		
		if (prevchrom_file1 == "" || chromfile1 != prevchrom_file1)
		{
			chromseq1.clear();
			cout << "read chromo\t"<<chromfile1<<endl;
			readchrom(chromfile1.c_str(), chromseq1);
			if (chromseq1 == "")
			{
				cout <<"empty chrom: "<<chromfile1<<endl;
				exit(1);
			}
			prevchrom_file1 = chromfile1;
		}

		if (prevchrom_file2 == "" || chromfile2 != prevchrom_file2)
		{
			chromseq2.clear();
			cout << "read chromo\t"<<chromfile2<<endl;
			readchrom(chromfile2.c_str(), chromseq2);
			if (chromseq2 == "")
			{
				cout <<"empty chrom: "<<chromfile2<<endl;
				exit(1);
			}
			prevchrom_file2 = chromfile2;
		}

		string prefixflankstr;
		string suffixflankstr;

		string seg1chromstr;
		string seg3chromstr;

		//cout << "prepare to fix"<<endl;

		if (strand1 == '+' && strand2 == '+')
		{
			if (offset1 + m_seg_len + m_seg_len + 2 + m_extend_bits > chromseq1.length() || offset2 < m_seg_len + 2 + m_extend_bits)
				continue;
			prefixflankstr = chromseq1.substr(offset1 + m_seg_len - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2) );
			suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2 - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));

			seg1chromstr = chromseq1.substr(offset1 + m_seg_len - seg1readstr.length(), seg1readstr.length());
			seg3chromstr = chromseq2.substr(offset2, seg3readstr.length());
		}
		else if (strand1 == '+' && strand2 == '-')
		{
			if (offset1 + m_seg_len + m_seg_len + 2 + m_extend_bits > chromseq1.length() || offset2 + m_seg_len + m_seg_len + 2 + m_extend_bits > chromseq2.length())
				continue;
			prefixflankstr = chromseq1.substr(offset1 + m_seg_len - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));
			suffixflankstr = chromseq2.substr(offset2 + m_seg_len - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));

			suffixflankstr = revcomp(suffixflankstr);

			seg1chromstr = chromseq1.substr(offset1 + m_seg_len - seg1readstr.length(), seg1readstr.length());
			seg3chromstr = chromseq2.substr(offset2 + m_seg_len - seg3readstr.length(), seg3readstr.length());

			seg3chromstr = revcomp(seg3chromstr);

		}
		else if (strand1 == '-' && strand2 == '-')
		{
			if (offset1 < m_seg_len + 2 + m_extend_bits || offset2 + m_seg_len + m_seg_len + 2 + m_extend_bits> chromseq2.length())
				continue;
			prefixflankstr = chromseq1.substr(offset1 - m_seg_len - 2 - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));
			suffixflankstr = chromseq2.substr(offset2 + m_seg_len - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));

			prefixflankstr = revcomp(prefixflankstr);

			suffixflankstr = revcomp(suffixflankstr);

			//seg1chromstr = chromseq1.substr(offset1, m_seg_len);
			//seg3chromstr = chromseq2.substr(offset2, m_seg_len);
			seg1chromstr = chromseq1.substr(offset1, seg1readstr.length());
			seg3chromstr = chromseq2.substr(offset2 + m_seg_len - seg3readstr.length(), seg3readstr.length());

			seg1chromstr = revcomp(seg1chromstr);
			seg3chromstr = revcomp(seg3chromstr);
		}
		else
		{
			if (offset1 < m_seg_len + 2 + m_extend_bits || offset2 < m_seg_len + 2 + m_extend_bits)
				continue;
			prefixflankstr = chromseq1.substr(offset1 - m_seg_len - 2 - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));
			suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2 - m_extend_bits, m_seg_len + 2 + (m_extend_bits * 2));

			prefixflankstr = revcomp(prefixflankstr);

			//seg1chromstr = chromseq1.substr(offset1, m_seg_len);
			//seg3chromstr = chromseq2.substr(offset2, m_seg_len);
			seg1chromstr = chromseq1.substr(offset1, seg1readstr.length());
			seg3chromstr = chromseq2.substr(offset2, seg3readstr.length());

			seg1chromstr = revcomp(seg1chromstr);
		}

		//cout << "extracted sequence "<<endl;

		//cout << m_reads[readID] << endl;
		//cout << prefixflankstr << endl;
		//cout << suffixflankstr << endl;
		//cout << seg1chromstr << endl;
		//cout << seg1readstr << endl;
		//cout << seg3chromstr << endl;
		//cout << seg3readstr << endl;

#ifdef DEBUG_CHECK
		if (seg1readstr.length() != seg1chromstr.length() || seg3readstr.length() != seg3chromstr.length())
		{
			cout << "read length does not match chrom length "<<endl;
			cout <<seg1readstr<<endl;
			cout <<seg1chromstr<<endl; 
			cout <<seg3readstr<<endl;
			cout <<seg3chromstr<<endl; 
		}

		for (size_t i = 0; i < seg1readstr.length(); ++i)
		{
			if (seg1readstr[i] >= 'a' && seg1readstr[i] <= 'z')
				seg1readstr[i] = seg1readstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg1chromstr.length(); ++i)
		{
			if (seg1chromstr[i] >= 'a' && seg1chromstr[i] <= 'z')
				seg1chromstr[i] = seg1chromstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg3readstr.length(); ++i)
		{
			if (seg3readstr[i] >= 'a' && seg3readstr[i] <= 'z')
				seg3readstr[i] = seg3readstr[i] + 'A' - 'a';
		}

		for (size_t i = 0; i < seg3chromstr.length(); ++i)
		{
			if (seg3chromstr[i] >= 'a' && seg3chromstr[i] <= 'z')
				seg3chromstr[i] = seg3chromstr[i] + 'A' - 'a';
		}

		size_t count1 = 0;
		for (size_t i = 0; i < seg3readstr.length(); ++i)
		{
			if (seg3readstr[i] != seg3chromstr[i])
				++count1;
		}

		if (count1 > m_max_mismatches * seg3readstr.length() / m_seg_len)
		{
			cout << "seg 3 mismatch exceed "<<endl;
			cout << seg3readstr<<endl;
			cout << seg3chromstr << endl;

			getchar();
		}

		size_t count2 = 0;

		for (size_t i = 0; i < seg1readstr.length(); ++i)
		{
			if (seg1readstr[i] != seg1chromstr[i])
				++count2;
		}

		if (count2 > m_max_mismatches * seg1readstr.length() / m_seg_len)
		{
			cout << "seg 1 mismatch exceed "<<endl;
			cout << seg1readstr<<endl;
			cout << seg1chromstr << endl;

			getchar();
		}
#endif
		//getchar();

		m_five_prim_suffix = WordPair(prefixflankstr);
		m_three_prim_prefix = WordPair(suffixflankstr);

		WordPair prefix_chrom = m_five_prim_suffix;

		m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

		prefix_chrom.left_shift(m_masks.left_shift_seg_width_ext - 2);//m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

		WordPair comb_chrom_seq;

		//cout <<"combine chrom"<<endl;
		prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on_ext, m_masks.suffix_seg_bits_on_ext, m_masks.comb_seg_bits_on_ext, m_three_prim_prefix, comb_chrom_seq);
		//                     (m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

		size_t max_loc, prim, left_mismatch = m_max_mismatches, rbits;
		m_matched_flank = 0;
		m_matched_bads = 0;

		//cout << "scoring"<<endl;

		//cout << "chromwp " << comb_chrom_seq.tostring2(m_masks.comb_seg_bits_on_ext)<<endl;
		//cout << "readswp " << dup_missed_read.tostring2(m_masks.comb_seg_bits_on_ext)<<endl;

		
		size_t score = Fixhole_score_selective_ext(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch, rbits);
					 //Fixfusion_score(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);

		//cout << "score:" << score<<endl;
		//cout << "max_loc:" << max_loc<<endl;
		//cout << "prim:" << prim<<endl;

		//getchar();
		size_t total_mismatch = other_mismatch + score;

		//cout <<  "total score:" <<total_mismatch << endl;

		//getchar();

		if (score <= m_max_mismatches && max_loc <= m_seg_len + (2 * m_extend_bits) && total_mismatch <= m_max_fusion_mismatch)
		{
			if (mark)
			{
				++m_cur_reads_count[tagidx];
			}
			int strandi = 0;
			if (strand == '-')
				strandi = 16;

			string prim5flankstr = prefixflankstr;
			string prim3flankstr = suffixflankstr;

			string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len + (2 * m_extend_bits) - max_loc), m_seg_len + (2 * m_extend_bits) - max_loc);

			comb_mapped_chromseq = comb_mapped_chromseq.substr(m_extend_bits, comb_mapped_chromseq.length() - m_extend_bits);

			//cout <<"comb_mapped_chromseq:"<<comb_mapped_chromseq<<endl; 

			string flankstr;

			FlankString(flankstr);			
			
			//<<strand1<<strand2<<'\t'<< chr1<< '_'<<chr2<<'\t';

			//revcomp flank string

			//string flankstr1 = flankstr.substr(0, 2);

			//if (strand1 == '-')
			//	flankstr1 = revcomp(flankstr1);

			//string flankstr2 = flankstr.substr(2, 2);

			//if (strand2 == '-')
			//	flankstr2 = revcomp(flankstr2);

			//flankstr = flankstr1 + flankstr2;

			size_t suffixst, prefixend;

			string prefixsyn, suffixsyn;

			size_t synlen = 60;

			if (strand1 == '+' && strand2 == '+')
			{
				prefixend = offset1 + m_seg_len + max_loc - m_extend_bits;
				suffixst = offset2 - m_seg_len + max_loc - m_extend_bits + 1;

				if (prefixend < synlen || suffixst + synlen >= chromseq2.length())
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend <<'\t' <<suffixst + synlen<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - synlen, synlen);
				suffixsyn = chromseq2.substr(suffixst - 1, synlen);

				ofs<<tagname<<'\t';

				//ofs << offset + m_seg_len + max_loc << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;
				ofs << chr1 << '\t' << strand1 << '\t' << prefixend - max_loc + m_extend_bits - seg1readstr.length() + 1 << '\t' <<prefixend 
					<<'\t'<<chr2 << '\t' << strand2 << '\t' << suffixst << '\t' << suffixst + m_seg_len - max_loc + m_extend_bits + seg3readstr.length() - 1
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr <<'\t' <<  "NM:i:"<<total_mismatch <<endl;

				
			}
			else if (strand1 == '+' && strand2 == '-')
			{
				prefixend = offset1 + m_seg_len + max_loc - m_extend_bits;
				suffixst = offset2 + m_seg_len + m_seg_len - max_loc + m_extend_bits;
				//ofs << offset + m_seg_len + max_loc << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

				if (prefixend < synlen || suffixst < synlen)
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend <<'\t' <<suffixst<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - synlen, synlen);
				suffixsyn = chromseq2.substr(suffixst - synlen, synlen);

				suffixsyn = revcomp(suffixsyn);

				ofs<<tagname<<'\t';

				ofs << chr1 << '\t' << strand1 << '\t' << prefixend - max_loc + m_extend_bits - seg1readstr.length() + 1 << '\t' <<prefixend 
					<<'\t'<<chr2 << '\t' << strand2 << '\t' << suffixst - (m_seg_len - max_loc + m_extend_bits + seg3readstr.length() - 1) << '\t' << suffixst
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/  << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;

				
			}
			else if (strand1 == '-' && strand2 == '-')
			{
				prefixend = offset1 - max_loc + m_extend_bits + 1;
				suffixst = offset2 + m_seg_len + m_seg_len - max_loc + m_extend_bits;

				if (prefixend + synlen >= chromseq1.length() || suffixst < synlen)
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend + synlen <<'\t' <<suffixst<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - 1, synlen);
				suffixsyn = chromseq2.substr(suffixst - synlen, synlen);

				prefixsyn = revcomp(prefixsyn);
				suffixsyn = revcomp(suffixsyn);

				ofs<<tagname<<'\t';

				//ofs << offset - max_loc + 1 << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;
				ofs << chr1 << '\t' << strand1 << '\t' << prefixend << '\t' <<prefixend + max_loc - m_extend_bits + seg1readstr.length() - 1 
					<<'\t'<<chr2 << '\t' << strand2 << '\t' << suffixst - (m_seg_len - max_loc + m_extend_bits + seg3readstr.length() - 1) << '\t' << suffixst
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;

				
			}
			else
			{
				prefixend  = offset1 - max_loc + m_extend_bits + 1;
				suffixst = offset2 - m_seg_len + max_loc - m_extend_bits + 1;

				if (prefixend + synlen >= chromseq1.length() || suffixst + synlen >= chromseq2.length())
				{
					cerr << "warning: fusion synlen exceed\t" << prefixend + synlen <<'\t' <<suffixst + synlen<<endl;
					continue;
				}

				prefixsyn = chromseq1.substr(prefixend - 1, synlen);
				suffixsyn = chromseq2.substr(suffixst - 1, synlen);

				prefixsyn = revcomp(prefixsyn);

				ofs<<tagname<<'\t';

				//ofs << offset - max_loc + 1 << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;
				ofs << chr1 << '\t' << strand1 << '\t' << prefixend << '\t' <<prefixend + max_loc - m_extend_bits + seg1readstr.length() - 1 
					<<'\t'<<chr2 << '\t' << strand2 << '\t' << suffixst << '\t' << suffixst + m_seg_len - max_loc + m_extend_bits + seg3readstr.length() - 1
					<< '\t' <<prefixsyn << suffixsyn/*seg1chromstr << comb_mapped_chromseq << seg3chromstr*/ << '\t' << flankstr << '\t' <<  "NM:i:"<<total_mismatch <<endl;

				
			}

			//ofs << '\t' <<flankstr <<endl;

			//cout <<"insert junction"<<endl;

			string chrom_ID = chr1;
			chrom_ID.append("_");
			chrom_ID.append(chr2);

			//cout << "find(suffixst) "<<endl;
			JUNCFusion_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

			if (iter_conj == conj_hash_map[chrom_ID].end())
			{
				//cout << "insert(JUNC_HASH: "<<endl;
				JUNCFusion_SUFFIX match_hash_temp;
				conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
				iter_conj = conj_hash_map[chrom_ID].find(suffixst);
			}

			JUNCFusion_SUFFIX_ITER iter_conj_pre;

			
			if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
			{
				//if (iter_conj_pre->second.flankstr != flankstr)
				//{
				//	cout <<"flank does not match"<<endl;
				//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
				//}
				//cout << "inc_hits(max_loc) "<<endl;
				//iter_conj_pre->second.inc_hits(max_loc + seg1chromstr.length() - m_extend_bits, (unsigned short)total_mismatch);	
			}	  
			else 
			{
				//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
				//cout << "insert(JUNC_SUFFIX "<<endl;
				//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

				string synstr = prefixsyn;
				synstr.append(suffixsyn);
				
				string strand12 = "++";
				strand12[0] = strand1;
				strand12[1] = strand2;
				//(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankstr, max_loc + seg1chromstr.length() - m_extend_bits, m_read_width, synstr, strand12, (unsigned short)total_mismatch)));
			}
		}

	}

	input_fs.close();
	ofs.close();

	//while(!input_fs.eof())
	//{
	//	//cout <<"read number"<<endl;
	//	//if (rec_num_str == "")
	//	//	continue;
	//	//char *endptr;
	//	//size_t rec_num = strtol(rec_num_str.c_str(), &endptr, 10);
	//	string line1, line2;

	//	getline(input_fs, line1);
	//	getline(input_fs, line2);

	//	if (line1 == "" || line2 == "")
	//		continue;

	//	vector<GroupedBwtRec> fix_hole_vec;

	//	char chromname[1000], tagname[1000], mapped_seq[1000], differs[1000], strand;
	//	int offset, seg_id, other_mapped;

	//	sscanf(line1.c_str(), "%s\t%d\t%c\t%s\t%d", tagname, &seg_id, &strand, chromname, &offset);
	//	fix_hole_vec.push_back(GroupedBwtRec(offset, seg_id, "0", 0, "0"));

	//	char chromname2[1000], tagname2[1000], mapped_seq2[1000], differs2[1000], strand2;
	//	int offset2, seg_id2, other_mapped2;
	//	sscanf(line2.c_str(), "%s\t%d\t%c\t%s\t%d", tagname2, &seg_id2, &strand2, chromname2, &offset2);
	//	fix_hole_vec.push_back(GroupedBwtRec(offset2, seg_id2, "0", 0, "0"));
	//	//}

	//	//string tagnamestr = tagname;
	//	//size_t last_idx = tagnamestr.find_last_of("_");
	//	//string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);

	//	if (m_reads_name_map_id.find(tagname) == m_reads_name_map_id.end())
	//	{
	//		cout << "can't find tagname "<<endl<< line1<<endl<<line2<<endl;
	//		continue;
	//	}

	//	string chr1 = chromname;

	//	string chr2 = chromname2;

	//	if (chr1 == chr2)
	//		return;

	//	size_t tagidx = m_reads_name_map_id[tagname] + 1;

	//	size_t missed_id = 2;

	//	//vector<GroupedBwtRec>::iterator bwt_rec_iter;

	//	//for (bwt_rec_iter = fix_hole_vec.begin(); bwt_rec_iter != fix_hole_vec.end() - 1; ++bwt_rec_iter)
	//	//{
	//	//	if (bwt_rec_iter->seg_id + 1 != (bwt_rec_iter + 1)->seg_id)
	//	//	{
	//	//		if (bwt_rec_iter->seg_id + 2 != (bwt_rec_iter + 1)->seg_id)
	//	//		{
	//	//			cout << "missed more than 1 segments"<<endl;
	//	//			exit(1);
	//	//		}
	//	//		else
	//	//		{
	//	//			missed_id = bwt_rec_iter->seg_id + 1;
	//	//			break;
	//	//		}
	//	//	}
	//	//}

	//	//if (bwt_rec_iter->offset + m_seg_len + m_seg_len + m_interal > (bwt_rec_iter + 1)->offset || bwt_rec_iter->offset + m_seg_len + m_seg_len + m_max_splice_dist < (bwt_rec_iter + 1)->offset)
	//	//	continue;

	//	size_t readID;

	//	WordPair dup_missed_read;

	//	readID = tagidx;
	//	m_read_words[readID].duplicate(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_masks.left_shift_seg_width, m_masks.right_shift_seg_width, dup_missed_read);

	//	string seg1readstr = m_reads[readID - 1];
	//	string seg3readstr = m_reads[readID + 1];

	//	string chromfile1 = m_chrom_files[m_chrom_names_map[chromname]];
	//	string chromfile2 = m_chrom_files[m_chrom_names_map[chromname2]];
	//	
	//	if (prevchrom_file1 == "" || chromfile1 != prevchrom_file1)
	//	{
	//		chromseq1.clear();
	//		//cout << "read chromo\t"<<chromfile<<endl;
	//		readchrom(chromfile1.c_str(), chromseq1);
	//		if (chromseq1 == "")
	//		{
	//			cout <<"empty chrom: "<<chromfile1<<endl;
	//			exit(1);
	//		}
	//		prevchrom_file1 = chromfile1;
	//	}

	//	if (prevchrom_file2 == "" || chromfile2 != prevchrom_file2)
	//	{
	//		chromseq2.clear();
	//		//cout << "read chromo\t"<<chromfile<<endl;
	//		readchrom(chromfile2.c_str(), chromseq2);
	//		if (chromseq2 == "")
	//		{
	//			cout <<"empty chrom: "<<chromfile2<<endl;
	//			exit(1);
	//		}
	//		prevchrom_file2 = chromfile2;
	//	}
	//	

	//	//cout << "get substring"<<endl;
	//	string prefixflankstr;// = chromseq.substr(bwt_rec_iter->offset + bwt_rec_iter->mapped_chrom_seq.length(), m_seg_len + 2);
	//	string suffixflankstr;// = chromseq.substr((bwt_rec_iter + 1)->offset - m_seg_len - 2, m_seg_len + 2);

	//	string seg1chromstr;
	//	string seg3chromstr;

	//	if (strand == '+' && strand2 == '+')
	//	{
	//		prefixflankstr = chromseq1.substr(offset + m_seg_len, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2, m_seg_len + 2);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);
	//	}
	//	else if (strand == '+' && strand2 == '-')
	//	{
	//		prefixflankstr = chromseq1.substr(offset + m_seg_len, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 + m_seg_len, m_seg_len + 2);

	//		suffixflankstr = revcomp(suffixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg3chromstr = revcomp(seg3chromstr);

	//	}
	//	else if (strand == '-' && strand2 == '-')
	//	{
	//		prefixflankstr = chromseq1.substr(offset - m_seg_len - 2, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 + m_seg_len, m_seg_len + 2);

	//		prefixflankstr = revcomp(prefixflankstr);

	//		suffixflankstr = revcomp(suffixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg1chromstr = revcomp(seg1chromstr);
	//		seg3chromstr = revcomp(seg3chromstr);
	//	}
	//	else
	//	{
	//		prefixflankstr = chromseq1.substr(offset - m_seg_len - 2, m_seg_len + 2);
	//		suffixflankstr = chromseq2.substr(offset2 - m_seg_len - 2, m_seg_len + 2);

	//		prefixflankstr = revcomp(prefixflankstr);

	//		seg1chromstr = chromseq1.substr(offset, m_seg_len);
	//		seg3chromstr = chromseq2.substr(offset2, m_seg_len);

	//		seg1chromstr = revcomp(seg1chromstr);
	//	}

	//	if (seg1readstr.length() != seg1chromstr.length() || seg3readstr.length() != seg3chromstr.length())
	//	{
	//		cout << "read length does not match chrom length "<<endl;
	//		cout <<seg1readstr<<endl;
	//		cout <<seg1chromstr<<endl; 
	//		cout <<seg3readstr<<endl;
	//		cout <<seg3chromstr<<endl; 
	//	}

	//	for (size_t i = 0; i < seg1readstr.length(); ++i)
	//	{
	//		if (seg1readstr[i] >= 'a' && seg1readstr[i] <= 'z')
	//			seg1readstr[i] = seg1readstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg1chromstr.length(); ++i)
	//	{
	//		if (seg1chromstr[i] >= 'a' && seg1chromstr[i] <= 'z')
	//			seg1chromstr[i] = seg1chromstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg3readstr.length(); ++i)
	//	{
	//		if (seg3readstr[i] >= 'a' && seg3readstr[i] <= 'z')
	//			seg3readstr[i] = seg3readstr[i] + 'A' - 'a';
	//	}

	//	for (size_t i = 0; i < seg3chromstr.length(); ++i)
	//	{
	//		if (seg3chromstr[i] >= 'a' && seg3chromstr[i] <= 'z')
	//			seg3chromstr[i] = seg3chromstr[i] + 'A' - 'a';
	//	}

	//	size_t count = 0;
	//	for (size_t i = 0; i < seg3readstr.length(); ++i)
	//	{
	//		if (seg3readstr[i] != seg3chromstr[i])
	//			++count;
	//	}

	//	if (count > m_max_mismatches)
	//	{
	//		cout << "seg 3 mismatch exceed "<<endl;
	//		cout << seg3readstr<<endl;
	//		cout << seg3chromstr << endl;

	//		getchar();
	//	}

	//	size_t count2 = 0;

	//	for (size_t i = 0; i < seg1readstr.length(); ++i)
	//	{
	//		if (seg1readstr[i] != seg1chromstr[i])
	//			++count2;
	//	}

	//	if (count2 > m_max_mismatches)
	//	{
	//		cout << "seg 1 mismatch exceed "<<endl;
	//		cout << seg1readstr<<endl;
	//		cout << seg1chromstr << endl;

	//		getchar();
	//	}

	//	m_five_prim_suffix = WordPair(prefixflankstr);
	//	m_three_prim_prefix = WordPair(suffixflankstr);

	//	WordPair prefix_chrom = m_five_prim_suffix;

	//	m_five_prim_suffix.left_shift(2); // make space for 2 bp of 3' flank string

	//	prefix_chrom.left_shift(m_masks.left_shift_seg_width - 2);//it's length is seg length + 2 basepair of flank string, ignore the last two base pairs

	//	WordPair comb_chrom_seq;

	//	//cout <<"combine chrom"<<endl;
	//	prefix_chrom.ps_combine(m_masks.prefix_seg_bits_on, m_masks.suffix_seg_bits_on, m_masks.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	//	size_t max_loc, prim, left_mismatch = m_max_mismatches;
	//	m_matched_flank = 0;
	//	m_matched_bads = 0;

	//	//cout << "scoring"<<endl;
	//	size_t score = Fixhole_score_selective(dup_missed_read, comb_chrom_seq, max_loc, prim, left_mismatch);

	//	//cout << "5 prim and 3 prim"<<endl;
	//	//cout << prefixflankstr << endl;
	//	//cout << suffixflankstr << endl;

	//	//cout <<"score: "<<score<<endl;
	//	//cout <<"max_loc: "<<max_loc<<endl;
	//	if (score <= m_max_mismatches && max_loc <= m_seg_len)
	//	{
	//		int strandi = 0;
	//		if (strand == '-')
	//			strandi = 16;

	//		string prim5flankstr = prefixflankstr;
	//		string prim3flankstr = suffixflankstr;
	//		string comb_mapped_chromseq = prim5flankstr.substr(0, max_loc) + prim3flankstr.substr(prim3flankstr.length() - (m_seg_len - max_loc), m_seg_len - max_loc);

	//		string flankstr;

	//		FlankString(flankstr);

	//		ofs<<tagname<<'\t'<<strand<<strand2<<'\t'<< chromname<< '_'<<chromname2<<'\t';//<<fix_hole_vec.front().offset + 1<<'\t'<<0<<'\t';

	//		size_t suffixst, prefixend;

	//		string prefixsyn, suffixsyn;

	//		if (strand == '+' && strand2 == '+')
	//		{
	//			prefixend = offset + m_seg_len + max_loc;
	//			suffixst = offset2 - m_seg_len + max_loc + 1;
	//			ofs << offset + m_seg_len + max_loc << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 30, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 1, 30);
	//		}
	//		else if (strand == '+' && strand2 == '-')
	//		{
	//			prefixend = offset + m_seg_len + max_loc;
	//			suffixst = offset2 + m_seg_len + m_seg_len - max_loc;
	//			ofs << offset + m_seg_len + max_loc << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 30, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 30, 30);

	//			suffixsyn = revcomp(suffixsyn);
	//		}
	//		else if (strand == '-' && strand2 == '-')
	//		{
	//			prefixend = offset - max_loc + 1;
	//			suffixst = offset2 + m_seg_len + m_seg_len - max_loc;
	//			ofs << offset - max_loc + 1 << '\t' << offset2 + m_seg_len + m_seg_len - max_loc<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 1, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 30, 30);

	//			prefixsyn = revcomp(prefixsyn);
	//			suffixsyn = revcomp(suffixsyn);
	//		}
	//		else
	//		{
	//			prefixend  = offset - max_loc + 1;
	//			suffixst = offset2 - m_seg_len + max_loc + 1;
	//			ofs << offset - max_loc + 1 << '\t' << offset2 - m_seg_len + max_loc + 1<<'\t'<<max_loc;

	//			prefixsyn = chromseq1.substr(prefixend - 1, 30);
	//			suffixsyn = chromseq2.substr(suffixst - 1, 30);

	//			prefixsyn = revcomp(prefixsyn);
	//		}

	//		ofs << '\t' <<flankstr <<endl;

	//		string chrom_ID = chromname;
	//		chrom_ID.append("_");
	//		chrom_ID.append(chromname2);

	//		//cout << "find(suffixst) "<<endl;
	//		JUNCFusion_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

	//		if (iter_conj == conj_hash_map[chrom_ID].end())
	//		{
	//			//cout << "insert(JUNC_HASH: "<<endl;
	//			JUNCFusion_SUFFIX match_hash_temp;
	//			conj_hash_map[chrom_ID].insert(JUNCFusion_HASH::value_type(suffixst, match_hash_temp));
	//			iter_conj = conj_hash_map[chrom_ID].find(suffixst);
	//		}

	//		JUNCFusion_SUFFIX_ITER iter_conj_pre;

	//		
	//		if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
	//		{
	//			//if (iter_conj_pre->second.flankstr != flankstr)
	//			//{
	//			//	cout <<"flank does not match"<<endl;
	//			//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
	//			//}
	//			//cout << "inc_hits(max_loc) "<<endl;
	//			iter_conj_pre->second.inc_hits(max_loc);	
	//		}	  
	//		else 
	//		{
	//			//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
	//			//cout << "insert(JUNC_SUFFIX "<<endl;
	//			//cout <<prefixend <<'\t'<< flankstr << '\t' <<max_loc <<'\t'<<m_read_width<<endl;

	//			string synstr = prefixsyn;
	//			synstr.append(suffixsyn);
	//			
	//			string strand12 = "++";
	//			strand12[0] = strand;
	//			strand12[1] = strand2;
	//			(iter_conj->second).insert(JUNCFusion_SUFFIX::value_type(prefixend, JuncInfoFusion(flankstr, max_loc, m_read_width, synstr, strand12)));
	//		}
	//	}
	//}
	
}

//void
//GenomeScan::TriggerMatchFixHole()
//{
//	    //prepare read for scoring, move down to avoid unnecessary operation
//	WordPair pre_wp;
//
//	size_t cur_prefix_offset = cur_seed.offset;
//
//	size_t slide_window_idx = cur_prefix_offset % m_window_size;
//
//	TagSeed* cur_TagSeed = &(m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index]);
//
//	if (mapped_TagSeed.find((size_t)cur_TagSeed) != mapped_TagSeed.end())
//		return;
//	else
//		mapped_TagSeed[(size_t)cur_TagSeed] = 1;
//
//	unsigned readID = cur_TagSeed->readID;
//
//	unsigned read_strand = cur_TagSeed->strand;
//
//	(read_strand)? 
//		m_read_words[readID].duplicate(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, m_masks.comp_left_shift_width, m_masks.comp_right_shift_width, pre_wp) :
//	m_read_words_rc[readID].duplicate(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, m_masks.comp_left_shift_width,m_masks.comp_right_shift_width, pre_wp);
//
//	if ((m_chrom_offset + 2 - cur_prefix_offset) < m_interal)
//		return;
//
//	size_t kmer_score = 0;
//
//	if (kmer_ID)
//		kmer_score = pre_wp.score_hmer(m_wp, m_masks.bad_kmer_mask);
//
//	m_left_mismatches = m_max_mismatches - m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index].score;
//
//	if (kmer_score > m_left_mismatches)
//		return;
//	else
//		m_left_mismatches -= kmer_score;
//
//	++m_candidate_count;
//
//	m_five_prim_suffix = m_slide_window[slide_window_idx].five_prim_suffix;
//
//	//the following makes the suffix
//	WordPair suffix_wp(m_wp);
//
//	////get the right buff_size + right suffix
//	//m_wp.combine(m_wp_reserve, m_read_width - m_seed_width, suffix_wp);
//
//	// - m_max_mismatch
//
//	//suffix_wp.right_shift(m_anchor_width);
//
//	m_three_prim_prefix = suffix_wp;
//
//	m_three_prim_prefix.right_shift(m_anchor_width);
//
//	//combine prefix and suffix for scoring  				
//	m_slide_window[slide_window_idx].prefix_wp.ps_combine(m_masks.comp_first_half_mask, m_masks.comp_second_half_mask, m_masks.comp_big_buff_mask, suffix_wp, m_shifted_wp);
//
//	size_t max_loc, prim;
//	m_matched_flank = 0;
//	m_matched_bads = 0;
//
//	//get match score
//	size_t score = (m_reverse_flag) ? Ps_score_selective_rev(pre_wp, max_loc, prim) : Ps_score_selective(pre_wp, max_loc, prim);
//
//	if (m_flank_case && !prim)
//		return;
//
//	score += m_slide_window[slide_window_idx].matchedTagsHalf[cur_seed.index].score + kmer_score;
//
//	//if (m_reverse_flag && m_read_width % 2 == 0 && max_loc == 0)
//	//	return;
//
//	if (score <= m_max_mismatches)
//	{
//		size_t lste_st  = cur_prefix_offset - m_read_width + 1;
//		size_t lste_end = cur_prefix_offset - m_read_width + (m_read_width + 1) / 2 + max_loc;
//		// - m_max_mismatch
//		size_t rste_st  = m_chrom_offset + 1 - m_anchor_width - (m_masks.score_buff_width - max_loc) + 1;
//		size_t rste_end = m_chrom_offset + 1;
//
//		if (m_islands_provide)
//		{
//			size_t prev_island_st= m_slide_window[slide_window_idx].cur_tac_pos->first;
//
//			size_t prev_island_ed = m_slide_window[slide_window_idx].cur_tac_pos->second;
//
//			if (lste_end < prev_island_st + 1 ||  rste_st < m_cur_island_start + 1)
//				return;
//
//			if (lste_end > prev_island_ed + 1) 
//			{
//				
//				PAIR_SIZE_T_VEC_ITER prefix_next_tac_pos = m_slide_window[slide_window_idx].cur_tac_pos + 1;
//				
//				if (!(prefix_next_tac_pos != m_tac_pos_end && prefix_next_tac_pos->first - prev_island_ed <= m_boundary && 
//					lste_end >= prefix_next_tac_pos->first + 1 && lste_end <= prefix_next_tac_pos->second + 1))
//				{
//					return;
//				}
//			}
//
//			if (rste_st > m_cur_island_end + 1)
//			{
//				PAIR_SIZE_T_VEC_ITER suffix_next_tac_pos = m_cur_tac_pos + 1;
//
//				if (!(suffix_next_tac_pos != m_tac_pos_end && suffix_next_tac_pos->first - m_cur_island_end <= m_boundary && 
//					rste_st >= suffix_next_tac_pos->first + 1 && rste_st <= suffix_next_tac_pos->second + 1))
//				{
//					return;
//				}
//			}
//		}
//
//		size_t prefix_width = max_loc + (m_read_width + 1) / 2;
//
//		m_score = score;
//
//		if (m_reverse_flag)
//		{
//			size_t temp_lste_st = m_chrom_sizes[chrom_ID] - lste_st + 1;
//			size_t temp_lste_end = m_chrom_sizes[chrom_ID] - lste_end + 1;
//			size_t temp_rste_st = m_chrom_sizes[chrom_ID] - rste_st + 1;
//			size_t temp_rste_end = m_chrom_sizes[chrom_ID] - rste_end + 1;
//			lste_end = temp_rste_st;
//			lste_st = temp_rste_end;
//			rste_st = temp_lste_end;
//			rste_end = temp_lste_st;
//		}
//
//		clock_t t1=clock();
//
//		WriteMapReads(prefix_width, m_slide_window[slide_window_idx].prefix_wp, prim, lste_end, rste_st, chrom_ID, readID, read_strand);
//
//		clock_t t2=clock();
//
//		m_output_time += (t2-t1)/(double)CLOCKS_PER_SEC;
//
//	}
//}

void
GenomeScan::StartFixHmer(int phase, size_t tissueID, size_t flankcase, bool avoid)
{
	//map<size_t, map<string, vector<SamFixtail> > >::iterator missedid_iter;

	map<string, vector<GroupedBwts*> >::iterator fixhmer_iter;

	map<string, vector<GroupedBwts*> >* fixhmer_ptr_ptr = 0;

	if (phase == 1)
		fixhmer_ptr_ptr = &m_fixhmer_ptr;
	else if (phase == 2)
		fixhmer_ptr_ptr = &m_fixhmer_ptr_rev;

	m_fixhmer_ptr;

	for (fixhmer_iter = fixhmer_ptr_ptr->begin(); fixhmer_iter != fixhmer_ptr_ptr->end(); ++fixhmer_iter)
	{
		if (m_chrom_names_map.find(fixhmer_iter->first) == m_chrom_names_map.end())
			continue;

		//test on chr1
		//if (fixhmer_iter->first != "chr1")
		//	continue;

		size_t chromID = m_chrom_names_map.find(fixhmer_iter->first)->second;

		cout << "fix hmer: "<<fixhmer_iter->first <<'\t'<<chromID<<endl;

		//CreateMapReadsFile(tissueID, chromID);

		FixHmer(chromID, fixhmer_iter->second, flankcase, avoid);

		Clear();

		//CloseCurMapReadsFile();
	}
	
	//for (missedid_iter = m_sam_fix_tail_map.begin(); missedid_iter != m_sam_fix_tail_map.end(); ++missedid_iter)
	//{
	//	if ((missedid_iter->first == 1 && phase == 1) || (phase == 2 && missedid_iter->first != 1))
	//		continue;

	//	map<string, vector<SamFixtail> >::iterator chrom_iter;

	//	for (chrom_iter = missedid_iter->second.begin(); chrom_iter != missedid_iter->second.end(); ++chrom_iter)
	//	{
	//		size_t chromID = m_chrom_names_map[chrom_iter->first];

	//		cout << chrom_iter->first<<endl;

	//		CreateMapReadsFile(tissueID, chromID);

	//		FixTail(chromID, chrom_iter->second);

			
	//	}
	//}
}

void
GenomeScan::StartMapReads(int phase, size_t tissueID)
{
	for (size_t chromID = 0; chromID < m_chrom_files.size(); ++chromID)
	{
		cout << "tissue : " <<m_tissue_names[tissueID]<<endl;
		cout << "chrom : " <<m_chrom_files[chromID]<<endl;

		CreateMapReadsFile(tissueID, chromID);

		clock_t t1=clock();

		MapReads(chromID);

		clock_t t2=clock();

		CloseCurMapReadsFile();

		for (size_t i = 0; i < m_open_reads.size(); ++i)
			m_open_reads[i].stats();

		m_total_time += (t2-t1)/(double)CLOCKS_PER_SEC;

		m_map_time += m_total_time - m_output_time;

		ofstream time_ofs;

		string fcounttime = m_path;
		fcounttime.append("log.txt");

		time_ofs.open(fcounttime.c_str(), ios_base::app);

		time_ofs << m_tissue_names[tissueID] << '\t';
		time_ofs << m_chrom_names[chromID]<<'\t';
		time_ofs << m_max_splice_dist<<'\t';
		time_ofs << m_interal<<'\t';
		time_ofs << m_seed_width<<'\t';
		time_ofs << m_max_mismatches<<'\t';
		time_ofs << m_boundary <<'\t';
		if (!m_prefix_match)
			time_ofs << "no_prefix_match\t";
		else
			time_ofs << "prefix_match\t";
		time_ofs << m_islands <<'\t';
		time_ofs << m_basecount<<'\t';
		time_ofs << m_candidate_count<<'\t';
		if (phase == 1)
			time_ofs << "First Phase\t";
		else if(phase == 2)
			time_ofs << "Second Phase\t";
		time_ofs << m_init_time<<'\t';
		time_ofs << m_map_time<<'\t';
		time_ofs << m_output_time <<'\t';

		#ifdef TIME
		time_ofs << m_enterhalf_count <<'\t';
		time_ofs << m_passhalf_count <<'\t';
		time_ofs << m_halftag_time <<'\t';
		time_ofs << m_checkseedtable_time <<'\t';
		time_ofs << m_matainslidewindow_time <<'\t';
		time_ofs << m_matainopenreads_time <<'\t';
		time_ofs << m_checkcandclearopenread_time <<'\t';
		time_ofs << m_triggermatch_time<<'\t';
		#endif

		time_ofs <<endl;

		time_ofs.close();


		m_map_time = 0;

		#ifdef TIME
		m_halftag_time = 0;

		m_checkseedtable_time = 0;

		m_matainslidewindow_time = 0;

		m_matainopenreads_time = 0;

		m_checkcandclearopenread_time = 0;
		
		m_triggermatch_time = 0;
		#endif

		m_output_time = 0;
		m_total_time = 0;

		Clear();
	}
}

void
GenomeScan::CallMapReads()
{
	if (VERBOSE)
		cerr << endl << "scanning chromosomes:" << endl;

	vector<string> mapped_reads;

	for (size_t tissueID = 0; tissueID < m_tissue_files.size(); ++tissueID)
	{
		///////InitReads(tissueID);

		cout <<"first phase\n";

		///////InitReadWords();

		m_reverse_flag = false;

		char fix_hole_file[1000];

		sprintf(fix_hole_file, "%s%sfixhole_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		m_ofs_fix_hole.open(fix_hole_file);

		char allmapreads_file[1000];

		sprintf(allmapreads_file, "%s%smapreads_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		mapped_reads.push_back(allmapreads_file);
		m_ofs_allmapreads.open(allmapreads_file);
		
		StartFixHmer(1, tissueID, m_flank_case, true);

		ClearReadWords();

		cout<<"second phase\n";

#ifndef DEBUG_CHECK
		InitReads(tissueID);
#endif

		InitReadWords_Rev();

#ifndef DEBUG_CHECK
		m_reads.clear();
#endif

		m_reverse_flag = true;

		StartFixHmer(2, tissueID, m_flank_case, true);

		char fix_head_tail_file[1000];

		sprintf(fix_head_tail_file, "%s%sfix_head_tail_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		mapped_reads.push_back(fix_head_tail_file);

		WriteFixedHeadAndTail(fix_head_tail_file, m_flank_case);

		m_ofs_fix_hole.close();

		m_ofs_allmapreads.close();

		ClearReadWords();

#ifndef DEBUG_CHECK
		InitReads(tissueID);
#endif

		InitReadWords();

#ifndef DEBUG_CHECK
		m_reads.clear();
#endif

		char fixed_hole_file[1000];

		sprintf(fixed_hole_file, "%s%sfixed_hole_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		mapped_reads.push_back(fixed_hole_file);

		cout << "fix hole after fix hmer"<<endl;

		//canonical
		FixHoleExt(fix_hole_file, fixed_hole_file, m_flank_case, false, true, false);

		cout << "fix hole"<<endl;
		
		char fixed_fixhole_file[1000];

		sprintf(fixed_fixhole_file, "%s%sfixed_fixhole_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		//canonical 
		FixHoleExt(m_fixhole_file.c_str(), fixed_fixhole_file, m_flank_case, false, true, false);

		UpdateReadCountAndClean();

		char fixed_fixhole_file_exceed[1000];

		sprintf(fixed_fixhole_file_exceed, "%s%sfixed_fixhole_exceed_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		cout << "fix hole exceed"<<endl;

		//canonical exceed
		FixHoleExt(m_fixhole_exceed_file.c_str(), fixed_fixhole_file_exceed, m_flank_case, true, true, true);

		UpdateReadCountAndClean();

		mapped_reads.push_back(fixed_fixhole_file);

		char fixed_fixhole_file_f0[1000];

		sprintf(fixed_fixhole_file_f0, "%s%sfixed_fixhole_f0_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		////fix hmer non-canonical
		//cout <<"first phase\n";

		///////InitReadWords();

		//m_reverse_flag = false;

		//char fix_hole_file_noncanon[1000];

		//sprintf(fix_hole_file_noncanon, "%s%sfixhole_noncanon_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		//m_ofs_fix_hole.open(fix_hole_file_noncanon);

//		char allmapreads_file_noncanon[1000];
//
//		sprintf(allmapreads_file_noncanon, "%s%smapreads_noncanon_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());
//
//		m_ofs_allmapreads.open(allmapreads_file_noncanon);
//		
//		StartFixHmer(1, tissueID, 0, true);
//
//		ClearReadWords();
//
//		cout<<"second phase\n";
//
//#ifndef DEBUG_CHECK
//		InitReads(tissueID);
//#endif
//
//		InitReadWords_Rev();
//
//#ifndef DEBUG_CHECK
//		m_reads.clear();
//#endif
//
//		m_reverse_flag = true;
//
//		StartFixHmer(2, tissueID, 0, true);
//
//		char fix_head_tail_file_noncanon[1000];
//
//		sprintf(fix_head_tail_file_noncanon, "%s%sfix_head_tail_noncanon_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());
//
//		WriteFixedHeadAndTail(fix_head_tail_file_noncanon, 0);
//
//		m_ofs_fix_hole.close();
//
//		m_ofs_allmapreads.close();

		/////fix hmer end


//		ClearReadWords();
//
//#ifndef DEBUG_CHECK
//		InitReads(tissueID);
//#endif
//
//		InitReadWords();
//
//#ifndef DEBUG_CHECK
//		m_reads.clear();
//#endif

		//non-canonical 

		char fixed_hole_file_f0[1000];

		sprintf(fixed_hole_file_f0, "%s%sfixed_hole_f0_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		char fixed_hole_file_f1_f0[1000];

		sprintf(fixed_hole_file_f1_f0, "%s%sfixed_hole_f1_f0_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		//FixHoleExt(fix_hole_file_noncanon, fixed_hole_file_f0, 0, true, true, false);

		FixHoleExt(fix_hole_file, fixed_hole_file_f1_f0, 0, true, true, false);
		
		FixHoleExt(m_fixhole_file.c_str(), fixed_fixhole_file_f0, 0, true, true, false);

		UpdateReadCountAndClean();
		
		char fixed_fixhole_file_exceed_f0[1000];

		sprintf(fixed_fixhole_file_exceed_f0, "%s%sfixed_fixhole_exceed_f0_%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[tissueID].c_str());

		//non-canonical exceed
		FixHoleExt(m_fixhole_exceed_file.c_str(), fixed_fixhole_file_exceed_f0, 0, true, true, true);

		UpdateReadCountAndClean();

		cout << "fix fusion"<<endl;

		string fixhole_fusion_file = m_bwtout;
		fixhole_fusion_file.append(".fusion");

		string fixhole_fusion_fixed_file = m_path;
		fixhole_fusion_fixed_file.append(m_original_dir);
		fixhole_fusion_fixed_file.append("fusion.mappedreads");

		string fusion_junction_file = m_path;
		fusion_junction_file.append(m_original_dir);
		fusion_junction_file.append("fusion_junctions.txt");
		ofstream fusion_junc_ofs(fusion_junction_file.c_str());

		CONJFusion_HASH_MAP conj_hash_map;

		FixHole_fusion(fixhole_fusion_file.c_str(), fixhole_fusion_fixed_file.c_str(), conj_hash_map, true, true);

		cout << "WriteJuncFusion"<<endl;
		WriteJuncFusion(conj_hash_map, fusion_junc_ofs);

		cout << "SortJuncFusion"<<endl;
		SortJuncFusion(fusion_junction_file.c_str());

		ClearReads();		

		m_init_time = 0;
	}

	cout << "remove duplicated mapped reads"<<endl;

	RemoveDupMapreads(mapped_reads);

	cout <<"seperate mapped reads"<<endl;
	SeprateMapreads(mapped_reads);

	cout <<"convert 2 junction"<<endl;
	Covert2Junc();

	//Sites_To_Regions();

	//if (m_ambiguous_file)
	//	Write_Non_Uniques();

	//// Output the results
	//ostream* out = (m_outfile) ? new ofstream(m_outfile) : &cout;
	//copy(m_hits.begin(), m_hits.end(), ostream_iterator<GenomicRegion>(*out, "\n"));
	//if (out != &cout) delete out;
}

void
GenomeScan::OutputJunction()
{
//{
//	char filename[1000];
//
//	string prefixmatch;
//	if (m_prefix_match)
//		prefixmatch = "Y";
//	else
//		prefixmatch = "N";
//
//	sprintf(filename, "%sjunctions_pm%s_x%d_h%d_n%d_m%d_i%d_ext%d.txt", m_path.c_str(), prefixmatch.c_str(), m_max_splice_dist, m_seed_width, m_anchor_width, m_max_mismatches, m_interal, m_boundary);
//	ofstream ofs(filename);
//
//	for (size_t chromID = 0; chromID < m_conj_hashs.size(); ++chromID)
//	{
//		JUNC_HASH_ITER iter_conj = m_conj_hashs[chromID].begin();
//
//		int j=0;
//		while (iter_conj!=m_conj_hashs[chromID].end())
//		{
//			hash_map<size_t, ConjInfo >::iterator iter_conj_pre = iter_conj->second.begin();
//			while (iter_conj_pre!=iter_conj->second.end())
//			{
//				string flankstr = FlankString(iter_conj_pre->second.flank, iter_conj_pre->second.bad);
//
//				if (!iter_conj_pre->second.bad && iter_conj_pre->second.rev)
//					reverse(flankstr.begin(), flankstr.end());
//
//				ofs <<m_chrom_names[chromID] <<'\t'<< iter_conj_pre->first/*iter_conj_pre->second.lsite*/ << "\t" << iter_conj->first/*iter_conj_pre->second.rsite*/ << '\t'<< 
//					iter_conj_pre->second.prim<<'\t'<< flankstr <<'\t';
//
//				for (size_t i = 0; i < iter_conj_pre->second.loc_count.size(); ++i)
//				{
//					if (iter_conj_pre->second.loc_count[i] > 0)
//						ofs << i + 1 << ';';
//				}
//				ofs << '\t';
//
//				double rank = 0;
//
//				for (size_t i = 0; i < iter_conj_pre->second.loc_count.size(); ++i)
//				{
//					if (iter_conj_pre->second.loc_count[i] > 0)
//					{
//						ofs << iter_conj_pre->second.loc_count[i] << ';';
//
//						double pi = (double)iter_conj_pre->second.loc_count[i] / (double)iter_conj_pre->second.hits;
//						rank += pi * log(pi);
//					}
//				}
//				ofs << '\t';
//
//				ofs << rank <<'\t' << iter_conj_pre->second.hits << '\t' << iter_conj_pre->second.kinds << endl;
//
//				j++;
//				iter_conj_pre++;
//			}
//
//			iter_conj++;
//		}
//
//		cout << "number of loops " << j << endl;
//	}
//
//	ofs.close();
}

void
GenomeScan::SeprateMapreads(vector<string>& comb_mapreads_files)
{
	vector<string>::iterator cm_iter;
	string prevchrom = "";
	ofstream ofs;

	string seprate_path = comb_mapreads_files.front();
	//seprate_path.append(comb_mapreads_files.front());

	seprate_path.append(".sepdir");

	string mkdir_cmdstr = "mkdir ";
	mkdir_cmdstr.append(seprate_path);
	system(mkdir_cmdstr.c_str());

	seprate_path.append("/");

	map<string, ofstream* > chrom_ofs_map; 

	for (cm_iter = comb_mapreads_files.begin(); cm_iter != comb_mapreads_files.end(); ++cm_iter)
	{
		ifstream ifs((*cm_iter).c_str());

		if (ifs.is_open())
		{
			string line;
			while (getline(ifs,line))
			{
				if (line == "")
					continue;
				char chromname[1000], readname[1000], chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				int prefixst,  strand_t, incorrect;

				sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s", 
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, chromseq, qualseq);

				char mapreads_file[1000];

				sprintf(mapreads_file, "%smapreads_%s.txt", seprate_path.c_str(), chromname);

				if (chrom_ofs_map.find(mapreads_file) == chrom_ofs_map.end())
				{
					ofstream* ofs = new ofstream;

					ofs->open(mapreads_file);

					chrom_ofs_map[mapreads_file] = ofs;
				}

				(*chrom_ofs_map[mapreads_file])<< line<<endl;
			}

			ifs.close();
		}
	}

	map<string, ofstream*>::iterator mf_iter;
	for (mf_iter = chrom_ofs_map.begin(); mf_iter != chrom_ofs_map.end(); ++mf_iter)
	{
		m_mapreads_files.push_back(mf_iter->first);

		mf_iter->second->close();

		delete mf_iter->second;
	}
}

void
GenomeScan::RemoveDupMapreads(vector<string>& comb_mapreads_files)
{
	map<string, int> mapped_reads;
	vector<string>::iterator cm_iter;
	for (cm_iter = comb_mapreads_files.begin(); cm_iter != comb_mapreads_files.end(); ++cm_iter)
	{
		ifstream ifs((*cm_iter).c_str());

		if (ifs.is_open())
		{
			string line;
			while (getline(ifs,line))
			{
				if (line == "")
					continue;

				if (mapped_reads.find(line) == mapped_reads.end())
					mapped_reads[line] = 1;
				else
					mapped_reads[line]++;
			}

			ifs.close();

		}
		else
		{
			cout << "can't open file "<< (*cm_iter) <<endl;
		}
	}

	char spliced_reads_file[1000];

	sprintf(spliced_reads_file, "%s%ssplicedreads_remdup%s.txt", m_path.c_str(), m_temp_dir.c_str(), m_tissue_names[0].c_str());

	ofstream ofs(spliced_reads_file);

	map<string, int>::iterator msi_iter;

	for (msi_iter = mapped_reads.begin(); msi_iter != mapped_reads.end(); ++msi_iter)
	{
		ofs<<msi_iter->first<<endl;
	}

	ofs.close();

	comb_mapreads_files.clear();

	comb_mapreads_files.push_back(spliced_reads_file);
}

void
GenomeScan::Covert2Junc()
{
	char junc_filename[1000];
	char wig_filename[1000];

	string prefixmatch;
	if (m_prefix_match)
		prefixmatch = "Y";
	else
		prefixmatch = "N";

	sprintf(junc_filename, "%s%scanonical_junctions.txt", m_path.c_str(), m_original_dir.c_str());
	ofstream ofs(junc_filename);

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	sprintf(wig_filename, "%scoverage.wig.txt", m_path.c_str());
	ofstream ofs_wig(wig_filename);

	string headline2 = "track type=bedGraph name=\"Mapsplice - read coverage\"";

	ofs_wig << headline2 << endl;

	for (size_t mi = 0; mi < m_mapreads_files.size(); ++mi)
	{
		ifstream ifs(m_mapreads_files[mi].c_str());

		CONJ_HASH_MAP conj_hash_map;

		VEC_COVER_BLOCK v_coverage_block;

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "")
					continue;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				int /*prim, */ prefixst, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/;

				unsigned short mismatch;
				//sscanf(line.c_str(), "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d", chromname, readname, &prim, flankseq, strand, readseq, prefixseq, suffixseq, &prefixlen, &prefixend, &suffixst);
				                                       //%dM\t%dN\t%dM

//61BEEAAXX:3:10:10644:16254#0/1_742255	0	chr1	4641	103	52M140N23M	*	0	0	GCAAAGGCTCCTCCGGGCCCCTCACCAGCCCCAGGTCCTTTCCCAGAGATGCCCTTGCGCCTCATGACCAGCTTG
//b`bbbbbbbbbbbbabbbbbbbbbbbbbbbbbabb_bbbbbbbbbbbbb_bbbbb_bbabbbbbbbbbbbbb`b_	NM:i:0	
				sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%hu", 
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway, /*&prefixlen, &spliceoutlen, &suffixlen,*/ chromseq, qualseq, &mismatch);

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_last_of("_");
				string tagidxstr = tagnamestr.substr(last_idx + 1, tagnamestr.length() - last_idx);
				size_t tagidx = atoi(tagidxstr.c_str()) - 1;

				//cout << line<<endl;
				vector<pair<size_t, size_t> > spliceway_vec;

				string splicewaystr = spliceway;
				size_t index = 0;
				while (true)
				{
					int maplen, intron;

					if (index == 0)
					{
						sscanf(splicewaystr.c_str() + index, "%dM", &maplen);
						spliceway_vec.push_back(make_pair(prefixst, maplen));

						//if (splicewaystr.find("M", index) == string::npos)
						//	break;

						index = splicewaystr.find("M", index) + 1;

						//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

						//getchar();

						
						
						//if (index == string::npos)
						//	break;
					}
					else
					{
						if (index == splicewaystr.length())
							break;
						sscanf(splicewaystr.c_str() + index, "%dN%dM", &intron, &maplen);
						spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + intron, maplen));
			
						if (splicewaystr.find("M", index) == string::npos || splicewaystr.find("M", index) == splicewaystr.length() - 1)
						{
							//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;
							//if (splicewaystr.find("M", index) == splicewaystr.length() - 1)
							//	getchar();
							break;
						}

						index = splicewaystr.find("M", index) + 1;

						//cout << spliceway_vec.back().first <<'\t'<<spliceway_vec.back().second<<'\t'<<index<<endl;

						//getchar();
						
					}
				}

				vector<pair<size_t, size_t> >::iterator vp_iter;
				for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end() - 1; ++vp_iter)
				{
					//suffixend = prefixst + spliceoutlen + (int)m_read_width;

					prefixst = vp_iter->first;
					size_t prefixend = vp_iter->first + vp_iter->second - 1;

					size_t suffixst = (vp_iter + 1)->first;// + spliceoutlen;

					size_t prefixlen = vp_iter->second;

					size_t suffixlen = (vp_iter + 1)->second;

					//string flankstr = flankseq;

					//coverage block
					v_coverage_block.push_back(CoverageBlock(prefixst + 1, prefixst + prefixlen, 1));
					v_coverage_block.push_back(CoverageBlock(suffixst, suffixlen, 1));

					string chrom_ID = chromname;

					//if (chrom_ID == "chr1" &&  prefixend == 4862 && suffixst == 5002)
					//{
					//	cout <<line<<endl;
					//	for (size_t i = 0; i < spliceway_vec.size(); ++i)
					//		cout << spliceway_vec[i].first<<'\t'<<spliceway_vec[i].second<<endl;
					//	getchar();
					//}
					JUNC_HASH_ITER iter_conj = conj_hash_map[chrom_ID].find(suffixst);

					if (iter_conj == conj_hash_map[chrom_ID].end())
					{
						JUNC_SUFFIX match_hash_temp;
						conj_hash_map[chrom_ID].insert(JUNC_HASH::value_type(suffixst, match_hash_temp));
						iter_conj = conj_hash_map[chrom_ID].find(suffixst);
					}

					JUNC_SUFFIX_ITER iter_conj_pre; 

					if ((iter_conj_pre = (iter_conj->second).find(prefixend))!= (iter_conj->second).end())
					{
						//if (iter_conj_pre->second.flankstr != flankstr)
						//{
						//	cout <<"flank does not match"<<endl;
						//	cout <<iter_conj_pre->second.flankstr << '\t'<<flankstr<<endl;
						//}
						iter_conj_pre->second.inc_hits(prefixlen, suffixlen, tagidx, mismatch);	
					}	  
					else 
					{
						//if (DEBUG) cerr << "new element in conj_hash another element again" << endl;
						(iter_conj->second).insert(JUNC_SUFFIX::value_type(prefixend, JuncInfo(/*prim, *//*flankseq, */prefixlen, suffixlen, m_read_width, tagidx, mismatch)));
					}
				}
			}
			ifs.close();
		}
		else cout << "Unable to open file";

		cout <<"write junc"<<endl;
		WriteJunc(conj_hash_map, ofs);

		cout <<"sort v_coverage_block"<<endl;
		sort(v_coverage_block.begin(), v_coverage_block.end(), compblk);

		cout <<"WriteCoverage"<<endl;
		//WriteCoverage(v_coverage_block, conj_hash_map.begin()->first, ofs_wig, mi);
	}

	ofs.close();

	cout <<"sort junc"<<endl;

	SortJunc(junc_filename);
}

void
GenomeScan::WriteJunc(const CONJ_HASH_MAP& conj_hash_map, ofstream& ofs)
{
	
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CONJ_HASH_MAP::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		size_t chrom_size = m_chrom_sizes[m_chrom_names_map[chm_iter->first]] - 1;

		JUNC_HASH_CITER iter_conj = chm_iter->second.begin();

		int j=0;
		while (iter_conj != chm_iter->second.end())
		{
			
			JUNC_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			while (iter_conj_pre!=iter_conj->second.end())
			{
				int kinds = 0;
				int hits = 0;
				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						++kinds;
						hits += iter_conj_pre->second.p[i];
					}
				}

				double rank = 0;

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						double pi = (double)iter_conj_pre->second.p[i] / (double)hits;
						rank += pi * log(pi);
					}
				}

				if (rank!=0)
					rank = -rank;

				if (rank < m_max_rank)
				{
					iter_conj_pre++;
					continue;
				}

				ofs <<chm_iter->first <<'\t'<< iter_conj_pre->first << '\t' << iter_conj->first<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;
				
				//int juncidlen = (int) log10((double)juncid);

				//string suffix0(8 - juncidlen

				//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';

				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
						//ofs << k << ',';
				//		break;
				//	}
				//}
				ofs << iter_conj_pre->second.max_prefix_len << ',';

				//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				//{
				//	if (iter_conj_pre->second.p[i] > 0)
				//	{
				//		ofs << m_read_width - i - 1 << ',';
				//		break;
				//	}
				//}
				ofs << iter_conj_pre->second.max_suffix_len << ',';
				ofs << '\t';

				//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				//{
				//	if (iter_conj_pre->second.p[k - 1] > 0)
				//	{
				//		ofs << 0 << ',';
				//		break;
				//	}
				//}
				ofs << 0 << ',';

				//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				//{
				//	if (iter_conj_pre->second.p[i] > 0)
				//	{
				//		ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
				//		break;
				//	}
				//}
				ofs << int (iter_conj->first + iter_conj_pre->second.max_suffix_len - iter_conj_pre->first + 1 )<<',';
				ofs << '\t';

				size_t intron_len = iter_conj->first - iter_conj_pre->first - 1;

				double ppower = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower = pow(1.0 - ppower, (double)chrom_size);

				double qpower = pow(0.25, double(iter_conj_pre->second.max_suffix_len));

				double pDpower = pow(1.0 - qpower, (double)intron_len);

				double lpq = 1.0 - (pNpower * pDpower);


				double ppower2 = pow(0.25, double(iter_conj_pre->second.max_prefix_len));

				double pNpower2 = pow(1.0 - ppower2, (double)intron_len );

				double qpower2 = pow(0.25, double(iter_conj_pre->second.max_suffix_len));

				double pDpower2 = pow(1.0 - qpower2, (double)chrom_size);

				double lpq2 = 1.0 - (pNpower2 * pDpower2);

				double lpqave = 1.0 - (lpq + lpq2) / 2;

				ofs << rank << '\t'<< lpqave/*iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr */
					<<'\t'<<iter_conj_pre->second.min_mismatch<<'\t' <<iter_conj_pre->second.max_mismatch<<'\t'
					<< (double)iter_conj_pre->second.sum_mismatch / (double)hits<<endl;

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}

struct JuncForSort{
	int juncst;
	int juncend;
	int hits;
	int kinds;
	string blocks;
	string blocksoffset;
	string rank;
	string lpq;
	unsigned short min_mismatch;
	unsigned short max_mismatch;
	double ave_mismatch;

	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& l, unsigned short min_mis, unsigned short max_mis, double ave_mis) : juncst(jst), juncend(jend), hits(hts), 
		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}
};

bool
compjunc(const JuncForSort& lhs, const JuncForSort& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

struct JuncForSortFusion{
	int juncst;
	int juncend;
	int hits;
	int kinds;
	string blocks;
	string blocksoffset;
	string rank;
	string flankcase;
	string flankchr;
	string strand;
	string synstr;
	unsigned short min_mismatch;
	unsigned short max_mismatch;
	double ave_mismatch;

	JuncForSortFusion(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& fc, const string& fcr, const string& sd, const string& st, unsigned short min_mis, unsigned short max_mis, double ave_mis) : 
	    juncst(jst), juncend(jend), hits(hts), kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), flankcase(fc), flankchr(fcr), strand(sd), synstr(st), 
		min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}
};

bool
compjuncfusion(const JuncForSortFusion& lhs, const JuncForSortFusion& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

void
GenomeScan::WriteJuncFusion(const CONJFusion_HASH_MAP& conj_hash_map, ofstream& ofs)
{
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CONJFusion_HASH_MAP::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		JUNCFusion_HASH_CITER iter_conj = chm_iter->second.begin();

		int j=0;
		while (iter_conj != chm_iter->second.end())
		{
			JUNCFusion_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			while (iter_conj_pre!=iter_conj->second.end())
			{
				int kinds = 0;
				int hits = 0;
				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						++kinds;
						hits += iter_conj_pre->second.p[i];
					}
				}

				double rank = 0;

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						double pi = (double)iter_conj_pre->second.p[i] / (double)hits;
						rank += pi * log(pi);
					}
				}

				if (rank!=0)
					rank = -rank;

				if (rank < m_max_rank)
				{
					iter_conj_pre++;
					continue;
				}

				ofs <<chm_iter->first <<'\t'<< iter_conj_pre->first << '\t' << iter_conj->first<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;
				
				//int juncidlen = (int) log10((double)juncid);

				//string suffix0(8 - juncidlen

				//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
				ofs << hits << "\t+\t"<<iter_conj_pre->first << '\t' << iter_conj->first<< "\t255,0,0\t"<<2<<'\t';

				for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				{
					if (iter_conj_pre->second.p[k - 1] > 0)
					{
						ofs << k << ',';
						break;
					}
				}

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						ofs << m_read_width - i - 1 << ',';
						break;
					}
				}
				ofs << '\t';

				for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
				{
					if (iter_conj_pre->second.p[k - 1] > 0)
					{
						ofs << 0 << ',';
						break;
					}
				}

				for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
				{
					if (iter_conj_pre->second.p[i] > 0)
					{
						ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
						break;
					}
				}
				ofs << '\t';

				ofs << rank << '\t'<< /*iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'<< iter_conj_pre->second.strand <<'\t'<< */iter_conj_pre->second.synstr 
					<<'\t'<<iter_conj_pre->second.min_mismatch<<'\t'<<iter_conj_pre->second.max_mismatch<<'\t'
					<<(double)iter_conj_pre->second.sum_mismatch / (double)hits<<endl;

				++juncid;
				j++;
				iter_conj_pre++;
			}

			iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}


void
GenomeScan::WriteJuncFusion(const CHROM_FUSION_JUNC_HASH_COMB& conj_hash_map, ofstream& ofs)
{
	string juncidstr = "JUNC_";
	size_t juncid = 1;
	
	CHROM_FUSION_JUNC_HASH_COMB::const_iterator chm_iter;
	
	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		FUSION_JUNC_HASH_COMB::const_iterator iter_conj = chm_iter->second.begin();

		int j=0;
		for (;iter_conj != chm_iter->second.end(); ++iter_conj)
		{
			size_t comb_offset = iter_conj->first;

			size_t prefix_end = comb_offset >> THIRTY_TWO;

			size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

			//JUNCFusion_SUFFIX_CITER iter_conj_pre = iter_conj->second.begin();
			//while (iter_conj_pre!=iter_conj->second.end())
			//{
			int kinds = 0;
			int hits = 0;
			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					++kinds;
					hits += iter_conj->second.p[i];
				}
			}

			double rank = 0;

			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					double pi = (double)iter_conj->second.p[i] / (double)hits;
					rank += pi * log(pi);
				}
			}

			if (rank!=0)
				rank = -rank;

			if (rank < m_max_rank)
			{
				//iter_conj_pre++;
				continue;
			}

			ofs <<chm_iter->first <<'\t'<< prefix_end << '\t' << suffix_st<< '\t'<<juncidstr<<juncid<<'\t'/*<< iter_conj_pre->second.prim<<'\t'<< iter_conj_pre->second.flankstr <<'\t'*/;

			//int juncidlen = (int) log10((double)juncid);

			//string suffix0(8 - juncidlen

			//ofs << rank <<'\t' << hits << '\t' << kinds << endl;
			ofs << hits << "\t"<<iter_conj->second.strand1 <<iter_conj->second.strand2 <</*"\t"<<prefix_end << '\t' << suffix_st<<*/ "\t255,0,0\t"<<2<<'\t' << 
				iter_conj->second.max_prefix_len << ',' << iter_conj->second.max_suffix_len  <<"\t0," << static_cast<int>(suffix_st - prefix_end)<<'\t';

			//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
			//{
			//	if (iter_conj_pre->second.p[k - 1] > 0)
			//	{
			//		ofs << k << ',';
			//		break;
			//	}
			//}

			//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
			//{
			//	if (iter_conj_pre->second.p[i] > 0)
			//	{
			//		ofs << m_read_width - i - 1 << ',';
			//		break;
			//	}
			//}
			//ofs << '\t';

			//for (size_t k = iter_conj_pre->second.p.size(); k > 0; --k)
			//{
			//	if (iter_conj_pre->second.p[k - 1] > 0)
			//	{
			//		ofs << 0 << ',';
			//		break;
			//	}
			//}

			//for (size_t i = 0; i < iter_conj_pre->second.p.size(); ++i)
			//{
			//	if (iter_conj_pre->second.p[i] > 0)
			//	{
			//		ofs << int (iter_conj->first - iter_conj_pre->first - m_read_width + i + 1 )<<',';
			//		break;
			//	}
			//}
			//ofs << '\t';

			string flankstr = FlankString(iter_conj->second.flank_bits, iter_conj->second.bads);

			//int flankcase = 0;

			//if (flankstr == "ATAC")
			//	flankcase = 1;
			//else if (flankstr == "CTAC")
			//	flankcase = 2;
			//else if (flankstr == "CTGC")
			//	 flankcase = 3;
			//else if (flankstr == "GCAG")
			//	 flankcase = 4;
			//else if (flankstr == "GTAG")
			//	 flankcase = 5;
			//else if (flankstr == "GTAT")
			//	 flankcase = 6;

			ofs << rank << '\t'<< iter_conj->second.prim<<'\t'<< flankstr <<'\t'/*<< iter_conj_pre->second.strand <<'\t'*/<< iter_conj->second.synstr 
				<<'\t'<<iter_conj->second.min_mismatch<<'\t'<<iter_conj->second.max_mismatch<<'\t'
				<<(double)iter_conj->second.sum_mismatch / (double)hits<<endl;

			++juncid;
			j++;
				//iter_conj_pre++;
			//}

			//iter_conj++;
		}

		cout << chm_iter->first<<'\t'<<"number of loops " << j << endl;
	}
}


void 
GenomeScan::SortJuncFusion(const char* juncfile)
{
	map<string, vector<JuncForSortFusion> > m_p;

    ifstream ifs(juncfile);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], flankcase[100], flankchr[100], strand12[10], synstr[1000];
			int juncst, juncend, prefixend, suffixst, kinds, hits;

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, flankcase, flankchr, strand12, synstr, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSortFusion(juncst, juncend, hits, kinds, blocks, blocksoffset, rank, flankcase, flankchr, strand12, synstr, min_mis, max_mis, ave_mis));
		}

		ifs.close();
	}
	map<string, vector<JuncForSortFusion> >::iterator mpit;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		sort(mpit->second.begin(), mpit->second.end(), compjuncfusion);
	}

	ofstream ofs(juncfile);

	string juncidstr = "JUNC_";
	size_t juncid = 1;

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		vector<JuncForSortFusion>::iterator vit;
		for (vit = mpit->second.begin(); vit != mpit->second.end(); ++vit)
		{
			ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<vit->strand <<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank <<'\t'<< vit->flankcase<<'\t'<< vit->flankchr<<'\t'<< vit->strand<<'\t'<< vit->synstr
				<<'\t'<<vit->min_mismatch << '\t'<< vit->max_mismatch<< '\t'<< vit->ave_mismatch<<endl;
			++juncid;
		}
	} 
}

void 
GenomeScan::SortJunc(const char* juncfile)
{
	map<string, vector<JuncForSort> > m_p;

    ifstream ifs(juncfile);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100], lpq[100];
			int juncst, juncend, prefixend, suffixst, kinds, hits;

			unsigned short min_mis, max_mis;

			double ave_mis;

			sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%hu\t%hu\t%lf", chromname,  &prefixend, &suffixst, 
				juncname, &hits, strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, lpq, &min_mis, &max_mis, &ave_mis);

			string chromstr = chromname;
			m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, kinds, blocks, blocksoffset,	rank, lpq, min_mis, max_mis, ave_mis));
		}

		ifs.close();
	}
	map<string, vector<JuncForSort> >::iterator mpit;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		sort(mpit->second.begin(), mpit->second.end(), compjunc);
	}

	ofstream ofs(juncfile);

	string juncidstr = "JUNC_";
	size_t juncid = 1;

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		string chromfile = m_chrom_files[m_chrom_names_map[mpit->first]];

		string chromseq;

		readchrom(chromfile.c_str(), chromseq);

		if (chromseq == "")
		{
			cout <<"empty chrom: "<<chromfile<<endl;
			exit(1);
		}

		vector<JuncForSort>::iterator vit;
		for (vit = mpit->second.begin(); vit != mpit->second.end(); ++vit)
		{
			string flankstr = chromseq.substr(vit->juncst, 2) + chromseq.substr(vit->juncend - 3, 2);

			for (size_t i = 0; i < flankstr.length(); ++i)
			{
				if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
					flankstr[i] = flankstr[i] + 'A' - 'a';
			}

			int flankcase = 0;

			if (flankstr == "ATAC")
				flankcase = 1;
			else if (flankstr == "CTAC")
				flankcase = 6;
			else if (flankstr == "CTGC")
				 flankcase = 3;
			else if (flankstr == "GCAG")
				 flankcase = 4;
			else if (flankstr == "GTAG")
				 flankcase = 5;
			else if (flankstr == "GTAT")
				 flankcase = 2;

			double il = 1.0 - (((double) (vit->juncend - vit->juncst - m_interal + 1)) / (double (m_max_splice_dist - m_interal + 2))); 

			ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<'+'<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< flankcase <<'\t'<<flankstr <<'\t'<<il<<'\t'<<vit->lpq
				<<'\t' << vit->min_mismatch<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<< endl;
			++juncid;
		}
	} 
}

void
GenomeScan::MergeCoverageBlock(const VEC_COVER_BLOCK& v_coverage_block, vector<unsigned short>& v_merged_coverage_block)
{
	VEC_COVER_BLOCK::const_iterator vcb_iter;
	
	// = v_coverage_block.begin();

	//v_merged_coverage_block.push_back(*vcb_iter);

	//++vcb_iter;
    
	for (vcb_iter = v_coverage_block.begin(); vcb_iter != v_coverage_block.end(); ++vcb_iter)
	{
		for (size_t idx = vcb_iter->blockst; idx <= vcb_iter->blockend; ++idx)
		{
			++v_merged_coverage_block[idx];
		}
		//if (vcb_iter->blockst >= v_merged_coverage_block.back().blockend)
		//	v_merged_coverage_block.push_back(*vcb_iter);
		//else
		//{
		//	if (vcb_iter->blockend < v_merged_coverage_block.back().blockend)
		//	{
		//		v_merged_coverage_block.back()
		//	}
		//	else
		//	{
		//	}
		//}
	}
}

void
GenomeScan::WriteCoverage(const VEC_COVER_BLOCK& v_coverage_block, const string& chromname, ofstream& ofs, const size_t& chromID)
{
	vector<unsigned short> v_merged_coverage_blcok(m_chrom_sizes[chromID] + 5, 0);

	size_t chrom_end = v_merged_coverage_blcok.size();

	MergeCoverageBlock(v_coverage_block, v_merged_coverage_blcok);

	size_t blockst = 1;
	for (size_t idx = 2; idx < chrom_end; ++idx)
	{
		if (v_merged_coverage_blcok[idx] != v_merged_coverage_blcok[blockst])
		{
			ofs<<blockst<<'\t'<<idx<<'\t'<< v_merged_coverage_blcok[blockst]<<endl;
			blockst = idx;
		}
	}	
}

void
GenomeScan::Sites_To_Regions()
{
	for (size_t i = 0; i < m_bests.size(); ++i)
	{
		if (m_bests[i].unique && m_bests[i].score <= m_max_mismatches)
		{
			m_hits.push_back(GenomicRegion(m_chrom_names[m_bests[i].chrom], m_bests[i].lsite_start, m_bests[i].lsite_end,
				m_bests[i].rsite_start, m_bests[i].rsite_end,m_read_names[i], m_bests[i].score, (m_bests[i].strand) ? '+' : '-'));
		}
	}
}

void
GenomeScan::Write_Non_Uniques()
{
	ofstream out(m_ambiguous_file);
	for (size_t i = 0; i < m_read_names.size(); ++i)
		if (!m_bests[i].unique)
			out << m_read_names[i] << endl;
	out.close();
}

int
main(int argc, const char **argv)
{
	try 
	{
		//read arguments
		cout << "init arguments"<<endl;
		CmdlineArgus argulist(argc, argv);

		//init GenomeScan

		cout << "init genomescan"<<endl;
		GenomeScan genomescan(argulist);

		genomescan.CallMapReads();
	}
	catch (const RMAPException &e) 
	{
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	catch (std::bad_alloc &ba) 
	{
		cerr << "ERROR: could not allocate memory" << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;

	return 0;
}
