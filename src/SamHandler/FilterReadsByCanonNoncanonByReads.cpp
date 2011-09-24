#include "SamHandler.h"

///compare to synthetic sam
int
main(int argc, char** argv)
{
	char* mps_sam = argv[1];

	char* junc_file = argv[2];

	char* filtered_canon = argv[3];

	char* filtered_noncanon = argv[4];

	char* fitlered_noncanon_canon = argv[5];

	char* filtered_ins = argv[6];

	char* stat_file = argv[7];

	//size_t maxhits = atoi(argv[4]);

	//char* unique_mapped = argv[5];

	SamHandler mps_samhandler;

	//cout << "read sam " << mps_sam << endl;

	//mps_samhandler.ReadAppendSamFile(mps_sam);

	//mps_samhandler.CollectStat();

	//cout << "sort sam " << mps_sam << endl;

	//mps_samhandler.SortSamByTagnameOffset();

	map<string, vector<JuncForSort> > junc_map;

	cout << "ReadSortedJunction"<<endl;

	mps_samhandler.ReadSortedJunction(junc_file, junc_map);

	//cout << "mark canon noncanon reads mapping"<<endl;

	//mps_samhandler.MarkCanonNoncanon(junc_map);

	cout << "filter canon noncanon " << mps_sam << endl;

	mps_samhandler.FilterCanonNoncaonByReads(mps_sam, filtered_canon, filtered_noncanon, fitlered_noncanon_canon, filtered_ins, stat_file, junc_map);
	
}