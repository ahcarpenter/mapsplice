all: check_input_files AddTagsToSam remove_index_sam index_reads parseCluster cluster FilterBWTByRegions newsam2junc_paired FilterSamByJunc SepSplicedUnspliced convert_to_abs_offset ReadRegions generate_bash_file parseSingle parsePER MapPERall MapPER pairing merge_paired_end_reads read_chromo_size mapsplice_search mapsplice_report dividereads dividereadsindexed mapsplice_segments newsam2junc RemDup FilterMultipleMappedByRead FilterReadsByCanonNoncanonByReads filterjuncbyROCarguNonCanonical junc2bed pairend remove_blankspace_perline juncdb_bwt2sam junc_db fusionsam2junc filterjuncbysmalldeletion filterjuncbyROCargu filterfusionjuncbyminmis filterbyjuncminmis filter_1hits countline separateuniquefusion merge_sam remove_blankspace count_canon_noncanon bowtie2sam

OPTFLAGS = -O3

#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)

check_input_files:
	g++ $(CFLAGS) -o bin/check_input_files src/check_input_files/check_input_files.cpp &>log

remove_index_sam:
	g++ $(CFLAGS) -o bin/remove_index_sam src/remove_index_sam/remove_index_sam.cpp &>log

AddTagsToSam:
	g++ $(CFLAGS) -o bin/AddTagsToSam  src/AddTagsToSam/AddTagsToSam.cpp &>log

index_reads:
	g++ $(CFLAGS) -o bin/index_reads src/index_reads/index_reads.cpp &>log

parseCluster:
	g++ $(CFLAGS) -o bin/parseCluster src/cluster/preparseFusion.cpp &>log

cluster:
	g++ $(CFLAGS) -o bin/cluster src/cluster/PE_match_junction.cpp &>log

FilterBWTByRegions:
	g++ $(CFLAGS) -o bin/FilterBWTByRegions src/FilterBWTByRegions/FilterBWTByRegions.cpp &>log

read_chromo_size:
	g++ $(CFLAGS) -o bin/read_chromo_size src/read_chromo_size/read_chromo_size.cpp &>log

mapsplice_search:
	g++ $(CFLAGS) -o bin/mapsplice_search src/mapsplice_search_pe/search_hole.cpp &>log

mapsplice_report:
	g++ $(CFLAGS) -o bin/mapsplice_report src/mapsplice_report_pe_fusion/mapsplice_report.cpp &>log

dividereadsindexed:
	g++ $(CFLAGS) -o bin/dividereadsindexed src/dividereadsindexed/dividereadsindexed.cpp &>log

dividereads:
	g++ $(CFLAGS) -o bin/dividereads src/dividereads/dividereads.cpp &>log
mapsplice_segments:
	g++ $(CFLAGS) -o bin/mapsplice_segments src/mapsplice_segments/mapsplice_segments.cpp &>log
newsam2junc:
	g++ $(CFLAGS) -o bin/newsam2junc src/newsam2junc/newsam2junc.cpp &>log
RemDup:
	g++ $(CFLAGS) -o bin/RemDup src/RemDup/RemDup.cpp &>log
FilterMultipleMappedByRead:
	g++ $(CFLAGS) -o bin/FilterMultipleMappedByRead src/SamHandler/FilterMultipleMappedByRead.cpp &>log
FilterReadsByCanonNoncanonByReads:
	g++ $(CFLAGS) -o bin/FilterReadsByCanonNoncanonByReads src/SamHandler/FilterReadsByCanonNoncanonByReads.cpp &>log
filterjuncbyROCarguNonCanonical:
	g++ $(CFLAGS) -o bin/filterjuncbyROCarguNonCanonical src/filterjuncbyROCarguNonCanonical/filterjuncbyROCarguNonCanonical.cpp &>log
junc2bed:
	g++ $(CFLAGS) -o bin/junc2bed src/junc2bed/junc2bed.cpp &>log
pairend:
	g++ $(CFLAGS) -o bin/pairend src/pairend/pairend.cpp &>log
remove_blankspace_perline:
	g++ $(CFLAGS) -o bin/remove_blankspace_perline src/remove_blankspace_perline/remove_blankspace_perline.cpp &>log
juncdb_bwt2sam:
	g++ $(CFLAGS) -o bin/juncdb_bwt2sam src/juncdb_bwt2sam/juncdb_bwt2sam.cpp &>log
junc_db:
	g++ $(CFLAGS) -o bin/junc_db src/junc_db_src/junction_seq_construction.cpp &>log
fusionsam2junc:
	g++ $(CFLAGS) -o bin/fusionsam2junc src/fusionsam2junc/fusionsam2junc.cpp &>log
filterjuncbysmalldeletion:
	g++ $(CFLAGS) -o bin/filterjuncbysmalldeletion src/filterjuncbysmalldeletion/filterjuncbysmalldeletion.cpp &>log
filterjuncbyROCargu:
	g++ $(CFLAGS) -o bin/filterjuncbyROCargu src/filterjuncbyROCargu/filterjuncbyROCargu.cpp &>log
filterfusionjuncbyminmis:
	g++ $(CFLAGS) -o bin/filterfusionjuncbyminmis src/filterfusionjuncbyminmis/filterfusionjuncbyminmis.cpp &>log
filterbyjuncminmis:
	g++ $(CFLAGS) -o bin/filterbyjuncminmis src/filterbyjuncminmis/filterbyjuncminmis.cpp &>log
filter_1hits:
	g++ $(CFLAGS) -o bin/filter_1hits src/filter_1hits/filter_1hits.cpp &>log
countline:
	g++ $(CFLAGS) -o bin/countline src/countline/countline.cpp &>log
separateuniquefusion:
	g++ $(CFLAGS) -o bin/separateuniquefusion src/separateuniquefusion/separateuniquefusion.cpp &>log
merge_sam:
	g++ $(CFLAGS) -o bin/merge_sam src/merge_sam/merge_sam.cpp &>log
remove_blankspace:
	g++ $(CFLAGS) -o bin/remove_blankspace src/remove_blankspace/remove_blankspace.cpp &>log
count_canon_noncanon:
	g++ $(CFLAGS) -o bin/count_canon_noncanon src/count_canon_noncanon/count_canon_noncanon.cpp &>log
bowtie2sam:
	g++ $(CFLAGS) -o bin/bowtie2sam src/bowtie2sam/bowtie2sam.cpp &>log
FilterSamByJunc:
	g++ $(CFLAGS) -o bin/FilterSamByJunc src/FilterSamByJunc/FilterSamByJunc.cpp &>log
SepSplicedUnspliced:
	g++ $(CFLAGS) -o bin/SepSplicedUnspliced src/SepSplicedUnspliced/SepSplicedUnspliced.cpp &>log
convert_to_abs_offset:
	g++ $(CFLAGS) -o bin/convert_to_abs_offset src/convert_to_abs_offset/convert_to_abs_offset.cpp &>log
ReadRegions:
	g++ $(CFLAGS) -o bin/ReadRegions src/ReadRegions/ReadRegions.cpp &>log
generate_bash_file:
	g++ $(CFLAGS) -o bin/generate_bash_file_and_run src/generate_bash_file/generate_bash_file_and_run.cpp &>log
parseSingle:
	g++ $(CFLAGS) -o bin/parseSingle src/MapPER/parseSingle.cpp &>log
parsePER:
	g++ $(CFLAGS) -o bin/parsePER src/MapPER/parsePER.cpp &>log
MapPERall:
	g++ $(CFLAGS) -o bin/MapPERall src/MapPER/MapPERall.cpp &>log
MapPER:
	g++ $(CFLAGS) -o bin/MapPER src/MapPER/MapPER.cpp &>log
pairing:
	g++ $(CFLAGS) -o bin/pairing src/pairing/pairing.cpp &>log
merge_paired_end_reads:
	g++ $(CFLAGS) -o bin/merge_paired_end_reads src/merge_paired_end_reads/merge_paired_end_reads.cpp &>log
newsam2junc_paired:
	g++ $(CFLAGS) -o bin/newsam2junc_paired src/newsam2junc_paired/newsam2junc_paired.cpp &>log
clean:
	rm bin/check_input_files bin/index_reads bin/cluster bin/parseCluster bin/FilterBWTByRegions bin/newsam2junc_paired bin/FilterSamByJunc bin/SepSplicedUnspliced bin/convert_to_abs_offset bin/ReadRegions bin/generate_bash_file_and_run bin/parseSingle bin/parsePER bin/MapPERall bin/MapPER bin/pairing bin/merge_paired_end_reads bin/read_chromo_size bin/mapsplice_search bin/mapsplice_report bin/dividereads bin/dividereadsindexed bin/mapsplice_segments bin/newsam2junc bin/RemDup bin/FilterMultipleMappedByRead bin/FilterReadsByCanonNoncanonByReads bin/filterjuncbyROCarguNonCanonical bin/junc2bed bin/pairend bin/remove_blankspace_perline bin/juncdb_bwt2sam bin/junc_db bin/fusionsam2junc bin/filterjuncbysmalldeletion bin/filterjuncbyROCargu bin/filterfusionjuncbyminmis bin/filterbyjuncminmis bin/filter_1hits bin/countline bin/separateuniquefusion bin/merge_sam bin/remove_blankspace bin/count_canon_noncanon bin/bowtie2sam

