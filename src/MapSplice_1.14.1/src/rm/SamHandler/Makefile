PROGS = FilterReadsByCanonNoncanonByReads FilterMultipleMappedByRead FilterReadsByCanonNoncanon filterbyunsplicedmapped merge_sort_sam compare2synthesam SamHandlerSam2junc filterbyanchor filterbyrepeated_reads filternotinunmapped filterbyisland filterbyintronlenhmer filterbysmallexon FilterMultipleMapped WriteStat Compare2Sam separate_spliced_unspliced
CXX = g++
CFLAGS = -Wall -lm
OPTFLAGS = -O3
DEBUGFLAGS = -g

#ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
#endif

ifndef NOPT
CFLAGS += $(OPTFLAGS)
endif

all:	$(PROGS)

%:	%.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	@-rm -f $(PROGS)

.PHONY: clean