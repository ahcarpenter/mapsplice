PROGS = newsam2junc
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