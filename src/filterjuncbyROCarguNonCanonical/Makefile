#  mapsplice.cpp
#  Mapsplice
#
#  Created by Kai Wang on 09/09/09.
#  Author: Kai Wang, Jinze Liu, Jan F. Prins
#  Copyright 2009 Kai Wang. All rights reserved.


PROGS = filterjuncbyROCarguNonCanonical

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
