OS := $(shell grep -E -s 'Ubuntu|Red Hat' /etc/issue)
HOSTNAME := $(shell hostname)

SVNREV := $(shell svn info | grep "Last Changed Rev:" | cut -d" " -f4)
CXX := g++


ifndef BOOST_PREFIX
BOOST_PREFIX := /usr
endif


override CXXFLAGS += -DSVN_VERSION="$(SVNREV)" -Wall -W -Wno-strict-aliasing -Wno-unknown-pragmas
override LDFLAGS += -lm -lz -lpthread

ifeq ($(HOSTNAME:flux-%=flux), flux)
override CXXFLAGS += -I$(BOOST_INC)
override LDFLAGS += -L$(BOOST_LINK) -lboost_program_options -lboost_iostreams -lboost_serialization -lboost_filesystem -lboost_system
else
override CXXFLAGS += -I$(BOOST_PREFIX)/include

BOOST_LIBS := boost_program_options boost_iostreams boost_serialization boost_filesystem boost_system
# So far the best way to eleviate the issues with different distros not properly symlinking boost libraries.
# If static building, might need to modify to search only .o or .a compiled libraries (kept generic to deal with osx & windows)
# <- Could completely solve this problem with cmake.... >_>
BOOST_FOUND_LIBS := $(foreach lib, $(BOOST_LIBS), $(if $(wildcard ${BOOST_PREFIX}/lib/*$(lib)*-mt*), $(lib)-mt, $(lib)))
override LDFLAGS += -L${BOOST_PREFIX}/lib $(addprefix -l, $(BOOST_FOUND_LIBS))
endif

all: release

.PHONY:	release debug profile valgrindO0 valgrindO1 valgrindO3

debug:	CXXFLAGS += -DDEBUG -DVERBOSE -O0 -g -Wextra
debug:	gr

release:	CXXFLAGS +=-O3
release:	gr

openmp:  CXXFLAGS += -DNDEBUG -O3 -fopenmp
openmp:  gr

profile:	CXXFLAGS += -pg -g
profile:	release

valgrindO0:	CXXFLAGS += -g -DNDEBUG -O0
valgrindO0:	gr

valgrindO1:	CXXFLAGS += -g -DNDEBUG -O1
valgrindO1:	gr

valgrindO3:	CXXFLAGS += -g -DNDEBUG -O3
valgrindO3:	gr

SRCS := params.cpp rand.cpp agent.cpp gr.cpp grsimulationgrid.cpp grgrid.cpp macrophage.cpp grsimulation.cpp onlinestat.cpp grdiffusion.cpp grdiffusionbtcs.cpp grdiffusionwrongbtcs.cpp grdiffusionadeswap.cpp grdiffusionftcsswap.cpp grvascular.cpp ttest.cpp areatest.cpp mtbtest.cpp tcell.cpp tcytotoxic.cpp tgamma.cpp tregulatory.cpp recruitmentbase.cpp recruitmentprob.cpp recruitmentlnode.cpp recruitmentlnodeproxy.cpp recruitmentlnodepure.cpp
OBJS := ${SRCS:.cpp=.o}

# If Makefile or gr.h changes, we need to recompile everything -- We can do this here because object files only depend on one file
$(OBJS):	gr.h Makefile

gr: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lhs:	params.o rand.o lhsmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

depend:	.depend
.depend:
	$(CXX) $(CXXFLAGS) -MM $(SRCS) > .depend

clean:
	rm -f ${OBJS} ${TINYXML_OBJS} grmain.o lhs.o gr lhs .depend

-include .depend
