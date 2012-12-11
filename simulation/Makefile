OS := $(shell grep -E 'Ubuntu|Red Hat' /etc/issue)

SVNREV := $(shell svn info | grep "Last Changed Rev:" | cut -d" " -f4)
CXX := g++


ifndef BOOST_PREFIX
BOOST_PREFIX := /usr
endif
override CXXFLAGS += -DSVN_VERSION="$(SVNREV)" -Wall -W -I$(BOOST_PREFIX)/include -Wno-strict-aliasing -Wno-unknown-pragmas

BOOST_LIBS := boost_program_options boost_iostreams boost_serialization boost_filesystem boost_system
# So far the best way to eleviate the issues with different distros not properly symlinking boost libraries.
# If static building, might need to modify to search only .o or .a compiled libraries (kept generic to deal with osx & windows)
# <- Could completely solve this problem with cmake.... >_>
BOOST_FOUND_LIBS := $(foreach lib, $(BOOST_LIBS), $(if $(wildcard ${BOOST_PREFIX}/lib/*$(lib)*-mt*), $(lib)-mt, $(lib)))
override LDFLAGS += -L${BOOST_PREFIX}/lib $(addprefix -l, $(BOOST_FOUND_LIBS)) -lm -lz -lpthread

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

# gcc compiler optimizations not included here...
gpu:	CXXFLAGS += -DGPU -fPIC
gpu:	LDFLAGS += -lcudart
gpu:	SRCS += cudautil.cu grdiffusionftcsswapgpu.cu grsimulationgpu.cu
gpu:	gr_gpu

SRCS := paramsbase.cpp params.cpp rand.cpp agent.cpp gr.cpp grsimulationgrid.cpp grgrid.cpp macrophage.cpp grsimulation.cpp onlinestat.cpp grdiffusion.cpp grdiffusionbtcs.cpp grdiffusionwrongbtcs.cpp grdiffusionadeswap.cpp grdiffusionftcsswap.cpp grvascular.cpp ttest.cpp areatest.cpp mtbtest.cpp tcell.cpp tcytotoxic.cpp tgamma.cpp tregulatory.cpp recruitmentbase.cpp recruitmentprob.cpp recruitmentlnode.cpp recruitmentlnodeproxy.cpp recruitmentlnodepure.cpp
OBJS := ${SRCS:.cpp=.o}
TINYXML_SRCS := tinyxml/tinystr.cpp tinyxml/tinyxml.cpp tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp
TINYXML_OBJS := ${TINYXML_SRCS:.cpp=.o}

# If Makefile or gr.h changes, we need to recompile everything -- We can do this here because object files only depend on one file
$(OBJS):	gr.h Makefile

tinyxml/%.o:	CXXFLAGS+= -DTIXML_USE_STL

%.cu.o:	%.cu gr.h Makefile
	nvcc -arch=sm_11 --compiler-options="$(CXXFLAGS)" -c -o $@ $<

gr_gpu: $(TINYXML_OBJS) $(OBJS) cudautil.cu grdiffusionftcsswapgpu.cu.o grsimulationgpu.cu.o grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr: $(TINYXML_OBJS) $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr1: $(TINYXML_OBJS) $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr2: $(TINYXML_OBJS) $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr3: $(TINYXML_OBJS) $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr4: $(TINYXML_OBJS) $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lhs:	$(TINYXML_OBJS) rand.o paramsbase.o params.o lhs.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

depend:	.depend
.depend:
	$(CXX) $(CXXFLAGS) -MM $(SRCS) > .depend

clean:
	rm -f ${OBJS} ${TINYXML_OBJS} grmain.o lhs.o gr gr1 gr2 gr3 gr4 lhs .depend

-include .depend
