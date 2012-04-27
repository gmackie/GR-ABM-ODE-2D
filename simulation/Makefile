OS := $(shell grep -E 'Ubuntu|Red Hat' /etc/issue)

ifndef BOOST_PREFIX
BOOST_PREFIX := /usr
endif

SVNREV := $(shell svn info | grep "Last Changed Rev:" | cut -d" " -f4)
CXX := g++

override CXXFLAGS += -DSVN_VERSION=\"$(SVNREV)\" -Wall -I$(BOOST_PREFIX)/include -Wno-strict-aliasing
override LDFLAGS += -L$(BOOST_PREFIX)/lib -lm -lboost_program_options -lboost_iostreams -lboost_serialization

all: release

.PHONY:	release debug profile valgrindO0 valgrindO1 valgrindO3

debug:	CXXFLAGS += -DDEBUG -DVERBOSE -O0 -g -Wextra
debug:	gr

release:	CXXFLAGS += -DNDEBUG -O3
release:	gr

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

SRCS := paramsbase.cpp params.cpp rand.cpp agent.cpp gr.cpp grsimulationgrid.cpp grgrid.cpp macrophage.cpp grsimulation.cpp onlinestat.cpp grdiffusion.cpp grdiffusionbtcs.cpp grdiffusionwrongbtcs.cpp grdiffusionadeswap.cpp grdiffusionftcsswap.cpp ttest.cpp areatest.cpp mtbtest.cpp tcell.cpp tcytotoxic.cpp tgamma.cpp tregulatory.cpp recruitmentbase.cpp recruitmentprob.cpp recruitmentlnode.cpp recruitmentlnodeproxy.cpp recruitmentlnodepure.cpp serialization.cpp 
OBJS := ${SRCS:.cpp=.o}
TINYXML_SRC := tinyxml/tinystr.cpp tinyxml/tinyxml.cpp tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp
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

clean:
	rm -f *.o tinyxml/*.o gr gr1 gr2 gr3 gr4 lhs

-include Makefile.depend
