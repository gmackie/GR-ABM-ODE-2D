OS := $(shell grep -E 'Ubuntu|Red Hat' /etc/issue)
DIM := 100

SVNREV := $(shell svn info | grep "Last Changed Rev:" | cut -d" " -f4)
CXX := g++
CXXFLAGS += -DSVN_VERSION=\"$(SVNREV)\" -msse -Wall -D__DIM__=$(DIM)


ifeq ($(OS), )
LDFLAGS += -lm -lboost_program_options
else
LDFLAGS += -lm -lboost_program_options-mt
endif

all: release

.PHONY:	release debug profile

debug:	CXXFLAGS += -DDEBUG -DVERBOSE -O0 -g -Wextra
debug:	gr

release:	CXXFLAGS += -DNDEBUG -O3
release:	gr

profile:	CXXFLAGS += -pg -g
profile:	release

# gcc compiler optimizations not included here...
gpu:	CXXFLAGS += -DGPU -fPIC
gpu:	LDFLAGS += -lcudart
gpu:	SRCS += cudautil.cu grdiffusionftcsswapgpu.cu grsimulationgpu.cu
gpu:	gr_gpu

SRCS := paramsbase.cpp params.cpp rand.cpp agent.cpp gr.cpp grsimulationgrid.cpp grgrid.cpp gridcell.cpp macrophage.cpp  grstat.cpp grsimulation.cpp onlinestat.cpp grdiffusion.cpp grdiffusionftcs.cpp grdiffusionbtcs.cpp grdiffusionwrongbtcs.cpp grdiffusionftcsswap.cpp ttest.cpp areatest.cpp mtbtest.cpp tcell.cpp tcytotoxic.cpp tgamma.cpp tregulatory.cpp recruitmentbase.cpp recruitmentprob.cpp recruitmentlnode.cpp recruitmentlnodepure.cpp serialization.cpp tinyxml/tinystr.cpp tinyxml/tinyxml.cpp tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp
OBJS := ${SRCS:.cpp=.o}

# If Makefile or gr.h changes, we need to recompile everything -- We can do this here because object files only depend on one file
$(OBJS):	gr.h Makefile

tinyxml/%.o:	CXXFLAGS+= -DTIXML_USE_STL

%.cu.o:	%.cu gr.h Makefile
	nvcc -arch=sm_11 --compiler-options="$(CXXFLAGS)" -c -o $@ $<

gr_gpu: $(OBJS) cudautil.cu grdiffusionftcsswapgpu.cu.o grsimulationgpu.cu.o grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr1: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr2: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr3: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr4: $(OBJS) grmain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lhs: $(OBJS) lhs.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o tinyxml/*.o gr gr1 gr2 gr3 gr4 lhs

-include Makefile.depend
