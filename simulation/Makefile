UBUNTU = $(shell uname -a | egrep "Ubuntu")

SVNREV = $(shell svn info | grep "Last Changed Rev:" | cut -d" " -f4)
CXX := g++
CXXFLAGS += $(GPUFLAG) -DSVN_VERSION="\"\\\"$(SVNREV)\\\"\"" -msse -Wall


ifeq ($(UBUNTU), )
LDFLAGS += -lm -lboost_program_options
else
LDFLAGS += -lm -lboost_program_options-mt
endif

all: release

.PHONY:	release debug profile

debug:
	@$(MAKE) CXXFLAGS="$(CXXFLAGS) $(PROFILE_FLAGS) -DDEBUG -DVERBOSE -O0 -g3" gr
release:
	@$(MAKE) CXXFLAGS="$(CXXFLAGS) -DNDEBUG -O3" gr
	
# PROFILE_FLAGS instead of using  CXXFLAGS="$(CXXFLAGS) -pg" avoids problems
# with the quotes surrounding the version number getting stripped on the recursive
# make invokation.
profile:
	@$(MAKE) PROFILE_FLAGS="-pg" debug
	
gr: paramsbase.o params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o serialization.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr1: paramsbase.o params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o serialization.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr2: paramsbase.o params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o serialization.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr3: paramsbase.o params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o serialization.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr4: paramsbase.o params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o serialization.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lhs: gr.o paramsbase.o params.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o lhs.o rand.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

#Makefile.depend: *.h *.cpp Makefile
#	$(CC) -M *.cpp > Makefile.depend

clean:
	rm -f *.o tinyxml/*.o gr gr1 gr2 gr3 gr4 lhs

-include Makefile.depend
