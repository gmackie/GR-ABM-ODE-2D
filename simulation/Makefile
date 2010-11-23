
ifeq ($(CXX),)
	CXX = g++
endif

ifeq ($(DEBUG),yes)
	LDFLAGS = -lm -lboost_program_options-mt
	CXXFLAGS = -Wall -g3 -DDEBUG -DVERBOSE
else
	LDFLAGS = -lm -lboost_program_options-mt
	CXXFLAGS = -Wall -DNDEBUG -ffast-math -fcaller-saves -finline-functions -O3
endif

all: gr

gr: params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr1: params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr2: params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr3: params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gr4: params.o rand.o agent.o gr.o grsimulationgrid.o grgrid.o gridcell.o macrophage.o grmain.o grstat.o grsimulation.o onlinestat.o grdiffusion.o grdiffusionftcs.o grdiffusionbtcs.o grdiffusionwrongbtcs.o grdiffusionftcsswap.o ttest.o areatest.o mtbtest.o tcell.o tcytotoxic.o tgamma.o tregulatory.o recruitmentbase.o recruitmentprob.o recruitmentlnode.o recruitmentlnodepure.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

lhs: gr.o params.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o lhs.o rand.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

par2csv: gr.o params.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o par2csv.o rand.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

#Makefile.depend: *.h *.cpp Makefile
#	$(CC) -M *.cpp > Makefile.depend

clean:
	rm -f *.o tinyxml/*.o gr gr1 gr2 gr3 gr4 lhs par2csv

-include Makefile.depend
