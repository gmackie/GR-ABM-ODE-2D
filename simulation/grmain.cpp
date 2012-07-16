/*
 * grmain.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "gr.h"
#include "params.h"
#include "grsimulation.h"
#include "stat.h"
#include "rand.h"
#include <sys/time.h>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include "recruitmentlnode.h"
#include "recruitmentlnodepure.h"
#include "recruitmentprob.h"
#include "recruitmentlnodeproxy.h"
#include <iostream>
#include <sstream>
#include <fstream>

namespace po = boost::program_options;
using namespace std;
using namespace boost::filesystem;

void printVersion()
{
	std::cout << "Version: " << GR_VERSION << std::endl;
}

void printUsage(const char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options] FILE\n" << desc << std::endl;
}

class oCSVStream {
private:
  std::ostream* f;
public:
  oCSVStream(std::ostream* s) : f(s) {}
  virtual ~oCSVStream() { delete f; }
  template<typename T> void operator<<(const T& data) { operator<<(std::string(data)); }

  template<typename T>
  void write(const T& data) { (*this)<<data; }
  void endRow() { (*f)<<std::endl; }
  virtual void saveRow(const GrSimulation& sim) = 0;
};
template<>
void oCSVStream::operator<<(const unsigned& data) { (*f)<<data<<','; }

template<>
void oCSVStream::operator<<(const int& data) { (*f)<<data<<','; }

template<>
void oCSVStream::operator<<(const long& data) { (*f)<<data<<','; }

template<>
void oCSVStream::operator<<(const unsigned long& data) { (*f)<<data<<','; }

template<>
void oCSVStream::operator<<(const float& data) {
  if(!isnan(data)) (*f)<<data<<',';
  else (*f)<<"NaN,";
}

template<>
void oCSVStream::operator<<(const double& data) {
  if(!isnan(data)) (*f)<<data<<',';
  else (*f)<<"NaN,";
}

template<>
void oCSVStream::operator<<(const char& data) { (*f)<<data<<','; }

//RFC 4180 Standard string handling
template<>
void oCSVStream::operator<<(const std::string& data) {
  std::string c(data);
  for(unsigned i=c.find('"'); i<c.size(); i=c.find('"', i+1))
    c.replace(i, 1, "\"\"");
  (*f)<<'"'<<c<<"\",";
}

class IntMtbStats : public oCSVStream {
public:
  IntMtbStats(std::ostream* s, GrSimulation& sim) : oCSVStream(s) { outputHeader(sim); }
  void outputHeader(const GrSimulation& sim) {
    write("time");
    size_t sz = sim.getStats().getIntMtbFreq().size();
    for(unsigned i=0;i<sz;i++){
      stringstream ss;
      ss<<(i+1);
      write(ss.str());
    }
    write("Mi Max"); write("Mi Min"); write("Mi Mean"); write("Mi Median"); write("Mi StdDev");
    write("Mci Max"); write("Mci Min"); write("Mci Mean"); write("Mci Median"); write("Mci StdDev");
    endRow();
  }
  void saveRow(const GrSimulation& sim) {
    const Stats& stat = sim.getStats();
    write(sim.getTime());
    {
      const std::vector<unsigned>& v = stat.getIntMtbFreq();
      for(unsigned i=0;i<v.size();i++) write(v[i]);
    }
    {
      namespace ba = boost::accumulators;
      const Stats::Stat& vi = stat.getMacIntMtbStats(Mac::MAC_INFECTED);
      if(ba::extract::count(vi) > 0){
        write(ba::extract::max(vi));
        write(ba::extract::min(vi));
        write(ba::extract::mean(vi));
        write(ba::extract::median(vi));
        write(sqrt(ba::extract::variance(vi)));
      }
      else { write(NAN); write(NAN); write(NAN); write(NAN); write(NAN); }
      const Stats::Stat& vic = stat.getMacIntMtbStats(Mac::MAC_INFECTED);
      if(ba::extract::count(vic) > 0) {
        write(ba::extract::max(vic));
        write(ba::extract::min(vic));
        write(ba::extract::mean(vic));
        write(ba::extract::median(vic));
        write(sqrt(ba::extract::variance(vic)));
      }
      else { write(NAN); write(NAN); write(NAN); write(NAN); write(NAN); }
    }
    endRow();
  }
};

class GeneralStats : public oCSVStream {
public:
  GeneralStats(std::ostream* s) : oCSVStream(s) { outputHeader(); }
  void outputHeader() {
    write("time");
    write("Mac"); write("Mr"); write("Mi"); write("Mci"); write("Ma"); write("Md");
    write("Tgam"); write("Tgam a"); write("Tgam dr"); write("Tgam ad"); write("Tgam reg");  write("Tgam d");
    write("Tcyt"); write("Tcyt a"); write("Tcyt reg"); write("Tcyt d");
    write("Treg"); write("Treg r"); write("Treg d");
    write("Int. Mtb."); write("Ext. Mtb."); write("repExtMtb"); write("NonRepl Ext. Mtb."); write("Tot Mtb.");
    write("TNF"); write("IntTNFR1"); write("TotMiMa kmRNA"); write("IL10"); write("CCL2"); write("CCL5"); write("CXCL9"); 
    write("AreaTNF"); write("AreaCellDensity"); write("LesionSize");
    write("MDC"); write("N4"); write("TH0"); write("TH1"); write("N8");
    write("T80"); write("T8"); write("TC"); write("TH0lung"); write("TH1lung");
    write("T80lung"); write("T8lung"); write("TClung");
    write("NrSourcesMac"); write("NrSourcesTgam"); write("NrSourcesTcyt"); write("NrSourcesTreg");
    write("NrSourcesMacActivatedFree"); write("NrSourcesTgamActivatedFree"); write("NrSourcesTcytActivatedFree"); write("NrSourcesTregActivatedFree");
    write("NrSourcesMacCrowded"); write("NrSourcesTgamCrowded"); write("NrSourcesTcytCrowded"); write("NrSourcesTregCrowded");
    write("NrCaseated");
    write("MacApoptosisTNF"); write("MrApoptTNF"); write("MiApoptTNF"); write("MciApoptTNF");  write("MaApoptTNF");
    write("TcellApoptTNF");
    write("Fas/FasL Killing");
    write("Cytotoxic Killing");
    write("NrOfMacsFullyInhibited");
    write("TgamQueued"); write("TcytQueued"); write("TregQueued");
    write("TgamQueuedDie"); write("TcytQueuedDie"); write("TregQueuedDie");
    write("TgamRecruited"); write("TcytRecruited"); write("TregRecruited");
    endRow();

  }
  void saveRow(const GrSimulation& sim) {
    const Stats& stats = sim.getStats();

    write(sim.getTime());

    write(stats.getNrOfMacs());
    for(size_t i=0;i<Mac::NSTATES;i++)
      write(stats.getNrOfMacs(Mac::State(i)));
    write(stats.getNrOfTgams());
    for(size_t i=0;i<Tgam::NSTATES;i++)
      write(stats.getNrOfTgams(Tgam::State(i)));
    write(stats.getNrOfTcyts());
    for(size_t i=0;i<Tcyt::NSTATES;i++)
      write(stats.getNrOfTcyts(Tcyt::State(i)));
    write(stats.getNrOfTregs());
    for(size_t i=0;i<Treg::NSTATES;i++)
      write(stats.getNrOfTregs(Treg::State(i)));

    Scalar repExtMtb = stats.getTotExtMtb() - stats.getTotNonRepExtMtb();
    write(stats.getTotIntMtb()); write(stats.getTotExtMtb()); write(repExtMtb); write(stats.getTotNonRepExtMtb()); write((stats.getTotIntMtb() + stats.getTotExtMtb()));

    write(stats.getTotTNF()); write(stats.getTotTNFR1int()); write(stats.getTotkmRNA()); write(stats.getTotIL10()); write(stats.getTotCCL2()); write(stats.getTotCCL5());  write(stats.getTotCXCL9());

    Scalar lesionSize = 2 * sqrt((0.0004 * stats.getAreaCellDensity()) / PI);
    write(stats.getAreaTNF()); write(stats.getAreaCellDensity()); write(lesionSize);

    write(stats.getMDC()); write(stats.getN4()); write(stats.getTH0()); write(stats.getTH1()); write(stats.getN8());
    write(stats.getT80()); write(stats.getT8()); write(stats.getTC());  write(stats.getTH0lung()); write(stats.getTH1lung());
    write(stats.getT80lung()); write(stats.getT8lung()); write(stats.getTClung());

    for(size_t i=0;i<NAGENTS;i++)
      write(stats.getNrSources((AgentType)i));
    for(size_t i=0;i<NAGENTS;i++)
      write(stats.getNrSourcesActive((AgentType)i));
    for(size_t i=0;i<NAGENTS;i++)
      write(stats.getNrSourcesCrowded((AgentType)i));

    write(stats.getNrCaseated());

    int totMacApoptosisTNF[Mac::NSTATES] = {0}; //Just temporary for Mohammed Fallahi
    int sumMacApoptosisTNF = 0;

    for(int i=0;i<Mac::MAC_DEAD;i++){    //Keep a running sum of deaths
      totMacApoptosisTNF[i]+=(stats.getMacApoptosisTNF((Mac::State)i));
      sumMacApoptosisTNF+=totMacApoptosisTNF[i];
    }

    write(sumMacApoptosisTNF);
    for(int i=0;i<Mac::MAC_DEAD;i++)
      write(totMacApoptosisTNF[i]);

    write(stats.getTcellApoptosisTNF());
    write(stats.getApoptosisFasFasL());
    write(stats.getKillCytotoxic());

    write(stats.getNrOfCellsTnfInhibited()/100);

    for(size_t i=TGAM;i<NAGENTS;i++)
      write(stats.getNrQueued((AgentType)i));
    for(size_t i=TGAM;i<NAGENTS;i++)
      write(stats.getNrQueuedDie((AgentType)i));
    for(size_t i=TGAM;i<NAGENTS;i++)
      write(stats.getNrRecruited((AgentType)i));

    oCSVStream::endRow();
  }
};

class MolecularTrackingStats : public oCSVStream {
public:
	MolecularTrackingStats(std::ostream* s) : oCSVStream(s) { outputHeader(); }

	void outputHeader() {
		write("time");
		write("cellID");
		write("cellType");
		write("cellState");

		write("intMtb"); // 0 for T cells

		// TNF associated attributes
		write("mTNF"); write("surfTNFR1"); write("surfTNFR2"); write("surfBoundTNFR1"); write("surfBoundTNFR2");
		write("intBoundTNFR1"); write("intBoundTNFR2"); write("mTNFRNA"); write("vTNFR1"); write("vTNFR2");
		write("kSynth"); write("kTACE"); write("kmRNA");

		// IL10 associated attributes
		write("surfIL10R"); write("vIL10R"); write("surfBoundIL10R"); write("kISynth");

		endRow();
	}

	// Actually save several rows, one for each agent to be tracked.
	void saveRow(const GrSimulation& sim)
	{

		MacList macList = sim.getMacList();
		for (MacList::iterator it = macList.begin(); it != macList.end(); it++)
		{
			if ((*it)->gettrackMolecularDynamics())
			{
				saveAgentRow(sim.getTime(), **it);
			}
		}

		// Don't track T cells for now.
		//TgamList tgamList = sim.getTgamList();
		//TcytList tcytList = sim.getTcytList();
		//TregList tregList = sim.getTregList();
	}

	void saveAgentRow(int time, Agent& agent)
	{
		write(time);
		write(agent.getID());
		write((int) agent.getAgentType());
		write((int) agent.getState());

		if (agent.getAgentType() == MAC)
		{
			Mac& m = dynamic_cast<Mac&>(agent);
			write(m.getIntMtb());
		}
		else
		{
			write(0);
			write(0.0);
		}

		// TNF associated attributes
		write(agent.getmTNF());
		write(agent.getsurfTNFR1());
		write(agent.getsurfTNFR2());
		write(agent.getsurfBoundTNFR1());
		write(agent.getsurfBoundTNFR2());
		write(agent.getintBoundTNFR1());
		write(agent.getintBoundTNFR2());
		write(agent.getmTNFRNA());
		write(agent.getvTNFR1());
		write(agent.getvTNFR2());
		write(agent.getkSynth());
		write(agent.getkTACE());
		write(agent.getkmRNA());

		// IL10 associated attributes
		write(agent.getsurfIL10R());
		write(agent.getvIL10R());
		write(agent.getsurfBoundIL10R());
		write(agent.getkISynth());

		oCSVStream::endRow();
  }
};


void saveState(const GrSimulation* pSim, int time, std::string dir=std::string("./")){
    int days, hours, minutes;
    char fname[20];
    assert(pSim);
    GrSimulation::convertSimTime(time, days, hours, minutes);
    sprintf(fname, "%03dd%02dh%02dm.state.gz", days, hours, minutes);
 
    namespace bio = boost::iostreams;
    bio::filtering_ostream out;
    out.push(bio::gzip_compressor());
    out.push(bio::file_sink((dir+std::string(fname)).c_str()));
    if(!out.good())
        std::cerr<<"Unable to open output file: "<<fname<<endl;
    else
    {
    	pSim->serialize(out);
    }
}

void printStats(const GrSimulation* pSim) {
  const Stats& stats = pSim->getStats();
  size_t sz = pSim->getGrid().getSize();
  printf("%-7d ", pSim->getTime());
  size_t i=0;
  printf("%3d - (", stats.getNrOfMacs());
  for(i=0;i<Mac::NSTATES-1;i++)
    printf("%d,", stats.getNrOfMacs(Mac::State(i)));
  printf("%d) ", stats.getNrOfMacs(Mac::State(i)));
  printf("%3d - (", stats.getNrOfTgams());
  for(i=0;i<Tgam::NSTATES-1;i++)
    printf("%d,", stats.getNrOfTgams(Tgam::State(i)));
  printf("%d) ", stats.getNrOfTgams(Tgam::State(i)));
  printf("%3d - (", stats.getNrOfTcyts());
  for(i=0;i<Tcyt::NSTATES-1;i++)
    printf("%d,", stats.getNrOfTcyts(Tcyt::State(i)));
  printf("%d) ", stats.getNrOfTcyts(Tcyt::State(i)));
  printf("%3d - (", stats.getNrOfTregs());
  for(i=0;i<Treg::NSTATES-1;i++)
    printf("%d,", stats.getNrOfTregs(Treg::State(i)));
  printf("%d) ", stats.getNrOfTregs(Treg::State(i)));

  printf("(%.5f, %.5f) ", stats.getTotExtMtb(), stats.getTotIntMtb());
  #define INV_SZ (1.0 / sz)
  printf("(%11.5f,%11.5f,%11.5f,%11.5f,%11.5f) ", stats.getTotTNF()*INV_SZ, stats.getTotIL10()*INV_SZ, stats.getTotCCL2()*INV_SZ, stats.getTotCCL5()*INV_SZ, stats.getTotCXCL9()*INV_SZ);
  printf("(");
  for(i=0;i<NAGENTS-1;i++)
    printf("%d,", stats.getNrSources((AgentType)i));
  printf("%d) ", stats.getNrSources((AgentType)i));
  printf("(%d, %.5f) ", stats.getNrCaseated(), stats.getTotNonRepExtMtb());
  printf("(%d)\n", stats.getNrOfCellsTnfInhibited()/100);
}

// Stopping criteria.
// Used to stop a simulation early when doing an lhs.
// Don't stop if the stopping time steps are 0, which should be the default,
// so old parameter files without these parameters will work without stopping prematurely.
bool shouldStop(int time, GrSimulation* pSim)
{
	const Stats& stats = pSim->getStats();
	Scalar totMtb = stats.getTotExtMtb() + stats.getTotIntMtb();
	int areaCellDensity = stats.getAreaCellDensity();

  // Stop if too large
	if ( (_PARAM(PARAM_MTB_STOPPING_TIME_STEP) > 0 && time == _PARAM(PARAM_MTB_STOPPING_TIME_STEP) && totMtb < _PARAM(PARAM_MTB_STOPPING_THRESHOLD) ) ||
		 (_PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_TIME_STEP) && time == _PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_TIME_STEP) && areaCellDensity < _PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_THRESHOLD)) )
    return true;
  // Stop if too small
	if ( (_PARAM(PARAM_MTB_STOPPING_TIME_STEP2) > 0 && time == _PARAM(PARAM_MTB_STOPPING_TIME_STEP2) && totMtb > _PARAM(PARAM_MTB_STOPPING_THRESHOLD2) ) ||
		 (_PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_TIME_STEP2) && time == _PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_TIME_STEP2) && areaCellDensity > _PARAM(PARAM_AREA_CELL_DENSITY_STOPPING_THRESHOLD2)) )
    return true;

	return false;
}

void performOutput(GrSimulation* pSim, int stateInterval, int csvInterval, bool screenDisplay, int timeToSimulate, std::string outputDir, std::vector<oCSVStream*> csvStreams, int time)
{
    // If not doing chunks then always save at the last time step.
    // If doing chunks only save at the last time step for the last chunk.
    // Except if stopping early then always save at the time step when we stop.
    bool saveAtLastTimeStep = (time == timeToSimulate || shouldStop(time, pSim));

    if (stateInterval > 0 && (time % stateInterval == 0  || saveAtLastTimeStep) )
    {
        saveState(pSim, time, outputDir);
    }

    if (time % csvInterval == 0 || saveAtLastTimeStep)  //Write out each csv file
      for(unsigned i=0;i<csvStreams.size();i++)
        csvStreams[i]->saveRow(*pSim);

    if (screenDisplay)
      printStats(pSim);
}

void run(GrSimulation* pSim, int stateInterval, int csvInterval, bool screenDisplay, int timeToSimulate, std::string outputDir, std::vector<oCSVStream*> csvStreams)
{
  if(screenDisplay)
  	cout << endl << "--seed " << g_Rand.getSeed() << endl;
  csvInterval = csvInterval < 1 ? 1 : csvInterval;

  // Needed for chunks. Any script that manages running chunks needs to know if the simulation stopped during a chunk
  // prior to the last chunk, so it knows not to run any subsequent chunks for a run.
  //bool earlyStop = false;

  // Needed when starting a simulation from a saved state.
  int time = pSim->getTime();

  // Only perform output before the simulation loop if this is time 0.
  // We don't want to repeat output for a time step when doing a simulation in chunks.
  if (time == 0)
  {
          performOutput(pSim, stateInterval, csvInterval, screenDisplay, timeToSimulate, outputDir, csvStreams, time);
  }

  // ++time so we don't repeat the simulation for the time step of a saved state - the last time step
  // performed for a prior chunk.
  for (++time; time <= timeToSimulate; time++)
  {
		pSim->solve();
		performOutput(pSim, stateInterval, csvInterval, screenDisplay, timeToSimulate, outputDir, csvStreams, time);

	//    gettimeofday(&end, NULL);
	//    acc( (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.0);
	//    printf("%7d,%10.4lf,%7d\n", pSim->getTime(), be::rolling_mean(acc), pSim->getStats().getNrOfAgents());

	    //Check stopping criteria.
	    if (shouldStop(time, pSim))
	    {
	      //earlyStop = true;
	      break;
	    }
  }
}
void buildSim(GrSimulation* pSim, DiffusionMethod diffMethod, RecruitmentMethod recrMethod, int odeMethod , bool tnfrDynamics, bool il10rDynamics, bool nfkbDynamics,
                bool adaptive, int tnfDepletionTimeStep, int il10DepletionTimeStep, bool tgammatransition, float areaTNFThreshold, float areaCellDensityThreshold) {
  
    cout << "\nODE SOLVER" << std::endl;
    cout << "----------" << std::endl;
    cout << "   " << odeMethod << "\n" << std::endl;
    
	pSim->setTnfrDynamics(tnfrDynamics || nfkbDynamics); // when NFkB is turned on, tnfr dynamics will be on automatically.
    
    cout << "Tunable Resolution" << std::endl;
    cout << "------------------" << std::endl;
    
    if (tnfrDynamics == 1 || nfkbDynamics == 1) {
        cout << "TNF  Dynamics   -  On" << std::endl;
    }
    else
    {
        cout << "TNF  Dynamics   -  Off" << std::endl;
    }
    pSim->setIl10rDynamics(il10rDynamics);
    if (il10rDynamics == 1) {
        cout << "IL10 Dynamics   -  On" << std::endl;
    }
    else
    {
        cout << "IL10 Dynamics   -  Off" << std::endl;
    }
    pSim->setNfkbDynamics(nfkbDynamics);
    if (nfkbDynamics == 1) {
        cout << "NFkB Dynamics   -  On" << std::endl;
    }
    else
    {
        cout << "NFkB Dynamics   -  Off" << std::endl;
    }
	
    pSim->setTgammaTransition(tgammatransition);
    
    if (tgammatransition == 1) {
        cout << "Treg Induction  -  On" << std::endl;
    }
    else
    {
        cout << "Treg Induction  -  Off" << std::endl;
    }

    pSim->setAdaptive(adaptive);
    if (adaptive == 1) {
        cout << "Adaptive Dynamics   -  On" << std::endl;
    }
    else
    {
        cout << "Adaptive Dynamics   -  Off" << std::endl;
    }
    
    pSim->setTnfDepletionTimeStep(tnfDepletionTimeStep);
    pSim->setIl10DepletionTimeStep(il10DepletionTimeStep);
	pSim->setRecruitmentMethod(recrMethod);
	pSim->setDiffusionMethod(diffMethod);

	// Parameters must be loaded, since since the base lymph ODE class, RecruitmentLnODE, uses parameters in its constructor.
	/* set recruitment method */
	pSim->setRecruitmentMethod(recrMethod);

    pSim->setODESolverMethod((ODESolvers::ODEMethod)odeMethod);
	
	//	Set area thresholds if specified on the command line.
	if (areaTNFThreshold >= 0)
		pSim->setAreaThreshold(areaTNFThreshold);
	
	if (areaCellDensityThreshold >= 0)
		pSim->setAreaThresholdCellDensity(areaCellDensityThreshold);
}

int main(int argc, char** argv)
{
  std::cout << std::endl;
  printVersion();

  unsigned long seed;
  size_t dim;
  std::string paramFile;
  std::string outputDir;

  double molecularTrackingRadius;
  std::string molecularTrackingFileName;
  std::string track_ids;
	
	/*
	 * Set seed to current time, in case not specified
	 * The seed adjustment is used for cluster runs.
	 * When a large number of jobs are submitted to a cluster
	 * some will start at essentially the same time - they will
	 * have the same value for curTime and so have the same seed.
	 * This is especially bad if they are using the same parameter file.
	 * For cluster jobs the seed adjustment command line argument
	 * is the cluster job number for that job which is xored with
	 * curTime to define a unique seed.
	 *
	 * */
	
	po::options_description general("General");
	general.add_options()
	("help,h", "Help message")
	("version,v", "Version number")
	("quiet,q", "Suppress printing statistics on the console.\n"
	            "This option is implied by --lhs.")
	("output-dir,o", po::value(&outputDir)->default_value("./"), "Output directory")
	("input-file,i", po::value(&paramFile), "Parameter file")
	("seed,s", po::value(&seed), "RNG seed, default is based on the current time");

  po::options_description stats("Statistics");
  stats.add_options()
  ("moi", "Generate csv file of internal MTB statistics, saved in -moi.csv")
  ("stats", "Generate csv file of general statistics, saved in .csv")
  ("csv-interval", po::value<unsigned>()->default_value(1), "CSV update interval (10 min timesteps)")
  ("molecular-track-radius", po::value(&molecularTrackingRadius)->default_value(0.0), "Radius from center of grid of initial cells to track molecular dynamics. 0 means don't track any cells.")
  ("molecular-track-ids", po::value(&track_ids), "Comma-seperated list of ids to track")
  ("molecular-track-file", po::value<std::string>(&molecularTrackingFileName), "File name to hold molecular dynamics cell tracking data");

  po::options_description sim_opts("Simulation");
  sim_opts.add_options()
  ("dim,d", po::value(&dim)->default_value(100), "Size of simulation grid")
  ("load,l", po::value<std::string>(), "Load from a saved state")
	("timesteps,t", po::value<unsigned>(), "Number of time steps to simulate\n"
									"Takes precedence over --days")
	("days", po::value<unsigned>()->default_value(200), "Number of days to simulate")
	("state-interval", po::value<unsigned>()->default_value(0), "State save interval (10 min timesteps)")
	("recr", po::value<unsigned>()->default_value(0), "recruitment:\n0 - probability\n1 - lymph node ode proxy\n2 - lymph node ode pure")
	("diffusion", po::value<unsigned>()->default_value(4),
	 "Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap\n4 - ADE Grid Swap")
	("area-tnf-threshold", po::value<float>()->default_value(0.5),"Threshold for granuloma area defined by TNF, in the range [0.0, 1.0]\n")
	("area-cell-density-threshold", po::value<float>()->default_value(0.5),"Threshold for granuloma area defined by cell density, in the range [0.0, 1.0]");

  po::options_description simdyn_opts("Simulation Dynamics");
  simdyn_opts.add_options()
    ("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
    ("il10r-dynamics", "Use molecular level IL10/IL10R dynamics in the model")
    ("NFkB-dynamics", "Use molecular level intracellular NFkB dynamics in the model")
    ("adaptive", "Use adaptive time-step ODE Solvers")
    ("odesolver", po::value<int>()->default_value(0),
     "ODE Solver Method:\n0 - Forward Euler\n1 - Euler Predictor-Corrector\n2 - Runge-Kutta 3rd Order\n3 - Runge-Kutta 4th Order\n4 - HuenEuler\n5 - Runge-Kutta Cache-Karp\n6 - Runge-Kutta Fehlberg\n7 - Bogacki-Shampine")
    ("Treg-induction", "Allow Tregs to be induced from Tgams in the model")
    ("tnf-depletion", po::value<int>()->default_value(-1), "The time step at which to stop secreting tnf, including by tnfr dynamics. -1: no depletion")
    ("il10-depletion", po::value<int>()->default_value(-1), "The time step at which to stop secreting il10, including by il10r dynamics. -1: no depletion");

	
  po::options_description lhs_opts("LHS");
  lhs_opts.add_options()
  ("lhs", "Running as part of an LHS run")
  ("seedadj", po::value<unsigned int>(), "Seed adjustment for cluster runs");

  po::options_description desc("Allowed Options");
  desc.add(general).add(stats).add(sim_opts).add(simdyn_opts).add(lhs_opts);
	
  po::variables_map vm;
	try
	{
		po::positional_options_description p;
		p.add("input-file", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		printUsage(argv[0], desc);
		return 1;
	}

  std::cout<<"GRID SIZE: NROWS: "<<dim<<" NCOLS: "<<dim<<std::endl;
  Pos pdim(dim, dim);
  
  if (vm.count("version"))
  {
	// Nothing to do - the version is always printed above.
    return 0;
  }
  
  if (vm.count("help"))
  {
    printUsage(argv[0], desc);
    return 0;
  }

  std::vector<size_t> ids;
  {
    std::stringstream ss(track_ids);
    string id;
    while(std::getline(ss, id, ','))
      ids.push_back(atoi(id.c_str()));
  }

  unsigned timeToSimulate;
  if (!vm.count("timesteps"))
    timeToSimulate = TIME_STEPS_PER_DAY * vm["days"].as<unsigned>();
  else 
    timeToSimulate = vm["timesteps"].as<unsigned>();
  
   if (!vm.count("seed"))
    seed = createTimeSeed();
  
  // Adjust the seed if a seed adjustment was specified.
  if (vm.count("seedadj"))
    seed = seed ^ (vm["seedadj"].as<unsigned int>());
  
  bool screenDisplay = !(vm.count("quiet") || vm.count("lhs"));

  // Must be done before making GrSimulation.
  // Also must be done before creating a lymph ODE recruitment object,
  // since the base lymph ODE class, RecruitmentLnODE, uses parameters in its constructor.
  if (!Params::getInstance(pdim)->fromXml(paramFile.c_str()))
    throw std::runtime_error("Unable to get parameters from file, cannot continue...");
  
  Params::getInstance()->setParam(PARAM_TNFODE_EN, vm.count("tnfr-dynamics") || vm.count("NFkB-dynamics"));
  Params::getInstance()->setParam(PARAM_IL10ODE_EN, vm.count("il10r-dynamics"));
  Params::getInstance()->setParam(PARAM_NFKBODE_EN, vm.count("NFkB-dynamics"));

  DiffusionMethod diffMethodEnum;
  switch (vm["diffusion"].as<unsigned>())
  {
    case 0:
      diffMethodEnum = DIFF_REC_EQ;
      break;
    case 3:
      diffMethodEnum = DIFF_REC_EQ_SWAP;
      break;
    case 4:
      diffMethodEnum = DIFF_ADE_SWAP;
      break;
    default:
      std::cerr<<"Unsupported Diffusion method: " << vm["diffusion"].as<unsigned>() << std::endl;
      printUsage(argv[0], desc);
      exit(1);
  }

  RecruitmentMethod recrMethod;

  switch (vm["recr"].as<unsigned>())
  {
    case 0:
    	recrMethod = RECR_PROB;
      break;
    case 1:
    	recrMethod = RECR_LN_ODE_PROXY;
      break;
    case 2:
    	recrMethod = RECR_LN_ODE_PURE;
      break;
    default:
      std::cerr<<"Unsupported recruitment method: " << vm["recr"].as<unsigned>() << std::endl;
      printUsage(argv[0], desc);
      exit(1);
  }


  GrSimulation* pSim = new GrSimulation(pdim);
  assert(pSim != NULL);
  if(vm.count("load")) {
    namespace bio = boost::iostreams;

    std::string savedf = vm["load"].as<std::string>();
	if (!exists(savedf))
	{
		cerr << "File to load, '" << savedf << "' does not exist." << endl;
		exit(1);
	}


    size_t found = string::npos;
    bio::filtering_istream in;
    if((found = savedf.rfind(".")) != string::npos && savedf.substr(found) == ".gz")
      in.push(bio::gzip_decompressor());   //Compressed?
    bio::file_source fileSource(savedf);
    in.push(fileSource);

    if(!in)
      throw std::runtime_error("Failed to open saved state file");
    pSim->deserialize(in);
  }

  buildSim(pSim, diffMethodEnum, recrMethod, vm["odesolver"].as<int>(), vm.count("tnfr-dynamics"), vm.count("il10r-dynamics"), vm.count("NFkB-dynamics"), vm.count("adaptive"), 
           vm["tnf-depletion"].as<int>(), vm["il10-depletion"].as<int>(),vm.count("Treg-induction"), vm["area-tnf-threshold"].as<float>(),
              vm["area-cell-density-threshold"].as<float>());

  if (!vm.count("load")){
    g_Rand.setSeed(seed);
    pSim->init();
  }

  pSim->initMolecularTracking(molecularTrackingRadius);
  pSim->initMolecularTracking(ids);

  if(outputDir[outputDir.size()-1] != '/')
    outputDir += '/';

  // Test if the output directory exists.
  if (!exists(outputDir))
  {
      std::cerr<<"Output directory '" << outputDir << "' doesn't exist."<<endl;
      exit(1);
  }


  // Write the seed to a file, so runs can be repeated.
  // This must come after any load of a saved state, so the correct seed is written,
  // since that loads a saved seed.
  std::string seedStreamName = outputDir + "seed";
  std::ofstream seedStream(seedStreamName.c_str());
  seedStream << g_Rand.getSeed() << std::endl;
  seedStream.close();

  std::vector<oCSVStream*> csvStreams;
  if(vm.count("moi")){
    stringstream ss;
    ss<<outputDir<<"moi"<<g_Rand.getSeed()<<".csv";
    ofstream* f = new ofstream(ss.str().c_str());
    csvStreams.push_back(new IntMtbStats(f, *pSim));
  }

  if(vm.count("stats")){
    stringstream ss;
    ss<<outputDir<<"seed"<<g_Rand.getSeed()<<".csv";
    ofstream* f = new ofstream(ss.str().c_str());
    csvStreams.push_back(new GeneralStats(f));
  }

  if(molecularTrackingRadius > 0 || ids.size() > 0)
  {
	  if (!molecularTrackingFileName.empty())
	  {
		stringstream ss;
		ss << outputDir << molecularTrackingFileName << g_Rand.getSeed() <<".csv";
		ofstream* f = new ofstream(ss.str().c_str());
		csvStreams.push_back(new MolecularTrackingStats(f));
	  }
  }

  run(pSim, vm["state-interval"].as<unsigned>(), vm["csv-interval"].as<unsigned>(),
      screenDisplay, timeToSimulate, outputDir, csvStreams);

  for(unsigned i=0;i<csvStreams.size();i++)
    delete csvStreams[i];
  delete pSim;
  return 0;
}
