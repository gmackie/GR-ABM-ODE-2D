/*
 * grmain.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "gr.h"
#include "params.h"
#include "grsimulation.h"
#include "grstat.h"
#include "rand.h"
#include <sys/time.h>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "recruitmentlnode.h"
#include "recruitmentlnodepure.h"
#include "recruitmentprob.h"
#include <iostream>
#include <sstream>
#include <fstream>

namespace po = boost::program_options;
using namespace std;

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
    size_t sz;
    sim.getStats().getIntMtbFreq(sz);
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
    const GrStat& stat = sim.getStats();
    write(sim.getTime());
    size_t sz;
    {
      const unsigned* v = stat.getIntMtbFreq(sz);
      for(unsigned i=0;i<sz;i++) write(v[i]);
    }
    {
      namespace ba = boost::accumulators;
      const GrStat::Stat* v = stat.getIntMtbStats(sz);
      if(ba::extract::count(v[MAC_INFECTED]) > 0){
        write(ba::extract::max(v[MAC_INFECTED]));
        write(ba::extract::min(v[MAC_INFECTED]));
        write(ba::extract::mean(v[MAC_INFECTED]));
        write(ba::extract::median(v[MAC_INFECTED]));
        write(sqrt(ba::extract::variance(v[MAC_INFECTED])));
      }
      else { write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); }
      if(ba::extract::count(v[MAC_CINFECTED]) > 0) {
        write(ba::extract::max(v[MAC_CINFECTED]));
        write(ba::extract::min(v[MAC_CINFECTED]));
        write(ba::extract::mean(v[MAC_CINFECTED]));
        write(ba::extract::median(v[MAC_INFECTED]));
        write(sqrt(ba::extract::variance(v[MAC_CINFECTED])));
      }
      else { write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); write(0.0/0.0); }
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
    write("Tgam"); write("Tgam a"); write("Tgam reg");  write("Tgam d");
    write("Tcyt"); write("Tcyt a"); write("Tcyt reg"); write("Tcyt d");
    write("Treg"); write("Treg r"); write("Treg d");
    write("Int. Mtb."); write("Ext. Mtb."); write("repExtMtb"); write("NonRepl Ext. Mtb."); write("Tot Mtb.");
	write("TNF"); write("IL10"); write("CCL2"); write("CCL5"); write("CXCL9"); 
    write("AreaTNF"); write("AreaCellDensity"); write("LesionSize");
    write("MDC"); write("N4"); write("TH0"); write("TH1"); write("N8");
    write("T80"); write("T8"); write("TC"); write("TH0lung"); write("TH1lung");
    write("T80lung"); write("T8lung"); write("TClung");
    write("Outcome (1)"); write("Outcome (2)");
    write("NrSourcesMac"); write("NrSourcesTgam"); write("NrSourcesTcyt"); write("NrSourcesTreg");
    write("NrCaseated");
    write("MacApoptosisTNF"); write("MrApoptTNF"); write("MiApoptTNF"); write("MciApoptTNF");  write("MaApoptTNF");
    write("TcellApoptTNF");
    write("MrActivationTNF");
    write("MiActivationTNF");
	write("NrOfMacsFullyInhibited");
    endRow();

  }
  void saveRow(const GrSimulation& sim) {
    const GrStat& stats = sim.getStats();

    write(sim.getTime());
    write(stats.getNrOfMac()); write(stats.getNrOfMacResting()); write(stats.getNrOfMacInfected()); write(stats.getNrOfMacCInfected()); write(stats.getNrOfMacActive()); write(stats.getNrOfMacDead());
    write(stats.getNrOfTgam()); write(stats.getNrOfTgamActive()); write(stats.getNrOfTgamDownRegulated()); write(stats.getNrOfTgamDead());
    write(stats.getNrOfTcyt()); write(stats.getNrOfTcytActive()); write(stats.getNrOfTcytDownRegulated()); write(stats.getNrOfTcytDead());
    write(stats.getNrOfTreg()); write(stats.getNrOfTregActive()); write(stats.getNrOfTregDead());

    FLOAT_TYPE repExtMtb = stats.getTotExtMtb() - stats.getTotNonRepExtMtb();
    write(stats.getTotIntMtb()); write(stats.getTotExtMtb()); write(repExtMtb); write(stats.getTotNonRepExtMtb()); write((stats.getTotIntMtb() + stats.getTotExtMtb()));

	  write(stats.getTotTNF()); write(stats.getTotIL10()); write(stats.getTotCCL2()); write(stats.getTotCCL5());  write(stats.getTotCXCL9());

    FLOAT_TYPE lesionSize = 2 * sqrt((0.0004 * stats.getAreaCellDensity()) / PI);
    write(stats.getAreaTNF()); write(stats.getAreaCellDensity()); write(lesionSize);

    write(stats.getMDC()); write(stats.getN4()); write(stats.getTH0()); write(stats.getTH1()); write(stats.getN8());
    write(stats.getT80()); write(stats.getT8()); write(stats.getTC());  write(stats.getTH0lung()); write(stats.getTH1lung());
    write(stats.getT80lung()); write(stats.getT8lung()); write(stats.getTClung());

    for (int i = 0; i < NOUTCOMES; i++)
    {
      switch (stats.getGrStatus(i))
      {
        case GR_CLEARANCE:
          write("Clearance");
          break;
        case GR_CONTAINMENT:
          write("Containment");
          break;
        case GR_CONTAINMENT_INCONSISTENT:
          write("Containment?");
          break;
        case GR_DISSEMINATION:
          write("Dissemination");
          break;
        case GR_DISSEMINATION_INCONSISTENT:
          write("Dissemination?");
          break;
        case GR_UNKNOWN:
          write("Unknown");
          break;
        case GR_NONE:
          write("None");
          break;
      }
    }
    write(stats.getNrSourcesMac()); write(stats.getNrSourcesTgam());  write(stats.getNrSourcesTcyt()); write(stats.getNrSourcesTreg());
    write(stats.getNrCaseated());

    int startState = 1; // Skip the dead state: apoptosis doesn't occur for an already dead mac.
    static int totMacApoptosisTNF[NMAC_STATES] = {0}; //Just temporary for Mohammed Fallahi
    int sumMacApoptosisTNF = 0;

    //
    for(int i=startState;i<NMAC_STATES;i++){    //Keep a running sum of deaths
      totMacApoptosisTNF[i]+=(stats.getNrMacApoptosisTNF((MacState)i));
      sumMacApoptosisTNF+=totMacApoptosisTNF[i];
    }

    write(sumMacApoptosisTNF);
    for(int i=startState;i<NMAC_STATES;i++)
      write(totMacApoptosisTNF[i]);

    write(stats.getNrTcellApoptosisTNF());
    write(stats.getNrRestingMacActivationTNF());
    write(stats.getNrInfMacActivationTNF());
	write(stats.getNrOfCellsInhibited()/100);
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
    else pSim->serialize(out);
}

void printStats(const GrSimulation* pSim) {
  const GrStat& stats = pSim->getStats();
  printf("%-7d ", pSim->getTime());
  printf("%3d - (%d,%d,%d,%d,%d) ", stats.getNrOfMac(), stats.getNrOfMacResting(), stats.getNrOfMacInfected(), stats.getNrOfMacCInfected(), stats.getNrOfMacActive(), stats.getNrOfMacDead());
  printf("%3d - (%d,%d,%d) ", stats.getNrOfTgam(), stats.getNrOfTgamActive(), stats.getNrOfTgamDownRegulated(), stats.getNrOfTgamDead());
  printf("%3d - (%d,%d,%d) ", stats.getNrOfTcyt(), stats.getNrOfTcytActive(), stats.getNrOfTcytDownRegulated(), stats.getNrOfTcytDead());
  printf("%3d - (%d,%d) ", stats.getNrOfTreg(), stats.getNrOfTregActive(), stats.getNrOfTregDead());
  printf("(%.5f, %.5f) ", stats.getTotExtMtb(), stats.getTotIntMtb());
  #define INV_SZ (1.0 / (NROWS*NCOLS))
  printf("(%11.5f,%11.5f,%11.5f,%11.5f,%11.5f) ", stats.getTotTNF()*INV_SZ, stats.getTotIL10()*INV_SZ, stats.getTotCCL2()*INV_SZ, stats.getTotCCL5()*INV_SZ, stats.getTotCXCL9()*INV_SZ);
  printf("(%d,%d,%d,%d) ", stats.getNrSourcesMac(), stats.getNrSourcesTgam(), stats.getNrSourcesTcyt(), stats.getNrSourcesTreg());
  printf("(%d, %.5f) ", stats.getNrCaseated(), stats.getTotNonRepExtMtb());
  printf("(%d)\n", stats.getNrOfCellsInhibited()/100);
}

// Stopping criteria.
// Used to stop a simulation early when doing an lhs.
// Don't stop if the stopping time steps are 0, which should be the default,
// so old parameter files without these parameters will work without stopping prematurely.
bool shouldStop(int time, GrSimulation* pSim)
{
	const GrStat& stats = pSim->getStats();
	FLOAT_TYPE totMtb = stats.getTotExtMtb() + stats.getTotIntMtb();
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

void run(GrSimulation* pSim, int stateInterval, int csvInterval, bool screenDisplay, int timeToSimulate, std::string outputDir, std::vector<oCSVStream*> csvStreams, bool lhs)
{
  if(screenDisplay)
  	cout << endl << "--seed " << g_Rand.getSeed() << endl;
  csvInterval = csvInterval < 1 ? 1 : csvInterval;
  for (int time = 0; time <= timeToSimulate; time += 1)
  {
    //if (!lhs) fprintf(stderr, "Done: %3.2f%%\r", 100.0*(time / float(timeToSimulate)));
    // Display and write output at the requested interval, and after the last time step.
    if (stateInterval > 0 && time % stateInterval == 0)
       saveState(pSim, time, outputDir);
    if (time % csvInterval == 0 || time == timeToSimulate)  //Write out each csv file
      for(unsigned i=0;i<csvStreams.size();i++)
        csvStreams[i]->saveRow(*pSim);
    if (screenDisplay)
      printStats(pSim);
    //Run the pSimulation one step
    if(time != timeToSimulate)
      pSim->solve();

    //Check stopping criteria.
    if (lhs && shouldStop(time, pSim))
    {
    	return;
    }
  }
}
void buildSim(GrSimulation* pSim, DiffusionMethod diffMethod, RecruitmentBase* pRecr, bool tnfrDynamics, bool il10rDynamics,
              bool nfkbDynamics, int tnfDepletionTimeStep, int il10DepletionTimeStep, float areaTNFThreshold, float areaCellDensityThreshold) {
  
	pSim->setTnfrDynamics(tnfrDynamics || nfkbDynamics); // when NFkB is turned on, tnfr dynamics will be on automatically.
    
    cout << "Tunable Resolution" << std::endl;
    cout << "------------------" << std::endl;
    
    if (tnfrDynamics == 1 || nfkbDynamics == 1) {
        cout << "TNF  Dynamics  -  On" << std::endl;
    }
    else
    {
        cout << "TNF  Dynamics  -  Off" << std::endl;
    }
    pSim->setIl10rDynamics(il10rDynamics);
    if (il10rDynamics == 1) {
        cout << "IL10 Dynamics  -  On" << std::endl;
    }
    else
    {
        cout << "IL10 Dynamics  -  Off" << std::endl;
    }
    pSim->setNfkbDynamics(nfkbDynamics);
    if (nfkbDynamics == 1) {
        cout << "NFkB Dynamics  -  On" << std::endl;
    }
    else
    {
        cout << "NFkB Dynamics  -  Off" << std::endl;
    }
	pSim->setTnfDepletionTimeStep(tnfDepletionTimeStep);
    pSim->setIl10DepletionTimeStep(il10DepletionTimeStep);
	pSim->setRecruitment(pRecr);
	pSim->setDiffusionMethod(diffMethod);
	
	//	Set area thresholds if specified on the command line.
	if (areaTNFThreshold >= 0)
		pSim->setAreaThreshold(areaTNFThreshold);
	
	if (areaCellDensityThreshold >= 0)
		pSim->setAreaThresholdCellDensity(areaCellDensityThreshold);
}

int main(int argc, char** argv)
{
  printVersion();
  std::cout << "GRID SIZE: NROWS: " << NROWS << " NCOLS: " << NCOLS << std::endl;
  unsigned long seed;
  std::string paramFile;
  std::string outputDir;
	
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
	struct timeval curTimeHiRes;
  if (gettimeofday(&curTimeHiRes, NULL) != 0)
		throw runtime_error("Error getting current high resolution time.");
	
	po::options_description general("General");
	general.add_options()
	("help,h", "Help message")
	("version,v", "Version number")
	("quiet,q", "Suppress printing statistics on the console.\n"
	            "This option is implied by --lhs.")
	("output-dir,o", po::value<std::string>(&outputDir)->default_value("./"), "Output directory")
	("input-file,i", po::value<std::string>(&paramFile), "Parameter file")
	("seed,s", po::value<unsigned long>(&seed), "RNG seed, default is based on the current time");

  po::options_description stats("Statistics");
  stats.add_options()
  ("moi", "Generate csv file of internal MTB statistics, saved in -moi.csv")
  ("stats", "Generate csv file of general statistics, saved in .csv")
	("csv-interval", po::value<unsigned>()->default_value(1), "CSV update interval (10 min timesteps)");

  po::options_description sim_opts("Simulation");
  sim_opts.add_options()
  ("load,l", po::value<std::string>(), "Load from a saved state")
	("timesteps,t", po::value<unsigned>(), "Number of time steps to simulate\n"
                                    "Takes precedence over --days")
	("days", po::value<unsigned>()->default_value(200), "Number of days to simulate")
	("state-interval", po::value<unsigned>()->default_value(0), "State save interval (10 min timesteps)")
	("diffusion,d", po::value<unsigned>()->default_value(3),
	 "Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap")
	("area-tnf-threshold", po::value<float>()->default_value(0.5),"Threshold for granuloma area defined by TNF, in the range [0.0, 1.0]\n")
	("area-cell-density-threshold", po::value<float>()->default_value(0.5),"Threshold for granuloma area defined by cell density, in the range [0.0, 1.0]");

  po::options_description simdyn_opts("Simulation Dynamics");
  simdyn_opts.add_options()
	("ode", "Use integrated lymph node ODE for recruitment")
	("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
    ("il10r-dynamics", "Use molecular level IL10/IL10R dynamics in the model")
	("NFkB-dynamics", "Use molecular level intracellular NFkB dynamics in the model")
	("tnf-depletion", po::value<int>()->default_value(-1), "The time step at which to stop secreting tnf, including by tnfr dynamics. -1: no depletion")
    ("il10-depletion", po::value<int>()->default_value(-1), "The time step at which to stop secreting il10, including by il10r dynamics. -1: no depletion");
	
  po::options_description lhs_opts("LHS");
  lhs_opts.add_options()
  ("lhs", "Running as part of an LHS run")
	("seedadj", po::value<unsigned long>(), "Seed adjustment for cluster runs");

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
  
  if (vm.count("version"))
  {
    printVersion();
    return 0;
  }
  
  if (vm.count("help"))
  {
    printUsage(argv[0], desc);
    return 0;
  }
  unsigned timeToSimulate;
  if (!vm.count("timesteps"))
    timeToSimulate = TIME_STEPS_PER_DAY * vm["days"].as<unsigned>();
  else 
    timeToSimulate = vm["timesteps"].as<unsigned>();
  
  // If we just use time in seconds we will get duplicate seeds when doing a large number of runs,
  // such as running a big LHS on a cluster.
  // If we use just the micro-second part of the high resolution time we only get seeds that vary
  // between 1 and 1,000,000.
  // So we XOR the two times, which gives a larger range of seeds, with the lower order part varying
  // with a higher resolution than we would get just by using seconds.
  if (!vm.count("seed"))
    seed = curTimeHiRes.tv_sec ^ curTimeHiRes.tv_usec;
  
  // Adjust the seed if a seed adjustment was specified.
  if (vm.count("seedadj"))
    seed = seed ^ (vm["seedadj"].as<unsigned long>());
  
  bool lhs = vm.count("lhs");

  // Write the seed to a file, so runs can be repeated, except for lhs runs.
  bool screenDisplay = !(vm.count("quiet") || vm.count("lhs"));
  if (!lhs)
  {
    std::ofstream seedStream("seed");
    seedStream << seed << std::endl;
    seedStream.close();
  }
  
  RecruitmentBase* pRecr = (vm.count("ode") ? (RecruitmentBase*)new RecruitmentLnODEPure() : (RecruitmentBase*)new RecruitmentProb());
  
  DiffusionMethod diffMethodEnum;
  switch (vm["diffusion"].as<unsigned>())
  {
    case 0:
      diffMethodEnum = DIFF_REC_EQ;
      break;
    case 3:
      diffMethodEnum = DIFF_REC_EQ_SWAP;
      break;
    default:
      std::cerr<<"Unsupported Diffusion method"<<std::endl;
      printUsage(argv[0], desc);
      exit(1);
  }

  if (!Params::getInstance(true)->fromXml(paramFile.c_str())) //Must be done before making GrSimulation
    throw std::runtime_error("Unable to get parameters from file, cannot continue...");

  GrSimulation* pSim = new GrSimulation();
  assert(pSim != NULL);
  if (vm.count("load")){
	std::string s = vm["load"].as<std::string>();
    std::ifstream f(s.c_str());

    if(!f)
    {
    	std::cerr << "Saved state " << s << " does not exist." << std::endl;
    	exit(1);
    }

    pSim->deserialize(f);
    f.close();
  }

  buildSim(pSim, diffMethodEnum, pRecr, vm.count("tnfr-dynamics"), vm.count("il10r-dynamics"), vm.count("NFkB-dynamics"),
           vm["tnf-depletion"].as<int>(), vm["il10-depletion"].as<int>(), vm["area-tnf-threshold"].as<float>(),
              vm["area-cell-density-threshold"].as<float>());
    
  if (!vm.count("load")){
    g_Rand.setSeed(seed);
    pSim->init();
  }
  if(outputDir[outputDir.size()-1] != '/')
    outputDir += '/';

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

  run(pSim, vm["state-interval"].as<unsigned>(), vm["csv-interval"].as<unsigned>(),
      screenDisplay, timeToSimulate, outputDir, csvStreams, lhs);

  for(unsigned i=0;i<csvStreams.size();i++)
    delete csvStreams[i];
  delete pSim;
  delete pRecr;
  return 0;
}
