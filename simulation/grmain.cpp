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

void printUsage(char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

void writeOutputHeaderOld(std::ofstream& outputFileStream, std::string inputFileName)
{
	outputFileStream <<inputFileName<<","<<"Mac" <<","<<"Mr"<<","<<"Mi"<<","<<"Mci"<<","<<"Ma"<<","<<"Tgam"<<
	","<<"Tgam a"<<","<<"Tgam reg"<<","<<"Tcyt"<<","<<"Tcyt a"<<","<<"Tcyt reg"<<","<<"Treg"<<
	","<<"Int. Mtb."<<","<<"Ext. Mtb."<<","<<"TNF"<<","<<"CCL2"<<","<<"CCL5"<<
	","<<"CXCL9"<<","<<"Mac Apop"<<","<< "MDC"<< ","<< "N4"<< ","
	<< "TH0"<< ","<< "TH1"<< ","<< "N8"<< ","<< "T80"<< ","<< "T8"<< ","
	<< "TC"<< "," << "TH0lung" << ","<< "TH1lung"<< ","<< "T80lung"<< ","
	<< "T8lung" << "," << "TClung"<< std::endl;
}

void writeOutputOld(std::ofstream& outputFileStream, GrSimulation& sim, int csvInterval)
{
	const GrStat& stats = sim.getStats();
	
	outputFileStream <<sim.getTime()/csvInterval <<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
	","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
	","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
	","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
	","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
	","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrMacApoptosisTNF() <<
	","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
	","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
	","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
	","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
	
}

void writeOutputHeader(std::ofstream& outputFileStream, std::string inputFileName)
{
	outputFileStream << "\"time\""
	<< ','
	<< "\"Mac\""
	<< ','
	<< "\"Mr\""
	<< ','
	<< "\"Mi\""
	<< ','
	<< "\"Mci\""
	<< ','
	<< "\"Ma\""
	<< ','
	<< "\"Md\""
	<< ','
	<< "\"Tgam\""
	<< ','
	<< "\"Tgam a\""
	<< ','
	<< "\"Tgam reg\""
	<< ','
	<< "\"Tgam d\""
	<< ','
	<< "\"Tcyt\""
	<< ','
	<< "\"Tcyt a\""
	<< ','
	<< "\"Tcyt reg\""
	<< ','
	<< "\"Tcyt d\""
	<< ','
	<< "\"Treg\""
	<< ','
	<< "\"Treg r\""
	<< ','
	<< "\"Treg d\""
	<< ','
	<< "\"Int. Mtb.\""
	<< ','
	<< "\"Ext. Mtb.\""
	<< ','
	<< "\"NonRepl Ext. Mtb.\""
	<< ','
	<< "\"Tot Mtb.\""
	<< ','
	<< "\"TNF\""
	<< ','
	<< "\"CCL2\""
	<< ','
	<< "\"CCL5\""
	<< ','
	<< "\"CXCL9\""
	<< ','
	<< "\"AreaTNF\""
	<< ','
	<< "\"AreaCellDensity\""
	<< ','
	<< "\"MDC\""
	<< ','
	<< "\"N4\""
	<< ','
	<< "\"TH0\""
	<< ','
	<< "\"TH1\""
	<< ','
	<< "\"N8\""
	<< ','
	<< "\"T80\""
	<< ','
	<< "\"T8\""
	<< ','
	<< "\"TC\""
	<< ','
	<< "\"TH0lung\""
	<< ','
	<< "\"TH1lung\""
	<< ','
	<< "\"T80lung\""
	<< ','
	<< "\"T8lung\""
	<< ','
	<< "\"TClung\""
	<< ','
	<< "\"Outcome (1)\""
	<< ','
	<< "\"Outcome (2)\""
	<< ','
	<< "\"NrSourcesMac\""
	<< ','
	<< "\"NrSourcesTgam\""
	<< ','
	<< "\"NrSourcesTcyt\""
	<< ','
	<< "\"NrSourcesTreg\""
	<< ','
	<< "\"NrCaseated\""
	<< ','
	<< "\"MacApoptosisTNF\""
	<< ','
	<< "\"MrApoptTNF\""
	<< ','
	<< "\"Mi&MciApoptTNF\""
	<< ','
	<< "\"MaApoptTNF\""
	<< ','
	<< "\"TcellApoptTNF\""
	<< ','
	<< "\"MrActivationTNF\""
	<< ','
	<< "\"MiActivationTNF\""
	<< "\n";
}

void writeOutput(std::ofstream& outputFileStream, GrSimulation& sim, int csvInterval)
{
	const GrStat& stats = sim.getStats();
	outputFileStream << sim.getTime()
	<< ','
	<< stats.getNrOfMac()
	<< ','
	<< stats.getNrOfMacResting()
	<< ','
	<< stats.getNrOfMacInfected()
	<< ','
	<< stats.getNrOfMacCInfected()
	<< ','
	<< stats.getNrOfMacActive()
	<< ','
	<< stats.getNrOfMacDead()
	<< ','
	<< stats.getNrOfTgam()
	<< ','
	<< stats.getNrOfTgamActive()
	<< ','
	<< stats.getNrOfTgamDownRegulated()
	<< ','
	<< stats.getNrOfTgamDead()
	<< ','
	<< stats.getNrOfTcyt()
	<< ','
	<< stats.getNrOfTcytActive()
	<< ','
	<< stats.getNrOfTcytDownRegulated()
	<< ','
	<< stats.getNrOfTcytDead()
	<< ','
	<< stats.getNrOfTreg()
	<< ','
	<< stats.getNrOfTregActive()
	<< ','
	<< stats.getNrOfTregDead()
	<< ','
	<< stats.getTotIntMtb()
	<< ','
	<< stats.getTotExtMtb()
	<< ','
	<< stats.getTotNonRepExtMtb()
	<< ','
	<< (stats.getTotIntMtb() + stats.getTotExtMtb())
	<< ','
	<< stats.getTotTNF()
	<< ','
	<< stats.getTotCCL2()
	<< ','
	<< stats.getTotCCL5()
	<< ','
	<< stats.getTotCXCL9()
	<< ','
	<< stats.getArea()
	<< ','
	<< stats.getAreaCellDensity()
	<< ','
	<< stats.getMDC()
	<< ','
	<< stats.getN4()
	<< ','
	<< stats.getTH0()
	<< ','
	<< stats.getTH1()
	<< ','
	<< stats.getN8()
	<< ','
	<< stats.getT80()
	<< ','
	<< stats.getT8()
	<< ','
	<< stats.getTC()
	<< ','
	<< stats.getTH0lung()
	<< ','
	<< stats.getTH1lung()
	<< ','
	<< stats.getT80lung()
	<< ','
	<< stats.getT8lung()
	<< ','
	<< stats.getTClung();
	
	for (int i = 0; i < NOUTCOMES; i++)
	{
		outputFileStream << ',';
		
		switch (stats.getGrStatus(i))
		{
			case GR_CLEARANCE:
				outputFileStream << "\"Clearance\"";
				break;
			case GR_CONTAINMENT:
				outputFileStream << "\"Containment\"";
				break;
			case GR_CONTAINMENT_INCONSISTENT:
				outputFileStream << "\"Containment?\"";
				break;
			case GR_DISSEMINATION:
				outputFileStream << "\"Dissemination\"";
				break;
			case GR_DISSEMINATION_INCONSISTENT:
				outputFileStream << "\"Dissemination?\"";
				break;
			case GR_UNKNOWN:
				outputFileStream << "\"Unknown\"";
				break;
			case GR_NONE:
				outputFileStream << "\"None\"";
				break;
		}
	}
	
	outputFileStream
	<< ','
	<< stats.getNrSourcesMac()
	<< ','
	<< stats.getNrSourcesTgam()
	<< ','
	<< stats.getNrSourcesTcyt()
	<< ','
	<< stats.getNrSourcesTreg()
	<< ','
	<< stats.getNrCaseated()
	<< ','
	<< stats.getNrMacApoptosisTNF()
	<< ','
	<< stats.getNrRestingMacApoptosisTNF()
	<< ','
	<< stats.getNrInfAndCinfMacApoptosisTNF()
	<< ','
	<< stats.getNrActivatedMacApoptosisTNF()
	<< ','
	<< stats.getNrTcellApoptosisTNF()
	<< ','
	<< stats.getNrRestingMacActivationTNF()
	<< ','
	<< stats.getNrInfMacActivationTNF()
	;
	
	outputFileStream << std::endl;
}
void saveState(const GrSimulation* pSim, int time, std::string dir=std::string(".")){
    int days, hours, minutes;
    char fname[17];
    assert(pSim);
    GrSimulation::convertSimTime(time, days, hours, minutes);
    sprintf(fname, "%03dd%02dh%02dm.state", days, hours, minutes);
    std::ofstream out(((dir+'/')+std::string(fname)).c_str(), std::ios_base::trunc);
    if(!out.good())
        std::cerr<<"Unable to open output file: "<<fname<<endl;
    else pSim->serialize(out);
    out.close();    //Don't want to flush here, would slow things down a lot
}


int run(unsigned long seed, const std::string& inputFileName, const std::string& outputFileName,
		int csvInterval, int stateInterval, bool screenDisplay, DiffusionMethod diffMethod, RecruitmentBase* pRecr, bool ode,
		bool tnfrDynamics, bool nfkbDynamics, int tnfDepletionTimeStep, int timeToSimulate, bool lhs,
		float areaTNFThreshold, float areaCellDensityThreshold)
{
	printVersion();
	std::cout << endl << "--seed " << seed << std::endl;
	
	if (!Params::getInstance(true)->fromXml(inputFileName.c_str()))
		return 1;
	
	g_Rand.setSeed(seed);
	
	// If an output file is requested, construct the complete output file name,
	// open the output file and write a header line.
	std::ofstream outputFileStream;
	
	if (outputFileName.size() > 0)
	{
    	std::ostringstream ofn;
        ofn << outputFileName << (lhs ? "/seed" : "-") << seed << ".csv";
        outputFileStream.open(ofn.str().c_str());
		
		writeOutputHeader(outputFileStream, inputFileName);
	}
	
	GrSimulation sim;     //Allocated on heap for larger grids, shouldn't effect runtimes
	
	sim.setTnfrDynamics(tnfrDynamics || nfkbDynamics); // when NFkB is turned on, tnfr dynamics will be on autamatically.
	sim.setNfkbDynamics(nfkbDynamics);
	sim.setTnfDepletionTimeStep(tnfDepletionTimeStep);
	sim.setRecruitment(pRecr);
	
	const GrStat& stats = sim.getStats();
	
	sim.setDiffusionMethod(diffMethod);
	
	//	Set area thresholds if specified on the command line.
	if (areaTNFThreshold >= 0)
	{
		sim.setAreaThreshold(areaTNFThreshold);
	}
	
	if (areaCellDensityThreshold >= 0)
	{
		sim.setAreaThresholdCellDensity(areaCellDensityThreshold);
	}
	
	sim.init();
	
	for (int time = 0; time <= timeToSimulate; time += 1)
	{
		
		// Display and write output at the requested interval, and after the last time step.
        if (stateInterval > 0 && time % stateInterval == 0)
            saveState(&sim, time, lhs ? outputFileName : ".");
		if (csvInterval <= 0 || time % csvInterval == 0 || time == timeToSimulate)
		{
			if (screenDisplay)
			{
				//			printf("%d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\n",
				printf("%d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\t(%d,%d,%d,%d)\t%d %f\n",
					   sim.getTime(),
					   stats.getNrOfMac(), stats.getNrOfMacResting(), stats.getNrOfMacInfected(),
					   stats.getNrOfMacCInfected(), stats.getNrOfMacActive(), stats.getNrOfMacDead(),
					   stats.getNrOfTgam(), stats.getNrOfTgamActive(),	stats.getNrOfTgamDownRegulated(),
					   stats.getNrOfTgamDead(), stats.getNrOfTcyt(), stats.getNrOfTcytActive(),
					   stats.getNrOfTcytDownRegulated(), stats.getNrOfTcytDead(),
					   stats.getNrOfTreg(), stats.getNrOfTregActive(), stats.getNrOfTregDead(),
					   stats.getTotExtMtb(), stats.getTotIntMtb(),
					   stats.getTotTNF() / (NROWS*NCOLS), stats.getTotCCL2() / (NROWS*NCOLS),
					   stats.getTotCCL5() / (NROWS*NCOLS), stats.getTotCXCL9() / (NROWS*NCOLS),
					   stats.getNrSourcesMac(), stats.getNrSourcesTgam(), stats.getNrSourcesTcyt(), stats.getNrSourcesTreg(),
					   stats.getNrCaseated(), stats.getTotNonRepExtMtb()
					   );
			}
			
        	if (outputFileStream.good())
			{
				writeOutput(outputFileStream, sim, csvInterval);
			}
		}
        //Run the simulation one step
		sim.solve();
	}
	
	return 0;
}

int main(int argc, char** argv)
{
	unsigned long seed;
	std::string inputFileName;
	std::string outputFileName;
	int csvInterval; // In time steps
	bool screenDisplay; // Whether or not to show statistics on the console.
	
	float areaTNFThreshold = -1;
	float areaCellDensityThreshold = -1;
	
	std::string lymphNodeODE;
	std::string lymphNodeTemp;
	int diffMethod;
    int stateInterval;
	int nDays;
    int timeToSimulate;
	bool ode;
	bool tnfrDynamics;
	bool nfkbDynamics;
	int tnfDepletionTimeStep;
	bool lhs;
	
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
	int hiResTimeResult = gettimeofday(&curTimeHiRes, NULL);
	if (hiResTimeResult == -1)
	{
		cerr <<" Error getting current high resolution time."<< endl;
		exit(1);
	}
	
	unsigned long seedadj;
	
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "Help message")
	("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
	("output-dir,o", po::value<std::string>(&outputFileName), "Output directory")
	("csv-interval", po::value<int>(&csvInterval)->default_value(1),
	 "CSV update interval (10 min timesteps)")
	("state-interval", po::value<int>(&stateInterval)->default_value(-1),
	 "State save interval (10 min timesteps)")
	("no-screen-display", "Suppress printing statistics on the console.\n"
	 "This option is implied by --lhs.")
	("seed,s", po::value<unsigned long>(&seed), "Seed")
	("seedadj", po::value<unsigned long>(&seedadj), "Seed adjustment for cluster runs")
	("diffusion,d", po::value<int>(&diffMethod)->default_value(3),
	 "Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap")
	("timesteps,t", po::value<int>(&timeToSimulate), "Number of time steps to simulate\nTakes precedence over --days")
	("days", po::value<int>(&nDays)->default_value(200), "Number of days to simulate")
	("area-tnf-threshold", po::value<float>(&areaTNFThreshold)->default_value(0.5),"Threshold for granuloma area defined by TNF, in the range [0.0, 1.0]\n")
	("area-cell-density-threshold", po::value<float>(&areaCellDensityThreshold)->default_value(0.5),"Threshold for granuloma area defined by cell density, in the range [0.0, 1.0]\n")
	("ode", "Use integrated lymph node ODE for recruitment")
	("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
	("NFkB-dynamics", "Use molecular level intracellular NFkB dynamics in the model")
	("tnf-depletion", po::value<int>(&tnfDepletionTimeStep)->default_value(-1), "The time step at which to stop secreting tnf, including by tnfr dynamics. -1: no depletion")
	("ln-ode,l", po::value<std::string>(&lymphNodeODE), "Lymph node application")
	("ln-ode-temp", po::value<std::string>(&lymphNodeTemp), "Lymph node temp file")
	("lhs", "Running as part of an LHS run")
	("version,v", "Version number");
	
	try
	{
		po::positional_options_description p;
		p.add("input-file", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
		
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
        if (!vm.count("timesteps"))
            timeToSimulate = 144 * nDays;
		
		if (csvInterval < 0)
		{
			printUsage(argv[0], desc);
			return 1;
		}
		
		// If we just use time in seconds we will get duplicate seeds when doing a large number of runs,
		// such as running a big LHS on a cluster.
		// If we use just the micro-second part of the high resolution time we only get seeds that vary
		// between 1 and 1,000,000.
		// So we XOR the two times, which gives a larger range of seeds, with the lower order part varying
		// with a higher resolution than we would get just by using seconds.
		if (!vm.count("seed"))
		{
			seed = curTimeHiRes.tv_sec ^ curTimeHiRes.tv_usec;
		}
		
		
		// Adjust the seed if a seed adjustment was specified.
		if (vm.count("seedadj"))
		{
			seed = seed ^ seedadj;
		}
		
		ode = vm.count("ode");
		tnfrDynamics = vm.count("tnfr-dynamics");
		nfkbDynamics = vm.count("NFkB-dynamics");
		
		lhs = vm.count("lhs");
		
		if (lhs && (outputFileName.size() == 0))
		{
			std::cerr << "LHS run but no output name specified." << std::endl;
			exit(1);
		}
		
		RecruitmentProb recr;
		RecruitmentLnODEPure pureOdeRecr;
		RecruitmentBase* pRecr;
		if (ode) pRecr = &pureOdeRecr; else pRecr = &recr;
		
		DiffusionMethod diffMethodEnum;
		switch (diffMethod)
		{
			case 0:
				diffMethodEnum = DIFF_REC_EQ;
				break;
			case 1:
				diffMethodEnum = DIFF_SOR_CORRECT;
				std::cerr << "The BTCS diffusion method is no longer supported." << std::endl;
				exit(1);
				break;
			case 2:
				diffMethodEnum = DIFF_SOR_WRONG;
				std::cerr << "The BTCS Wrong diffusion method is no longer supported." << std::endl;
				exit(1);
				break;
			case 3:
				diffMethodEnum = DIFF_REC_EQ_SWAP;
				break;
			default:
				printUsage(argv[0], desc);
				return 1;
		}
		
	    screenDisplay = !vm.count("no-screen-display");
		
		// Write the seed to a file, so runs can be repeated, except for lhs runs.
		if (!lhs)
		{
			std::ofstream seedStream;
			seedStream.open("seed");
			seedStream << seed << std::endl;
		}
	    else {
			// For lhs runs, we don't want to print to the screen unless explicitly told otherwise
			screenDisplay = false;
	    }
		
		return run(seed, inputFileName, outputFileName, csvInterval, stateInterval, screenDisplay, diffMethodEnum, pRecr, ode,
				   tnfrDynamics, nfkbDynamics, tnfDepletionTimeStep, timeToSimulate, lhs, areaTNFThreshold, areaCellDensityThreshold );
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}
	
	return 0;
}
