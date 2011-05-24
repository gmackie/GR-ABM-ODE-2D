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
#include <time.h>
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
                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
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
			<< "\"TNF\""
			<< ','
			<< "\"CCL2\""
			<< ','
			<< "\"CCL5\""
			<< ','
			<< "\"CXCL9\""
			<< ','
			<< "\"Area\""
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
			<< ","
			<< "\"NonRepl Ext. Mtb.\""
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
			<< stats.getTotNonRepExtMtb()
			;

		outputFileStream << std::endl;
}
void saveState(const GrSimulation* pSim, int time) {
    int days, hours, minutes;
    char fname[14];
    assert(pSim);
    GrSimulation::convertSimTime(time, days, hours, minutes);
    sprintf(fname, "%03d%02d%02d.state", days, hours, minutes);
    std::ofstream out(fname, std::ios_base::trunc);
    if(!out.good())
        std::cerr<<"Unable to open output file: "<<fname<<std::endl;
    else pSim->serialize(out);
    out.close();    //Don't want to flush here, would slow things down a lot
}


int run(unsigned long seed, const std::string& inputFileName, const std::string& outputFileName, int csvInterval, int stateInterval, DiffusionMethod diffMethod, RecruitmentBase* pRecr, bool ode, bool tnfrDynamics, int timeToSimulate)
{
	std::cout << endl << "--seed " << seed << std::endl;

	if (!Params::getInstance(true)->fromXml(inputFileName.c_str()))
		return 1;

	g_Rand.setSeed(seed);

	// If an output file is requested, construct the complete output file name,
	// open the output file and write a header line.
	std::string fullOutputFileName = "";
	std::ofstream outputFileStream;

	if (outputFileName.size() > 0)
	{
		std::ostringstream ofn;
//DBG	ofn << timeToSimulate / 144 << "-" << outputFileName << "-" << seed << ".csv";
		ofn << outputFileName << "-" << seed << ".csv";
		fullOutputFileName = ofn.str();
	    outputFileStream.open(fullOutputFileName.c_str());

		writeOutputHeader(outputFileStream, inputFileName);
	}

	GrSimulation sim;

	sim.setTnfrDynamics(tnfrDynamics);
	sim.setRecruitment(pRecr);

	const GrStat& stats = sim.getStats();

	sim.setDiffusionMethod(diffMethod);
	sim.init();

	// Write output at time 0, the initial simulation state.
	if (fullOutputFileName.size() > 0)
	{
		writeOutput(outputFileStream, sim, csvInterval);
	}

	for (int time = 1; time <= timeToSimulate; time += 1)
	{
		sim.solve();

		// Display and write output at the requested interval, and after the last time step.
        if (stateInterval > 0 && time % stateInterval == 0)
            saveState(&sim, time);
		if (time % csvInterval == 0 || time == timeToSimulate)
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

			if (fullOutputFileName.size() > 0)
			{
				writeOutput(outputFileStream, sim, csvInterval);
			}
		}
	}

	return 0;
}

int main(int argc, char** argv)
{
	unsigned long seed;
	std::string inputFileName;
	std::string outputFileName;
	int csvInterval; // In time steps
	std::string lymphNodeODE;
	std::string lymphNodeTemp;
	int diffMethod;
    int stateInterval;
	int nDays;
    int timeToSimulate;
	bool ode;
	bool tnfrDynamics;

	/* set seed to current time, in case not specified */
	time_t curTime;
	time(&curTime);

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
		("output-file,o", po::value<std::string>(&outputFileName), "Output file name")
		("csv-interval", po::value<int>(&csvInterval)->default_value(1),
						"CSV update interval (10 min timesteps)")
		("state-interval", po::value<int>(&stateInterval)->default_value(-1),
						"State save interval (10 min timesteps)")
		("seed,s", po::value<unsigned long>(&seed)->default_value((unsigned long) curTime), "Seed")
		("diffusion,d", po::value<int>(&diffMethod)->default_value(3),
				"Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap")
		("days", po::value<int>(&nDays)->default_value(200), "Number of days to simulate")
		("ode", "Use integrated lymph node ODE for recruitment")
		("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
		("ln-ode,l", po::value<std::string>(&lymphNodeODE), "Lymph node application")
		("ln-ode-temp", po::value<std::string>(&lymphNodeTemp), "Lymph node temp file")
		("timesteps,t", po::value<int>(&timeToSimulate), "Number of time steps to simulate\nTakes precedence over --days")
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
		
		ode = vm.count("ode");
		tnfrDynamics = vm.count("tnfr-dynamics");

		//Recruitment recrLN(lymphNodeODE, lymphNodeTemp);
		//Recruitment
		//Recruitment* pRecr = (vm.count("ln-ode") && vm.count("ln-ode-temp")) ? &recrLN : NULL;

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
			break;
		case 2:
			diffMethodEnum = DIFF_SOR_WRONG;
			break;
		case 3:
			diffMethodEnum = DIFF_REC_EQ_SWAP;
			break;
		default:
			printUsage(argv[0], desc);
			return 1;
		}

		return run(seed, inputFileName, outputFileName, csvInterval, stateInterval, diffMethodEnum, pRecr, ode, tnfrDynamics, timeToSimulate);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
