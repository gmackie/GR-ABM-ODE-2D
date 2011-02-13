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
	outputFileStream
		<<inputFileName
		<< "," << "Mac"
		<< "," << "Mr"
		<< "," << "Mi"
		<< "," << "Mci"
		<< "," << "Ma"
		<< "," << "Mac dead"
		<< "," << "Tgam"
		<< "," << "Tgam a"
		<< "," << "Tgam reg"
		<< "," << "Tgam dead"
		<< "," << "Tcyt"
		<< "," << "Tcyt a"
		<< "," << "Tcyt reg"
		<< "," << "Tcyt dead"
		<< "," << "Treg"
		<< "," << "Treg a"
		<< "," << "Treg dead"
		<< "," << "Int. Mtb."
		<< "," << "Ext. Mtb."
		<< "," << "TNF"
		<< "," << "CCL2"
		<< "," << "CCL5"
		<< "," << "CXCL9"
		<< "," << "NrSourcesMac"
		<< "," << "NrSourcesTgam"
		<< "," << "NrSourcesTcyt"
		<< "," << "NrSourcesTreg"
		<< std::endl;
}

// The same output as written to the screen.
void writeOutput(std::ofstream& outputFileStream, GrSimulation& sim, int csvInterval)
{
	const GrStat& stats = sim.getStats();
		outputFileStream
			<< sim.getTime()/csvInterval
			<< "," << stats.getNrOfMac()
			<< "," << stats.getNrOfMacResting()
			<< "," << stats.getNrOfMacInfected()
			<< "," << stats.getNrOfMacCInfected()
			<< "," << stats.getNrOfMacActive()
			<< "," << stats.getNrOfMacDead()
			<< "," << stats.getNrOfTgam()
			<< "," << stats.getNrOfTgamActive()
			<< "," << stats.getNrOfTgamDownRegulated()
			<< "," << stats.getNrOfTgamDead()
			<< "," << stats.getNrOfTcyt()
			<< "," << stats.getNrOfTcytActive()
			<< "," << stats.getNrOfTcytDownRegulated()
			<< "," << stats.getNrOfTcytDead()
			<< "," << stats.getNrOfTreg()
			<< "," << stats.getNrOfTregActive()
			<< "," << stats.getNrOfTregDead()
			<< "," << stats.getTotExtMtb()
			<< "," << stats.getTotIntMtb()
			<< "," << stats.getTotTNF() / (NROWS*NCOLS)
			<< "," << stats.getTotCCL2() / (NROWS*NCOLS)
			<< "," << stats.getTotCCL5() / (NROWS*NCOLS)
			<< "," << stats.getTotCXCL9() / (NROWS*NCOLS)
			<< "," << stats.getNrSourcesMac()
			<< "," << stats.getNrSourcesTgam()
			<< "," << stats.getNrSourcesTcyt()
			<< "," << stats.getNrSourcesTreg()
			<< std::endl;
		;
}

int run(unsigned long seed, const std::string& inputFileName, const std::string& outputFileName, int csvInterval, DiffusionMethod diffMethod, int simulationDays, RecruitmentBase* pRecr, bool ode, bool tnfrDynamics)
{
	if (!Params::getInstance(true)->fromXml(inputFileName.c_str()))
		return 1;

	g_Rand.setSeed(seed);

	const int timeToSimulate = 144 * simulationDays;

	// If an output file is requested, construct the complete output file name,
	// open the output file and write a header line.
	std::string fullOutputFileName = "";
	std::ofstream outputFileStream;

	if (outputFileName.size() > 0)
	{
		std::ostringstream ofn;
		ofn << simulationDays << "-" << outputFileName << "-" << seed << ".csv";
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
		if (time % csvInterval == 0 || time == timeToSimulate)
		{
//			printf("%d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\n",
			printf("%d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\t(%d,%d,%d,%d)\n",
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
				stats.getNrSourcesMac(), stats.getNrSourcesTgam(), stats.getNrSourcesTcyt(), stats.getNrSourcesTreg()
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
	int nDays;
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
		("seed,s", po::value<unsigned long>(&seed)->default_value((unsigned long) curTime), "Seed")
		("diffusion,d", po::value<int>(&diffMethod)->default_value(3),
				"Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap")
		("days", po::value<int>(&nDays)->default_value(200), "Number of days to simulate")
		("ode", "Use integrated lymph node ODE for recruitment")
		("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
		("ln-ode,l", po::value<std::string>(&lymphNodeODE), "Lymph node application")
		("ln-ode-temp,t", po::value<std::string>(&lymphNodeTemp), "Lymph node temp file")
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

		return run(seed, inputFileName, outputFileName, csvInterval, diffMethodEnum, nDays, pRecr, ode, tnfrDynamics);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
