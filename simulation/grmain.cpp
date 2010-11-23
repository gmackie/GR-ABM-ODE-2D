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

namespace po = boost::program_options;

void printVersion()
{
	std::cout << "Version: " << GR_VERSION << std::endl;
}

void printUsage(char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

int run(unsigned long seed, const std::string& inputFileName, DiffusionMethod diffMethod, int simulationDays, RecruitmentBase* pRecr, bool ode, bool tnfrDynamics)
{
	if (!Params::getInstance(true)->fromXml(inputFileName.c_str()))
		return 1;

	g_Rand.setSeed(seed);

	const int timeToSimulate = 144 * simulationDays;
	GrSimulation sim;

	sim.setTnfrDynamics(tnfrDynamics);
	sim.setRecruitment(pRecr);

	const GrStat& stats = sim.getStats();

	sim.setDiffusionMethod(diffMethod);
	sim.init();

	for (int time = 0; time <= timeToSimulate; time += 1)
	{
		printf("%d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\n",
			sim.getTime(), 
			stats.getNrOfMac(), stats.getNrOfMacResting(), stats.getNrOfMacInfected(),
			stats.getNrOfMacCInfected(), stats.getNrOfMacActive(), stats.getNrOfMacDead(),
			stats.getNrOfTgam(), stats.getNrOfTgamActive(),	stats.getNrOfTgamDownRegulated(),
			stats.getNrOfTgamDead(), stats.getNrOfTcyt(), stats.getNrOfTcytActive(),
			stats.getNrOfTcytDownRegulated(), stats.getNrOfTcytDead(),
			stats.getNrOfTreg(), stats.getNrOfTregActive(), stats.getNrOfTregDead(),
			stats.getTotExtMtb(), stats.getTotIntMtb(),
			stats.getTotTNF() / (NROWS*NCOLS), stats.getTotCCL2() / (NROWS*NCOLS), 
			stats.getTotCCL5() / (NROWS*NCOLS), stats.getTotCXCL9() / (NROWS*NCOLS));
		sim.solve();
	}

	return 0;
}

int main(int argc, char** argv)
{
	unsigned long seed;
	std::string inputFileName;
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
		("seed,s", po::value<unsigned long>(&seed)->default_value((unsigned long) curTime), "Seed")
		("diffusion,d", po::value<int>(&diffMethod)->default_value(0),
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

		switch (diffMethod)
		{
		case 0:
			return run(seed, inputFileName, DIFF_REC_EQ, nDays, pRecr, ode, tnfrDynamics);
		case 1:
			return run(seed, inputFileName, DIFF_SOR_CORRECT, nDays, pRecr, ode, tnfrDynamics);
		case 2:
			return run(seed, inputFileName, DIFF_SOR_WRONG, nDays, pRecr, ode, tnfrDynamics);
		case 3:
			return run(seed, inputFileName, DIFF_REC_EQ_SWAP, nDays, pRecr, ode, tnfrDynamics);
		default:
			printUsage(argv[0], desc);
			return 1;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
