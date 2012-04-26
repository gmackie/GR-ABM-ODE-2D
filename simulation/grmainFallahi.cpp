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
#include <time.h>
#include <boost/program_options.hpp>
#include "recruitmentlnode.h"
#include "recruitmentlnodepure.h"
#include "recruitmentprob.h"
#include <iostream>
#include <fstream>

using namespace std;

namespace po = boost::program_options;

void printVersion()
{
	std::cout << "Version: " << GR_VERSION << std::endl;
}

void printUsage(char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

bool stopCondition(int timeToSimulate, int time, Stats stats)
{
	return (stats.getTotExtMtb() == 0 &&
			stats.getTotIntMtb() == 0 &&
			stats.getTotTNF() < DBL_EPSILON * 10.0);
}

int run(unsigned long seed, const std::string& inputFileName, DiffusionMethod diffMethod, int simulationDays, RecruitmentBase* pRecr, bool ode, bool tnfrDynamics, bool useStopCondition)
{
	if (!Params::getInstance(true)->fromXml(inputFileName.c_str()))
		return 1;

	g_Rand.setSeed(seed);

	const int timeToSimulate = 144 * simulationDays;

	GrSimulation sim;

	sim.setTnfrDynamics(tnfrDynamics);
	sim.setRecruitment(pRecr);

	const Stats& stats = sim.getStats();

	sim.setDiffusionMethod(diffMethod);
	sim.init();

    std::ofstream ExcelOutputFile7, ExcelOutputFile20, ExcelOutputFile50, ExcelOutputFile100, ExcelOutputFile200, ExcelOutputFileAllDays, ExcelMolInputFile;
    std::ofstream ExcelOutputFileAllTimeSteps;

//    ExcelOutputFileAllTimeSteps.open("144-outputAllTimeSteps.csv");
    ExcelOutputFileAllDays.open("200-outputAllDays.csv");

//        ExcelOutputFile7.open("/home/simeone/immunology/Test/output7.csv",std::ios::app);
//        ExcelOutputFile20.open("/home/simeone/immunology/Test/output20.csv",std::ios::app);
//        ExcelOutputFile50.open("/home/simeone/immunology/Test/output50.csv",std::ios::app);
//        ExcelOutputFile100.open("/home/simeone/immunology/Test/output100.csv",std::ios::app);
//        ExcelOutputFile200.open("/home/simeone/immunology/Test/output200.csv",std::ios::app);
//        ExcelOutputFileAllDays.open("10-outputAllTimeSteps.csv",std::ios::app);
//        ExcelMolInputFile.open("/home/simeone/immunology/Test/input.csv",std::ios::app);

  /*    ExcelOutputFile7.open("/home/simeone/immunology/MDC_test/output7.csv",std::ios::app);
        ExcelOutputFile20.open("/home/simeone/immunology/MDC_test/output20.csv",std::ios::app);
        ExcelOutputFile50.open("/home/simeone/immunology/MDC_test/output50.csv",std::ios::app);
        ExcelOutputFile100.open("/home/simeone/immunology/MDC_test/output100.csv",std::ios::app);
        ExcelOutputFile200.open("/home/simeone/immunology/MDC_test/output200.csv",std::ios::app);
        ExcelOutputFileAllDays.open("/home/simeone/immunology/MDC_test/outputAllDays.csv",std::ios::app);
        ExcelMolInputFile.open("/home/simeone/immunology/MDC_test/input.csv",std::ios::app);

        ExcelOutputFile7.open("/home/simeone/immunology/Dissemination_test/output7.csv",std::ios::app);
        ExcelOutputFile20.open("/home/simeone/immunology/Dissemination_test/output20.csv",std::ios::app);
        ExcelOutputFile50.open("/home/simeone/immunology/Dissemination_test/output50.csv",std::ios::app);
        ExcelOutputFile100.open("/home/simeone/immunology/Dissemination_test/output100.csv",std::ios::app);
        ExcelOutputFile200.open("/home/simeone/immunology/Dissemination_test/output200.csv",std::ios::app);
        ExcelOutputFileAllDays.open("/home/simeone/immunology/Dissemination_test/outputAllDays.csv",std::ios::app);
        ExcelMolInputFile.open("/home/simeone/immunology/Dissemination_test/input.csv",std::ios::app);

        ExcelOutputFile7.open("/home/simeone/immunology/Containment_test/output7.csv",std::ios::app);
        ExcelOutputFile20.open("/home/simeone/immunology/Containment_test/output20.csv",std::ios::app);
        ExcelOutputFile50.open("/home/simeone/immunology/Containment_test/output50.csv",std::ios::app);
        ExcelOutputFile100.open("/home/simeone/immunology/Containment_test/output100.csv",std::ios::app);
        ExcelOutputFile200.open("/home/simeone/immunology/Containment_test/output200.csv",std::ios::app);
        ExcelOutputFileAllDays.open("/home/simeone/immunology/Containment_test/outputAllDays_1.csv",std::ios::app);
        ExcelMolInputFile.open("/home/simeone/immunology/Containment_test/input.csv",std::ios::app);

        ExcelOutputFile7.open("/home/simeone/immunology/Test/output7_1.csv",std::ios::app);
        ExcelOutputFile20.open("/home/simeone/immunology/Test/output20_1.csv",std::ios::app);
        ExcelOutputFile50.open("/home/simeone/immunology/Test/output50_1.csv",std::ios::app);
        ExcelOutputFile100.open("/home/simeone/immunology/Test/output100_1.csv",std::ios::app);
        ExcelOutputFile200.open("/home/simeone/immunology/Test/output200_1.csv",std::ios::app);
        ExcelOutputFileAllDays.open("/home/simeone/immunology/Test/outputAllDays_1.csv",std::ios::app);
        ExcelMolInputFile.open("/home/simeone/immunology/Test/input.csv",std::ios::app);

                ExcelOutputFile7.open("/home/simeone/immunology/Test/output7.csv",std::ios::app);
                ExcelOutputFile20.open("/home/simeone/immunology/Test/output20.csv",std::ios::app);
                ExcelOutputFile50.open("/home/simeone/immunology/Test/output50.csv",std::ios::app);
                ExcelOutputFile100.open("/home/simeone/immunology/Test/output100.csv",std::ios::app);
                ExcelOutputFile200.open("/home/simeone/immunology/Test/output200.csv",std::ios::app);
                ExcelOutputFileAllDays.open("/home/simeone/immunology/Test/outputAllDays.csv",std::ios::app);
                ExcelMolInputFile.open("/home/simeone/immunology/Test/input.csv",std::ios::app);



                ExcelOutputFile7.open("/home/simeone/immunology/Test/output7_3.csv",std::ios::app);
                ExcelOutputFile20.open("/home/simeone/immunology/Test/output20_3.csv",std::ios::app);
                ExcelOutputFile50.open("/home/simeone/immunology/Test/output50_3.csv",std::ios::app);
                ExcelOutputFile100.open("/home/simeone/immunology/Test/output100_3.csv",std::ios::app);
                ExcelOutputFile200.open("/home/simeone/immunology/Test/output200_3.csv",std::ios::app);
                ExcelOutputFileAllDays.open("/home/simeone/immunology/Test/outputAllDays_3.csv",std::ios::app);
                ExcelMolInputFile.open("/home/simeone/immunology/Test/input.csv",std::ios::app);

                ExcelOutputFile7.open("/home/simeone/immunology/Test/output7_2.csv",std::ios::app);
                ExcelOutputFile20.open("/home/simeone/immunology/Test/output20_2.csv",std::ios::app);
                ExcelOutputFile50.open("/home/simeone/immunology/Test/output50_2.csv",std::ios::app);
                ExcelOutputFile100.open("/home/simeone/immunology/Test/output100_2.csv",std::ios::app);
                ExcelOutputFile200.open("/home/simeone/immunology/Test/output200_2.csv",std::ios::app);
                ExcelOutputFileAllDays.open("/home/simeone/immunology/Test/outputAllDays_2.csv",std::ios::app);
                ExcelMolInputFile.open("/home/simeone/immunology/Test/input.csv",std::ios::app);
*/
                /* ExcelOutputFile20.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/output20.csv",std::ios::app);
                ExcelOutputFile40.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/output40.csv",std::ios::app);
                ExcelOutputFile100.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/output100.csv",std::ios::app);
                ExcelOutputFile200.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/output200.csv",std::ios::app);
                ExcelOutputFileAllDays.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/outputAllDays.csv",std::ios::app);
                ExcelMolInputFile.open("/Users/lindelab/Desktop/QT-GR ABM/LHS/results1/input.csv",std::ios::app);*/

	ExcelOutputFileAllDays <<inputFileName<<","<<"Mac" <<","<<"Mr"<<","<<"Mi"<<","<<"Mci"<<","<<"Ma"<<","<<"Tgam"<<
        ","<<"Tgam a"<<","<<"Tgam reg"<<","<<"Tcyt"<<","<<"Tcyt a"<<","<<"Tcyt reg"<<","<<"Treg"<<
        ","<<"Int. Mtb."<<","<<"Ext. Mtb."<<","<<"TNF"<<","<<"CCL2"<<","<<"CCL5"<<
        ","<<"CXCL9"<<","<<"Mac Apop"<<","<< "MDC"<< ","<< "N4"<< ","
        << "TH0"<< ","<< "TH1"<< ","<< "N8"<< ","<< "T80"<< ","<< "T8"<< ","
        << "TC"<< "," << "TH0lung" << ","<< "TH1lung"<< ","<< "T80lung"<< ","
        << "T8lung" << "," << "TClung"<< std::endl;
        /*<<"Mr Apop"<<","<<"Mi&Mci Apop"<<","<<"Ma Apop"<<","<<"T Apop"<<","<<"Mr NFkB"<<
	","<<"Mi NFkB"<<","<<"R1 Mr"<<","<<"R1 Mi&Mci"<<","<<"R1 Ma"<<","<<"R1 Tcell"<<","<<"i/R1 Mr"<<","<<"i/R1 Mi&Mci"<<
	","<<"i/R1 Ma"<<","<<"i/R1 Tcell"<< std::endl;
*/
	/*
	 ExcelMolInputFile <<"INPUT"<<","<<"muMDC_LN"<<","<<"sn4"<<","<<"muN4"<<","<<"k13"<<","<<"hs13"<<","<<"k14"<<","<<"k15"<<","<<"rho2"<<
	 ","<<"k20a"<<","<<"hs20a"<<","<<"csi1"<<","<<"csi1a"<<","<<"sn8"<<","<<"muN8"<<","<<"wT80"<<","<<"k16"<<","<<"hs16"<<
         ","<<"k17"<<","<<"hs17"<<","<<"k18"<<","<<"rho3"<<","<<"k24a"<<","<<"hs24a"<<","<<"csi2"<<","<<"csi2a"<<","<<"csi2b"<<","<<"scaling"<<
         ","<<"m"<<","<<"scalingMDC"<<","<<"D_Chemokine"<< ","<<"degChemokine" <<","<<"minChemotaxis"<<","<<"maxChemotaxis"<<","<<"pMr-killMtb"<<","<<"NFkB-threshMtb"<<
	 ","<<"MtbThreshBurst"<<","<<"pSTAT1Tgam"<<","<<"MacThreshRec"<<","<<"pMacRec"<<","<<"TmoveMac"<<","<<"TmoveT"<<","<<"pTcellRec"<<","<<"pApopFas"<<","<<"TgamThreshRec"<<
	 ","<<"pKillMacTcyt"<<","<<"TcytThreshRec"<<","<<"TregThreshRec"<<","<<"intMtbGrowth"<<","<<"extMtbGrowth"<< std::endl;
	 */
//    ExcelMolInputFile <<inputFileName<<","<<_PARAM(PARAM_muMDC_LN)<<","<<_PARAM(PARAM_initN4)<<","<<_PARAM(PARAM_muN4)<<
//	","<<_PARAM(PARAM_k13)<<","<<_PARAM(PARAM_hs13)<<","<<_PARAM(PARAM_k14)<<","<<_PARAM(PARAM_k15)<<","<<_PARAM(PARAM_rho2)<<
//	","<<_PARAM(PARAM_k20a)<<","<<_PARAM(PARAM_hs20a)<<","<<_PARAM(PARAM_csi1)<<","<<_PARAM(PARAM_csi1a)<<
//	","<<_PARAM(PARAM_initN8)<<","<<_PARAM(PARAM_muN8)<<","<<_PARAM(PARAM_wT80)<<","<<_PARAM(PARAM_k16)<<","<<_PARAM(PARAM_hs16)<<
//	","<<_PARAM(PARAM_k17)<<","<<_PARAM(PARAM_hs17)<<","<<_PARAM(PARAM_k18)<<","<<_PARAM(PARAM_rho3)<<","<<_PARAM(PARAM_k24a)<<
//	","<<_PARAM(PARAM_hs24a)<<","<<_PARAM(PARAM_csi2)<<","<<_PARAM(PARAM_csi2a)<<","<<_PARAM(PARAM_csi2b)<<","<<_PARAM(PARAM_scaling)<<
//    ","<<_PARAM(PARAM_m)<<","<<_PARAM(PARAM_scaling_MDC)<<","<<_PARAM(PARAM_scaling_LUNG)<<","<<_PARAM(PARAM_scaling_LN)<<
//    ","<<_PARAM(PARAM_GR_D_CHEMOKINES)<<","<<_PARAM(PARAM_GR_DEG_CHEMOKINES)<<","<<_PARAM(PARAM_GR_MIN_CHEMOTAXIS)<<
//	","<<_PARAM(PARAM_GR_MAX_CHEMOTAXIS)<<","<<_PARAM(PARAM_MAC_PROB_KILL_R_EXTMTB)<<","<<_PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB)<<
//	","<<_PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB)<<","<<_PARAM(PARAM_MAC_PROB_STAT1_TGAM)<<","<<_PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT)<<
//	","<<_PARAM(PARAM_MAC_PROB_RECRUITMENT)<<","<<_PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC)<<","<<_PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL)<<
//	","<<_PARAM(PARAM_TCELL_PROB_RECRUITMENT)<<","<<_PARAM(PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL)<<","<<_PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT)<<
//	","<<_PARAM(PARAM_TCYT_PROB_KILL_MAC)<<","<<_PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT)<<","<<_PARAM(PARAM_TREG_THRESHOLD_RECRUITMENT)<<
//	","<<_PARAM(PARAM_INTMTB_GROWTH_RATE)<<","<<_PARAM(PARAM_EXTMTB_GROWTH_RATE)<<","<< std::endl;

	for (int time = 0; time <= timeToSimulate; time += 1)
	{
		if (time % 144 == 0)
//		if (time % 1 == 0)
		{
			printf("%lu %d\t %d - (%d,%d,%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d,%d)\t%d - (%d,%d)\t(%f,%f)\t(%f,%f,%f,%f)\n - ",
				   seed,
	  			   sim.getTime()/144, 
				   stats.getNrOfMac(), stats.getNrOfMacResting(), stats.getNrOfMacInfected(),
				   stats.getNrOfMacCInfected(), stats.getNrOfMacActive(), stats.getNrOfMacDead(),
				   stats.getNrOfTgam(), stats.getNrOfTgamActive(),	stats.getNrOfTgamDownRegulated(),
				   stats.getNrOfTgamDead(), stats.getNrOfTcyt(), stats.getNrOfTcytActive(),
				   stats.getNrOfTcytDownRegulated(), stats.getNrOfTcytDead(),
				   stats.getNrOfTreg(), stats.getNrOfTregActive(), stats.getNrOfTregDead(),
				   stats.getTotExtMtb(), stats.getTotIntMtb(),
				   stats.getTotTNF() / (NROWS*NCOLS), stats.getTotCCL2() / (NROWS*NCOLS), 
				   stats.getTotCCL5() / (NROWS*NCOLS), stats.getTotCXCL9() / (NROWS*NCOLS)),
				   stats.getTH0(), stats.getTH1(),stats.getN8(), stats.getT80(),stats.getT8(),
				   stats.getTC(),stats.getTH0lung(), stats.getTH1lung(),stats.getT80lung(),
				   stats.getT8lung(), stats.getTClung();

//			ExcelOutputFileAllTimeSteps <<sim.getTime()<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//			","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//			","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                        ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//			","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                        ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                        ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                        ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                        ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                        ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;


			ExcelOutputFileAllDays <<sim.getTime()/144<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
			","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
			","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
                        ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
			","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
                        ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
                        ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
                        ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
                        ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
                        ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
                        /*<<","<<stats.getNrRestingMacApoptosisTNF()<<
			","<<stats.getNrInfAndCinfMacApoptosisTNF()<<","<<stats.getNrActivatedMacApoptosisTNF()<<","<<stats.getNrTcellApoptosisTNF()<<
			","<<stats.getNrRestingMacNFkBActivationTNF()<<","<<stats.getNrInfMacNFkBActivationTNF()<<
			","<<stats.getTotSurfBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
			","<<stats.getTotSurfBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
			","<<stats.getTotSurfBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
			","<<stats.getTotSurfBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<
			","<<stats.getTotIntBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
			","<<stats.getTotIntBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
			","<<stats.getTotIntBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
                        ","<<stats.getTotIntBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<std::endl;*/
//                        if (time == 1008)
//                                ExcelOutputFile7 <<inputFileName<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//                                ","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//                                ","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                                ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//                                ","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                                ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                                ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                                ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                                ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
//                        if (time == 2880)
//				ExcelOutputFile20 <<inputFileName<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//				","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//				","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                                ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//				","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                                ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                                ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                                ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                                ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
                        /*<<","<<stats.getNrRestingMacApoptosisTNF()<<
				","<<stats.getNrInfAndCinfMacApoptosisTNF()<<","<<stats.getNrActivatedMacApoptosisTNF()<<","<<stats.getNrTcellApoptosisTNF()<<
				","<<stats.getNrRestingMacNFkBActivationTNF()<<","<<stats.getNrInfMacNFkBActivationTNF()<<
				","<<stats.getTotSurfBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotSurfBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotSurfBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
				","<<stats.getTotSurfBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<
				","<<stats.getTotIntBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotIntBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotIntBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
                                ","<<stats.getTotIntBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<< std::endl;*/
//			if (time == 5760)
//                                ExcelOutputFile50 <<inputFileName<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//				","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//				","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                                ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//				","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                                ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                                ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                                ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                                ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
                        /*<<","<<stats.getNrRestingMacApoptosisTNF()<<
				","<<stats.getNrInfAndCinfMacApoptosisTNF()<<","<<stats.getNrActivatedMacApoptosisTNF()<<","<<stats.getNrTcellApoptosisTNF()<<
				","<<stats.getNrRestingMacNFkBActivationTNF()<<","<<stats.getNrInfMacNFkBActivationTNF()<<
				","<<stats.getTotSurfBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotSurfBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotSurfBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
				","<<stats.getTotSurfBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<
				","<<stats.getTotIntBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotIntBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotIntBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
                                ","<<stats.getTotIntBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<std::endl;*/
//			if (time == 14400)
//				ExcelOutputFile100 <<inputFileName<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//				","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//				","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                                ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//				","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                                ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                                ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                                ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                                ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
                        /*<<","<<stats.getNrRestingMacApoptosisTNF()<<
				","<<stats.getNrInfAndCinfMacApoptosisTNF()<<","<<stats.getNrActivatedMacApoptosisTNF()<<","<<stats.getNrTcellApoptosisTNF()<<
				","<<stats.getNrRestingMacNFkBActivationTNF()<<","<<stats.getNrInfMacNFkBActivationTNF()<<
				","<<stats.getTotSurfBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotSurfBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotSurfBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
				","<<stats.getTotSurfBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<
				","<<stats.getTotIntBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotIntBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotIntBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
                                ","<<stats.getTotIntBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<std::endl;*/
//			if (time == 28800)
//				ExcelOutputFile200 <<inputFileName<<","<<stats.getNrOfMac()<<","<<stats.getNrOfMacResting()<<","<<stats.getNrOfMacInfected()<<
//				","<<stats.getNrOfMacCInfected()<<","<<stats.getNrOfMacActive()<<","<<stats.getNrOfTgam()<<","<<stats.getNrOfTgamActive()<<
//				","<<stats.getNrOfTgamDownRegulated()<<","<<stats.getNrOfTcyt()<<","<<stats.getNrOfTcytActive()<<","<<stats.getNrOfTcytDownRegulated()<<
//                                ","<<stats.getNrOfTregActive()<<","<<stats.getTotIntMtb()<<","<<stats.getTotExtMtb()<<
//				","<<stats.getTotTNF() / (NROWS*NCOLS)<<","<<stats.getTotCCL2() / (NROWS*NCOLS)<<","<<stats.getTotCCL5() / (NROWS*NCOLS)<<
//                                ","<<stats.getTotCXCL9() / (NROWS*NCOLS)<<","<<stats.getNrApoptosisTNF() <<
//                                ","<<stats.getMDC()<< ","<< stats.getN4()<<","<< stats.getTH0()<<","<< stats.getTH1()<<
//                                ","<< stats.getN8()<< ","<< stats.getT80()<< ","<< stats.getT8()<< ","<< stats.getTC()<<
//                                ","<<stats.getTH0lung()<< ","<< stats.getTH1lung()<< ","<< stats.getT80lung()<<
//                                ","<< stats.getT8lung()<< ","<< stats.getTClung()<<std::endl;
                        /*<<","<<stats.getNrRestingMacApoptosisTNF()<<
				","<<stats.getNrInfAndCinfMacApoptosisTNF()<<","<<stats.getNrActivatedMacApoptosisTNF()<<","<<stats.getNrTcellApoptosisTNF()<<
				","<<stats.getNrRestingMacNFkBActivationTNF()<<","<<stats.getNrInfMacNFkBActivationTNF()<<
				","<<stats.getTotSurfBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotSurfBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotSurfBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
				","<<stats.getTotSurfBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<
				","<<stats.getTotIntBoundTNFR1RestingMac() / stats.getNrOfMacResting()<<
				","<<stats.getTotIntBoundTNFR1InfAndCinfMac() / (stats.getNrOfMacInfected() + stats.getNrOfMacCInfected())<<
				","<<stats.getTotIntBoundTNFR1ActiveMac() / stats.getNrOfMacActive()<<
                                ","<<stats.getTotIntBoundTNFR1Tcell() / (stats.getNrOfTgam() + stats.getNrOfTcyt() + stats.getNrOfTreg())<<std::endl;*/
		}
		sim.solve();
		if (useStopCondition && stopCondition(timeToSimulate, time, stats)) {
			break;
		}
	}
//	ExcelOutputFile7.close();
//	ExcelOutputFile20.close();
//	ExcelOutputFile50.close();
//	ExcelOutputFile100.close();
//	ExcelOutputFile200.close();
//	ExcelMolInputFile.close();

	ExcelOutputFileAllDays.close();


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
	bool useStopCondition = false;

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
		("use-stop-condition,c", "Stop simulation when stop condition occurs (usually clearance)")
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
		useStopCondition = vm.count("use-stop-condition");

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
			return run(seed, inputFileName, DIFF_REC_EQ, nDays, pRecr, ode, tnfrDynamics, useStopCondition);
		case 1:
			return run(seed, inputFileName, DIFF_SOR_CORRECT, nDays, pRecr, ode, tnfrDynamics, useStopCondition);
		case 2:
			return run(seed, inputFileName, DIFF_SOR_WRONG, nDays, pRecr, ode, tnfrDynamics, useStopCondition);
		case 3:
			return run(seed, inputFileName, DIFF_REC_EQ_SWAP, nDays, pRecr, ode, tnfrDynamics, useStopCondition);
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

