/*
 * grvascular.cpp
 *
 * Created by ncilfone
 * 09/25/2012
 */


#include "grvascular.h"
#include "grsimulation.h"
#include "gr.h"

Vascular::Vascular()
    : //_bloodConcentrationIL10(0.0), // blood concentration of IL10 in mol/L
      //_bloodConcentrationTNF(0.0), // blood concentration of TNF in mol/L
      _bloodConcentrationINH(0.0), //blood concentration of INH in mol/L
      _bloodConcentrationRIF(0.0)       //blood concentration of RIF in mol/L
{ 
}

Vascular::~Vascular()
{
}

void Vascular::solveVascularSources(GrGrid& grid, double dt, const int time, int DiffStep)
{
    int realTime = static_cast<int>(time*600 + DiffStep*_PARAM(_timestepDiffusion));

    if (isDoseStartTime(realTime))
    {
        //addDose(realTime, _PARAM(PARAM_GR_VASCULAR_IL10_DOSE), _PARAM(PARAM_GR_VASCULAR_IL10_DOSE_INTERVAL), _bloodConcentrationIL10);
        //addDose(realTime, _PARAM(PARAM_GR_VASCULAR_TNF_DOSE), _PARAM(PARAM_GR_VASCULAR_TNF_DOSE_INTERVAL), _bloodConcentrationTNF);
        addDose(realTime, _PARAM(_bloodINHdose), _PARAM(_dosageIntervalINH), _bloodConcentrationINH);
        addDose(realTime, _PARAM(_bloodRIFdose), _PARAM(_dosageIntervalRIF), _bloodConcentrationRIF);

        const std::vector<Pos>& sources = grid.getSources();

        for (PosVector::const_iterator it = sources.begin(); it != sources.end(); it++)
        {
            // TNF
            /*double localTNF = getAverageOfNeighbours(grid.TNF(), *it, grid.getRange());
            double tnfChange = calculateFluxChange(localTNF, _bloodConcentrationTNF, dt);
            grid.incTNF(*it, tnfChange);
            modulateBlood(_bloodConcentrationTNF, tnfChange);*/


            // IL10
            /*double localIL10 = getAverageOfNeighbours(grid.il10(), *it, grid.getRange());
            double il10Change = calculateFluxChange(localIL10, _bloodConcentrationIL10, dt);
            grid.incil10(*it, il10Change);
            modulateBlood(_bloodConcentrationIL10, il10Change);*/


            // INH
            double oldINH = grid.INH(*it);
            double localINH = getAverageOfNeighbours(grid.INH(), *it, grid.getRange(), MW_INH);
            double INHChange = calculateFluxChange(localINH, _bloodConcentrationINH, dt, _PARAM(_vascularPermeabilityINH), MW_INH);
            grid.incINH(*it, INHChange);
            double newINH = grid.INH(*it);
            modulateBlood(_bloodConcentrationINH, INHChange, MW_INH);

            // RIF
            double localRIF = getAverageOfNeighbours(grid.RIF(), *it, grid.getRange(), MW_RIF);
            double RIFChange = calculateFluxChange(localRIF, _bloodConcentrationRIF, dt, _PARAM(_vascularPermeabilityRIF), MW_RIF);
            grid.incRIF(*it, RIFChange);
            modulateBlood(_bloodConcentrationRIF, RIFChange, MW_RIF);


        }
        //updateBloodODE(_bloodConcentrationIL10, _PARAM(PARAM_GR_VASCULAR_IL10_CLEARANCE), dt);
        //updateBloodODE(_bloodConcentrationTNF, _PARAM(PARAM_GR_VASCULAR_TNF_CLEARANCE), dt);
        updateBloodODE(_bloodConcentrationINH, _PARAM(_degRateINH_Plasma_Fast_Acet), dt);
        updateBloodODE(_bloodConcentrationRIF, _PARAM(_degRateRIF_Plasma), dt);
    }
}


double Vascular::getAverageOfNeighbours(Scalar* grid, const Pos& pos, const Pos& dim, const Scalar MW)
{
    double sumC = 0.0;

    Scalar up, dn, lt, rt, ct;

    up = grid[Indexer::padInd(dim,((pos.x)-1),pos.y)];
    lt = grid[Indexer::padInd(dim,((pos.x)),((pos.y)-1))];
    dn = grid[Indexer::padInd(dim,((pos.x)+1),pos.y)];
    rt = grid[Indexer::padInd(dim,((pos.x)),((pos.y)+1))];
    ct = grid[Indexer::padInd(dim,((pos.x)),((pos.y)))];

//    std::cout << "Up: " << up << "  Lt: " << lt << "  Dn: " << dn << "  Rt: " << rt << std::endl;

    sumC = (up + dn + lt + rt + ct)/(5.0)/MW;

    return (sumC); // Returns average concentration mol/L
}

void Vascular::updateBloodODE(double& bloodConcentration, double clearance, double dt)
{
    // Analytical solution to the one compartment blood ODE with first order clearance
    Scalar degRate = exp((-1.0 *  clearance * dt));
    bloodConcentration = bloodConcentration * degRate;
    //bloodConcentration = bloodConcentration * clearance;

}
