/*
 * GrDiffusionFTCS_Swap.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#include "grdiffusionftcsswap.h"

GrDiffusionFTCS_Swap::GrDiffusionFTCS_Swap()
: _cutOffValue(0.000001)
{
}

GrDiffusionFTCS_Swap::~GrDiffusionFTCS_Swap()
{
}

void GrDiffusionFTCS_Swap::diffuse(GrSimulationGrid& grSim) const
{
	GrGrid& currentGrid = grSim.getCurrentGrid();
	GrGrid& nextGrid = grSim.getNextGrid();

	const double muTNF = _PARAM(PARAM_GR_D_TNF) * 6 / (4e-6);
	const double muShedTNFR2 = _PARAM(PARAM_GR_D_SHED_TNFR2) * 6 / (4e-6);
	const double muChemokines = _PARAM(PARAM_GR_D_CHEMOKINES) * 6 / (4e-6);
    const double muIL10 = _PARAM(PARAM_GR_D_IL10) * 6 / (4e-6);
	const double degTNF = _PARAM(PARAM_GR_DEG_TNF);
	const double degChemokines = _PARAM(PARAM_GR_DEG_CHEMOKINES);
    const double degIL10 = _PARAM(PARAM_GR_DEG_IL10);
	const double ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
	const double ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);


	const double dt = 6; // time-step (sec)

	// TNF and CCs are zero
	GridCell zeroCell;

	// dt = 6s, solve for 1 timestep, 6 seconds
	for (int t = 0; t < 1; t++)
	{
		for (int i = 0; i < NROWS; i++)
		{
			for (int j = 0; j < NCOLS; j++)
			{
				GridCell& newCell = nextGrid(i, j);

				GridCell& cell = currentGrid(i, j);
				GridCell& cell_i_min_1_j = (i > 0) ? currentGrid(i - 1, j) : zeroCell;
				GridCell& cell_i_plus_1_j =  (i < NROWS - 1) ? currentGrid(i + 1, j) : zeroCell;
				GridCell& cell_i_j_min_1 = (j > 0) ? currentGrid(i, j - 1) : zeroCell;
				GridCell& cell_i_j_plus_1 = (j < NCOLS - 1) ? currentGrid(i, j + 1) : zeroCell;

				double tnf_i_j_old = cell.getTNF();
				double tnf_i_min_1_j = cell_i_min_1_j.getTNF();
				double tnf_i_plus_1_j = cell_i_plus_1_j.getTNF();
				double tnf_i_j_min_1 = cell_i_j_min_1.getTNF();
				double tnf_i_j_plus_1 = cell_i_j_plus_1.getTNF();

				double res = tnf_i_j_old +
				muTNF * (tnf_i_min_1_j + tnf_i_plus_1_j + tnf_i_j_min_1 + tnf_i_j_plus_1 - 4 * tnf_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

				res *= degTNF;
				newCell.setTNF(res);

				double shedtnfr2_i_j_old = cell.getShedTNFR2();
				double shedtnfr2_i_min_1_j = cell_i_min_1_j.getShedTNFR2();
				double shedtnfr2_i_plus_1_j = cell_i_plus_1_j.getShedTNFR2();
				double shedtnfr2_i_j_min_1 = cell_i_j_min_1.getShedTNFR2();
				double shedtnfr2_i_j_plus_1 = cell_i_j_plus_1.getShedTNFR2();

				res = shedtnfr2_i_j_old +
				muShedTNFR2 * (shedtnfr2_i_min_1_j + shedtnfr2_i_plus_1_j + shedtnfr2_i_j_min_1 + shedtnfr2_i_j_plus_1 - 4 * shedtnfr2_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

				// degradation of shedTNFR2 is neglected because unbinding reaction is much faster
				res = res * (1 - _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * dt);
				newCell.setShedTNFR2(res);
				double temp = newCell.getTNF();
				newCell.setTNF(temp + _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * res * dt);

                
                double il10_i_j_old = cell.getIL10();
                double il10_i_min_1_j = cell_i_min_1_j.getIL10();
                double il10_i_plus_1_j = cell_i_plus_1_j.getIL10();
                double il10_i_j_min_1 = cell_i_j_min_1.getIL10();
                double il10_i_j_plus_1 = cell_i_j_plus_1.getIL10();
                
                
                res = il10_i_j_old + muIL10 * (il10_i_min_1_j + il10_i_plus_1_j + il10_i_j_min_1 + il10_i_j_plus_1 - 4 * il10_i_j_old);
                
                if (res <= _cutOffValue) {
                    res = 0;
                }
                
                res *= degIL10;
                newCell.setIL10(res);
                
                
				double ccl2_i_j_old = cell.getCCL2();
				double ccl2_i_min_1_j = cell_i_min_1_j.getCCL2();
				double ccl2_i_plus_1_j = cell_i_plus_1_j.getCCL2();
				double ccl2_i_j_min_1 = cell_i_j_min_1.getCCL2();
				double ccl2_i_j_plus_1 = cell_i_j_plus_1.getCCL2();

				res = ccl2_i_j_old +
				muChemokines * (ccl2_i_min_1_j + ccl2_i_plus_1_j + ccl2_i_j_min_1 + ccl2_i_j_plus_1 - 4 * ccl2_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;

				newCell.setCCL2(res);
				newCell.setCCL5(res * ratioCCL5toCCL2);
				newCell.setCXCL9(res * ratioCXCL9toCCL2);

				double macAttractant_i_j_old = cell.getMacAttractant();
				double macAttractant_i_min_1_j = cell_i_min_1_j.getMacAttractant();
				double macAttractant_i_plus_1_j = cell_i_plus_1_j.getMacAttractant();
				double macAttractant_i_j_min_1 = cell_i_j_min_1.getMacAttractant();
				double macAttractant_i_j_plus_1 = cell_i_j_plus_1.getMacAttractant();

				res = macAttractant_i_j_old +
				muChemokines * (macAttractant_i_min_1_j + macAttractant_i_plus_1_j + macAttractant_i_j_min_1 + macAttractant_i_j_plus_1 - 4 * macAttractant_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;
				newCell.setMacAttractant(res);
			}	/* j, columns */
		}	/*i, rows */

		grSim.swap();

	}	/* t, timesteps */
}
