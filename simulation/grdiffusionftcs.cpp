/*
 * grdiffusionftcs.cpp
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#include "grdiffusionftcs.h"
#include "params.h"

GrDiffusionFTCS::GrDiffusionFTCS()
	: _cutOffValue(0.000001)
	, _cpyGrid()
{
}

GrDiffusionFTCS::~GrDiffusionFTCS()
{
}

void GrDiffusionFTCS::diffuse(GrGrid& grid) const
{
	// dt = 6, (dx)^2 = 4e-6, see thesis
	const double muTNF = _PARAM(PARAM_GR_D_TNF) * 6 / (4e-6);
	const double muChemokines = _PARAM(PARAM_GR_D_CHEMOKINES) * 6 / (4e-6);
	const double degTNF = _PARAM(PARAM_GR_DEG_TNF);
	const double degChemokines = _PARAM(PARAM_GR_DEG_CHEMOKINES);
	const double dAttractant = _PARAM(PARAM_GR_SEC_RATE_ATTRACTANT);

	// TNF and CCs are zero in zeroCell
	GridCell zeroCell;

	GrGrid* pNewGrid = &_cpyGrid;
	GrGrid* pOldGrid = &grid;

	// We need to copy the grid before doing diffusion, so the copy has
	// all the pointers to the agents. Otherwise secretion won't work
	// when the copy is the "old" grid.
	*pNewGrid = *pOldGrid;

	// dt = 6s, solve for 10 minutes = 100 * 6 seconds
	for (int t = 0; t < 100; t++)
	{
		for (int i = 0; i < NROWS; i++)
		{
			for (int j = 0; j < NCOLS; j++)
			{
				GridCell& newCell = (*pNewGrid)(i, j);
				GridCell& cell = (*pOldGrid)(i, j);
				GridCell& cell_i_min_1_j = (i > 0) ? (*pOldGrid)(i - 1, j) : zeroCell;
				GridCell& cell_i_plus_1_j =  (i < NROWS - 1) ? (*pOldGrid)(i + 1, j) : zeroCell;
				GridCell& cell_i_j_min_1 = (j > 0) ? (*pOldGrid)(i, j - 1) : zeroCell;
				GridCell& cell_i_j_plus_1 = (j < NCOLS - 1) ? (*pOldGrid)(i, j + 1) : zeroCell;

				// secrete
				Agent* pAgent = cell.getAgent(0);
				if (pAgent)
					pAgent->secrete(*pOldGrid);

				pAgent = cell.getAgent(1);
				if (pAgent)
					pAgent->secrete(*pOldGrid);

				if (cell.isCaseated())
					cell.incMacAttractant(dAttractant);

				/* TNF */
				double tnf_i_j = cell.getTNF();
				double res = tnf_i_j +
					muTNF * (cell_i_min_1_j.getTNF() + cell_i_plus_1_j.getTNF() +
							cell_i_j_min_1.getTNF() + cell_i_j_plus_1.getTNF() - 4 * tnf_i_j);

				if (res <= _cutOffValue)
					res = 0;

				res *= degTNF;
				newCell.setTNF(res);

				/* CCL2 */
				double ccl2_i_j = cell.getCCL2();
				res = ccl2_i_j +
					muChemokines * (cell_i_min_1_j.getCCL2() + cell_i_plus_1_j.getCCL2() +
							cell_i_j_min_1.getCCL2() + cell_i_j_plus_1.getCCL2() - 4 * ccl2_i_j);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;
				newCell.setCCL2(res);

				/* CCL5 */
				double ccl5_i_j = cell.getCCL5();
				res = ccl5_i_j +
					muChemokines * (cell_i_min_1_j.getCCL5() + cell_i_plus_1_j.getCCL5() +
							cell_i_j_min_1.getCCL5() + cell_i_j_plus_1.getCCL5() - 4 * ccl5_i_j);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;
				newCell.setCCL5(res);

				/* CXCL9 */
				double cxcl9_i_j = cell.getCXCL9();
				res = cxcl9_i_j +
					muChemokines * (cell_i_min_1_j.getCXCL9() + cell_i_plus_1_j.getCXCL9() +
							cell_i_j_min_1.getCXCL9() + cell_i_j_plus_1.getCXCL9() - 4 * cxcl9_i_j);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;
				newCell.setCXCL9(res);

				/* Mac attractant */
				double mac_i_j = cell.getMacAttractant();
				res = mac_i_j +
					muChemokines * (cell_i_min_1_j.getMacAttractant() + cell_i_plus_1_j.getMacAttractant() +
							cell_i_j_min_1.getMacAttractant() + cell_i_j_plus_1.getMacAttractant() - 4 * mac_i_j);

				if (res <= _cutOffValue)
					res = 0;

				res *= degChemokines;
				newCell.setMacAttractant(res);
			}	/* j, columns */
		}	/*i, rows */

		std::swap(pOldGrid, pNewGrid);
	}	/* t, timesteps */
}
