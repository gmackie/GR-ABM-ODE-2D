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

//Assumes padding on borders
static void diffuse(const Scalar* grid, Scalar* newgrid, const Pos& dim, const Scalar diffuse, const Scalar degrade, const Scalar cutOff)
{
  Scalar up, dn, lt, rt, ct;
#define DT (6.0)
  for(int i=0; i<GETROW(dim); i++)
    for(int j=0; j<GETCOL(dim); j++)
      {
        up  = grid[Indexer::padInd(dim, i-1,j)];
        lt  = grid[Indexer::padInd(dim, i,j-1)];
        ct  = grid[Indexer::padInd(dim, i,j)];
        rt  = grid[Indexer::padInd(dim, i,j+1)];
        dn  = grid[Indexer::padInd(dim, i+1,j)];
        Scalar& nc = newgrid[Indexer::padInd(dim, i,j)];
        nc = ct + diffuse * (up + dn + lt + rt - Scalar(4.0) * ct);
        //Degradation equation here is incorrect, the effect is degradation lags
        //behind diffusion.  We make the assumption diffuse << degrade to
        //show degradation is much slower and incurs no stability issues
        nc *= (1.0 - degrade * DT);
        nc = nc <= cutOff ? Scalar(0.0) : nc;
      }
}
/*static void diffuse_recon(const Scalar* u, const Scalar* v, Scalar* grid) {
  #define DT (6.0)
  for(int i=0;i<GETROW(dim);i++)
    for(int j=0;j<GETCOL(dim);j++)
    {
	grid[Indexer] = (u[Indexer::padInd] + v[Indexer])/2;
    }
}*/

void GrDiffusionFTCS_Swap::diffuse(GrSimulationGrid& grSim, const int time) const
{
  GrGrid& grid = grSim.getCurrentGrid();
  GrGrid& nextGrid = grSim.getNextGrid();
  #pragma omp parallel sections //if(_threaded)   //Spawns a new team of threads
  {
    #pragma omp section
    {
      ::diffuse(grid.TNF(), nextGrid.TNF(), grid.getRange(), _PARAM(_diffusivityTNF) * 6 / (4e-6), _PARAM(_kDeg), _cutOffValue);
    }
    #pragma omp section
    {
      ::diffuse(grid.shedTNFR2(), nextGrid.shedTNFR2(), grid.getRange(), _PARAM(_diffusivityShedTNFR2) * 6 / (4e-6), 1.0, _cutOffValue);
    }
    #pragma omp section
    {
      ::diffuse(grid.macAttractant(), nextGrid.macAttractant(), grid.getRange(), _PARAM(_diffusivityChemokines) * 6 / (4e-6), _PARAM(_ChemokinekDeg), _cutOffValue);
    }
    #pragma omp section
    {
      ::diffuse(grid.CCL2(), nextGrid.CCL2(), grid.getRange(), _PARAM(_diffusivityChemokines) * 6 / (4e-6), _PARAM(_ChemokinekDeg), _cutOffValue);
    }
    #pragma omp section
    {
      ::diffuse(grid.il10(), nextGrid.il10(), grid.getRange(), _PARAM(_diffusivityIL10) * 6 / (4e-6), _PARAM(_Ikdeg), _cutOffValue);
    }
  } //omp parallel

  const Scalar dt = 6.0;
  const Scalar ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);
  const Scalar ratioCXCL9toCCL2 = _PARAM(Mac_dCXCL9) / _PARAM(Mac_dCCL2);
  Pos p;
  for(p.x = 0; p.x < grid.getRange().x; ++p.x)
    for(p.y = 0; p.y < grid.getRange().y; ++p.y)
      {
        nextGrid.setshedTNFR2(p, nextGrid.shedTNFR2(p) * (1 - _PARAM(_KD2) * _PARAM(_kOn2) * dt));
        nextGrid.incTNF(p, (_PARAM(_KD2) * _PARAM(_kOn2) * nextGrid.shedTNFR2(p) * dt));
        nextGrid.setCXCL9(p, nextGrid.CCL2(p) * ratioCXCL9toCCL2);
        nextGrid.setCCL5(p, nextGrid.CCL2(p) * ratioCCL5toCCL2);
      }
  grSim.swap();
}
/*
void GrDiffusionFTCS_Swap::diffuse(GrSimulationGrid& grSim) const
{
	GrGrid& currentGrid = grSim.getCurrentGrid();
	GrGrid& nextGrid = grSim.getNextGrid();

	const double muTNF = _PARAM(_diffusivityTNF) * 6 / (4e-6);
	const double muShedTNFR2 = _PARAM(_diffusivityShedTNFR2) * 6 / (4e-6);
	const double muChemokines = _PARAM(_diffusivityChemokines) * 6 / (4e-6);
    const double muIL10 = _PARAM(_diffusivityIL10) * 6 / (4e-6);
    const double degTNFode = _PARAM(_kDeg);
    const double degIL10ode = _PARAM(_Ikdeg);
    const double degChemokinesode = _PARAM(_ChemokinekDeg);
	const double ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);
	const double ratioCXCL9toCCL2 = _PARAM(Mac_dCXCL9) / _PARAM(Mac_dCCL2);


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

                double resdeg = res - (res * degTNFode * dt);

				newCell.setTNF(resdeg);

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
				res = res * (1 - _PARAM(_KD2) * _PARAM(_kOn2) * dt);
				newCell.setShedTNFR2(res);

                double temp = newCell.getTNF();
				newCell.setTNF(temp + _PARAM(_KD2) * _PARAM(_kOn2) * res * dt);


                double il10_i_j_old = cell.getIL10();
                double il10_i_min_1_j = cell_i_min_1_j.getIL10();
                double il10_i_plus_1_j = cell_i_plus_1_j.getIL10();
                double il10_i_j_min_1 = cell_i_j_min_1.getIL10();
                double il10_i_j_plus_1 = cell_i_j_plus_1.getIL10();


                res = il10_i_j_old + muIL10 * (il10_i_min_1_j + il10_i_plus_1_j + il10_i_j_min_1 + il10_i_j_plus_1 - 4 * il10_i_j_old);

                if (res <= _cutOffValue) {
                    res = 0;
                }

                resdeg = res - (res * degIL10ode * dt);

//
//                if (res > 0 || resdeg > 0) {
//                    double diff = (res - resdeg);
//                    cout << diff << std::endl;
//                }

//                newCell.setIL10(res);


//                if ( i == 50 && j == 50) {
//                    cout << il10_i_j_old << "    " << resdeg << "     " << res << std::endl;
//                }


                newCell.setIL10(resdeg);

				double ccl2_i_j_old = cell.getCCL2();
				double ccl2_i_min_1_j = cell_i_min_1_j.getCCL2();
				double ccl2_i_plus_1_j = cell_i_plus_1_j.getCCL2();
				double ccl2_i_j_min_1 = cell_i_j_min_1.getCCL2();
				double ccl2_i_j_plus_1 = cell_i_j_plus_1.getCCL2();

				res = ccl2_i_j_old +
				muChemokines * (ccl2_i_min_1_j + ccl2_i_plus_1_j + ccl2_i_j_min_1 + ccl2_i_j_plus_1 - 4 * ccl2_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

                resdeg = res - (res * degChemokinesode * dt);

				newCell.setCCL2(resdeg);
				newCell.setCCL5(resdeg * ratioCCL5toCCL2);
				newCell.setCXCL9(resdeg * ratioCXCL9toCCL2);

				double macAttractant_i_j_old = cell.getMacAttractant();
				double macAttractant_i_min_1_j = cell_i_min_1_j.getMacAttractant();
				double macAttractant_i_plus_1_j = cell_i_plus_1_j.getMacAttractant();
				double macAttractant_i_j_min_1 = cell_i_j_min_1.getMacAttractant();
				double macAttractant_i_j_plus_1 = cell_i_j_plus_1.getMacAttractant();

				res = macAttractant_i_j_old +
				muChemokines * (macAttractant_i_min_1_j + macAttractant_i_plus_1_j + macAttractant_i_j_min_1 + macAttractant_i_j_plus_1 - 4 * macAttractant_i_j_old);

				if (res <= _cutOffValue)
					res = 0;

                resdeg = res - (res * degChemokinesode * dt);

				newCell.setMacAttractant(resdeg);
			}	// j, columns
		}	// i, rows

		grSim.swap();

	}	// t, timesteps
}
*/
