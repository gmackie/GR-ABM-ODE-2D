/*
 * grdiffusionwrongbtcs.cpp
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#include "grdiffusionwrongbtcs.h"
#include "grgrid.h"
#include "params.h"

GrDiffusionWrongBTCS::GrDiffusionWrongBTCS()
{
}

GrDiffusionWrongBTCS::~GrDiffusionWrongBTCS()
{
}

void GrDiffusionWrongBTCS::diffuse(GrSimulationGrid&, const int time) const
{
  std::cerr << " GrDiffusionWrongBTCS::diffuse no longer supported: not updated for grid swapping or mac attractant" << std::endl;
  exit(1);

#if 0
  GrGrid& grid = simGrid.getGrid();

  const double muTNF = _PARAM(_diffusivityTNF) * 6 / (4e-6);
  const double muChemokines = _PARAM(_diffusivityChemokines) * 6 / (4e-6);
  const double degTNF = _PARAM(_degRateTNF);
  const double degChemokines = _PARAM(_degRateChemokines);
  const double dAttractant = _PARAM(_dAttractant);
  const double ok = 0.001; // allowable error for good enough convergence

  // we have to use the heap here, since oldGrid does not fit in the stack
  GrGrid* pOldGrid = (GrGrid*) malloc(sizeof(GrGrid));

  double errTNF = 2 * ok;
  double errCCL2 = 2 * ok;
  double errCCL5 = 2 * ok;
  double errCXCL9 = 2 * ok;

  // dt = 6s, solve for 10 minutes = 100 * 6 seconds
  for (int t = 0; t < 100; t++)
    {
      for (int i = 0; i < NROWS; i++)
        {
          for (int j = 0; j < NCOLS; j++)
            {
              GridCell& cell = grid(i, j);

              // secrete
              Agent* pAgent = cell.getAgent(0);
              if (pAgent)
                {
                  printf("SEC (%d,%d)\n", i, j);
                  pAgent->secrete(grid, false);
                }

              pAgent = cell.getAgent(1);
              if (pAgent)
                {
                  printf("SEC (%d,%d)\n", i, j);
                  pAgent->secrete(grid, false);
                }

              if (cell.isCaseated())
                cell.incMacAttractant(dAttractant);
            }
        }

      // update pOldGrid
      memcpy(pOldGrid, &grid, sizeof(GrGrid));

      for (int k = 0; k < 100 && (errTNF > ok || errCCL2 > ok || errCCL5 > ok || errCXCL9 > ok);  k++)
        {
          double localErrTNF = 0, localErrCCL2 = 0, localErrCCL5 = 0, localErrCXCL9 = 0;
          for (int i = NROWS - 1; i >= 0; i--)
            {
              for (int j = NCOLS - 1; j >= 0; j--)
                {
                  GridCell& oldCell = (*pOldGrid)(i, j);
                  GridCell& cell = grid(i, j);
                  GridCell& cell_i_min_1_j = (i > 0) ? (*pOldGrid)(i - 1, j) : cell;
                  GridCell& cell_i_plus_1_j =  (i < NROWS - 1) ? (*pOldGrid)(i + 1, j) : cell;
                  GridCell& cell_i_j_min_1 = (j > 0) ? (*pOldGrid)(i, j - 1) : cell;
                  GridCell& cell_i_j_plus_1 = (j < NCOLS - 1) ? (*pOldGrid)(i, j + 1) : cell;

                  if (errTNF > ok)
                    {
                      double tnf_i_j_old = oldCell.getTNF();
                      double tnf_i_j = cell.getTNF();
                      double tnf_i_min_1_j = (i > 0) ? cell_i_min_1_j.getTNF() : 0;
                      double tnf_i_plus_1_j = (i < NROWS - 1) ? cell_i_plus_1_j.getTNF() : 0;
                      double tnf_i_j_min_1 = (j > 0) ? cell_i_j_min_1.getTNF() : 0;
                      double tnf_i_j_plus_1 = (j < NCOLS - 1) ? cell_i_j_plus_1.getTNF() : 0;

                      double res = tnf_i_j_old - (muTNF + 1) * tnf_i_j +
                                   0.25 * muTNF * (tnf_i_min_1_j + tnf_i_plus_1_j + tnf_i_j_min_1 + tnf_i_j_plus_1);
                      grid(i, j).incTNF(res);

                      localErrTNF += fabs(res);
                    }
                  if (errCCL2 > ok)
                    {
                      double ccl2_i_j_old = oldCell.getCCL2();
                      double ccl2_i_j = cell.getCCL2();
                      double ccl2_i_min_1_j = (i > 0) ? cell_i_min_1_j.getCCL2() : 0;
                      double ccl2_i_plus_1_j = (i < NROWS - 1) ? cell_i_plus_1_j.getCCL2() : 0;
                      double ccl2_i_j_min_1 = (j > 0) ? cell_i_j_min_1.getCCL2() : 0;
                      double ccl2_i_j_plus_1 = (j < NCOLS - 1) ? cell_i_j_plus_1.getCCL2() : 0;

                      double res = ccl2_i_j_old - (muChemokines + 1) * ccl2_i_j +
                                   0.25 * muChemokines * (ccl2_i_min_1_j + ccl2_i_plus_1_j + ccl2_i_j_min_1 + ccl2_i_j_plus_1);
                      grid(i, j).incCCL2(res);

                      localErrCCL2 += fabs(res);
                    }
                  if (errCCL5 > ok)
                    {
                      double ccl5_i_j_old = oldCell.getCCL5();
                      double ccl5_i_j = cell.getCCL5();
                      double ccl5_i_min_1_j = (i > 0) ? cell_i_min_1_j.getCCL5() : 0;
                      double ccl5_i_plus_1_j = (i < NROWS - 1) ? cell_i_plus_1_j.getCCL5() : 0;
                      double ccl5_i_j_min_1 = (j > 0) ? cell_i_j_min_1.getCCL5() : 0;
                      double ccl5_i_j_plus_1 = (j < NCOLS - 1) ? cell_i_j_plus_1.getCCL5() : 0;

                      double res = ccl5_i_j_old - (muChemokines + 1) * ccl5_i_j +
                                   0.25 * muChemokines * (ccl5_i_min_1_j + ccl5_i_plus_1_j + ccl5_i_j_min_1 + ccl5_i_j_plus_1);
                      grid(i, j).incCCL5(res);

                      localErrCCL5 += fabs(res);
                    }
                  if (errCXCL9 > ok)
                    {
                      double cxcl9_i_j_old = oldCell.getCXCL9();
                      double cxcl9_i_j = cell.getCXCL9();
                      double cxcl9_i_min_1_j = (i > 0) ? cell_i_min_1_j.getCXCL9() : 0;
                      double cxcl9_i_plus_1_j = (i < NROWS - 1) ? cell_i_plus_1_j.getCXCL9() : 0;
                      double cxcl9_i_j_min_1 = (j > 0) ? cell_i_j_min_1.getCXCL9() : 0;
                      double cxcl9_i_j_plus_1 = (j < NCOLS - 1) ? cell_i_j_plus_1.getCXCL9() : 0;

                      double res = cxcl9_i_j_old - (muChemokines + 1) * cxcl9_i_j +
                                   0.25 * muChemokines * (cxcl9_i_min_1_j + cxcl9_i_plus_1_j + cxcl9_i_j_min_1 + cxcl9_i_j_plus_1);
                      grid(i, j).incCXCL9(res);

                      localErrCXCL9 += fabs(res);
                    }
                }	/* j, columns */
            }	/*i, rows */

          errTNF = localErrTNF;
          errCCL2 = localErrCCL2;
          errCCL5 = localErrCCL5;
          errCXCL9 = localErrCXCL9;
        }	/* k, #iterations */

      // update _oldGrid, and compute degradation
      for (int i = 0; i < NROWS; i++)
        {
          for (int j = 0; j < NCOLS; j++)
            {
              GridCell& cell = grid(i, j);
              cell.setTNF(cell.getTNF() * degTNF);
              cell.setCCL2(cell.getCCL2() * degChemokines);
              cell.setCCL5(cell.getCCL5() * degChemokines);
              cell.setCXCL9(cell.getCXCL9() * degChemokines);
            }
        }
    }	/* t, timesteps */

  //delete pOldGrid;
  free(pOldGrid);
#endif
}
