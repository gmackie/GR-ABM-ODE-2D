/*
 * GrDiffusionADE_Swap.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: ncilfone
 */

#include "grdiffusionadeswap.h"

GrDiffusionADE_Swap::GrDiffusionADE_Swap()
  : _cutOffValue(0.00000000000001)
{
}

GrDiffusionADE_Swap::~GrDiffusionADE_Swap()
{
}

using namespace std;

// 2-D Alternating Direction Explicit Algorithm

// On the Solution of the Diffusion Equations by Numerical Methods
// H.Z. BARAKAT and J.A. Clark
// Journal of Heat Transfer November, 1966 pp 421-427

// Note that the derivation in this paper uses i=cols and j=rows
// We use i=rows and j=cols thus our coefficients b and c are different
// than what is shown in the paper (they are switched)

// Please see the pdf included in the docs folder of GR-ABM-ODE for more information

// This function is the main ADE diffusion algorithm.
//Assumes padding on borders
static void diffuse_u(const Scalar* __restrict__ grid_u, Scalar* __restrict__ newgrid_u, const Pos& dim, const Scalar u_diffuse, const Scalar cutOff)
{

  register Scalar up_u, dn_u, lt_u, rt_u, ct_u;

  // Diffusion Constants for Grid and Time
  const Scalar dt = _PARAM(PARAM_GR_DT_DIFFUSION); //time step (sec)
  const Scalar dx2 = 20e-4*20e-4;
  const Scalar dy2 = 20e-4*20e-4;

  // Diffusion Constants for Any Molecule
  const Scalar lambda = u_diffuse * dt * ((1.0/dx2) + (1.0/dy2));
  const Scalar coefficient_a = (1-lambda)/(1+lambda);
  const Scalar coefficient_b = ((u_diffuse * dt)/(dy2))/(1+lambda);
  const Scalar coefficient_c = ((u_diffuse * dt)/(dx2))/(1+lambda);


  // ADE SCHEME FOR U
  // up: i-1, j n+1
  // dn: i+1, j n
  // lt: i, j-1 n+1
  // rt: i, j+1 n
  // ct: i, j   n

  for(int i=0; i<GETROW(dim); i++)
    {
      for(int j=0; j<GETCOL(dim); j++)
        {
          up_u  = newgrid_u[Indexer::padInd(dim, i-1,j)];
          lt_u  = newgrid_u[Indexer::padInd(dim, i,j-1)];
          ct_u  = grid_u[Indexer::padInd(dim, i,j)];
          rt_u  = grid_u[Indexer::padInd(dim, i,j+1)];
          dn_u  = grid_u[Indexer::padInd(dim, i+1,j)];

          Scalar& new_u = newgrid_u[Indexer::padInd(dim, i,j)];

          new_u = (coefficient_a * ct_u) + (coefficient_b * (up_u + dn_u)) + (coefficient_c * (lt_u + rt_u));
          new_u = new_u <= cutOff ? Scalar(0.0) : new_u;

        }
    }
}

////Assumes padding on borders
static void diffuse_v(const Scalar* __restrict__ grid_v, Scalar* __restrict__ newgrid_v, const Pos& dim, const Scalar v_diffuse, const Scalar cutOff)
{

  register Scalar up_v, dn_v, lt_v, rt_v, ct_v;

  // Diffusion Constants for Grid and Time
  const Scalar dt = _PARAM(PARAM_GR_DT_DIFFUSION); //time step (sec)
  const Scalar dx2 = 20e-4*20e-4;
  const Scalar dy2 = 20e-4*20e-4;

  // Diffusion Constants for Any Molecule
  const Scalar lambda = v_diffuse * dt * ((1.0/dx2) + (1.0/dy2));
  const Scalar coefficient_a = (1-lambda)/(1+lambda);
  const Scalar coefficient_b = ((v_diffuse * dt)/(dy2))/(1+lambda);
  const Scalar coefficient_c = ((v_diffuse * dt)/(dx2))/(1+lambda);

  // ADE SCHEME FOR V
  // up: i-1, j n
  // dn: i+1, j n+1
  // lt: i, j-1 n
  // rt: i, j+1 n+1
  // ct: i, j   n
  for(int i=GETROW(dim)-1; i >=0; i--)
    {
      for(int j=GETCOL(dim)-1; j >= 0; j--)
        {
          up_v  = grid_v[Indexer::padInd(dim, i-1,j)];
          lt_v  = grid_v[Indexer::padInd(dim, i,j-1)];
          ct_v  = grid_v[Indexer::padInd(dim, i,j)];
          rt_v  = newgrid_v[Indexer::padInd(dim, i,j+1)];
          dn_v  = newgrid_v[Indexer::padInd(dim, i+1,j)];

          Scalar& new_v = newgrid_v[Indexer::padInd(dim, i,j)];
#if 0   //TODO: make loops faster with manual prefetching
          __builtin_prefetch(grid_v + Indexer::padInd(dim, i-1, j-1), 0);
          __builtin_prefetch(grid_v + Indexer::padInd(dim, i, j-2), 0);
          __builtin_prefetch(newgrid_v + Indexer::padInd(dim, i, j-1), 1);
#endif  //PREFETCH

          new_v = (coefficient_a * ct_v) + (coefficient_b * (up_v + dn_v)) + (coefficient_c * (lt_v + rt_v));
          new_v = new_v <= cutOff ? Scalar(0.0) : new_v;

        }
    }
}

// This function is used to reconcile the u and v values to the actual concentration for the cytokines.
static void diffuse_avg(Scalar* __restrict__ newgrid_u, Scalar* __restrict__ newgrid_v, Scalar* __restrict__ newgrid, const Pos& dim)
{
  for(int i=0; i<GETROW(dim); i++)
    {
      for(int j=0; j<GETCOL(dim); j++)
        {
          const Scalar _u  = newgrid_u[Indexer::padInd(dim, i,j)];
          const Scalar _v  = newgrid_v[Indexer::padInd(dim, i,j)];

          Scalar& newconc = newgrid[Indexer::padInd(dim,i,j)];

          newconc = (_u + _v)/2.0;
        }
    }

}

// This function is used to reconcile the u and v values to the actual concentration for the chemokines.
// It takes in extra grids for CCL5 and CXCL9 since they are ratios of CCL2
static void diffuse_avg_ratio(Scalar* newgrid_u, Scalar* newgrid_v, Scalar* newgridccl2, Scalar* newgridccl5, Scalar* newgridcxcl9, const Pos& dim)
{
  // NOTE: THE _SEC_RATE_CCL should both be multiplied by the MOLECULAR_DT but since these terms will cancel
  // in the division they are not used to keep the error of the floating point division down.

  const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
  const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

  for(int i=0; i<GETROW(dim); i++)
    {
      for(int j=0; j<GETCOL(dim); j++)
        {
          Scalar _u  = newgrid_u[Indexer::padInd(dim, i,j)];
          Scalar _v  = newgrid_v[Indexer::padInd(dim, i,j)];

          Scalar& newconc = newgridccl2[Indexer::padInd(dim,i,j)];

          newconc = (_u + _v)/2;

          Scalar& ccl5 = newgridccl5[Indexer::padInd(dim,i,j)];
          Scalar& cxcl9 = newgridcxcl9[Indexer::padInd(dim,i,j)];

          ccl5 = newconc * ratioCCL5toCCL2;
          cxcl9 = newconc * ratioCXCL9toCCL2;
        }
    }

}

// This function is used to degrade the molecules on the grid after they have diffused.
static void diffuse_degrade(GrGrid& nextGrid)
{
  const Scalar dt = _PARAM(PARAM_GR_DT_DIFFUSION);

  const Scalar TNFdegRate = pow(2.7183, (-1.0 *  _PARAM(PARAM_GR_K_DEG) * dt));
  const Scalar IL10degRate = pow(2.7183, (-1.0 *  _PARAM(PARAM_GR_I_K_DEG) * dt));
  const Scalar CHEMOKINEdegRate = pow(2.7183, (-1.0 *  _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));

  //Degradation equation here is separated based on operator splitting causing it to lag
  //behind diffusion.  We make the assumption diffuse << degrade to
  //show degradation is much slower and incurs no stability issues. This argument should hold for
  //most degradation rates we use in the code. Please refer to the concept of Thiele Modulus for a mathematical
  //derivation of this argument.

  // Degradation is solved by the analytical solution to the linear degradation ODE
  // instead of by a numerical solver. This should vastly increase the accuracy of the
  // degradation in the ABM code. This is allowed due to the concept of operator splitting

  Pos p;

//    double Total=0.0;

  for(p.x = 0; p.x < nextGrid.getRange().x; ++p.x)
    for(p.y = 0; p.y < nextGrid.getRange().y; ++p.y)
      {
        // Degradation of sTNF
        nextGrid.setTNF(p, nextGrid.TNF(p) * TNFdegRate);
//            Total += nextGrid.TNF(p);
//            nextGrid.incTNF(p, (-1.0 * nextGrid.TNF(p) * _PARAM(PARAM_GR_K_DEG) * dt));

        // Degradation of il10
        nextGrid.setil10(p, nextGrid.il10(p) * IL10degRate);
//            nextGrid.setil10(p, (-1.0 * nextGrid.il10(p) * _PARAM(PARAM_GR_I_K_DEG) * dt));

        // Degradation of CCL2
        nextGrid.setCCL2(p, nextGrid.CCL2(p) * CHEMOKINEdegRate);
//            nextGrid.incCCL2(p, (-1.0 * nextGrid.CCL2(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));

        // Degradation of CCL5
        nextGrid.setCCL5(p, nextGrid.CCL5(p) * CHEMOKINEdegRate);
//            nextGrid.incCCL5(p, (-1.0 * nextGrid.CCL5(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));

        // Degradation of CXCL9
        nextGrid.setCXCL9(p, nextGrid.CXCL9(p) * CHEMOKINEdegRate);
//            nextGrid.incCXCL9(p, (-1.0 * nextGrid.CXCL9(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));

        if (_PARAM(PARAM_GR_SEC_RATE_ATTRACTANT) != 0)
          {
            // Degradation of MacAttractant
            nextGrid.setmacAttractant(p, nextGrid.macAttractant(p) * CHEMOKINEdegRate);
//                nextGrid.incmacAttractant(p, (-1.0 * nextGrid.macAttractant(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));
          }

        // Unbinding of sTNF from ShedTNFR2
        nextGrid.incshedTNFR2(p, (-1.0 * nextGrid.shedTNFR2(p) * _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * dt));

        // Addition of sTNF from unbinding from ShedTNFR2
        nextGrid.incTNF(p, (1.0 * nextGrid.shedTNFR2(p) * _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * dt));

      }

//    std::cout << "Total Conc: " << Total << std::endl;
}


void GrDiffusionADE_Swap::diffuse(GrSimulationGrid& grSim) const
{
  GrGrid& grid = grSim.getCurrentGrid();
  GrGrid& nextGrid = grSim.getNextGrid();

  // Note that the same function is used to diffuse both u and v.
  // This is due to the fact that the v matrix is stored backwards thus allowing us to iterate in the forward direction

  #pragma omp parallel sections //if(_threaded)   //Spawns a new team of threads
  {
    #pragma omp section
    {
      ::diffuse_u(grid.u_TNF(), nextGrid.u_TNF(), grid.getRange(), _PARAM(PARAM_GR_D_TNF), _cutOffValue);
      ::diffuse_v(grid.v_TNF(), nextGrid.v_TNF(), grid.getRange(), _PARAM(PARAM_GR_D_TNF), _cutOffValue);
      ::diffuse_avg(nextGrid.u_TNF(), nextGrid.v_TNF(), nextGrid.TNF(), grid.getRange());
    }
    #pragma omp section
    {
      ::diffuse_u(grid.u_shedTNFR2(), nextGrid.u_shedTNFR2(), grid.getRange(), _PARAM(PARAM_GR_D_SHED_TNFR2), _cutOffValue);
      ::diffuse_v(grid.v_shedTNFR2(), nextGrid.v_shedTNFR2(), grid.getRange(), _PARAM(PARAM_GR_D_SHED_TNFR2), _cutOffValue);
      ::diffuse_avg(nextGrid.u_shedTNFR2(), nextGrid.v_shedTNFR2(), nextGrid.shedTNFR2(), grid.getRange());
    }
    #pragma omp section
    {
      if (_PARAM(PARAM_GR_SEC_RATE_ATTRACTANT) != 0)
        {
          std::cout << "Why is this running?!?" << std::endl;
          ::diffuse_u(grid.u_macAttractant(), nextGrid.u_macAttractant(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
          ::diffuse_v(grid.v_macAttractant(), nextGrid.v_macAttractant(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
          ::diffuse_avg(nextGrid.u_macAttractant(), nextGrid.v_macAttractant(), nextGrid.macAttractant(), grid.getRange());
        }

    }
    #pragma omp section
    {
      ::diffuse_u(grid.u_CCL2(), nextGrid.u_CCL2(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
      ::diffuse_v(grid.v_CCL2(), nextGrid.v_CCL2(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
      ::diffuse_avg_ratio(nextGrid.u_CCL2(), nextGrid.v_CCL2(), nextGrid.CCL2(), nextGrid.CCL5(), nextGrid.CXCL9(), grid.getRange());

    }
    #pragma omp section
    {
      ::diffuse_u(grid.u_il10(), nextGrid.u_il10(), grid.getRange(), _PARAM(PARAM_GR_D_IL10), _cutOffValue);
      ::diffuse_v(grid.v_il10(), nextGrid.v_il10(), grid.getRange(), _PARAM(PARAM_GR_D_IL10), _cutOffValue);
      ::diffuse_avg(nextGrid.u_il10(), nextGrid.v_il10(), nextGrid.il10(), grid.getRange());
    }
  } //omp parallel

  ::diffuse_degrade(nextGrid);

  grSim.swap();
}

