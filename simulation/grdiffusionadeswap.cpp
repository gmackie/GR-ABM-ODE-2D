/*
 * GrDiffusionADE_Swap.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: ncilfone
 */

#include "grdiffusionadeswap.h"

GrDiffusionADE_Swap::GrDiffusionADE_Swap()
: _cutOffValue(0.000001)
{
}

GrDiffusionADE_Swap::~GrDiffusionADE_Swap()
{
}


//Assumes padding on borders
static void diffuse_u(const Scalar* __restrict__ grid_u, Scalar* __restrict__ newgrid_u, const Pos& dim, const Scalar u_diffuse, const Scalar cutOff) {

    register Scalar up_u, dn_u, lt_u, rt_u, ct_u;

    // Diffusion Constants for Grid and Time
    const Scalar dt = _PARAM(PARAM_GR_DT_DIFFUSION); //time step (sec)
    const Scalar dx = 20e-4;
    const Scalar dy = 20e-4;
    
    // Diffusion Constants for Any Molecule
    const Scalar lambda = u_diffuse * dt * ((1/pow(dx,2)) + (1/pow(dy,2)));
    const Scalar coefficient_a = (1-lambda)/(1+lambda);
    const Scalar coefficient_b = ((u_diffuse * dt)/(pow(dy,2)))/(1+lambda);
    const Scalar coefficient_c = ((u_diffuse * dt)/(pow(dx,2)))/(1+lambda);
    
    // ADE SCHEME FOR U
    // up: i-1, j n+1
    // dn: i+1, j n
    // lt: i, j-1 n+1
    // rt: i, j+1 n
    // ct: i, j   n
    
    for(int i=0;i<GETROW(dim);i++)
    {
        for(int j=0;j<GETCOL(dim);j++)
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

//Assumes padding on borders
static void diffuse_v(const Scalar* __restrict__ grid_v, Scalar* __restrict__ newgrid_v, const Pos& dim, const Scalar v_diffuse, const Scalar cutOff) {
    
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
    for(int i=GETROW(dim)-1;i >=0; i--)
    {
        for(int j=GETCOL(dim)-1;j >= 0; j--)
        {
            up_v  = grid_v[Indexer::padInd(dim, i-1,j)];
            lt_v  = grid_v[Indexer::padInd(dim, i,j-1)];
            ct_v  = grid_v[Indexer::padInd(dim, i,j)];
            rt_v  = newgrid_v[Indexer::padInd(dim, i,j+1)];
            dn_v  = newgrid_v[Indexer::padInd(dim, i+1,j)];

            Scalar& new_v = newgrid_v[Indexer::padInd(dim, i,j)];
#if 1
            __builtin_prefetch(grid_v + Indexer::padInd(dim, i-1, j-1), 0);
            __builtin_prefetch(grid_v + Indexer::padInd(dim, i, j-2), 0);
            __builtin_prefetch(newgrid_v + Indexer::padInd(dim, i, j-1), 1);
#endif  //PREFETCH
            
            new_v = (coefficient_a * ct_v) + (coefficient_b * (up_v + dn_v)) + (coefficient_c * (lt_v + rt_v));
            new_v = new_v <= cutOff ? Scalar(0.0) : new_v;
            
        }
    }
}

static void diffuse_avg(Scalar* __restrict__ newgrid_u, Scalar* __restrict__ newgrid_v, Scalar* __restrict__ newgrid, const Pos& dim)
{
    for(int i=0;i<GETROW(dim);i++)
    {
        for(int j=0;j<GETCOL(dim);j++)
        {
            const Scalar _u  = newgrid_u[Indexer::padInd(dim, i,j)];
            const Scalar _v  = newgrid_v[Indexer::padInd(dim, dim.x - i-1, dim.y - j-1)];
            
            Scalar& newconc = newgrid[Indexer::padInd(dim,i,j)];
            
            newconc = (_u + _v)/2.0;
            
        }
    }
    
}

static void diffuse_avg_ratio(Scalar* newgrid_u, Scalar* newgrid_v, Scalar* newgridccl2, Scalar* newgridccl5, Scalar* newgridcxcl9, const Pos& dim)
{
    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
	const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
    
    for(int i=0;i<GETROW(dim);i++)
    {
        for(int j=0;j<GETCOL(dim);j++)
        {
            Scalar _u  = newgrid_u[Indexer::padInd(dim, i,j)];
            Scalar _v  = newgrid_v[Indexer::padInd(dim, dim.x - i - 1, dim.y - j - 1)];
            
            Scalar& newconc = newgridccl2[Indexer::padInd(dim,i,j)];

            newconc = (_u + _v)/2;
            
            Scalar& ccl5 = newgridccl5[Indexer::padInd(dim,i,j)];
            Scalar& cxcl9 = newgridcxcl9[Indexer::padInd(dim,i,j)];
            
            ccl5 = newconc * ratioCCL5toCCL2;
            cxcl9 = newconc * ratioCXCL9toCCL2;
        }
    }
    
}


static void diffuse_degrade(GrGrid& nextGrid)
{
    const Scalar dt = _PARAM(PARAM_GR_DT_DIFFUSION);

    Pos p;
    for(p.x = 0; p.x < nextGrid.getRange().x; ++p.x)
        for(p.y = 0; p.y < nextGrid.getRange().y; ++p.y) {
            
            // Degradation of sTNF
            nextGrid.incTNF(p, (-1.0 * nextGrid.TNF(p) * _PARAM(PARAM_GR_K_DEG) * dt));
            
            // Degradation of sTNF
            nextGrid.incil10(p, (-1.0 * nextGrid.il10(p) * _PARAM(PARAM_GR_I_K_DEG) * dt));
            
            // Degradation of CCL2
            nextGrid.incCCL2(p, (-1.0 * nextGrid.CCL2(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));
            
            // Degradation of CCL5
            nextGrid.incCCL5(p, (-1.0 * nextGrid.CCL5(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));
            
            // Degradation of CXCL9
            nextGrid.incCXCL9(p, (-1.0 * nextGrid.CXCL9(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));
            
            // Degredation of MacAttractant
            nextGrid.incmacAttractant(p, (-1.0 * nextGrid.macAttractant(p) * _PARAM(PARAM_GR_CHEMOKINE_K_DEG) * dt));
            
            // Unbinding of sTNF from ShedTNFR2
            nextGrid.incshedTNFR2(p, (-1.0 * nextGrid.shedTNFR2(p) * _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * dt));
            
            // Addition of sTNF from unbinding from ShedTNFR2
            nextGrid.incTNF(p, (1.0 * nextGrid.shedTNFR2(p) * _PARAM(PARAM_GR_KD2) * _PARAM(PARAM_GR_K_ON2) * dt));

        }
}


void GrDiffusionADE_Swap::diffuse(GrSimulationGrid& grSim) const {
	GrGrid& grid = grSim.getCurrentGrid();
	GrGrid& nextGrid = grSim.getNextGrid();
    
    
#pragma omp parallel sections //if(_threaded)   //Spawns a new team of threads
    {
#pragma omp section
        {
            ::diffuse_u(grid.u_TNF(), nextGrid.u_TNF(), grid.getRange(), _PARAM(PARAM_GR_D_TNF), _cutOffValue);
            ::diffuse_u(grid.v_TNF(), nextGrid.v_TNF(), grid.getRange(), _PARAM(PARAM_GR_D_TNF), _cutOffValue);
            ::diffuse_avg(nextGrid.u_TNF(), nextGrid.v_TNF(), nextGrid.TNF(), grid.getRange());
        }
#pragma omp section
        {
            ::diffuse_u(grid.u_shedTNFR2(), nextGrid.u_shedTNFR2(), grid.getRange(), _PARAM(PARAM_GR_D_SHED_TNFR2), _cutOffValue);
            ::diffuse_u(grid.v_shedTNFR2(), nextGrid.v_shedTNFR2(), grid.getRange(), _PARAM(PARAM_GR_D_SHED_TNFR2), _cutOffValue);
            ::diffuse_avg(nextGrid.u_shedTNFR2(), nextGrid.v_shedTNFR2(), nextGrid.shedTNFR2(), grid.getRange());
        }
#pragma omp section
        {
            ::diffuse_u(grid.u_macAttractant(), nextGrid.u_macAttractant(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
            ::diffuse_u(grid.v_macAttractant(), nextGrid.v_macAttractant(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
            ::diffuse_avg(nextGrid.u_macAttractant(), nextGrid.v_macAttractant(), nextGrid.macAttractant(), grid.getRange());
        }
#pragma omp section
        {
            ::diffuse_u(grid.u_CCL2(), nextGrid.u_CCL2(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
            ::diffuse_u(grid.v_CCL2(), nextGrid.v_CCL2(), grid.getRange(), _PARAM(PARAM_GR_D_CHEMOKINES), _cutOffValue);
            ::diffuse_avg_ratio(nextGrid.u_CCL2(), nextGrid.v_CCL2(), nextGrid.CCL2(), nextGrid.CCL5(), nextGrid.CXCL9(), grid.getRange());

        }
#pragma omp section
        {
            ::diffuse_u(grid.u_il10(), nextGrid.u_il10(), grid.getRange(), _PARAM(PARAM_GR_D_IL10), _cutOffValue);
            ::diffuse_u(grid.v_il10(), nextGrid.v_il10(), grid.getRange(), _PARAM(PARAM_GR_D_IL10), _cutOffValue);
            ::diffuse_avg(nextGrid.u_il10(), nextGrid.v_il10(), nextGrid.il10(), grid.getRange());
        }
    } //omp parallel
    
    ::diffuse_degrade(nextGrid);
    
    grSim.swap();
}

