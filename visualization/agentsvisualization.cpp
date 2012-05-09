/*
 * agentsvisualization.cpp
 *
 *  Created on: 25-nov-2009
 *      Author: El-Kebir, M.
 */

#include "agentsvisualization.h"
#include <vector>
#include <assert.h>

AgentsVisualization::AgentsVisualization(int DIM, const ScalarAgentGrid* pScalarAgentGrid)
	: Visualization(DIM)
	, _pScalarAgentGrid(pScalarAgentGrid)
	, _drawGrid(true)
	, _drawTgam(true)
	, _drawTcyt(true)
	, _drawTreg(true)
	, _drawCas(true)
	, _drawSrc(true)
	, _drawExtMtb(true)
    , _drawm1m2(false)
	, _gridAlpha(0.7f)
	, _gridHeight(0.0f)
	, _drawSrcMac(false)
	, _drawSrcTgam(false)
	, _drawSrcTcyt(false)
	, _drawSrcTreg(false)
	, _selRow(-1)
	, _selCol(-1)
{
    memset(_macFilter, 1, sizeof(_macFilter));
	_deltaX = _MAX_X / (_DIM);
	_deltaY = _MAX_Y / (_DIM);
}

AgentsVisualization::~AgentsVisualization()
{
}

void AgentsVisualization::visualize(bool blend, const Simulation*, const ColorMap*) const
{
	if (blend)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		glDisable(GL_DEPTH_TEST);
	}
	else
	{
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
	}

	if (_drawGrid)
	{
		drawGrid();
	}

	if (_selRow != -1 && _selCol != -1)
	{
		glLineWidth(2.0f);
		glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
		drawMark(_selRow, _selCol);
		glLineWidth(1.0f);
	}

	const std::vector<ScalarAgentItem>& grid = _pScalarAgentGrid->getGrid();
  double minRatio = DBL_MAX, maxRatio = 0;
  if(_drawm1m2 == 2) {
    for(int i=0;i<_DIM;i++)
      for(int j=0;j<_DIM;j++)
        for(int k=0;k<2;k++)
        {
          if(!grid[i*_DIM+j]._pAgent[k]) continue;
          switch(grid[i*_DIM+j]._pAgent[k]->getAgentType()) {
          case MAC:
          {
              const Mac* pMac = static_cast<const Mac*>(grid[i*_DIM+j]._pAgent[k]);
              if(!_macFilter[pMac->getState()][ENBL]) continue;
              char state = (pMac->getNFkB() << NFKB) | (pMac->getStat1() << STAT1) | (pMac->isDeactivated() << DEACT);
              state |= (state == 0) << OTHER;
              state &= (_macFilter[pMac->getState()][NFKB] << NFKB) | (_macFilter[pMac->getState()][STAT1] << STAT1) | (_macFilter[pMac->getState()][DEACT] << DEACT) | (_macFilter[pMac->getState()][OTHER] << OTHER);
              if(!state) continue;
          } break;
          case TGAM: { if(!_drawTreg) continue; } break;
          case TREG: { if(!_drawTgam) continue; } break;
          case TCYT: { if(!_drawTcyt) continue; } break;
          }
          const double stnfr = std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundTNFR1());
          const double sil10r = std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundIL10R());
          minRatio = std::min(minRatio, stnfr / sil10r);
          maxRatio = std::max(maxRatio, stnfr / sil10r);
        }
  }

	for (int i = 0; i < _DIM; i++)
	{
		for (int j = 0; j < _DIM; j++)
		{
			int val = grid[i * _DIM + j]._bitMask;
			if (GET_BIT(val, ScalarAgentGrid::_bitSrc) && _drawSrc &&
					(!_drawSrcMac || GET_BIT(val, ScalarAgentGrid::_bitSrcMac)) &&
					(!_drawSrcTgam || GET_BIT(val, ScalarAgentGrid::_bitSrcTgam)) &&
					(!_drawSrcTcyt || GET_BIT(val, ScalarAgentGrid::_bitSrcTcyt)) &&
					(!_drawSrcTreg || GET_BIT(val, ScalarAgentGrid::_bitSrcTreg)))
			{
				glColor4f(0.8f, 0.8f, 0.8f, _gridAlpha);
				drawQuad(i, j);
			}
			else if (GET_BIT(val, ScalarAgentGrid::_bitCas) && _drawCas)
			{
				glColor4f(1.0f, 1.0f, 1.0f, _gridAlpha);
				drawCross(i, j);
			}
            else if(_drawm1m2)
            {
              for(int k=0;k<2;k++)
              {
                if(grid[i*_DIM+j]._pAgent[k])
                {
                  if(_drawm1m2 == 1){
                      if((grid[i*_DIM+j]._pAgent[k]->getsurfBoundTNFR1() / grid[i*_DIM+j]._pAgent[k]->getsurfBoundIL10R()) > _m1m2thres)
                          glColor4f(0,1,0, _gridAlpha);   //Green if over the threshold
                      else
                          glColor4f(1,0,0, _gridAlpha);    //Red otherwise
                  }
                  else {
                      if(minRatio != maxRatio) {
                          const double stnfr =  std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundTNFR1());
                          const double sil10r = std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundIL10R());
                          const double ratio = ( (stnfr / sil10r) - minRatio) / (maxRatio - minRatio);
                          glColor4f(1-ratio, ratio, 0, _gridAlpha);
                      }
                      else
                          glColor4f(1, 0, 0, _gridAlpha);
                  }
                  switch(grid[i*_DIM+j]._pAgent[k]->getAgentType()) {
                  case MAC:
                  {
                      const Mac* pMac = static_cast<const Mac*>(grid[i*_DIM+j]._pAgent[k]);
                      if(!_macFilter[pMac->getState()][ENBL]) continue;
                      char state = (pMac->getNFkB() << NFKB) | (pMac->getStat1() << STAT1) | (pMac->isDeactivated() << DEACT);
                      state |= (state == 0) << OTHER;
                      state &= (_macFilter[pMac->getState()][NFKB] << NFKB) | (_macFilter[pMac->getState()][STAT1] << STAT1) | (_macFilter[pMac->getState()][DEACT] << DEACT) | (_macFilter[pMac->getState()][OTHER] << OTHER);
                      if(!state) continue;
                  } break;
                  case TGAM: { if(!_drawTreg) continue; } break;
                  case TREG: { if(!_drawTgam) continue; } break;
                  case TCYT: { if(!_drawTcyt) continue; } break;
                  }
                  drawQuad(i, j);
                }
              }
            }
			else if (GET_BIT(val, ScalarAgentGrid::_bitTgam) && _drawTgam)
			{
				glColor4f(1.0f, 0.71f, 0.76f, _gridAlpha);
				drawQuad(i, j);
			}
			else if (GET_BIT(val, ScalarAgentGrid::_bitTcyt) && _drawTcyt)
			{
				glColor4f(0.5f, 0.0f, 0.5f, _gridAlpha);
				drawQuad(i, j);
			}
			else if (GET_BIT(val, ScalarAgentGrid::_bitTreg) && _drawTreg)
			{
				glColor4f(0.0f, 1.0f, 1.0f, _gridAlpha);
				drawQuad(i, j);
			}
			else if (GET_BIT(val, ScalarAgentGrid::_bitMac))
			{
                const Mac* pMac = NULL;
                for(int k=0;k<2;k++)
                    if(grid[i*_DIM+j]._pAgent[k]) {
                        pMac = static_cast<const Mac*>(grid[i*_DIM+j]._pAgent[k]->getAgentType() == MAC ? grid[i*_DIM+j]._pAgent[k] : pMac);
                    }
                Q_CHECK_PTR(pMac);
                if(!_macFilter[pMac->getState()][ENBL]) continue;
                char state = (pMac->getNFkB() << NFKB) | (pMac->getStat1() << STAT1) | (pMac->isDeactivated() << DEACT);
                state |= (state == 0) << OTHER;
                state &= (_macFilter[pMac->getState()][NFKB] << NFKB) | (_macFilter[pMac->getState()][STAT1] << STAT1) | (_macFilter[pMac->getState()][DEACT] << DEACT) | (_macFilter[pMac->getState()][OTHER] << OTHER);
                if(!state) continue;

                if (GET_BIT(val, ScalarAgentGrid::_bitMacResting) && _macFilter[Mac::MAC_RESTING][ENBL])
                  glColor4f(0.0f, 1.0f, 0.0f, _gridAlpha);
                else if (GET_BIT(val, ScalarAgentGrid::_bitMacInfected) && _macFilter[Mac::MAC_INFECTED][ENBL])
                  glColor4f(1.0f, 0.65f, 0.0f, _gridAlpha);
                else if (GET_BIT(val, ScalarAgentGrid::_bitMacCInfected) && _macFilter[Mac::MAC_CINFECTED][ENBL])
                  glColor4f(1.0f, 0.0f, 0.0f, _gridAlpha);
                else if (GET_BIT(val, ScalarAgentGrid::_bitMacActive) && _macFilter[Mac::MAC_ACTIVE][ENBL])
                  glColor4f(0.0f, 0.0f, 1.0f, _gridAlpha);
                else
                  continue;
                drawQuad(i, j);
			}
			else if (GET_BIT(val, ScalarAgentGrid::_bitExtMtb) && _drawExtMtb)
			{
				glColor4f(0.67f, 0.67f, 0.0f, _gridAlpha);
				drawQuad(i, j);
			}
		}
	}

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
}

void AgentsVisualization::drawGrid() const
{
	glColor4f(0.1f, 0.1f, 0.1f, _gridAlpha);
	glLineWidth(0.2f);

	glBegin(GL_LINES);
	for (int i = 0; i <= _DIM; i++)
	{
		glVertex3f(i * _deltaX, _MIN_Y, _gridHeight);
		glVertex3f(i * _deltaX, _MAX_Y, _gridHeight);
	}
	for (int i = 0; i <= _DIM; i++)
	{
		glVertex3f(_MIN_X, i * _deltaY, _gridHeight);
		glVertex3f(_MAX_X, i * _deltaY, _gridHeight);
	}
	glEnd();

	glLineWidth(1.0f);
}
