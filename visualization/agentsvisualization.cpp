/*
 * agentsvisualization.cpp
 *
 *  Created on: 25-nov-2009
 *      Author: El-Kebir, M.
 */

#include "agentsvisualization.h"
#include "scalardatasets/scalaragentgrid.h"
#include <vector>
#include <assert.h>

float rand_float()
  { return rand() * 2.0f / RAND_MAX - 1.0f; }

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
  , _drawSquares(true)
  , _gridAlpha(0.7f)
  , _gridHeight(0.0f)
  , _drawSrcMac(false)
  , _drawSrcTgam(false)
  , _drawSrcTcyt(false)
  , _drawSrcTreg(false)
  , _drawm1m2(false)
  , _m1m2thres(0.0)
  , _selRow(-1)
  , _selCol(-1)
  , stratify(boost::extents[DIM][DIM][2])
{
  std::generate(stratify.origin(), stratify.origin() + stratify.num_elements(), rand_float);
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
  if(_drawm1m2 == 2)
    {
      for(int i=0; i<_DIM; i++)
        for(int j=0; j<_DIM; j++)
          for(int k=0; k<2; k++)
            {
              if(!grid[i*_DIM+j]._pAgent[k]) continue;
              switch(grid[i*_DIM+j]._pAgent[k]->getAgentType())
                {
                case MAC:
                {
                  const Mac* pMac = static_cast<const Mac*>(grid[i*_DIM+j]._pAgent[k]);
                  if(!_macFilter[pMac->getState()][AgentsWidget::ENBL]) continue;
                  char state = (pMac->getNFkB() << AgentsWidget::NFKB) | (pMac->getStat1() << AgentsWidget::STAT1) | (pMac->isDeactivated() << AgentsWidget::DEACT);
                  state |= (state == 0) << AgentsWidget::OTHER;
                  state &= (_macFilter[pMac->getState()][AgentsWidget::NFKB] << AgentsWidget::NFKB)
                           | (_macFilter[pMac->getState()][AgentsWidget::STAT1] << AgentsWidget::STAT1)
                           | (_macFilter[pMac->getState()][AgentsWidget::DEACT] << AgentsWidget::DEACT)
                           | (_macFilter[pMac->getState()][AgentsWidget::OTHER] << AgentsWidget::OTHER);
                  if(!state) continue;
                }
                break;
                default: continue;
                }
              const double stnfr = std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundTNFR1());
              const double sil10r = std::max(1e-5, grid[i*_DIM+j]._pAgent[k]->getsurfBoundIL10R());
              minRatio = std::min(minRatio, stnfr / sil10r);
              maxRatio = std::max(maxRatio, stnfr / sil10r);
            }
    }

  for (int i = 0; i < _DIM; i++)
    for (int j = 0; j < _DIM; j++)
      {
        int val = grid[i * _DIM + j]._bitMask;
        if (GET_BIT(val, ScalarAgentGrid::_bitCas) && _drawCas)
        {
          glColor4f(1.0f, 1.0f, 1.0f, _gridAlpha);
          drawCross(i, j);
        }
        if (GET_BIT(val, ScalarAgentGrid::_bitSrc) && _drawSrc &&
            (!_drawSrcMac || GET_BIT(val, ScalarAgentGrid::_bitSrcMac)) &&
            (!_drawSrcTgam || GET_BIT(val, ScalarAgentGrid::_bitSrcTgam)) &&
            (!_drawSrcTcyt || GET_BIT(val, ScalarAgentGrid::_bitSrcTcyt)) &&
            (!_drawSrcTreg || GET_BIT(val, ScalarAgentGrid::_bitSrcTreg)))
          {
            glColor4f(0.8f, 0.8f, 0.8f, _gridAlpha);
            drawQuad(i, j);
          }
        if(_drawSquares) {  //Setup a priority for each agents drawn in case of drawing squares.
                 if(_drawTgam && GET_BIT(val, ScalarAgentGrid::_bitTgam)) drawTcell(TGAM, i, j);
            else if(_drawTcyt && GET_BIT(val, ScalarAgentGrid::_bitTcyt)) drawTcell(TCYT, i, j);
            else if(_drawTreg && GET_BIT(val, ScalarAgentGrid::_bitTreg)) drawTcell(TREG, i, j);
            else {
              int k=0;
              for(k=0;k<2;k++) {
                if(!(grid[i*_DIM+j]._pAgent[k] && (grid[i*_DIM+j]._pAgent[k])->getAgentType() == MAC)) continue;
                else if(drawMac(static_cast<const Mac*>(grid[i*_DIM+j]._pAgent[k]), i, j, k, minRatio, maxRatio))
                  k=3;  //Finish the loop and don't try to draw any more.
              }
              if(k>3)
                continue; //Drew *a* mac, skip extmtb later.
            }
        }
        else {
          for(int k=0;k<2;k++)
          {
            const Agent* agent = grid[i*_DIM+j]._pAgent[k];
            if(agent == NULL) continue;
            switch(agent->getAgentType())
            {
              case MAC:
               drawMac(static_cast<const Mac*>(agent), i, j, k, minRatio, maxRatio);
               break;
              default:
               drawTcell(agent->getAgentType(), i, j, k);
               break;
            }
          }
        }
        if (GET_BIT(val, ScalarAgentGrid::_bitExtMtb) && _drawExtMtb)
        {
          glColor4f(0.67f, 0.67f, 0.0f, _gridAlpha);
          drawQuad(i, j);
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
void AgentsVisualization::drawTcell(const AgentType a_t, int x, int y, int p) const
{
  switch(a_t) {
    case TGAM:
    {
      if(!_drawTgam) return;
      glColor4f(1.0f, 0.71f, 0.76f, _gridAlpha);
      break;
    }
    case TREG:
    {
      if(!_drawTreg) return;
      glColor4f(0.0f, 1.0f, 1.0f, _gridAlpha);
      break;
    }
    case TCYT:
    {
      if(!_drawTcyt) return;
      glColor4f(0.5f, 0.0f, 0.5f, _gridAlpha);
      break;
    }
    default:
    {
      qWarning("Not a valid tcell type!"); return;
    }
  }
  if(_drawSquares)
    drawQuad(x, y);
  else
  {
    float rngx = stratify[x][y][p]*0.25f+0.75f;
    float rngy = stratify[x][y][1-p]*0.25f+0.75f;
    drawCell((x+rngx*(2*p-1)*(0.5f-0.3f)), (y+rngy*(2*p-1)*(0.5f-0.3f)), 0.005f, 0.3f);
  }
}
bool AgentsVisualization::drawMac(const Mac* pMac, int x, int y, int p, const double minRatio, const double maxRatio) const
{
  if(!pMac) return false;
  if(pMac->getAgentType() != MAC) return false;
  char state = (pMac->getNFkB() << AgentsWidget::NFKB) | (pMac->getStat1() << AgentsWidget::STAT1) | (pMac->isDeactivated() << AgentsWidget::DEACT);
  state |= (state == 0) << AgentsWidget::OTHER;
  state &= (_macFilter[pMac->getState()][AgentsWidget::NFKB] << AgentsWidget::NFKB)
           | (_macFilter[pMac->getState()][AgentsWidget::STAT1] << AgentsWidget::STAT1)
           | (_macFilter[pMac->getState()][AgentsWidget::DEACT] << AgentsWidget::DEACT)
           | (_macFilter[pMac->getState()][AgentsWidget::OTHER] << AgentsWidget::OTHER);
  if(!(state && _macFilter[pMac->getState()][AgentsWidget::ENBL]))
    return false; //Don't draw macs in a disabled state

  if(_drawm1m2 == 1)
  {
    if(pMac->getM1M2Ratio() > _m1m2thres)
      glColor4f(0,1,0, _gridAlpha);   //Green if over the threshold
    else
      glColor4f(1,0,0, _gridAlpha);    //Red otherwise
  }
  else if(_drawm1m2 != 0)
  {
    if(minRatio != maxRatio)
    {
      const double ratio = (pMac->getM1M2Ratio() - minRatio) / (maxRatio - minRatio);
      glColor4f(1-ratio, ratio, 0, _gridAlpha);
    }
    else
      glColor4f(1, 0, 0, _gridAlpha);
  }
  else {
    switch(pMac->getState())
    {
      case Mac::MAC_RESTING:
        glColor4f(0.0f, 1.0f, 0.0f, _gridAlpha); break;
      case Mac::MAC_INFECTED:
        glColor4f(1.0f, 0.65f, 0.0f, _gridAlpha); break;
      case Mac::MAC_CINFECTED:
        glColor4f(1.0f, 0.0f, 0.0f, _gridAlpha); break;
      case Mac::MAC_ACTIVE:
        glColor4f(0.0f, 0.0f, 1.0f, _gridAlpha); break;
      default:
        qWarning("Unknown state '%d' for mac", pMac->getState()); return false;
    }
  }
  if(_drawSquares)
    drawQuad(x, y);
  else
  {
    float rngx = stratify[x][y][p]*0.25f+0.75f;
    float rngy = stratify[x][y][1-p]*0.25f+0.75f;
    drawCell((x+rngx*(2*p-1)*(0.5f-0.375f)), (y+rngy*(2*p-1)*(0.5f-0.375f)), 0.00f, 0.375f);
  }
  return true;
}
