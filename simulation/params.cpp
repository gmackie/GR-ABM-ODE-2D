#include "params.h"
Params::Params(const Pos& dim) :
#define P(type, n, p, defau, d, u, min, max)  \
  _##p##_##n(type(defau)), \
  _##p##_##n##Desc(type(min), type(max), #n, boost::optional<type>(defau), std::string(d), std::string(u), make_path(#p)),
#include "params.def"
  _dim(dim)
{

}

// Define any computed parameters including any that are computed based on other parameters.
void Params::computeParams()
{

  if(getMac_initDensityDesc().wasRead)
  {
    setMac_initNumber((int)(getMac_initDensity() * double(_dim.x*_dim.y)));
  }
  else if(getMac_initNumberDesc().wasRead)
  {
    setMac_initDensity((int)(getMac_initNumber() / double(_dim.x*_dim.y)));
  }
  else
    throw std::runtime_error("Initial resting macs not specified in parameter file");

  if(get_sourceDensityDesc().wasRead)
  {
    set_nrSources((int)(get_sourceDensity() * double(_dim.x*_dim.y)));
  }
  else if(get_nrSourcesDesc().wasRead)
  {
    set_sourceDensity(get_nrSources() / double(_dim.x*_dim.y));
  }
  else
    throw std::runtime_error("Sources not specified in parameter file");                                       

  if (!get_effectRecCCL2Desc().wasRead)
  {
    set_effectRecCCL2(getMac_dTNF() / getMac_dCCL2());
  }

  if (!get_effectRecCCL5Desc().wasRead)
  {
    set_effectRecCCL5(getMac_dTNF() / getMac_dCCL5());
  }

  if (!get_effectRecCXCL9Desc().wasRead)
  {
    set_effectRecCXCL9(getMac_dTNF() / getMac_dCXCL9());
  }

  //Scalar fluxFactorTimeDelta = getParam(PARAM_TCELL_LYMPH_PROXY_NONLINEAR_TIME_FULL) -  getParam(PARAM_TCELL_LYMPH_PROXY_NONLINEAR_TIME_START);

  //Scalar fluxFactorM = 1.0/fluxFactorTimeDelta;                                                              
  //setParam(PARAM_TCELL_LYMPH_PROXY_NONLINEAR_M, fluxFactorM);                                                
  
  //Scalar fluxFactorB = 1.0 - (getParam(PARAM_TCELL_LYMPH_PROXY_NONLINEAR_TIME_FULL)/fluxFactorTimeDelta);    
  //setParam(PARAM_TCELL_LYMPH_PROXY_NONLINEAR_B, fluxFactorB);                                                
}

