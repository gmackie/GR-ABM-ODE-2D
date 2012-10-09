/*
 * RecruitmentLnODEPure.cpp
 *
 *  Created on: Jul 4, 2010
 *      Author: mohammed
 */

#include "recruitmentlnodepure.h"
#include "stat.h"

BOOST_CLASS_EXPORT_IMPLEMENT(RecruitmentLnODEPure)

RecruitmentLnODEPure::RecruitmentLnODEPure()
  : RecruitmentLnODE("", "")
{
}

RecruitmentLnODEPure::~RecruitmentLnODEPure()
{
}

void RecruitmentLnODEPure::solveODE(const int, const Stats&, Stats& stats)
{
  const double dt = 1.0/14400.0;

  static const double muMDC_LN= _PARAM(PARAM_muMDC_LN);
  static const double init_n4 = _PARAM(PARAM_initN4  );
  static const double init_n8 = _PARAM(PARAM_initN8  );
  static const double muN4    = _PARAM(PARAM_muN4    );
  static const double k13     = _PARAM(PARAM_k13     );
  static const double hs13    = _PARAM(PARAM_hs13    );
  static const double k14     = _PARAM(PARAM_k14     );
  static const double k15     = _PARAM(PARAM_k15     );
  static const double rho2    = _PARAM(PARAM_rho2    );
  static const double k20a    = _PARAM(PARAM_k20a    );
  static const double hs20a   = _PARAM(PARAM_hs20a   );
  static const double csi1    = _PARAM(PARAM_csi1    );
  static const double csi1a   = _PARAM(PARAM_csi1a   );
  static const double muN8    = _PARAM(PARAM_muN8    );
  static const double wT80    = _PARAM(PARAM_wT80    );
  static const double k16     = _PARAM(PARAM_k16     );
  static const double hs16    = _PARAM(PARAM_hs16    );
  static const double k17     = _PARAM(PARAM_k17     );
  static const double hs17    = _PARAM(PARAM_hs17    );
  static const double k18     = _PARAM(PARAM_k18     );
  static const double rho3    = _PARAM(PARAM_rho3    );
  static const double k24a    = _PARAM(PARAM_k24a    );
  static const double hs24a   = _PARAM(PARAM_hs24a   );
  static const double csi2    = _PARAM(PARAM_csi2    );
  static const double csi2a   = _PARAM(PARAM_csi2a   );
  static const double csi2b   = _PARAM(PARAM_csi2b   );
  static const double scaling = _PARAM(PARAM_scaling );
  static const double m       = _PARAM(PARAM_m       );
  static const double sn4     = init_n4 * muN4;
  static const double sn8     = init_n8 * muN8;

  double& MDC = _odeInitialConditions[_idxMDC];
  double& N4 = _odeInitialConditions[_idxNaiveCD4];
  double& TH0 = _odeInitialConditions[_idxTH0];
  double& TH1 = _odeInitialConditions[_idxTH1];
  double& N8 = _odeInitialConditions[_idxNaiveCD8];
  double& T80 = _odeInitialConditions[_idxT80];
  double& T8 = _odeInitialConditions[_idxT8];
  double& TC = _odeInitialConditions[_idxPrecursorCTL];
  double& TH0lung = _odeInitialConditions[_idxEffectorTH0];
  double& TH1lung = _odeInitialConditions[_idxEffectorTH1];
  double& T80lung = _odeInitialConditions[_idxEffectorT80];
  double& T8lung = _odeInitialConditions[_idxEffectorT8];
  double& TClung = _odeInitialConditions[_idxCTL];

  for (int i = 0; i < 100; i++)
    {
      // y(1), MDC LN [MDC] - Mature/Licensed Dendritic Cell (in the Lymph Node)
      double dMDC = - muMDC_LN*MDC*dt;
      // y(2), Naive CD4+ T cells in the Lymph Node  [N4]
      double dN4 = (sn4 + k13*N4*(MDC/(MDC+hs13)) - muN4*N4 - k14*N4*MDC)*dt;
      // y(3), Precursor Th1 LN  [TH0]
      double dTH0 = (k14*N4*MDC + k15*TH0*(1-TH0/rho2) - k20a*TH0*(MDC/(MDC+hs20a)) - csi1*TH0)*dt;
      // y(4), Th1 LN  [TH1]
      double dTH1 = (k20a*TH0*(MDC/(MDC+hs20a)) - csi1a*TH1)*dt;
      // y(5), Naive CD8+ T cells in the Lymph Node  [N8]
      double dN8 = (sn8 + k16*N8*(MDC/(MDC+hs16)) - muN8*N8 - k17*MDC*N8*((TH1+wT80*TH0)/(TH1+wT80*TH0+hs17)))*dt;
      // y(6), CTL precursor cells in the Lymph Node  T80LN  [T80]
      double dT80 = (k17*MDC*N8*((TH1+wT80*TH0)/(TH1+wT80*TH0+hs17)) + k18*T80*(1-T80/rho3) - k24a*T80*(MDC/(MDC+hs24a)) - csi2*T80)*dt;
      // y(7), T8 cells in the Lymph Node,  [T8]
      double dT8 = (m*k24a*T80*(MDC/(MDC+hs24a)) - csi2a*T8)*dt;
      // CTL cells in the Lymph Node, [TC]
      double dTC = (m*k24a*T80*(MDC/(MDC+hs24a)) - csi2b*TC)*dt;
      // y(9), precursor TH1 migration from the LN [TH0lung]
      double dTH0lung=scaling*csi1*TH0*dt;
      // y(10), TH1 migration from the LN [TH1lung]
      double dTH1lung=scaling*csi1a*TH1*dt;
      // y(11), precursor T8 migration from the LN [T80lung]
      double dT80lung = scaling*csi2*T80*dt;
      //  y(12), T8 migration from the LN [T8lung]
      double dT8lung = scaling*csi2a*T8*dt;
      // y(13), CTL migration from the LN[TClung]
      double dTClung = scaling*csi2b*TC*dt;

      MDC+=dMDC; // MDC - LN [1]
      N4+=dN4; // Naive CD4 LN (N4) - [2]
      TH0+=dTH0; // precursor TH1 LN - [3]
      TH1+=dTH1; // Th1 LN LN [4]
      N8+=dN8; // Naive CD8 LN (N8) - [5]
      T80+=dT80;// T80LN - [6]
      T8+=dT8; // T8 LN [7]
      TC+=dTC; // CTL LN [8]
      TH0lung+=dTH0lung; // precursor TH1 migration from the LN [9]
      TH1lung+=dTH1lung; // TH1 migration from the LN [10]
      T80lung+=dT80lung; // precursor T8 migration from the LN [11]
      T8lung+=dT8lung; // T8 migration from the LN [12]
      TClung+=dTClung; // CTL migration from the LN[13]
    }

  stats.setMDC(MDC);
  stats.setN4(N4);
  stats.setTH0(TH0);
  stats.setTH1(TH1);
  stats.setN8(N8);
  stats.setT80(T80);
  stats.setT8(T8);
  stats.setTC(TC);
  stats.setTH0lung(TH0lung);
  stats.setTH1lung(TH1lung);
  stats.setT80lung(T80lung);
  stats.setT8lung(T8lung);
  stats.setTClung(TClung);


  /*std::cout << "MDC\t" << _odeInitialConditions[_idxMDC] << std::endl;
  std::cout << "N4\t" << _odeInitialConditions[_idxNaiveCD4] << std::endl;
  std::cout << "TH0\t" << _odeInitialConditions[_idxTH0] << std::endl;
  std::cout << "TH1\t" << _odeInitialConditions[_idxTH1] << std::endl;
  std::cout << "N8\t" << _odeInitialConditions[_idxNaiveCD8] << std::endl;
  std::cout << "T80\t" << _odeInitialConditions[_idxT80] << std::endl;
  std::cout << "T8\t" << _odeInitialConditions[_idxT8] << std::endl;
  std::cout << "TC\t" << _odeInitialConditions[_idxPrecursorCTL] << std::endl;
  std::cout << "TH0l\t" << _odeInitialConditions[_idxEffectorTH0] << std::endl;
  std::cout << "TH1l\t" << _odeInitialConditions[_idxEffectorTH1] << std::endl;
  std::cout << "T80l\t" << _odeInitialConditions[_idxEffectorT80] << std::endl;
  std::cout << "T8l\t" << _odeInitialConditions[_idxEffectorT8] << std::endl;
  std::cout << "TCl\t" << _odeInitialConditions[_idxCTL] << std::endl;*/
}
