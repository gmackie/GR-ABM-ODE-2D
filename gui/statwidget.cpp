#include "statwidget.h"

StatWidget::StatWidget(QWidget* pParent)
    : QWidget(pParent)
{
	_ui.setupUi(this);
}

StatWidget::~StatWidget()
{
}

void StatWidget::updateLabels(const GrStat& stats)
{
	QString str = QString("%1 - (%2,%3,%4,%5,%6)").arg(stats.getNrOfMac()).
			arg(stats.getNrOfMacResting()).arg(stats.getNrOfMacInfected()).
			arg(stats.getNrOfMacCInfected()).arg(stats.getNrOfMacActive()).
			arg(stats.getNrOfMacDead());
	_ui.labelMac->setText(str);

	str = QString("%1 - (%2,%3,%4)").arg(stats.getNrOfTgam()).
			arg(stats.getNrOfTgamActive()).arg(stats.getNrOfTgamDownRegulated()).
			arg(stats.getNrOfTgamDead());
	_ui.labelTgam->setText(str);

	str = QString("%1 - (%2,%3,%4)").arg(stats.getNrOfTcyt()).
			arg(stats.getNrOfTcytActive()).arg(stats.getNrOfTcytDownRegulated()).
			arg(stats.getNrOfTcytDead());
	_ui.labelTcyt->setText(str);

	str = QString("%1 - (%2,%3)").arg(stats.getNrOfTreg()).
			arg(stats.getNrOfTregActive()).arg(stats.getNrOfTregDead());
	_ui.labelTreg->setText(str);

	str = QString("(%1,%2)").arg(stats.getTotExtMtb(), 0, 'f', 2).
			arg(stats.getTotIntMtb(), 0, 'f', 2);
	_ui.labelMtb->setText(str);

	str = QString("(%1,%2)").arg(stats.getNrCaseated()).
			arg(stats.getTotNonRepExtMtb(), 0, 'f', 2);
	_ui.labelMtbNonRep->setText(str);

	str = QString("%1 - (%2,%3,%4,%5,%6)").arg(stats.getNrOfMacNFkB()).
			arg(stats.getNrOfMacNFkBResting()).arg(stats.getNrOfMacNFkBInfected()).
			arg(stats.getNrOfMacNFkBCInfected()).arg(stats.getNrOfMacNFkBActive()).
			arg(stats.getNrOfMacNFkBDead());
	_ui.labelMacNFkB->setText(str);

	str = QString("%1 - (%2,%3,%4,%5,%6)").arg(stats.getNrOfMacStat1()).
			arg(stats.getNrOfMacStat1Resting()).arg(stats.getNrOfMacStat1Infected()).
			arg(stats.getNrOfMacStat1CInfected()).arg(stats.getNrOfMacStat1Active()).
			arg(stats.getNrOfMacStat1Dead());
	_ui.labelMacStat1->setText(str);

	str = QString("%1 - (%2,%3,%4,%5,%6)").arg(stats.getNrOfMacDeact()).
			arg(stats.getNrOfMacDeactResting()).arg(stats.getNrOfMacDeactInfected()).
			arg(stats.getNrOfMacDeactCInfected()).arg(stats.getNrOfMacDeactActive()).
			arg(stats.getNrOfMacDeactDead());
	_ui.labelMacDeact->setText(str);

	str = QString("(%1,%2)").arg(stats.getNrMacApoptosisTNF()).arg(stats.getNrApoptosisFasFasL());
	_ui.labelApoptosis->setText(str);

	str = QString("(%1,%2,%3,%4)").
		arg(stats.getNrSourcesMac()).
		arg(stats.getNrSourcesTgam()).
		arg(stats.getNrSourcesTcyt()).
		arg(stats.getNrSourcesTreg());
	_ui.labelSources->setText(str);

	str = QString("(%1,%2,%3,%4)").
		arg(stats.getNrSourcesActiveMac()).
		arg(stats.getNrSourcesActiveTgam()).
		arg(stats.getNrSourcesActiveTcyt()).
		arg(stats.getNrSourcesActiveTreg());
	_ui.labelSourcesAct->setText(str);

	str = QString("%1").arg(stats.getTotMacAttractant(), 0, 'f', 2);
	_ui.labelMacAttractant->setText(str);

	str = QString("%1").arg(stats.getTotTNF(), 0, 'f', 2);
	_ui.labelTNF->setText(str);

	str = QString("%1").arg(stats.getTotCCL2(), 0, 'f', 2);
	_ui.labelCCL2->setText(str);

	str = QString("%1").arg(stats.getTotCCL5(), 0, 'f', 2);
	_ui.labelCCL5->setText(str);

	str = QString("%1").arg(stats.getTotCXCL9(), 0, 'f', 2);
	_ui.labelCXCL9->setText(str);

	// update granuloma area
	str = QString("%1").arg(stats.getAreaTNF());
	_ui.labelGranulomaAreaVal->setText(str);

	// update ODE stats
	str = QString("%1").arg(stats.getMDC());
	_ui.labeldMDC->setText(str);

	str = QString("(%1,%2,%3)").arg(stats.getNrTgamQueued()).
		arg(stats.getNrTcytQueued()).arg(stats.getNrTregQueued());
	_ui.labelTcellQueue->setText(str);

	str = QString("%1").arg(stats.getFluxTgam());
	_ui.labelTgamFlux->setText(str);

	str = QString("%1").arg(stats.getFluxTcyt());
	_ui.labelTcytFlux->setText(str);

	str = QString("%1").arg(stats.getFluxTreg());
	_ui.labelTregFlux->setText(str);
}
