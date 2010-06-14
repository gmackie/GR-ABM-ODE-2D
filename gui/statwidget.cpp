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

	str = QString("(%1,%2)").arg(stats.getNrApoptosisTNF()).arg(stats.getNrApoptosisFasFasL());
	_ui.labelApoptosis->setText(str);

	int nSrcMac = stats.getNrSourcesMac();
	int nSrcTgam = stats.getNrSourcesTgam();
	int nSrcTcyt = stats.getNrSourcesTcyt();
	int nSrcTreg = stats.getNrSourcesTreg();

	str = QString("(%1,%2,%3,%4)").arg(nSrcMac).arg(nSrcTgam).arg(nSrcTcyt).arg(nSrcTreg);
	_ui.labelSources->setText(str);

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
	str = QString("%1").arg(stats.getArea());
	_ui.labelGranulomaAreaVal->setText(str);
}
