#include "statwidget.h"

StatWidget::StatWidget(QWidget* pParent)
    : QWidget(pParent)
{
	_ui.setupUi(this);
}

StatWidget::~StatWidget()
{
}

static QString toString(const Scalar& v) {
  return QString::number(v, 'f', 2);
}
template<typename T>
static QString toString(const T& v) {
  return QString("%1").arg(v);
}
template<typename InputIterator>
static void toString(InputIterator start, InputIterator end, QString& s) {
  s = QString('(');
  while(start != end)
    s += toString(*(start++)) + ',';
  s[s.size()-1] = ')';
}
void StatWidget::updateLabels(const Stats& stats)
{
  QString str;
  toString(stats.getNrOfMacsArray().begin(), stats.getNrOfMacsArray().end(), str);
	str = QString("%1 - %2").arg(stats.getNrOfMacs()).arg(str);
	_ui.labelMac->setText(str);

  toString(stats.getNrOfTgamsArray().begin(), stats.getNrOfTgamsArray().end(), str);
	str = QString("%1 - %2").arg(stats.getNrOfTgams()).arg(str);
	_ui.labelTgam->setText(str);

  toString(stats.getNrOfTcytsArray().begin(), stats.getNrOfTcytsArray().end(), str);
	str = QString("%1 - %2").arg(stats.getNrOfTcyts()).arg(str);
	_ui.labelTcyt->setText(str);

  toString(stats.getNrOfTregsArray().begin(), stats.getNrOfTregsArray().end(), str);
	str = QString("%1 - %2").arg(stats.getNrOfTregs()).arg(str);
	_ui.labelTreg->setText(str);

	str = QString("(%1,%2)").arg(stats.getTotExtMtb(), 0, 'f', 2).
			arg(stats.getTotIntMtb(), 0, 'f', 2);
	_ui.labelMtb->setText(str);

	str = QString("(%1,%2)").arg(stats.getNrCaseated()).
			arg(stats.getTotNonRepExtMtb(), 0, 'f', 2);
	_ui.labelMtbNonRep->setText(str);

  toString(stats.getMacNFkBArray().begin(), stats.getMacNFkBArray().end(), str);
	str = QString("%1 - %2").arg(stats.getMacNFkB()).arg(str);
	_ui.labelMacNFkB->setText(str);

  toString(stats.getMacStat1Array().begin(), stats.getMacStat1Array().end(), str);
	str = QString("%1 - %2").arg(stats.getMacStat1()).arg(str);
	_ui.labelMacStat1->setText(str);

  toString(stats.getMacDeactArray().begin(), stats.getMacDeactArray().end(), str);
	str = QString("%1 - %2").arg(stats.getMacDeact()).arg(str);
	_ui.labelMacDeact->setText(str);

	str = QString("(%1,%2)").arg(stats.getMacApoptosisTNF()).arg(stats.getApoptosisFasFasL());
	_ui.labelApoptosis->setText(str);

  toString(stats.getNrSourcesArray().begin(), stats.getNrSourcesArray().end(), str);
	str = QString("%1 - %2").arg(stats.getTotNrSources()).arg(str);
	_ui.labelSources->setText(str);

  toString(stats.getNrSourcesActiveArray().begin(), stats.getNrSourcesActiveArray().end(), str);
	str = QString("%1 - %2").arg(stats.getTotNrSourcesActive()).arg(str);
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
	
	str = QString("%1").arg(stats.getTotIL10(), 0, 'f', 2);
	_ui.labelIL10->setText(str);

	// update granuloma area
	str = QString("%1").arg(stats.getAreaTNF());
	_ui.labelGranulomaAreaVal->setText(str);

	// update ODE stats
	str = QString("%1").arg(stats.getMDC());
	_ui.labeldMDC->setText(str);

	str = QString("(%1,%2,%3)").arg(stats.getNrQueued(TGAM)).
		arg(stats.getNrQueued(TCYT)).arg(stats.getNrQueued(TREG));
	_ui.labelTcellQueue->setText(str);

	str = QString("%1").arg(stats.getFlux(TGAM));
	_ui.labelTgamFlux->setText(str);

	str = QString("%1").arg(stats.getFlux(TCYT));
	_ui.labelTcytFlux->setText(str);

	str = QString("%1").arg(stats.getFlux(TREG));
	_ui.labelTregFlux->setText(str);
}
