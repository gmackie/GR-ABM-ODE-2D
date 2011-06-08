#include "paramwindow.h"
#include <QFileDialog>
#include <QMessageBox>

ParamWindow::ParamWindow(MainInterface* pItfc, QWidget* parent)
    : QWidget(parent)
    , _pItfc(pItfc)
    , _pGR(NULL)
	, _pMac(NULL)
	, _pTcell(NULL)
	, _pTgam(NULL)
	, _pTcyt(NULL)
	, _pTreg(NULL)
	, _pMtb(NULL)
{
	_ui.setupUi(this);
	_ui.treeWidget->setColumnWidth(0, 200);
	_ui.treeWidget->setColumnWidth(1, 50);
	_ui.treeWidget->setColumnWidth(2, 100);
	update();

	connect(_ui.pushButtonLoad, SIGNAL(clicked(bool)), this, SLOT(loadParams(void)));
	connect(_ui.pushButtonSave, SIGNAL(clicked(bool)), this, SLOT(saveParams(void)));
}

ParamWindow::~ParamWindow()
{
}

void ParamWindow::newItem(QTreeWidgetItem* pParentItem, ParamDoubleType param)
{
	const Params* params = Params::getInstance();

	QTreeWidgetItem* pItem = new QTreeWidgetItem(pParentItem, param);
	pItem->setText(0, QString("%1").arg(params->getName(param).c_str()));
	pItem->setText(1, QString("%1").arg(params->getParam(param)));
	pItem->setText(2, QString("%1").arg(params->getUnit(param).c_str()));
	pItem->setText(3, QString("%1").arg(params->getDescription(param).c_str()));
}

void ParamWindow::newItem(QTreeWidgetItem* pParentItem, ParamIntType param)
{
	const Params* params = Params::getInstance();

	QTreeWidgetItem* pItem = new QTreeWidgetItem(pParentItem);
	pItem->setText(0, QString("%1").arg(params->getName(param).c_str()));
	pItem->setText(1, QString("%1").arg(params->getParam(param)));
	pItem->setText(2, QString("%1").arg(params->getUnit(param).c_str()));
	pItem->setText(3, QString("%1").arg(params->getDescription(param).c_str()));
}

void ParamWindow::update()
{
	QTreeWidget* pTreeWidget = _ui.treeWidget;
	pTreeWidget->clear();

	_pGR = new QTreeWidgetItem(pTreeWidget);
	_pGR->setExpanded(true);
	_pGR->setText(0, "GR");

	newItem(_pGR, PARAM_GR_NR_SOURCES);
	newItem(_pGR, PARAM_GR_NR_KILLINGS_FOR_CASEATION);
	newItem(_pGR, PARAM_GR_D_TNF);
	newItem(_pGR, PARAM_GR_D_CHEMOKINES);
	newItem(_pGR, PARAM_GR_DEG_TNF);
	newItem(_pGR, PARAM_GR_DEG_CHEMOKINES);
	newItem(_pGR, PARAM_GR_THRESHOLD_APOPTOSIS_TNF);
	newItem(_pGR, PARAM_GR_K_APOPTOSIS);
	newItem(_pGR, PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR);
	newItem(_pGR, PARAM_GR_K_APOPTOSIS_MOLECULAR);
	newItem(_pGR, PARAM_GR_PROB_APOPTOSIS_TNF);
	newItem(_pGR, PARAM_GR_MIN_CHEMOTAXIS);
	newItem(_pGR, PARAM_GR_MAX_CHEMOTAXIS);
	newItem(_pGR, PARAM_GR_SEC_RATE_ATTRACTANT);

	newItem(_pGR, PARAM_muMDC_LN);
	newItem(_pGR, PARAM_initN4);
	newItem(_pGR, PARAM_initN8);
	newItem(_pGR, PARAM_muN4);
	newItem(_pGR, PARAM_scaling_LN);
	newItem(_pGR, PARAM_scaling_LUNG);
	newItem(_pGR, PARAM_k13);
	newItem(_pGR, PARAM_hs13);
	newItem(_pGR, PARAM_k14);
	newItem(_pGR, PARAM_k15);
	newItem(_pGR, PARAM_rho2);
	newItem(_pGR, PARAM_k20a);
	newItem(_pGR, PARAM_hs20a);
	newItem(_pGR, PARAM_csi1);
	newItem(_pGR, PARAM_csi1a);
	newItem(_pGR, PARAM_muN8);
	newItem(_pGR, PARAM_wT80);
	newItem(_pGR, PARAM_k16);
	newItem(_pGR, PARAM_hs16);
	newItem(_pGR, PARAM_k17);
	newItem(_pGR, PARAM_hs17);
	newItem(_pGR, PARAM_k18);
	newItem(_pGR, PARAM_rho3);
	newItem(_pGR, PARAM_k24a);
	newItem(_pGR, PARAM_hs24a);
	newItem(_pGR, PARAM_csi2);
	newItem(_pGR, PARAM_csi2a);
	newItem(_pGR, PARAM_csi2b);
	newItem(_pGR, PARAM_scaling);
	newItem(_pGR, PARAM_m);
	newItem(_pGR, PARAM_scaling_MDC);

	if (Params::getInstance()->getUseRecruitmentWeights())
	{
		newItem(_pGR, PARAM_GR_WEIGHT_TNF_RECRUITMENT);
		newItem(_pGR, PARAM_GR_WEIGHT_CCL2_RECRUITMENT);
		newItem(_pGR, PARAM_GR_WEIGHT_CCL5_RECRUITMENT);
		newItem(_pGR, PARAM_GR_WEIGHT_CXCL9_RECRUITMENT);
	}

	_pMac = new QTreeWidgetItem(_pGR);
	_pMac->setExpanded(true);
	_pMac->setText(0, "Mac");
	_pMac->setText(3, "Macrophage specific parameters");

	newItem(_pMac, PARAM_MAC_AGE);
	newItem(_pMac, PARAM_MAC_A_AGE);
	newItem(_pMac, PARAM_MAC_INIT_NUMBER);
	newItem(_pMac, PARAM_MAC_TIMESPAN_REGULATED);
	newItem(_pMac, PARAM_MAC_MOVEMENT_RESTING);
	newItem(_pMac, PARAM_MAC_MOVEMENT_ACTIVE);
	newItem(_pMac, PARAM_MAC_MOVEMENT_INFECTED);
	newItem(_pMac, PARAM_MAC_SEC_RATE_TNF);
	newItem(_pMac, PARAM_MAC_SEC_RATE_CCL2);
	newItem(_pMac, PARAM_MAC_SEC_RATE_CCL5);
	newItem(_pMac, PARAM_MAC_SEC_RATE_CXCL9);
	newItem(_pMac, PARAM_MAC_THRESHOLD_NFKB_TNF);
	newItem(_pMac, PARAM_MAC_K_NFKB);
	newItem(_pMac, PARAM_MAC_THRESHOLD_NFKB_TNF_MOLECULAR);
	newItem(_pMac, PARAM_MAC_K_NFKB_MOLECULAR);
	newItem(_pMac, PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB);
	newItem(_pMac, PARAM_MAC_THRESHOLD_BURST_CI_INTMTB);
	newItem(_pMac, PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
	newItem(_pMac, PARAM_MAC_NR_UPTAKE_A_EXTMTB);
	newItem(_pMac, PARAM_MAC_NR_UPTAKE_RI_EXTMTB);
	newItem(_pMac, PARAM_MAC_PROB_KILL_R_EXTMTB);
	newItem(_pMac, PARAM_MAC_PROB_STAT1_TGAM);
	newItem(_pMac, PARAM_MAC_PROB_RECRUITMENT);
	newItem(_pMac, PARAM_MAC_THRESHOLD_RECRUITMENT);

	_pTcell = new QTreeWidgetItem(_pGR);
	_pTcell->setExpanded(true);
	_pTcell->setText(0, "Tcell");
	_pTcell->setText(3, "Tcell specific parameters");

	newItem(_pTcell, PARAM_TCELL_AGE);
	newItem(_pTcell, PARAM_TCELL_TIME_RECRUITMENT_ENABLED);
	newItem(_pTcell, PARAM_TCELL_PROB_MOVE_TO_MAC);
	newItem(_pTcell, PARAM_TCELL_PROB_MOVE_TO_TCELL);
	newItem(_pTcell, PARAM_TCELL_PROB_RECRUITMENT);

	_pTgam = new QTreeWidgetItem(_pTcell);
	_pTgam->setExpanded(true);
	_pTgam->setText(0, "Tgam");
	_pTgam->setText(3, "Tgam specific parameters");

	newItem(_pTgam, PARAM_TGAM_TIMESPAN_REGULATED);
	newItem(_pTgam, PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL);
	newItem(_pTgam, PARAM_TGAM_PROB_RECRUITMENT);
	newItem(_pTgam, PARAM_TGAM_THRESHOLD_RECRUITMENT);

	_pTcyt = new QTreeWidgetItem(_pTcell);
	_pTcyt->setExpanded(true);
	_pTcyt ->setText(0, "Tcyt");
	_pTcyt->setText(3, "Tcyt specific parameters");

	newItem(_pTcyt, PARAM_TCYT_TIMESPAN_REGULATED);
	newItem(_pTcyt, PARAM_TCYT_PROB_KILL_MAC);
	newItem(_pTcyt, PARAM_TCYT_PROB_KILL_MAC_CLEANLY);
	newItem(_pTcyt, PARAM_TCYT_PROB_RECRUITMENT);
	newItem(_pTcyt, PARAM_TCYT_THRESHOLD_RECRUITMENT);

	_pTreg = new QTreeWidgetItem(_pTcell);
	_pTreg->setExpanded(true);
	_pTreg->setText(0, "Treg");
	_pTreg->setText(3, "Treg specific parameters");

	newItem(_pTreg, PARAM_TREG_PROB_RECRUITMENT);
	newItem(_pTreg, PARAM_TREG_THRESHOLD_RECRUITMENT);

	_pMtb = new QTreeWidgetItem(_pGR);
	_pMtb->setExpanded(true);
	_pMtb->setText(0, "Mtb");
	_pMtb->setText(3, "Mtb specific parameters");

	newItem(_pMtb, PARAM_INTMTB_GROWTH_RATE);
	newItem(_pMtb, PARAM_EXTMTB_GROWTH_RATE);
	newItem(_pMtb, PARAM_EXTMTB_UPPER_BOUND);

	_ui.treeWidget->insertTopLevelItem(0, _pGR);
}

void ParamWindow::loadParams()
{
	Simulation& sim = _pItfc->getSimulation();

	QString fileName = QFileDialog::getOpenFileName(this, "Load parameters", "", "*.xml");
	if (fileName != QString::null)
	{
		sim.lock();
		if (!Params::reinit(fileName.toLatin1().data()))
		{
			QMessageBox::critical(this, "Lung ABM",
					"Failed to open file '" + fileName + "'.",
					QMessageBox::Ok, QMessageBox::Ok);
		}
		sim.unlock();
	}
}

void ParamWindow::saveParams()
{
	QString fileName = QFileDialog::getSaveFileName(this, "Save parameters", "", "*.xml");
	if (fileName != QString::null)
	{
		if (!Params::getInstance()->toXml(fileName.toLatin1().data()))
		{
			QMessageBox::critical(this, "Lung ABM",
					"Failed to save file '" + fileName + "'.",
					QMessageBox::Ok, QMessageBox::Ok);
		}
	}
}
