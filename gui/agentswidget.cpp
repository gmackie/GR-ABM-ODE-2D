#include "agentswidget.h"

AgentsWidget::AgentsWidget(AgentsVisualization* pAgentsVisualization, QWidget *parent)
    : QWidget(parent)
    , _pAgentsVisualization(pAgentsVisualization)
{
	assert(pAgentsVisualization);
    _ui.setupUi(this);

    connect(_ui.checkBoxDrawAgentMa, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMci, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMi, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMr, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentTgam, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentTcyt, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentTreg, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentCaseation, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentSources, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentExtMtb, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxAgentDrawGrid, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.horizontalSliderAgentsAlpha, SIGNAL(valueChanged(int)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.doubleSpinBoxAgentsGridHeight, SIGNAL(valueChanged(double)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.doubleSpinBoxAgentsGridHeight, SIGNAL(valueChanged(double)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentSourcesMac, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentSourcesTgam, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentSourcesTcyt, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentSourcesTreg, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.comboBoxM1M2, SIGNAL(activated(int)), this, SLOT(updateM1M2Settings(void)));
    connect(_ui.doubleSpinBoxM1M2Threshold, SIGNAL(valueChanged(double)), this, SLOT(updateAgentsSettings(void)));

    connect(_ui.checkBoxMaDeact, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMaStat1, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMaNFkB, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMaOther, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMa, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));

    connect(_ui.checkBoxMrDeact, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMrStat1, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMrNFkB, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMr, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMrOther, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));

    connect(_ui.checkBoxMiDeact, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMiStat1, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMiNFkB, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMi, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMiOther, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));

    connect(_ui.checkBoxMciDeact, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMciStat1, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMciNFkB, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentMci, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxMciOther, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));

    updateAgentsSettings();
}

AgentsWidget::~AgentsWidget()
{
}

void AgentsWidget::setAgentSelection(int row, int col)
{
	_pAgentsVisualization->setSelection(row, col);
}

void AgentsWidget::updateM1M2Settings()
{
    _ui.checkBoxDrawAgentTcyt->setChecked(false);
    _ui.checkBoxDrawAgentTgam->setChecked(false);
    _ui.checkBoxDrawAgentTreg->setChecked(false);
    _ui.checkBoxDrawAgentExtMtb->setChecked(false);
    _ui.doubleSpinBoxM1M2Threshold->setEnabled(_ui.comboBoxM1M2->currentIndex() == 1);
    updateAgentsSettings();
}

void AgentsWidget::updateAgentsSettings()
{
    /// Temporary
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, AgentsVisualization::DEACT, _ui.checkBoxMaDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, AgentsVisualization::STAT1, _ui.checkBoxMaStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, AgentsVisualization::NFKB, _ui.checkBoxMaNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, AgentsVisualization::OTHER, _ui.checkBoxMaOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, AgentsVisualization::ENBL, _ui.checkBoxDrawAgentMa->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, AgentsVisualization::DEACT, _ui.checkBoxMrDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, AgentsVisualization::STAT1, _ui.checkBoxMrStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, AgentsVisualization::NFKB, _ui.checkBoxMrNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, AgentsVisualization::OTHER, _ui.checkBoxMrOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, AgentsVisualization::ENBL, _ui.checkBoxDrawAgentMr->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, AgentsVisualization::DEACT, _ui.checkBoxMiDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, AgentsVisualization::STAT1, _ui.checkBoxMiStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, AgentsVisualization::NFKB, _ui.checkBoxMiNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, AgentsVisualization::OTHER, _ui.checkBoxMiOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, AgentsVisualization::ENBL, _ui.checkBoxDrawAgentMi->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, AgentsVisualization::DEACT, _ui.checkBoxMciDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, AgentsVisualization::STAT1, _ui.checkBoxMciStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, AgentsVisualization::NFKB, _ui.checkBoxMciNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, AgentsVisualization::OTHER, _ui.checkBoxMciOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, AgentsVisualization::ENBL, _ui.checkBoxDrawAgentMci->isChecked());

    _pAgentsVisualization->setGridAlpha(_ui.horizontalSliderAgentsAlpha->value() / _SLIDER_AGENTS_ALPHA);
    _pAgentsVisualization->setDrawGrid(_ui.checkBoxAgentDrawGrid->isChecked());
    _pAgentsVisualization->setDrawTgam(_ui.checkBoxDrawAgentTgam->isChecked());
    _pAgentsVisualization->setDrawTreg(_ui.checkBoxDrawAgentTreg->isChecked());
    _pAgentsVisualization->setDrawTcyt(_ui.checkBoxDrawAgentTcyt->isChecked());
    _pAgentsVisualization->setDrawCas(_ui.checkBoxDrawAgentCaseation->isChecked());
    _pAgentsVisualization->setDrawSrc(_ui.checkBoxDrawAgentSources->isChecked());
    _pAgentsVisualization->setDrawExtMtb(_ui.checkBoxDrawAgentExtMtb->isChecked());
    _pAgentsVisualization->setGridHeight((float) _ui.doubleSpinBoxAgentsGridHeight->value());
    _pAgentsVisualization->setPredicates(_ui.checkBoxDrawAgentSourcesMac->isChecked(),
                                         _ui.checkBoxDrawAgentSourcesTgam->isChecked(),
                                         _ui.checkBoxDrawAgentSourcesTcyt->isChecked(),
                                         _ui.checkBoxDrawAgentSourcesTreg->isChecked());
    _pAgentsVisualization->setDrawM1M2(_ui.comboBoxM1M2->currentIndex(),
                                       _ui.doubleSpinBoxM1M2Threshold->value());

	emit updateGL();
}
