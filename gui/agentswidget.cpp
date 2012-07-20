#include "agentswidget.h"
#include "simulation/gr.h"
#include "simulation/macrophage.h"
#include "visualization/agentsvisualization.h"

AgentsWidget::AgentsWidget(AgentsVisualization* pAgentsVisualization, QWidget *parent)
    : QWidget(parent)
    , _pAgentsVisualization(pAgentsVisualization)
{
    Q_CHECK_PTR(pAgentsVisualization);
    _ui.setupUi(this);

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

    QButtonGroup* group = new QButtonGroup(this);
    group->setExclusive(false);
    group->addButton(_ui.checkBoxDrawAgentMa, ENBL);
    group->addButton(_ui.checkBoxMaStat1, STAT1);
    group->addButton(_ui.checkBoxMaNFkB, NFKB);
    group->addButton(_ui.checkBoxMaDeact, DEACT);
    group->addButton(_ui.checkBoxMaOther, OTHER);
    connect(group, SIGNAL(buttonClicked(QAbstractButton*)), SLOT(maGroupActivated(QAbstractButton*)));
    connect(group, SIGNAL(buttonClicked(int)), SLOT(updateAgentsSettings(void)));

    group = new QButtonGroup(this);
    group->setExclusive(false);
    group->addButton(_ui.checkBoxDrawAgentMr, ENBL);
    group->addButton(_ui.checkBoxMrStat1, STAT1);
    group->addButton(_ui.checkBoxMrNFkB, NFKB);
    group->addButton(_ui.checkBoxMrDeact, DEACT);
    group->addButton(_ui.checkBoxMrOther, OTHER);
    connect(group, SIGNAL(buttonClicked(QAbstractButton*)), SLOT(mrGroupActivated(QAbstractButton*)));
    connect(group, SIGNAL(buttonClicked(int)), SLOT(updateAgentsSettings(void)));

    group = new QButtonGroup(this);
    group->setExclusive(false);
    group->addButton(_ui.checkBoxDrawAgentMi, ENBL);
    group->addButton(_ui.checkBoxMiStat1, STAT1);
    group->addButton(_ui.checkBoxMiNFkB, NFKB);
    group->addButton(_ui.checkBoxMiDeact, DEACT);
    group->addButton(_ui.checkBoxMiOther, OTHER);
    connect(group, SIGNAL(buttonClicked(QAbstractButton*)), SLOT(miGroupActivated(QAbstractButton*)));
    connect(group, SIGNAL(buttonClicked(int)), SLOT(updateAgentsSettings(void)));

    group = new QButtonGroup(this);
    group->setExclusive(false);
    group->addButton(_ui.checkBoxDrawAgentMci, ENBL);
    group->addButton(_ui.checkBoxMciStat1, STAT1);
    group->addButton(_ui.checkBoxMciNFkB, NFKB);
    group->addButton(_ui.checkBoxMciDeact, DEACT);
    group->addButton(_ui.checkBoxMciOther, OTHER);
    connect(group, SIGNAL(buttonClicked(QAbstractButton*)), SLOT(mciGroupActivated(QAbstractButton*)));
    connect(group, SIGNAL(buttonClicked(int)), SLOT(updateAgentsSettings(void)));

    connect(_ui.checkBoxDrawAgentTgam, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentTcyt, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));
    connect(_ui.checkBoxDrawAgentTreg, SIGNAL(toggled(bool)), this, SLOT(updateAgentsSettings(void)));

    updateAgentsSettings();
}

AgentsWidget::~AgentsWidget()
{
}

void AgentsWidget::maGroupActivated(QAbstractButton* b) {
    emit agentFilterChanged(MAC, Mac::MAC_ACTIVE, b->group()->id(b), b->isChecked());
}

void AgentsWidget::mrGroupActivated(QAbstractButton* b) {
    QButtonGroup* group = b->group();
    Q_CHECK_PTR(group);
    emit agentFilterChanged(MAC, Mac::MAC_RESTING, group->id(b), b->isChecked());
}

void AgentsWidget::miGroupActivated(QAbstractButton* b) {
    emit agentFilterChanged(MAC, Mac::MAC_INFECTED, b->group()->id(b), b->isChecked());
}

void AgentsWidget::mciGroupActivated(QAbstractButton* b) {
    emit agentFilterChanged(MAC, Mac::MAC_CINFECTED, b->group()->id(b), b->isChecked());
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
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, DEACT, _ui.checkBoxMaDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, STAT1, _ui.checkBoxMaStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, NFKB, _ui.checkBoxMaNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, OTHER, _ui.checkBoxMaOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_ACTIVE, ENBL, _ui.checkBoxDrawAgentMa->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, DEACT, _ui.checkBoxMrDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, STAT1, _ui.checkBoxMrStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, NFKB, _ui.checkBoxMrNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, OTHER, _ui.checkBoxMrOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_RESTING, ENBL, _ui.checkBoxDrawAgentMr->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, DEACT, _ui.checkBoxMiDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, STAT1, _ui.checkBoxMiStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, NFKB, _ui.checkBoxMiNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, OTHER, _ui.checkBoxMiOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_INFECTED, ENBL, _ui.checkBoxDrawAgentMi->isChecked());

    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, DEACT, _ui.checkBoxMciDeact->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, STAT1, _ui.checkBoxMciStat1->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, NFKB, _ui.checkBoxMciNFkB->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, OTHER, _ui.checkBoxMciOther->isChecked());
    _pAgentsVisualization->setDrawMac(Mac::MAC_CINFECTED, ENBL, _ui.checkBoxDrawAgentMci->isChecked());

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
