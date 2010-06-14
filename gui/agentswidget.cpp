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
}

AgentsWidget::~AgentsWidget()
{
}

void AgentsWidget::setAgentSelection(int row, int col)
{
	_pAgentsVisualization->setSelection(row, col);
}

void AgentsWidget::updateAgentsSettings()
{
	_pAgentsVisualization->setGridAlpha(_ui.horizontalSliderAgentsAlpha->value() / _SLIDER_AGENTS_ALPHA);
	_pAgentsVisualization->setDrawGrid(_ui.checkBoxAgentDrawGrid->isChecked());
	_pAgentsVisualization->setDrawMacResting(_ui.checkBoxDrawAgentMr->isChecked());
	_pAgentsVisualization->setDrawMacInfected(_ui.checkBoxDrawAgentMi->isChecked());
	_pAgentsVisualization->setDrawMacCInfected(_ui.checkBoxDrawAgentMci->isChecked());
	_pAgentsVisualization->setDrawMacActive(_ui.checkBoxDrawAgentMa->isChecked());
	_pAgentsVisualization->setDrawTgam(_ui.checkBoxDrawAgentTgam->isChecked());
	_pAgentsVisualization->setDrawTreg(_ui.checkBoxDrawAgentTreg->isChecked());
	_pAgentsVisualization->setDrawTcyt(_ui.checkBoxDrawAgentTcyt->isChecked());
	_pAgentsVisualization->setDrawCas(_ui.checkBoxDrawAgentCaseation->isChecked());
	_pAgentsVisualization->setDrawSrc(_ui.checkBoxDrawAgentSources->isChecked());
	_pAgentsVisualization->setDrawExtMtb(_ui.checkBoxDrawAgentExtMtb->isChecked());
	_pAgentsVisualization->setGridHeight((float) _ui.doubleSpinBoxAgentsGridHeight->value());
	_pAgentsVisualization->setPredicates(_ui.checkBoxDrawAgentSourcesMac->isChecked(),
			_ui.checkBoxDrawAgentSourcesTgam->isChecked(), _ui.checkBoxDrawAgentSourcesTcyt->isChecked(),
			_ui.checkBoxDrawAgentSourcesTreg->isChecked());

	emit updateGL();
}
