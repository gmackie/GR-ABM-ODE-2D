#include "agentswidget.h"
#include "simulation/gr.h"
#include "simulation/macrophage.h"
#include "visualization/agentsvisualization.h"
#include <QSettings>

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

  connect(_ui.checkBoxDrawAgentTgam, SIGNAL(toggled(bool)), this, SLOT(tgamFilterChanged(bool)));
  connect(_ui.checkBoxDrawAgentTcyt, SIGNAL(toggled(bool)), this, SLOT(tcytFilterChanged(bool)));
  connect(_ui.checkBoxDrawAgentTreg, SIGNAL(toggled(bool)), this, SLOT(tregFilterChanged(bool)));

  updateAgentsSettings();
}

AgentsWidget::~AgentsWidget()
{
}

void AgentsWidget::maGroupActivated(QAbstractButton* b)
{
  emit agentFilterChanged(MAC, Mac::MAC_ACTIVE, b->group()->id(b), b->isChecked());
}

void AgentsWidget::mrGroupActivated(QAbstractButton* b)
{
  QButtonGroup* group = b->group();
  Q_CHECK_PTR(group);
  emit agentFilterChanged(MAC, Mac::MAC_RESTING, group->id(b), b->isChecked());
}

void AgentsWidget::miGroupActivated(QAbstractButton* b)
{
  emit agentFilterChanged(MAC, Mac::MAC_INFECTED, b->group()->id(b), b->isChecked());
}

void AgentsWidget::mciGroupActivated(QAbstractButton* b)
{
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

  _pAgentsVisualization->setGridAlpha((1.0*_ui.horizontalSliderAgentsAlpha->value()) / _ui.horizontalSliderAgentsAlpha->maximum());
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

void AgentsWidget::tgamFilterChanged(bool en)
{
  emit agentFilterChanged(TGAM, 0, ENBL, en);
}
void AgentsWidget::tregFilterChanged(bool en)
{
  emit agentFilterChanged(TREG, 0, ENBL, en);
}
void AgentsWidget::tcytFilterChanged(bool en)
{
  emit agentFilterChanged(TCYT, 0, ENBL, en);
}

void AgentsWidget::saveSettings() const
{
  QSettings settings(QApplication::applicationDirPath() + "/.config", QSettings::IniFormat);
  settings.beginGroup("AgentSettings");
  settings.setValue("caseation", _ui.checkBoxDrawAgentCaseation->isChecked());
  settings.setValue("sources", _ui.checkBoxDrawAgentSources->isChecked() );
  settings.setValue("extmtb", _ui.checkBoxDrawAgentExtMtb->isChecked() );
  settings.setValue("grid", _ui.checkBoxAgentDrawGrid->isChecked() );
  settings.setValue("alpha", _ui.horizontalSliderAgentsAlpha->value() );
  settings.setValue("gridHeight", _ui.doubleSpinBoxAgentsGridHeight->value() );
  settings.setValue("MacSrcs", _ui.checkBoxDrawAgentSourcesMac->isChecked() );
  settings.setValue("TgamSrcs", _ui.checkBoxDrawAgentSourcesTgam->isChecked() );
  settings.setValue("TcytSrcs", _ui.checkBoxDrawAgentSourcesTcyt->isChecked() );
  settings.setValue("TregSrcs", _ui.checkBoxDrawAgentSourcesTreg->isChecked() );
  settings.setValue("SquareAgents", _ui.checkboxSquareAgents->isChecked() );
  settings.setValue("M1M2", _ui.comboBoxM1M2->currentIndex() );
  settings.setValue("M1M2Thres", _ui.doubleSpinBoxM1M2Threshold->value() );

  settings.setValue("MaENBL", _ui.checkBoxDrawAgentMa->isChecked());
  settings.setValue("MaSTAT1", _ui.checkBoxMaStat1->isChecked());
  settings.setValue("MaNFkB", _ui.checkBoxMaNFkB->isChecked());
  settings.setValue("MaDEACT", _ui.checkBoxMaDeact->isChecked());
  settings.setValue("MaOTHER", _ui.checkBoxMaOther->isChecked());

  settings.setValue("MrENBL", _ui.checkBoxDrawAgentMr->isChecked());
  settings.setValue("MrSTAT1", _ui.checkBoxMrStat1->isChecked());
  settings.setValue("MrNFkB", _ui.checkBoxMrNFkB->isChecked());
  settings.setValue("MrDEACT", _ui.checkBoxMrDeact->isChecked());
  settings.setValue("MrOTHER", _ui.checkBoxMrOther->isChecked());

  settings.setValue("MiENBL", _ui.checkBoxDrawAgentMi->isChecked());
  settings.setValue("MiSTAT1", _ui.checkBoxMiStat1->isChecked());
  settings.setValue("MiNFkB", _ui.checkBoxMiNFkB->isChecked());
  settings.setValue("MiDEACT", _ui.checkBoxMiDeact->isChecked());
  settings.setValue("MiOTHER", _ui.checkBoxMiOther->isChecked());

  settings.setValue("MciENBL", _ui.checkBoxDrawAgentMci->isChecked());
  settings.setValue("MciSTAT1", _ui.checkBoxMciStat1->isChecked());
  settings.setValue("MciNFkB", _ui.checkBoxMciNFkB->isChecked());
  settings.setValue("MciDEACT", _ui.checkBoxMciDeact->isChecked());
  settings.setValue("MciOTHER", _ui.checkBoxMciOther->isChecked());

  settings.setValue("TgamENBL", _ui.checkBoxDrawAgentTgam->isChecked());
  settings.setValue("TcytENBL", _ui.checkBoxDrawAgentTcyt->isChecked());
  settings.setValue("TregENBL", _ui.checkBoxDrawAgentTreg->isChecked());
  settings.endGroup();
}
void AgentsWidget::loadSettings()
{
  blockSignals(true); //Prevent updates to glwindow until all is done
  QSettings settings(QApplication::applicationDirPath() + "/.config", QSettings::IniFormat);
  settings.beginGroup("AgentSettings");
  _ui.checkBoxDrawAgentCaseation->setChecked(  settings.value("caseation", _ui.checkBoxDrawAgentCaseation->isChecked()).toBool());
  _ui.checkBoxDrawAgentSources->setChecked(    settings.value("sources",   _ui.checkBoxDrawAgentSources->isChecked() ).toBool());
  _ui.checkBoxDrawAgentExtMtb->setChecked(     settings.value("extmtb",    _ui.checkBoxDrawAgentExtMtb->isChecked() ).toBool());
  _ui.checkBoxAgentDrawGrid->setChecked(       settings.value("grid",      _ui.checkBoxAgentDrawGrid->isChecked() ).toBool());
  _ui.horizontalSliderAgentsAlpha->setValue(   settings.value("alpha",     _ui.horizontalSliderAgentsAlpha->value() ).toInt());
  _ui.doubleSpinBoxAgentsGridHeight->setValue( settings.value("gridHeight", _ui.doubleSpinBoxAgentsGridHeight->value() ).toBool());
  _ui.checkBoxDrawAgentSourcesMac->setChecked( settings.value("MacSrcs",   _ui.checkBoxDrawAgentSourcesMac->isChecked() ).toBool());
  _ui.checkBoxDrawAgentSourcesTgam->setChecked(settings.value("TgamSrcs",  _ui.checkBoxDrawAgentSourcesTgam->isChecked() ).toBool());
  _ui.checkBoxDrawAgentSourcesTcyt->setChecked(settings.value("TcytSrcs",  _ui.checkBoxDrawAgentSourcesTcyt->isChecked() ).toBool());
  _ui.checkBoxDrawAgentSourcesTreg->setChecked(settings.value("TregSrcs",  _ui.checkBoxDrawAgentSourcesTreg->isChecked() ).toBool());
  _ui.checkboxSquareAgents->setChecked(        settings.value("SquareAgents", _ui.checkboxSquareAgents->isChecked() ).toBool());
  _ui.comboBoxM1M2->setCurrentIndex(           settings.value("M1M2",      _ui.comboBoxM1M2->currentIndex() ).toInt());
  _ui.doubleSpinBoxM1M2Threshold->setValue(    settings.value("M1M2Thres", _ui.doubleSpinBoxM1M2Threshold->value() ).toDouble());

  _ui.checkBoxDrawAgentMa->setChecked(settings.value("MaENBL",  _ui.checkBoxDrawAgentMa->isChecked()).toBool());
  _ui.checkBoxMaStat1->setChecked(    settings.value("MaSTAT1", _ui.checkBoxMaStat1->isChecked()).toBool());
  _ui.checkBoxMaNFkB->setChecked(     settings.value("MaNFkB",  _ui.checkBoxMaNFkB->isChecked()).toBool());
  _ui.checkBoxMaDeact->setChecked(    settings.value("MaDEACT", _ui.checkBoxMaDeact->isChecked()).toBool());
  _ui.checkBoxMaOther->setChecked(    settings.value("MaOTHER", _ui.checkBoxMaOther->isChecked()).toBool());

  _ui.checkBoxDrawAgentMr->setChecked(settings.value("MrENBL",  _ui.checkBoxDrawAgentMr->isChecked()).toBool());
  _ui.checkBoxMrStat1->setChecked(    settings.value("MrSTAT1", _ui.checkBoxMrStat1->isChecked()).toBool());
  _ui.checkBoxMrNFkB->setChecked(     settings.value("MrNFkB",  _ui.checkBoxMrNFkB->isChecked()).toBool());
  _ui.checkBoxMrDeact->setChecked(    settings.value("MrDEACT", _ui.checkBoxMrDeact->isChecked()).toBool());
  _ui.checkBoxMrOther->setChecked(    settings.value("MrOTHER", _ui.checkBoxMrOther->isChecked()).toBool());

  _ui.checkBoxDrawAgentMi->setChecked(settings.value("MiENBL",  _ui.checkBoxDrawAgentMi->isChecked()).toBool());
  _ui.checkBoxMiStat1->setChecked(    settings.value("MiSTAT1", _ui.checkBoxMiStat1->isChecked()).toBool());
  _ui.checkBoxMiNFkB->setChecked(     settings.value("MiNFkB",  _ui.checkBoxMiNFkB->isChecked()).toBool());
  _ui.checkBoxMiDeact->setChecked(    settings.value("MiDEACT", _ui.checkBoxMiDeact->isChecked()).toBool());
  _ui.checkBoxMiOther->setChecked(    settings.value("MiOTHER", _ui.checkBoxMiOther->isChecked()).toBool());

  _ui.checkBoxDrawAgentMci->setChecked(settings.value("MciENBL",  _ui.checkBoxDrawAgentMci->isChecked()).toBool());
  _ui.checkBoxMciStat1->setChecked(    settings.value("MciSTAT1", _ui.checkBoxMciStat1->isChecked()).toBool());
  _ui.checkBoxMciNFkB->setChecked(     settings.value("MciNFkB",  _ui.checkBoxMciNFkB->isChecked()).toBool());
  _ui.checkBoxMciDeact->setChecked(    settings.value("MciDEACT", _ui.checkBoxMciDeact->isChecked()).toBool());
  _ui.checkBoxMciOther->setChecked(    settings.value("MciOTHER", _ui.checkBoxMciOther->isChecked()).toBool());

  _ui.checkBoxDrawAgentTgam->setChecked(settings.value("TgamENBL", _ui.checkBoxDrawAgentTgam->isChecked()).toBool());
  _ui.checkBoxDrawAgentTcyt->setChecked(settings.value("TcytENBL", _ui.checkBoxDrawAgentTcyt->isChecked()).toBool());
  _ui.checkBoxDrawAgentTreg->setChecked(settings.value("TregENBL", _ui.checkBoxDrawAgentTreg->isChecked()).toBool());
  settings.endGroup();
  blockSignals(false);
  emit updateGL();
}

void AgentsWidget::on_checkboxSquareAgents_toggled(bool checked)
{
  _pAgentsVisualization->setDrawSquares(checked);
  emit updateGL();
}
