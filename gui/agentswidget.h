#ifndef AGENTSWIDGET_H
#define AGENTSWIDGET_H

#include <QtGui/QWidget>
#include "ui_agentswidget.h"

class AgentsVisualization;

class AgentsWidget : public QWidget
{
  Q_OBJECT
public:
  enum MacSecondStates { ENBL, STAT1, NFKB, DEACT, OTHER, NUM};
  AgentsWidget(AgentsVisualization* pAgentsVisualization, QWidget* parent = 0);
  ~AgentsWidget();

public slots:
  void updateM1M2Settings();
  void updateAgentsSettings();
  void setAgentSelection(int row, int col);
  void saveSettings() const;
  void loadSettings();

private slots:
  void maGroupActivated(QAbstractButton* b);
  void mrGroupActivated(QAbstractButton* b);
  void miGroupActivated(QAbstractButton* b);
  void mciGroupActivated(QAbstractButton* b);
  void tgamFilterChanged(bool en);
  void tregFilterChanged(bool en);
  void tcytFilterChanged(bool en);

  void on_checkboxSquareAgents_toggled(bool checked);

signals:
  void updateGL();
  void agentFilterChanged(int, int, int, bool);

private:
  Ui::AgentsWidgetClass _ui;
  AgentsVisualization* _pAgentsVisualization;
};

#endif // AGENTSWIDGET_H
