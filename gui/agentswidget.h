#ifndef AGENTSWIDGET_H
#define AGENTSWIDGET_H

#include <QtGui/QWidget>
#include "visualization/agentsvisualization.h"
#include "ui_agentswidget.h"

class AgentsWidget : public QWidget
{
    Q_OBJECT

public:
    AgentsWidget(AgentsVisualization* pAgentsVisualization, QWidget* parent = 0);
    ~AgentsWidget();

public slots:
	void updateAgentsSettings();
	void setAgentSelection(int row, int col);

signals:
	void updateGL();

private:
    Ui::AgentsWidgetClass _ui;
    AgentsVisualization* _pAgentsVisualization;
};

#endif // AGENTSWIDGET_H
