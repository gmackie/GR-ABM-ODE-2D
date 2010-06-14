#ifndef STATWIDGET_H
#define STATWIDGET_H

#include <QtGui/QWidget>
#include "ui_statwidget.h"
#include "maininterface.h"

class StatWidget : public QWidget
{
    Q_OBJECT

public:
    StatWidget(QWidget* pParent = 0);
    ~StatWidget();
    void updateLabels(const GrStat& stats);

private:
    Ui::StatWidgetClass _ui;
};

#endif // STATWIDGET_H
