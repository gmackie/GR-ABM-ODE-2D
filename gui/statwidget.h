#ifndef STATWIDGET_H
#define STATWIDGET_H

#include <QtGui/QWidget>
#include "ui_statwidget.h"
#include "simulation/gr.h"

class StatWidget : public QWidget
{
  Q_OBJECT
  friend struct StatWidgetNameVisitor;
  friend struct StatWidgetUpdateVisitor;
public:
  StatWidget(const Stats& stats, QWidget* pParent = 0);
  ~StatWidget();
  void updateLabels(const Stats& stats);

private:
  Ui::StatWidgetClass _ui;
};

#endif // STATWIDGET_H
