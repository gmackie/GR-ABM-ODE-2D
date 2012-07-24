#ifndef PARAMWINDOW_H
#define PARAMWINDOW_H

#include <QtGui/QWidget>
#include <QTreeWidgetItem>
#include "ui_paramwindow.h"
#include "simulation/params.h"

class MainInterface;

class ParamWindow : public QWidget
{
  Q_OBJECT

public:
  ParamWindow(MainInterface* pItfc, QWidget* parent = 0);
  ~ParamWindow();
  void update();

public slots:
  void loadParams();
  void saveParams();

private:
  MainInterface* _pItfc;
  Ui::ParamWindowClass _ui;

  void newItem(QTreeWidgetItem* pParentItem, ParamDoubleType param);
  void newItem(QTreeWidgetItem* pParentItem, ParamIntType param);
  void processElement(const TiXmlElement* pElement, QTreeWidgetItem* pParentTreeWidgetItem);
};

#endif // PARAMWINDOW_H
