#ifndef PARAMWINDOW_H
#define PARAMWINDOW_H

#include <QtGui/QWidget>
#include <QTreeWidgetItem>
#include "ui_paramwindow.h"
#include "simulation/params.h"
#include "maininterface.h"

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
    QTreeWidgetItem* _pGR;
    QTreeWidgetItem* _pMac;
    QTreeWidgetItem* _pTcell;
    QTreeWidgetItem* _pTgam;
    QTreeWidgetItem* _pTcyt;
    QTreeWidgetItem* _pTreg;
    QTreeWidgetItem* _pMtb;
    void newItem(QTreeWidgetItem* pParentItem, ParamDoubleType param);
    void newItem(QTreeWidgetItem* pParentItem, ParamIntType param);
};

#endif // PARAMWINDOW_H
