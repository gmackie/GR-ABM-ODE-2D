#include "paramwindow.h"
#include "simulation/tinyxml/tinyxml.h"
#include <QFileDialog>
#include <QMessageBox>

ParamWindow::ParamWindow(MainInterface* pItfc, QWidget* parent)
    : QWidget(parent)
    , _pItfc(pItfc)
 {
	_ui.setupUi(this);
	_ui.treeWidget->setColumnWidth(0, 200);
	_ui.treeWidget->setColumnWidth(1, 50);
	_ui.treeWidget->setColumnWidth(2, 100);
	update();

	connect(_ui.pushButtonLoad, SIGNAL(clicked(bool)), this, SLOT(loadParams(void)));
	connect(_ui.pushButtonSave, SIGNAL(clicked(bool)), this, SLOT(saveParams(void)));
}

ParamWindow::~ParamWindow()
{
}

void ParamWindow::newItem(QTreeWidgetItem* pParentItem, ParamDoubleType param)
{
	const Params* params = Params::getInstance();

	QTreeWidgetItem* pItem = new QTreeWidgetItem(pParentItem, param);
	pItem->setText(0, QString("%1").arg(params->getName(param).c_str()));
	pItem->setText(1, QString("%1").arg(params->getParam(param)));
	pItem->setText(2, QString("%1").arg(params->getUnit(param).c_str()));
	pItem->setText(3, QString("%1").arg(params->getDescription(param).c_str()));
}

void ParamWindow::newItem(QTreeWidgetItem* pParentItem, ParamIntType param)
{
	const Params* params = Params::getInstance();

	QTreeWidgetItem* pItem = new QTreeWidgetItem(pParentItem);
	pItem->setText(0, QString("%1").arg(params->getName(param).c_str()));
	pItem->setText(1, QString("%1").arg(params->getParam(param)));
	pItem->setText(2, QString("%1").arg(params->getUnit(param).c_str()));
	pItem->setText(3, QString("%1").arg(params->getDescription(param).c_str()));
}

// Add each XML element and its parameters (XML attributes) to the tree view.
// Also do this for each of its child elements.
void ParamWindow::processElement(const TiXmlElement* pElement, QTreeWidgetItem* pParentTreeWidgetItem)
{
	const Params* params = Params::getInstance();

	pParentTreeWidgetItem->setExpanded(true);
	pParentTreeWidgetItem->setText(0, pElement->Value());

	XmlElement xmlElement = params->findElementDescription(pElement->Value());
	if (xmlElement == NULL_NODE)
	{
		return;
	}

	pParentTreeWidgetItem->setText(3,  ParamsBase::getXmlElementDescription(xmlElement).c_str());

	// Process the parameters for the element attributes.
	// If an attribute name does not match any parameter name then ignore it.
	// If an attribute name does match a parameter name then display it.
	const TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	while (pAttrib)
	{
		std::string attributeName = std::string(pAttrib->Name());

		int parameterIndex = params->findParameterDescription(attributeName, pElement);

		if (parameterIndex >= 0)
		{
			if (params->isDouble(parameterIndex))
			{
				newItem(pParentTreeWidgetItem, (ParamDoubleType) parameterIndex);
			}
			else
			{
				newItem(pParentTreeWidgetItem, (ParamIntType) (params->intIndex(parameterIndex)));
			}

		}

		pAttrib=pAttrib->Next();
	}

	// Process the parameters from the children of this element that are themselves XML elements.
	const TiXmlNode* child = 0;
	while( (child = pElement->IterateChildren( child )) )
	{
		// Ignore children that are not XML elements. Also ignore the Init element and its children.
		if (child->ToElement() && strcmp(child->Value(), ParamsBase::getXmlElementName(INIT_NODE).c_str()))
		{
			QTreeWidgetItem* pChildTreeWidgetItem = new QTreeWidgetItem(pParentTreeWidgetItem);
			processElement(child->ToElement(), pChildTreeWidgetItem);
		}
	}
}

// Place all the parameters in a tree view window.
// Each XML element appears as a node in the tree view, with its parameters (XML attributes) and child
// elements as sub-nodes in the tree view. The items to show in the tree view are determined by a
// traversal of the tiny XML document object created from reading the parameter file. This means that only
// parameters that are present in the parameter file and are defined in the ParamsBase class are shown.
// Any parameters not present in the parameter file are not shown, such as those that are calculated in
// the code or have default values used if not present in the parameter file. Also any parameter in a
// parameter file that is not defined in the ParamsBase class is not shown.
void ParamWindow::update()
{
	QTreeWidget* pTreeWidget = _ui.treeWidget;
	pTreeWidget->clear();

	const TiXmlElement* pRootElement = Params::getInstance()->getXmlDoc().RootElement();

	// Process the parameters starting at the root XML element of the parameter file.
	// This only processes parameters that were present in the parameter file that was read.
	// We have to define the QTreeWidgetItem here, not in processElement, because there
	// are different QTreeWidgetItem constructors with a constructor argument of a parent
	// widget of QTreeWidget and for a parent widget of QTreeWidgetItem.
	QTreeWidgetItem* pTreeWidgetItem = new QTreeWidgetItem(pTreeWidget);
	processElement(pRootElement, pTreeWidgetItem);
	_ui.treeWidget->insertTopLevelItem(0, pTreeWidgetItem);
}

void ParamWindow::loadParams()
{
	Simulation& sim = _pItfc->getSimulation();

	QString fileName = QFileDialog::getOpenFileName(this, "Load parameters", "", "*.xml");
	if (fileName != QString::null)
	{
		sim.lock();
		if (!Params::reinit(fileName.toLatin1().data()))
		{
			QMessageBox::critical(this, "Lung ABM",
					"Failed to open file '" + fileName + "'.",
					QMessageBox::Ok, QMessageBox::Ok);
		}
		sim.unlock();
	}
}

void ParamWindow::saveParams()
{
	QString fileName = QFileDialog::getSaveFileName(this, "Save parameters", "", "*.xml");
	if (fileName != QString::null)
	{
		if (!Params::getInstance()->toXml(fileName.toLatin1().data()))
		{
			QMessageBox::critical(this, "Lung ABM",
					"Failed to save file '" + fileName + "'.",
					QMessageBox::Ok, QMessageBox::Ok);
		}
	}
}
