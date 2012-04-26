#include "glwindow.h"

GLWindow::GLWindow(MainInterface* pItfc, QWidget* parent)
    : QWidget(parent)
    , _ui()
    , _pItfc(pItfc)
	, _selRow(-1)
	, _selCol(-1)
	, _printTime(_PRINT_TIME)
	, _printOutcome(_PRINT_OUTCOME)
{
	for (int i = 0; i < NOUTCOMES; i++)
	{
		_status[i] = GR_NONE;
	}

	_ui.setupUi(this);

    _ui.glWidget->dim = pItfc->getSimulation().getSize();
	connect(_ui.glWidget, SIGNAL(updateSelection(int, int)), this, SIGNAL(updateSelection(int, int)));
	connect(_ui.glWidget, SIGNAL(visualize(void)), this, SLOT(visualize(void)));
	connect(_ui.glWidget, SIGNAL(printText(void)), this, SLOT(printText(void)));
	connect(this, SIGNAL(set2DView(void)), _ui.glWidget, SLOT(set2DView(void)));
	connect(this, SIGNAL(set3DView(void)), _ui.glWidget, SLOT(set3DViewHeight(void)));
}

GLWindow::~GLWindow()
{
}

void GLWindow::setColorMap(ColorMap* pCurrentColorMap)
{
	_ui.colorMapWidget->setColorMap(pCurrentColorMap);
}

void GLWindow::updateWindow()
{
	_ui.glWidget->updateGL();
}

void GLWindow::toggleFullScreen()
{
	if (isFullScreen())
	{
		_ui.colorMapWidget->show();
		_ui.labelLegend1->show();
		_ui.labelLegend2->show();
		_ui.labelLegend3->show();
		_ui.labelLegend4->show();
		_ui.labelLegend5->show();

		showNormal();
	}
	else
	{
		_ui.colorMapWidget->hide();
		_ui.labelLegend1->hide();
		_ui.labelLegend2->hide();
		_ui.labelLegend3->hide();
		_ui.labelLegend4->hide();
		_ui.labelLegend5->hide();

		showFullScreen();
	}
}

void GLWindow::updateSelectedCellStats()
{
	if (_selRow != -1 && _selCol != -1)
	{
		const ScalarAgentGrid* pAgentGrid = static_cast<ScalarAgentGrid*>(_pItfc->getScalarAgentGrid());
		const ScalarAgentItem& item
			= pAgentGrid->getGrid()[_selRow * _ui.glWidget->dim.x + _selCol];

		QString str = QString("pos = (%1,%2)").arg(_selCol).arg(_selRow);

		if (GET_BIT(item._bitMask, ScalarAgentGrid::_bitCas))
		{
			str += ", caseous";
		}
		else
		{
			str += QString(", killings = %1").arg(item._nKillings);
		}

		if (GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrc))
		{
			str += ", source";

			bool srcMac = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcMac);
			bool srcTgam= GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTgam);
			bool srcTcyt = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTcyt);
			bool srcTreg = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTreg);

			str += " (";

			bool first = false;
			if (srcMac)
			{
				str += "M";
				first = true;
			}

			if (srcTgam)
			{
				if (first) str += ",";
				else first = true;
				str += "Tg";
			}

			if (srcTcyt)
			{
				if (first) str += ",";
				else first = true;
				str += "Tc";
			}

			if (srcTreg)
			{
				if (first) str += ",";
				else first = true;
				str += "Tr";
			}
			str += ")";
		}

		str += QString(", Be = %1").arg(item._extMtb, 0, 'f', 2);

		glColor3f(1.0f, 0.0f, 0.0f);
		_ui.glWidget->renderText(10, 20, str);

		str = QString("(%1, %2, %3, %4, %5)").arg(item._attractant, 0, 'f', 2).
			arg(item._TNF, 0, 'f', 2).arg(item._CCL2, 0, 'f', 2).
			arg(item._CCL5, 0, 'f', 2).arg(item._CXCL9, 0, 'f', 2);
		_ui.glWidget->renderText(10, 40, str);

#if 0
		// For debugging granuloma boundary definition.
		GridCell cell = _pItfc->getSimulation().getGrGrid()(_selRow, _selCol);
		bool occupied = cell.isOccupied();
		str = QString("(%1, %2, %3").arg(item._TNF, 0, 'f', 2).
									  arg(item._attractant, 0, 'f', 2).
									  arg(item._extMtb, 0, 'f', 2);
		str += occupied ? ", O )" : ", NO)";

		_ui.glWidget->renderText(10, 60, str);
#endif

		str = getAgentStr(item._pAgent[0]);
		_ui.glWidget->renderText(10, _ui.glWidget->height() - 10, str);

		str = getAgentStr(item._pAgent[1]);
		_ui.glWidget->renderText(10, _ui.glWidget->height() - 30, str);
	}
}

QString GLWindow::getAgentStr(const Agent* pAgent)
{
	QString res;

	if (Mac::isMac(pAgent))
	{
		const Mac* pMac = dynamic_cast<const Mac*>(pAgent);
		switch (pMac->getState())
		{
		case Mac::MAC_RESTING:
			res = "Mr";
			break;
		case Mac::MAC_INFECTED:
			res = "Mi";
			break;
		case Mac::MAC_CINFECTED:
			res = "Mci";
			break;
		case Mac::MAC_ACTIVE:
		{
			res = "Ma";

			int days, hours, minutes;
			GrSimulation::convertSimTime(pMac->getActivationTime(), days, hours, minutes);
			res += QString(", act = %1d%2h%3m").arg(days, 3, 10, QChar('0')).
				arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));

			break;
		}
		case Mac::MAC_DEAD:
			return res;
		}

		if (pMac->getNFkB())
		{
			res += ", NFkB";
		}

		if (pMac->getStat1())
		{
			res += ", stat1";
		}

		if (pMac->isDeactivated())
		{
			res += ", deact";
		}

		res += QString(", Bi = %1").arg(pMac->getIntMtb(), 0, 'f', 2);
		
		res += QString(", [sTNF/TNFR1] = %1").arg(pMac->getSurfBoundTNFR1(), 0, 'f', 2);
		
		res += QString(", norACT = %1").arg(pMac->getNormalizedACT(), 0, 'f', 2);
	}
	else if (Tgam::isTgam(pAgent))
	{
		const Tgam* pTgam = dynamic_cast<const Tgam*>(pAgent);
		switch (pTgam->getState())
		{
		case Tgam::TGAM_ACTIVE:
			res = "Tg,a";
			break;
		case Tgam::TGAM_DOWN_REGULATED:
			res = "Tg,reg";
			break;
		case Tgam::TGAM_DEAD:
			return res;
		}
	}
	else if (Tcyt::isTcyt(pAgent))
	{
		const Tcyt* pTcyt = dynamic_cast<const Tcyt*>(pAgent);
		switch (pTcyt->getState())
		{
		case Tcyt::TCYT_ACTIVE:
			res = "Tc";
			break;
		case Tcyt::TCYT_DOWN_REGULATED:
			res = "Tc,reg";
			break;
		case Tcyt::TCYT_DEAD:
			return res;
		}
	}
	else if (Treg::isTreg(pAgent))
	{
		const Treg* pTreg = dynamic_cast<const Treg*>(pAgent);
		switch (pTreg->getState())
		{
		case Treg::TREG_ACTIVE:
			res = "Tr";
			break;
		case Treg::TREG_DEAD:
			return res;
		}
	}

	if (!res.isEmpty())
	{
		int days, hours, minutes;
/*		GrSimulation::convertSimTime(pAgent->getBirthTime(), days, hours, minutes);

		if (pAgent->getBirthTime() < 0)
		{
			res += QString(", birth = -%1d%2h%3m").arg(-1 * days, 3, 10, QChar('0')).
				arg(-1 * hours, 2, 10, QChar('0')).arg(-1 * minutes, 2, 10, QChar('0'));
		}
		else
		{
			res += QString(", birth = %1d%2h%3m").arg(days, 3, 10, QChar('0')).
				arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));
		}
*/
		GrSimulation::convertSimTime(pAgent->getDeathTime(), days, hours, minutes);

		if (pAgent->getDeathTime() < 0)
		{
			res += QString(", death = -%1d%2h%3m").arg(-1 * days, 3, 10, QChar('0')).
				arg(-1 * hours, 2, 10, QChar('0')).arg(-1 * minutes, 2, 10, QChar('0'));
		}
		else
		{
			res += QString(", death = %1d%2h%3m").arg(days, 3, 10, QChar('0')).
				arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));
		}
	}

	return res;
}

QImage GLWindow::grabFrameBuffer()
{
	return _ui.glWidget->renderPixmap().toImage();
}

void GLWindow::updateColorMapLabels(float min, float max)
{
	float f2 = (max - min) * .75 + min;
	float f3 = (max - min) * 0.5 + min;
	float f4 = (max - min) * .25 + min;

	_ui.labelLegend1->setText(QString::number(max, '.', _COLORMAP_LEGEND_PRECISION));
	_ui.labelLegend2->setText(QString::number(f2, '.', _COLORMAP_LEGEND_PRECISION));
	_ui.labelLegend3->setText(QString::number(f3, '.', _COLORMAP_LEGEND_PRECISION));
	_ui.labelLegend4->setText(QString::number(f4, '.', _COLORMAP_LEGEND_PRECISION));
	_ui.labelLegend5->setText(QString::number(min, '.', _COLORMAP_LEGEND_PRECISION));
}

void GLWindow::selectCell(int row, int col)
{
	if (_selRow == row && _selCol == col)
	{
		// unselect
		_selRow = _selCol = -1;
	}
	else
	{
		_selRow = row;
		_selCol = col;
	}

	updateSelectedCellStats();
	updateWindow();
}

void GLWindow::moveSelectionLeft()
{
	if (_selRow != -1 && _selCol != -1)
	{
		emit updateSelection(_selRow, (_selCol - 1 + _ui.glWidget->dim.x) % _ui.glWidget->dim.x);
	}
}

void GLWindow::moveSelectionRight()
{
	if (_selRow != -1 && _selCol != -1)
	{
		emit updateSelection(_selRow, (_selCol + 1 + _ui.glWidget->dim.x) % _ui.glWidget->dim.x);
	}
}

void GLWindow::moveSelectionUp()
{
	if (_selRow != -1 && _selCol != -1)
	{
		emit updateSelection((_selCol + 1 + _ui.glWidget->dim.y) % _ui.glWidget->dim.y, _selCol);
	}
}

void GLWindow::moveSelectionDown()
{
	if (_selRow != -1 && _selCol != -1)
	{
		emit updateSelection((_selCol - 1 + _ui.glWidget->dim.y) % _ui.glWidget->dim.y, _selCol);
	}
}

void GLWindow::visualize()
{
	_pItfc->visualize();
}

void GLWindow::printText()
{
    const Simulation* pSim = &(_pItfc->getSimulation());
	updateSelectedCellStats();

	glColor3f(1.0f, 0.0f, 0.0f);

	int xOffset = 200;

	if (_printTime)
	{
		_ui.glWidget->renderText(_ui.glWidget->width() - xOffset, 20,
			Simulation::getTimeStr(pSim->getTime(), _pItfc->getTime()));
	}

	if (_printOutcome)
	{
		QString text;
		for (int i = 0; i < NOUTCOMES; i++)
		{
			if (i != 0)
				text += ", ";

			switch (_status[i])
			{
			case GR_NONE:
				text += "None";
				break;
			case GR_UNKNOWN:
				text += "Unknown";
				break;
			case GR_CLEARANCE:
				text += "Clearance";
				break;
			case GR_CONTAINMENT:
				text += "Containment";
				break;
			case GR_CONTAINMENT_INCONSISTENT:
				text += "Containment?";
				break;
			case GR_DISSEMINATION:
				text += "Dissemination";
				break;
			case GR_DISSEMINATION_INCONSISTENT:
				text += "Dissemination?";
				break;

			}
		}

		_ui.glWidget->renderText(_ui.glWidget->width() - xOffset, 40, text);
	}
}

void GLWindow::resizeGLWidget(int width, int height)
{
	int dWidth = this->width() - _ui.glWidget->width();
	int dHeight = this->height() - _ui.glWidget->height();

	resize(width + dWidth, height + dHeight);
}

void GLWindow::closeEvent(QCloseEvent* pEvent)
{
	pEvent->accept();
	qApp->quit();
}
