/*
 * colorMapWidget.cpp
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#include "colormapwidget.h"
#include "grviz.h"
#include <QPainterPath>
#include <QPainter>
#include <QColor>

ColorMapWidget::ColorMapWidget(QWidget* parent)
	: QWidget(parent)
	, _pColorMap(NULL)
{
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(false);
}

ColorMapWidget::~ColorMapWidget()
{
}

void ColorMapWidget::paintEvent(QPaintEvent*)
{
    if (_pColorMap)
    {
		int nColor = _pColorMap->getNrBands();
	    assert(nColor >= 1);

		/* Draw color map */
		QPainter colorMapPainter(this);
		colorMapPainter.scale(width() / (double)nColor, 1);

		float R, G, B;
		for (int i = 0; i < nColor; i++)
		{
			float value = i;
			value = value / nColor + FLT_EPSILON;

			_pColorMap->map(value, R, G, B);

			QColor color = QColor::fromRgbF(R, G, B);

			colorMapPainter.setPen(color);
			colorMapPainter.setBrush(color);
			colorMapPainter.drawRect(i, 0, 1, height());
		}
    }
    else
    {
    	QPainter painter(this);
    	painter.setPen(Qt::black);
    	painter.setBrush(Qt::black);
    	painter.drawRect(0, 0, width(), height());
    }
}
