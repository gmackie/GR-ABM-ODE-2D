/*
 * colorMapWidget.h
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#ifndef COLORMAPWIDGET_H_
#define COLORMAPWIDGET_H_

#include <QWidget>
#include <QColor>
#include <QPair>

class ColorMap;

typedef std::vector<QColor> ColorVector;

class ColorMapWidget : public QWidget
{
	Q_OBJECT

public:
    ColorMapWidget(QWidget* parent = 0);
    ~ColorMapWidget();
    void setColorMap(ColorMap* pColorMap);
    const ColorMap* getColorMap();

protected:
    void paintEvent(QPaintEvent* event);

private:
    ColorMap* _pColorMap;
};

inline const ColorMap* ColorMapWidget::getColorMap()
{
    return _pColorMap;
}

inline void ColorMapWidget::setColorMap(ColorMap* pColorMap)
{
    _pColorMap = pColorMap;
    update();
}

#endif /* COLORMAPWIDGET_H_ */
