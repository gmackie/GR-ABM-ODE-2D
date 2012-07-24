#ifndef GLWINDOW_H
#define GLWINDOW_H

#include <QtGui/QWidget>
#include "ui_glwindow.h"
#include "simulation/stat.h"
#include <QString>

class ColorMap;
class MainInterface;
class QImage;
class QCloseEvent;
class QTreeWidget;
class QTreeWidgetItem;

class GLWindow : public QWidget
{
    Q_OBJECT
    int _trackid;

public:
    GLWindow(MainInterface* pItfc, QWidget* parent = 0);
    ~GLWindow();
    QImage grabFrameBuffer();
    void resizeGLWidget(int width, int height);
    void setGrStatus(int index, GrStatus status);
    bool getPrintTime() const;
    bool getPrintOutcome() const;

public slots:
	void setColorMap(ColorMap* pCurrentColorMap);
	void updateWindow();
	void updateColorMapLabels(float min, float max);
	void toggleFullScreen();
    void updateSelectedCellStats();
	void selectCell(int row, int col);
	void moveSelectionLeft();
	void moveSelectionRight();
	void moveSelectionUp();
	void moveSelectionDown();
	void visualize();
	void printText();
	void setPrintTime(bool value);
	void setPrintOutcome(bool value);

protected slots:
    void updateTracking(QTreeWidgetItem*);
signals:
	void updateSelection(int row, int col);
	void set2DView();
	void set3DView();

protected:
	void closeEvent(QCloseEvent* pEvent);

private:
    Ui::GLWindowClass _ui;
    MainInterface* _pItfc;
    int _selRow;
    int _selCol;
    GrStatus _status[NOUTCOMES];
    bool _printTime;
    bool _printOutcome;
    QTreeWidget* agentInfoWindow;

    QString getAgentStr(const Agent* pAgent);
};

inline void GLWindow::setPrintTime(bool value)
{
	_printTime = value;
	updateWindow();
}

inline void GLWindow::setPrintOutcome(bool value)
{
	_printOutcome = value;
	updateWindow();
}

inline bool GLWindow::getPrintTime() const
{
	return _printTime;
}

inline bool GLWindow::getPrintOutcome() const
{
	return _printOutcome;
}

inline void GLWindow::setGrStatus(int index, GrStatus status)
{
	assert(0 <= index && index < NOUTCOMES);
	_status[index] = status;
}

#endif // GLWINDOW_H
