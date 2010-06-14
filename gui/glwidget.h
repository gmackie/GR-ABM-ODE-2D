#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include "maininterface.h"
#include "vectordatasets/vector.h"

#define PI 3.141592653

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget* parent = 0);
    ~GLWidget();
    void setSimulation(Simulation* pSimulation);
    void setVisualization(Visualization* pVisualization);
    void setDepth(float depth);

public slots:
	void updateWindow();
	void set2DView();
	void set3DViewHeight();
	void setAllowMoveCamera(bool allow);

signals:
	void updateSelection(int row, int col);
	void visualize();
	void printText();

protected:
	void paintGL();
	void resizeGL(int width, int height);
	void mouseMoveEvent(QMouseEvent* pEvent);
	void mousePressEvent(QMouseEvent* pEvent);
	void mouseReleaseEvent(QMouseEvent* pEvent);
	void wheelEvent(QWheelEvent* pEvent);

private:
	double _radius;
	double _phi;
	double _theta;
	double _deltaX;
	double _deltaY;
	int _lmx;
	int _lmy;
	bool _leftMouseButtonPressed;
	GLdouble _modelViewMatrix[16];
	GLdouble _projectionMatrix[16];
	GLint _viewport[4];
	float _depth;
	bool _allowMoveCamera;

	void updateProjectionMatrix();
	void updateModelViewMatrix();
	void updateViewport(int width, int height);
	vec3f mapWindowCoordinatesToModelCoordinates(int mx, int my);
	void handleAddSeed(int mx, int my);
};

inline void GLWidget::setDepth(float depth)
{
	_depth = depth;
}

inline void GLWidget::setAllowMoveCamera(bool allow)
{
	_allowMoveCamera = allow;
}

#endif // GLWIDGET_H
