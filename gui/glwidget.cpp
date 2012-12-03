#include <QtGui>
#include <QtOpenGL>
#include "glwidget.h"

#include <iostream>
#include <QGLWidget>
#include <QFileDialog>
#include "mainwindow.h"
#include <math.h>

GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent)
  , _radius(100.0)
  , _phi(1.5 * PI)
  , _theta(0)
  , _deltaX(0)
  , _deltaY(0)
  , _lmx(0)
  , _lmy(0)
  , _leftMouseButtonPressed(false)
  , _modelViewMatrix()
  , _projectionMatrix()
  , _viewport()
  , _depth(0.0f)
  , _allowMoveCamera(true)
{
  this->setFocusPolicy(Qt::ClickFocus);
  this->setAutoBufferSwap(true);
}

GLWidget::~GLWidget()
{

}

void GLWidget::initializeGL()
{
  glEnable(GL_MULTISAMPLE);
  glMatrixMode(GL_MODELVIEW);
  glShadeModel(GL_SMOOTH);

  // Create light components
  float ambientLight[] = { 0.9f, 0.9f, 0.9f, 1.0f };
  float diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
  float specularLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };

  // Assign created components to GL_LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);

  float specReflection[] = { 0.1f, 0.1f, 0.1f, 1.0f };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specReflection);
  glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 40);

  updateViewport(width(), height());
  updateProjectionMatrix();
  updateModelViewMatrix();
}
void GLWidget::updateWindow()
{
  if(!_update_posted)
    QMetaObject::invokeMethod(this, "update", Qt::QueuedConnection);
  _update_posted = true;  //Don't queue any more screen updates right now
}

void GLWidget::updateProjectionMatrix()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(60.0, width()/(height()+1.0f), 2.0, 300.0);  //for perspective viewing
  glGetDoublev(GL_PROJECTION_MATRIX, _projectionMatrix);
}

void GLWidget::updateModelViewMatrix()
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  float position1[] = { 0.0f, 0.0f, 60.0f, 0.0f };
  //float position1[] = { 0.0f, 0.0f, 160.0f, 0.0f };

  double x = _deltaX + _radius * cos(_phi) * cos(_theta);
  double y = _deltaY + _radius * sin(_theta);
  double z = -_radius * sin(_phi) * cos(_theta);

  gluLookAt(x, y, z, _deltaX, _deltaY, 0, 0, 1, 0);
  glTranslatef(-(_MAX_X / 2.0f), -(_MAX_Y / 2.0f), _depth);
  glLightfv(GL_LIGHT0, GL_POSITION, position1);

  glGetDoublev(GL_MODELVIEW_MATRIX, _modelViewMatrix);
}

void GLWidget::updateViewport(int width, int height)
{
  glViewport(0, 0, (GLsizei)width, (GLsizei)height);
  glGetIntegerv(GL_VIEWPORT, _viewport);
}

void GLWidget::paintGL()
{
  _update_posted = false;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_LINE_SMOOTH);

  emit visualize();
  emit printText();

  glFlush();
}

void GLWidget::resizeGL(int width, int height)
{
  updateViewport(width, height);
  updateProjectionMatrix();
  updateModelViewMatrix();
  updateGL();
}

void GLWidget::mousePressEvent(QMouseEvent* pEvent)
{
  _leftMouseButtonPressed = pEvent->buttons() & Qt::LeftButton;
  _lmy = pEvent->y();
  _lmx = pEvent->x();
}

void GLWidget::mouseReleaseEvent(QMouseEvent* pEvent)
{
  if (_leftMouseButtonPressed)
    {
      _leftMouseButtonPressed = false;

      int mx = pEvent->x();
      int my = pEvent->y();
      int winHeight = size().height();
      int winWidth = size().width();

      if (mx < 0 || mx > winWidth || my < 0 || my > winHeight)
        return;

      // Compute the array index that corresponds to the cursor location
      vec3f coordinates = mapWindowCoordinatesToModelCoordinates(mx, my);

      int X = (int) floorf(dim.x / _MAX_X * coordinates[0]);
      int Y = (int) floorf(dim.y / _MAX_Y * coordinates[1]);

      if (X > (dim.x - 1))
        return;
      if (Y > (dim.y - 1))
        return;
      if (X < 0)
        return;
      if (Y < 0)
        return;

      emit updateSelection(Y, X);
    }
}

void GLWidget::wheelEvent(QWheelEvent* pEvent)
{
  if (!_allowMoveCamera)
    return;

  int numDegrees = pEvent->delta() / 8;
  int numSteps = -(numDegrees / 15);

  _radius += (0.5 * numSteps * 4);

  updateModelViewMatrix();
  updateGL();
}

void GLWidget::mouseMoveEvent(QMouseEvent* pEvent)
{
  if ((pEvent->buttons() & Qt::MidButton) && _allowMoveCamera)
    {
      int mx = pEvent->x();
      int my = pEvent->y();
      int winHeight = size().height();
      int winWidth = size().width();

      double deltaPhi = (((mx - _lmx) / (double) winWidth) * 2 * PI);
      double deltaTheta = (((my - _lmy) / (double) winHeight) * 2 * PI);

      _phi += deltaPhi;

      if (-PI / 2 <= (deltaTheta + _theta) && (deltaTheta + _theta) <= PI / 2)
        _theta += deltaTheta;

      _lmy = my;
      _lmx = mx;

      updateModelViewMatrix();

      updateGL();
    }
  else if ((pEvent->buttons() & Qt::RightButton) && _allowMoveCamera)
    {
      int mx = pEvent->x();
      int my = pEvent->y();

      _deltaX -= (mx - _lmx) / 2;
      _deltaY += (my - _lmy) / 2;

      _lmy = my;
      _lmx = mx;

      updateModelViewMatrix();

      updateGL();
    }
}

vec3f GLWidget::mapWindowCoordinatesToModelCoordinates(int mx, int my)
{
  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixd(_modelViewMatrix);

  float z;
  glReadPixels(mx, _viewport[3] - my, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);

  double position[3];
  gluUnProject(mx, _viewport[3] - my, z, _modelViewMatrix, _projectionMatrix, _viewport,
               &position[0], &position[1], &position[2]);

  vec3f result;
  result[0] = position[0];
  result[1] = position[1];
  result[2] = position[2];

  return result;
}

void GLWidget::set3DViewHeight()
{
  _deltaX = _deltaY = 0;
  _radius = 128.0;
  _phi = 2 * PI - PI / 2.0;
  _theta = -PI / 4.0;

  updateModelViewMatrix();

  updateGL();
}

void GLWidget::set2DView()
{
  _deltaX = _deltaY = 0;
  _radius = 100.0;
  _phi = 1.5 * PI;
  _theta = 0.0;

  updateModelViewMatrix();

  updateGL();
}
