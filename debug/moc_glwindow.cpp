/****************************************************************************
** Meta object code from reading C++ file 'glwindow.h'
**
** Created: Tue Jun 7 11:46:13 2011
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../gui/glwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'glwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GLWindow[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
      17,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // signals: signature, parameters, type, tag, flags
      18,   10,    9,    9, 0x05,
      43,    9,    9,    9, 0x05,
      55,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      84,   67,    9,    9, 0x0a,
     107,    9,    9,    9, 0x0a,
     130,  122,    9,    9, 0x0a,
     164,    9,    9,    9, 0x0a,
     183,    9,    9,    9, 0x0a,
     209,   10,    9,    9, 0x0a,
     229,    9,    9,    9, 0x0a,
     249,    9,    9,    9, 0x0a,
     270,    9,    9,    9, 0x0a,
     288,    9,    9,    9, 0x0a,
     308,    9,    9,    9, 0x0a,
     320,    9,    9,    9, 0x0a,
     338,  332,    9,    9, 0x0a,
     357,  332,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GLWindow[] = {
    "GLWindow\0\0row,col\0updateSelection(int,int)\0"
    "set2DView()\0set3DView()\0pCurrentColorMap\0"
    "setColorMap(ColorMap*)\0updateWindow()\0"
    "min,max\0updateColorMapLabels(float,float)\0"
    "toggleFullScreen()\0updateSelectedCellStats()\0"
    "selectCell(int,int)\0moveSelectionLeft()\0"
    "moveSelectionRight()\0moveSelectionUp()\0"
    "moveSelectionDown()\0visualize()\0"
    "printText()\0value\0setPrintTime(bool)\0"
    "setPrintOutcome(bool)\0"
};

const QMetaObject GLWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_GLWindow,
      qt_meta_data_GLWindow, 0 }
};

const QMetaObject *GLWindow::metaObject() const
{
    return &staticMetaObject;
}

void *GLWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GLWindow))
        return static_cast<void*>(const_cast< GLWindow*>(this));
    return QWidget::qt_metacast(_clname);
}

int GLWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateSelection((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 1: set2DView(); break;
        case 2: set3DView(); break;
        case 3: setColorMap((*reinterpret_cast< ColorMap*(*)>(_a[1]))); break;
        case 4: updateWindow(); break;
        case 5: updateColorMapLabels((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        case 6: toggleFullScreen(); break;
        case 7: updateSelectedCellStats(); break;
        case 8: selectCell((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 9: moveSelectionLeft(); break;
        case 10: moveSelectionRight(); break;
        case 11: moveSelectionUp(); break;
        case 12: moveSelectionDown(); break;
        case 13: visualize(); break;
        case 14: printText(); break;
        case 15: setPrintTime((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: setPrintOutcome((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 17;
    }
    return _id;
}

// SIGNAL 0
void GLWindow::updateSelection(int _t1, int _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GLWindow::set2DView()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void GLWindow::set3DView()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}
QT_END_MOC_NAMESPACE
