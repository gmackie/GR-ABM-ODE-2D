/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon Nov 22 16:26:05 2010
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../gui/mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
      50,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,
      23,   11,   11,   11, 0x05,
      59,   49,   11,   11, 0x05,
      93,   85,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
     127,   11,   11,   11, 0x0a,
     145,   11,   11,   11, 0x0a,
     164,  156,   11,   11, 0x0a,
     184,  156,   11,   11, 0x0a,
     203,  156,   11,   11, 0x0a,
     223,  156,   11,   11, 0x0a,
     245,  156,   11,   11, 0x0a,
     275,  269,   11,   11, 0x0a,
     296,  269,   11,   11, 0x0a,
     317,  269,   11,   11, 0x0a,
     343,  269,   11,   11, 0x0a,
     368,  269,   11,   11, 0x0a,
     395,  269,   11,   11, 0x0a,
     412,  269,   11,   11, 0x0a,
     428,   11,   11,   11, 0x0a,
     445,   11,   11,   11, 0x0a,
     468,   11,   11,   11, 0x0a,
     489,  269,   11,   11, 0x0a,
     520,   11,   11,   11, 0x0a,
     541,   11,   11,   11, 0x0a,
     558,   11,   11,   11, 0x0a,
     581,   11,   11,   11, 0x0a,
     606,   11,   11,   11, 0x0a,
     633,  269,   11,   11, 0x0a,
     664,  269,   11,   11, 0x0a,
     695,  269,   11,   11, 0x0a,
     713,  156,   11,   11, 0x0a,
     736,   11,   11,   11, 0x0a,
     748,   11,   11,   11, 0x0a,
     761,   11,   11,   11, 0x0a,
     783,   11,   11,   11, 0x0a,
     806,   11,   11,   11, 0x0a,
     831,   11,   11,   11, 0x0a,
     858,   11,   11,   11, 0x0a,
     876,   11,   11,   11, 0x0a,
     894,   11,   11,   11, 0x0a,
     912,   11,   11,   11, 0x0a,
     930,   11,   11,   11, 0x0a,
     948,  156,   11,   11, 0x0a,
     970,   11,   11,   11, 0x0a,
     996,   11,   11,   11, 0x0a,
    1020,   11,   11,   11, 0x0a,
    1046,   11,   11,   11, 0x0a,
    1053,   11,   11,   11, 0x0a,
    1065,   11,   11,   11, 0x0a,
    1077,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0updateGL()\0updateSelectedCellStats()\0"
    "pColorMap\0updateColorMap(ColorMap*)\0"
    "min,max\0updateColorMapLabels(float,float)\0"
    "toggleAnimation()\0setBlend()\0checked\0"
    "setDrawAgents(bool)\0setDrawSmoke(bool)\0"
    "setDrawGlyphs(bool)\0setDrawIsolines(bool)\0"
    "setDrawHeightPlot(bool)\0value\0"
    "setGlyphScaling(int)\0setColorMap(QString)\0"
    "setMappingMethod(QString)\0"
    "setSmokeDataset(QString)\0"
    "setColorMapSource(QString)\0updateDelay(int)\0"
    "updateSeed(int)\0updateColorMap()\0"
    "updateColorMapLabels()\0updateHeightMinMax()\0"
    "updateDiffusionMethod(QString)\0"
    "updateStopCriteria()\0rescaleRequest()\0"
    "rescaleHeightRequest()\0updateIsolinesSettings()\0"
    "updateHeightPlotSettings()\0"
    "setGlyphVectorDataset(QString)\0"
    "setGlyphScalarDataset(QString)\0"
    "setGlyph(QString)\0setGlyphClamping(bool)\0"
    "dumpGrids()\0showParams()\0selectSmokeColorMap()\0"
    "selectGlyphsColorMap()\0selectIsolinesColorMap()\0"
    "selectHeightPlotColorMap()\0selectColorMap0()\0"
    "selectColorMap1()\0selectColorMap2()\0"
    "selectColorMap3()\0selectColorMap4()\0"
    "setEnableOutput(bool)\0updateGranulomaSettings()\0"
    "updateOutcomeSettings()\0"
    "updateOutcomeParameters()\0stop()\0"
    "loadState()\0saveState()\0takePictureSnapshot()\0"
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, 0 }
};

const QMetaObject *MainWindow::metaObject() const
{
    return &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateGL(); break;
        case 1: updateSelectedCellStats(); break;
        case 2: updateColorMap((*reinterpret_cast< ColorMap*(*)>(_a[1]))); break;
        case 3: updateColorMapLabels((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        case 4: toggleAnimation(); break;
        case 5: setBlend(); break;
        case 6: setDrawAgents((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: setDrawSmoke((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: setDrawGlyphs((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: setDrawIsolines((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: setDrawHeightPlot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: setGlyphScaling((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: setColorMap((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 13: setMappingMethod((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: setSmokeDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: setColorMapSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 16: updateDelay((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: updateSeed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: updateColorMap(); break;
        case 19: updateColorMapLabels(); break;
        case 20: updateHeightMinMax(); break;
        case 21: updateDiffusionMethod((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 22: updateStopCriteria(); break;
        case 23: rescaleRequest(); break;
        case 24: rescaleHeightRequest(); break;
        case 25: updateIsolinesSettings(); break;
        case 26: updateHeightPlotSettings(); break;
        case 27: setGlyphVectorDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 28: setGlyphScalarDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 29: setGlyph((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 30: setGlyphClamping((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 31: dumpGrids(); break;
        case 32: showParams(); break;
        case 33: selectSmokeColorMap(); break;
        case 34: selectGlyphsColorMap(); break;
        case 35: selectIsolinesColorMap(); break;
        case 36: selectHeightPlotColorMap(); break;
        case 37: selectColorMap0(); break;
        case 38: selectColorMap1(); break;
        case 39: selectColorMap2(); break;
        case 40: selectColorMap3(); break;
        case 41: selectColorMap4(); break;
        case 42: setEnableOutput((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 43: updateGranulomaSettings(); break;
        case 44: updateOutcomeSettings(); break;
        case 45: updateOutcomeParameters(); break;
        case 46: stop(); break;
        case 47: loadState(); break;
        case 48: saveState(); break;
        case 49: takePictureSnapshot(); break;
        default: ;
        }
        _id -= 50;
    }
    return _id;
}

// SIGNAL 0
void MainWindow::updateGL()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void MainWindow::updateSelectedCellStats()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void MainWindow::updateColorMap(ColorMap * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void MainWindow::updateColorMapLabels(float _t1, float _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
QT_END_MOC_NAMESPACE
