/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Tue Jun 7 11:46:14 2011
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
      51,   12, // methods
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
     146,   11,   11,   11, 0x0a,
     164,   11,   11,   11, 0x0a,
     183,  175,   11,   11, 0x0a,
     203,  175,   11,   11, 0x0a,
     222,  175,   11,   11, 0x0a,
     242,  175,   11,   11, 0x0a,
     264,  175,   11,   11, 0x0a,
     294,  288,   11,   11, 0x0a,
     315,  288,   11,   11, 0x0a,
     336,  288,   11,   11, 0x0a,
     362,  288,   11,   11, 0x0a,
     387,  288,   11,   11, 0x0a,
     414,  288,   11,   11, 0x0a,
     431,  288,   11,   11, 0x0a,
     447,   11,   11,   11, 0x0a,
     464,   11,   11,   11, 0x0a,
     487,   11,   11,   11, 0x0a,
     508,  288,   11,   11, 0x0a,
     539,   11,   11,   11, 0x0a,
     560,   11,   11,   11, 0x0a,
     577,   11,   11,   11, 0x0a,
     600,   11,   11,   11, 0x0a,
     625,   11,   11,   11, 0x0a,
     652,  288,   11,   11, 0x0a,
     683,  288,   11,   11, 0x0a,
     714,  288,   11,   11, 0x0a,
     732,  175,   11,   11, 0x0a,
     755,   11,   11,   11, 0x0a,
     767,   11,   11,   11, 0x0a,
     780,   11,   11,   11, 0x0a,
     802,   11,   11,   11, 0x0a,
     825,   11,   11,   11, 0x0a,
     850,   11,   11,   11, 0x0a,
     877,   11,   11,   11, 0x0a,
     895,   11,   11,   11, 0x0a,
     913,   11,   11,   11, 0x0a,
     931,   11,   11,   11, 0x0a,
     949,   11,   11,   11, 0x0a,
     967,  175,   11,   11, 0x0a,
     989,   11,   11,   11, 0x0a,
    1015,   11,   11,   11, 0x0a,
    1039,   11,   11,   11, 0x0a,
    1065,   11,   11,   11, 0x0a,
    1072,   11,   11,   11, 0x0a,
    1084,   11,   11,   11, 0x0a,
    1096,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0updateGL()\0updateSelectedCellStats()\0"
    "pColorMap\0updateColorMap(ColorMap*)\0"
    "min,max\0updateColorMapLabels(float,float)\0"
    "updateTimeBox(int)\0toggleAnimation()\0"
    "setBlend()\0checked\0setDrawAgents(bool)\0"
    "setDrawSmoke(bool)\0setDrawGlyphs(bool)\0"
    "setDrawIsolines(bool)\0setDrawHeightPlot(bool)\0"
    "value\0setGlyphScaling(int)\0"
    "setColorMap(QString)\0setMappingMethod(QString)\0"
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
        case 4: updateTimeBox((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: toggleAnimation(); break;
        case 6: setBlend(); break;
        case 7: setDrawAgents((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: setDrawSmoke((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: setDrawGlyphs((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: setDrawIsolines((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: setDrawHeightPlot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: setGlyphScaling((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: setColorMap((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: setMappingMethod((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: setSmokeDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 16: setColorMapSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 17: updateDelay((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: updateSeed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 19: updateColorMap(); break;
        case 20: updateColorMapLabels(); break;
        case 21: updateHeightMinMax(); break;
        case 22: updateDiffusionMethod((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 23: updateStopCriteria(); break;
        case 24: rescaleRequest(); break;
        case 25: rescaleHeightRequest(); break;
        case 26: updateIsolinesSettings(); break;
        case 27: updateHeightPlotSettings(); break;
        case 28: setGlyphVectorDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 29: setGlyphScalarDataset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 30: setGlyph((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 31: setGlyphClamping((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 32: dumpGrids(); break;
        case 33: showParams(); break;
        case 34: selectSmokeColorMap(); break;
        case 35: selectGlyphsColorMap(); break;
        case 36: selectIsolinesColorMap(); break;
        case 37: selectHeightPlotColorMap(); break;
        case 38: selectColorMap0(); break;
        case 39: selectColorMap1(); break;
        case 40: selectColorMap2(); break;
        case 41: selectColorMap3(); break;
        case 42: selectColorMap4(); break;
        case 43: setEnableOutput((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 44: updateGranulomaSettings(); break;
        case 45: updateOutcomeSettings(); break;
        case 46: updateOutcomeParameters(); break;
        case 47: stop(); break;
        case 48: loadState(); break;
        case 49: saveState(); break;
        case 50: takePictureSnapshot(); break;
        default: ;
        }
        _id -= 51;
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
