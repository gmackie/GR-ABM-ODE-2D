TEMPLATE = app
TARGET = grviz-lung

# CONFIG += debug
CONFIG += debug_and_release
QT += core \
    gui \
    opengl
HEADERS += simulation/recruitmentlnode.h \
    simulation/recruitmentprob.h \
    simulation/recruitmentbase.h \
    scalardatasets/scalarcelldensitydataset.h \
    scalardatasets/scalartnfattrextmtb.h \
    scalardatasets/scalarattractantdataset.h \
    scalardatasets/scalaragentgridbase.h \
    gui/agentswidget.h \
    gui/statwidget.h \
    colormaps/blackwhite.h \
    colormaps/colormap.h \
    colormaps/coolwarm.h \
    colormaps/fire.h \
    colormaps/fixed.h \
    colormaps/greenred.h \
    colormaps/rainbow.h \
    datasets/dataset.h \
    datasets/grid.h \
    glyphs/glyph.h \
    glyphs/glypharrow.h \
    glyphs/glyphcone.h \
    glyphs/glyphhedgehog.h \
    glyphs/glyphtexture.h \
    gui/colormapwidget.h \
    gui/glwidget.h \
    gui/glwindow.h \
    gui/mainwindow.h \
    gui/paramwindow.h \
    scalardatasets/scalaragentgrid.h \
    scalardatasets/scalarccl2dataset.h \
    scalardatasets/scalarccl5dataset.h \
    scalardatasets/scalarcxcl9dataset.h \
    scalardatasets/scalardataset.h \
    scalardatasets/scalardivergencedataset.h \
    scalardatasets/scalarextmtbdataset.h \
    scalardatasets/scalargrid.h \
    scalardatasets/scalarintmtbdataset.h \
    scalardatasets/scalarnormalizer.h \
    scalardatasets/scalartnfdataset.h \
    simulation/agent.h \
    simulation/areatest.h \
    simulation/gr.h \
    simulation/grdiffusion.h \
    simulation/grdiffusionbtcs.h \
    simulation/grdiffusionftcs.h \
    simulation/grdiffusionwrongbtcs.h \
    simulation/grgrid.h \
    simulation/gridcell.h \
    simulation/grsimulation.h \
    simulation/grstat.h \
    simulation/macrophage.h \
    simulation/mtbtest.h \
    simulation/onlinestat.h \
    simulation/params.h \
    simulation/rand.h \
    simulation/tcell.h \
    simulation/tcytotoxic.h \
    simulation/tgamma.h \
    simulation/tregulatory.h \
    simulation/ttest.h \
    simulation/tinyxml/tinystr.h \
    simulation/tinyxml/tinyxml.h \
    vectordatasets/vector.h \
    vectordatasets/vectordataset.h \
    vectordatasets/vectorgradientdataset.h \
    vectordatasets/vectorgrid.h \
    visualization/agentsvisualization.h \
    visualization/glyphvisualization.h \
    visualization/heightplotvisualization.h \
    visualization/invisiblequadvisualisation.h \
    visualization/isolinesvisualization.h \
    visualization/smokevisualization.h \
    visualization/visualization.h \
    grviz.h \
    maininterface.h \
    simulation.h \
    snapshot.h
SOURCES += simulation/recruitmentlnode.cpp \
    simulation/recruitmentprob.cpp \
    simulation/recruitmentbase.cpp \
    gui/agentswidget.cpp \
    gui/statwidget.cpp \
    colormaps/blackwhite.cpp \
    colormaps/colormap.cpp \
    colormaps/coolwarm.cpp \
    colormaps/fire.cpp \
    colormaps/fixed.cpp \
    colormaps/greenred.cpp \
    colormaps/rainbow.cpp \
    gui/colormapwidget.cpp \
    gui/glwidget.cpp \
    gui/glwindow.cpp \
    gui/mainwindow.cpp \
    gui/paramwindow.cpp \
    scalardatasets/scalaragentgrid.cpp \
    scalardatasets/scalargrid.cpp \
    scalardatasets/scalarnormalizer.cpp \
    simulation/rand.cpp \
    simulation/mtbtest.cpp \
    simulation/ttest.cpp \
    simulation/grdiffusionbtcs.cpp \
    simulation/grdiffusionwrongbtcs.cpp \
    simulation/grdiffusionftcs.cpp \
    simulation/grdiffusion.cpp \
    simulation/areatest.cpp \
    simulation/onlinestat.cpp \
    simulation/grstat.cpp \
    simulation/gridcell.cpp \
    simulation/agent.cpp \
    simulation/grgrid.cpp \
    simulation/grsimulation.cpp \
    simulation/gr.cpp \
    simulation/macrophage.cpp \
    simulation/params.cpp \
    simulation/tcell.cpp \
    simulation/tcytotoxic.cpp \
    simulation/tgamma.cpp \
    simulation/tinyxml/tinystr.cpp \
    simulation/tinyxml/tinyxml.cpp \
    simulation/tinyxml/tinyxmlerror.cpp \
    simulation/tinyxml/tinyxmlparser.cpp \
    simulation/tregulatory.cpp \
    vectordatasets/vectorgrid.cpp \
    visualization/agentsvisualization.cpp \
    visualization/glyphvisualization.cpp \
    visualization/heightplotvisualization.cpp \
    visualization/invisiblequadvisualisation.cpp \
    visualization/isolinesvisualization.cpp \
    visualization/smokevisualization.cpp \
    main.cpp \
    maininterface.cpp \
    simulation.cpp \
    snapshot.cpp
FORMS += gui/agentswidget.ui \
    gui/statwidget.ui \
    gui/paramwindow.ui \
    gui/glwindow.ui \
    gui/mainwindow.ui
unix:LIBS += -lboost_program_options-mt
macx { 
    LIBS -= -lboost_program_options-mt
    LIBS += -lboost_program_options
}
win32 { 
    CONFIG -= flat
    LIBS += -Lc:\Qt\boost\lib \
        -llibboost_program_options-mgw44-mt-s
    INCLUDEPATH += c:\Qt\boost
}
debug { 
    DEFINES -= NDEBUG
    OBJECTS_DIR = debug
}
release { 
    DEFINES += NDEBUG
    OBJECTS_DIR = release
}
