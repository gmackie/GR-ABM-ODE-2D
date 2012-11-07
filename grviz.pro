TEMPLATE = app
TARGET = grviz-lung

# CONFIG += debug

macx { INSTALLBASE = /opt/local }
else:unix { INSTALLBASE = /usr }

CONFIG *= debug_and_release qt opengl qwt
QT += core gui opengl

HEADERS += scalardatasets/scalartotmtbdataset.h \
    simulation/params.h \
    simulation/paramsbase.h \
    simulation/grdiffusionadeswap.h \
    simulation/grdiffusionftcsswap.h \
    simulation/grsimulationgrid.h \
    simulation/tinyxml/tinyxml.h \
    simulation/tinyxml/tinystr.h \
    simulation/recruitmentlnodepure.h \
    simulation/recruitmentlnodeproxy.h \    
    simulation/recruitmentlnode.h \
    simulation/recruitmentprob.h \
    simulation/recruitmentbase.h \
    simulation/stat.h \
    simulation/stat.def \
    scalardatasets/scalarcelldensitydataset.h \
    scalardatasets/scalartnfattrextmtb.h \
    scalardatasets/scalarattractantdataset.h \
    scalardatasets/scalaragentgridbase.h \
    gui/agentswidget.h \
    gui/graphcontroller.h \
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
    scalardatasets/scalaril10dataset.h \
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
    simulation/grdiffusionwrongbtcs.h \
    simulation/grgrid.h \
    simulation/grsimulation.h \
    simulation/macrophage.h \
    simulation/mtbtest.h \
    simulation/onlinestat.h \
    simulation/rand.h \
    simulation/tcell.h \
    simulation/tcytotoxic.h \
    simulation/tgamma.h \
    simulation/tregulatory.h \
    simulation/ttest.h \
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
    snapshot.h \
    gui/graphviewer.h \
    gui/agenthistogram.h \
    simulation/numericalMethods.h \
    scalardatasets/scalarkillingsdataset.h \
    scalardatasets/scalarindexeddataset.h
SOURCES += \
    simulation/params.cpp \
    simulation/paramsbase.cpp \
    simulation/grdiffusionadeswap.cpp \
    simulation/grdiffusionftcsswap.cpp \
    simulation/grsimulationgrid.cpp \
    simulation/tinyxml/tinyxmlparser.cpp \
    simulation/tinyxml/tinyxmlerror.cpp \
    simulation/tinyxml/tinyxml.cpp \
    simulation/tinyxml/tinystr.cpp \
    simulation/recruitmentlnodepure.cpp \
    simulation/recruitmentlnodeproxy.cpp \
    simulation/recruitmentlnode.cpp \
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
    gui/graphcontroller.cpp \
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
    simulation/grdiffusion.cpp \
    simulation/areatest.cpp \
    simulation/onlinestat.cpp \
    simulation/agent.cpp \
    simulation/grgrid.cpp \
    simulation/grsimulation.cpp \
    simulation/gr.cpp \
    simulation/macrophage.cpp \
    simulation/tcell.cpp \
    simulation/tcytotoxic.cpp \
    simulation/tgamma.cpp \
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
    snapshot.cpp \
    gui/graphviewer.cpp \
    gui/agenthistogram.cpp
FORMS += gui/agentswidget.ui \
    gui/statwidget.ui \
    gui/paramwindow.ui \
    gui/glwindow.ui \
    gui/mainwindow.ui \
    gui/graphviewer.ui \
    gui/agenthistogram.ui
QMAKE_CXXFLAGS_RELEASE -= -O2 -Wno-strict-aliasing
QMAKE_CXXFLAGS_RELEASE += -O3 -Wno-strict-aliasing

BOOST_PREFIX=$$(BOOST_PREFIX) #Grab the environment variable
isEmpty(BOOST_PREFIX) {
  warning("BOOST_PREFIX not defined, defaulting to /usr")
  BOOST_PREFIX=/usr
}
INCLUDEPATH *= $$quote($(BOOST_PREFIX)/include)

BOOST_LIBS=boost_program_options boost_iostreams boost_serialization
LIBS *= -L$${BOOST_PREFIX}/lib
for(lib, BOOST_LIBS) {
  !exists($${BOOST_PREFIX}/lib/*$${lib}*):error("Unable to find boost library:" $${lib})
  exists($${BOOST_PREFIX}/lib/*$${lib}*-mt*):lib=$${lib}-mt
  LIBS *= -l$${lib}
}

unix:!macx:LIBS *= -lGLU

exists( .git/ ) {
      VERSION = $$quote($$system(git svn find-rev HEAD))
} else : exists( .svn/ ) {
      VERSION = $$quote($$system(svn info | awk \'/^Last Changed Rev:/ {print $4}\'))
} else {
      VERSION = "Unknown"
}

DEFINES += SVN_VERSION=\"$$VERSION\"
DEFINES += TIXML_USE_STL

macx { 
    # Have qmake create make files that put the executable
    # in a file in the build directory, rather than in a Mac
    # application bundle. Ex. put the executable in grviz-lung
    # rather than in grviz-lung.app/Contents/MacOS/grviz-lung.
    CONFIG -= app_bundle
    LIBS *= -framework OpenGL
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

# This is for g++ code only.  Use this if you are compiling the executable for ONE MACHINE only
# DO NOT move the generated executable to another machine with this option on
