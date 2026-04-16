# xat.pro
TEMPLATE = app
macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 14.0
}
CONFIG -= qt
CONFIG += console c++23
CONFIG += sdk_no_version_check
# Nom et emplacement de l'exécutable
#DESTDIR = ~/bin
TARGET = xat
#message($$(HOME))
DESTDIR = $$(HOME)/bin

# FORCE le chemin de sortie (clé)
#QMAKE_LFLAGS += -Wl,-o,$$DESTDIR/$$TARGET

CONFIG -= app_bundle

# Fichiers sources
SOURCES += \
    main.cpp \
    IO.cpp \
    PGM.cpp \
    coding.cpp \
    Thinning.cpp \
    fib.cpp \
    bitstr.cpp \
    arithcoder.cpp \
    qsmodel.cpp \
    rangecod.cpp \
    matrix.cpp \
    Edge.cpp \
    Point2D.cpp \
    PointGrid.cpp \
    Triangle.cpp \
    Triangulation.cpp \
    VoronoiCell.cpp \
    Test.cpp

# Fichiers headers
HEADERS += \
    IO.h \
    PGM.h \
    coding.h \
    Thinning.h \
    fib.h \
    bitstr.h \
    arithcoder.h \
    qsmodel.h \
    rangecod.h \
    matrix.h \
    Edge.h \
    Point2D.h \
    PointGrid.h \
    Triangle.h \
    Triangulation.h \
    VoronoiCell.h \
    Test.h

# Options de compilation
QMAKE_CXXFLAGS += -Wall -W -Wno-sign-compare -Wno-unused-label

# Options debug / release
CONFIG(debug, debug|release) {
    QMAKE_CXXFLAGS += -g
} else {
    QMAKE_CXXFLAGS += -O3 -DNDEBUG
}

# Librairies
LIBS += -lm
