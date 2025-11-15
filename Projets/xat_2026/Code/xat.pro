TEMPLATE = app
CONFIG += console c++23
CONFIG -= qt
CONFIG -= app_bundle

# Nom de l'exécutable
#TARGET = xat
DESTDIR = ~/bin
#TARGET = xat

# Dossier de build et dossier où installer l'exécutable (optionnel)
#DESTDIR = $$HOME/bin

# Liste des fichiers sources
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
    Test.cpp

# Liste des fichiers header
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
    Test.h

# Options de compilation
QMAKE_CXXFLAGS += -Wall -W -Wno-sign-compare -Wno-unused-label

# Options de release vs debug
CONFIG(debug, debug|release) {
    QMAKE_CXXFLAGS += -g
} else {
    QMAKE_CXXFLAGS += -O3 -DNDEBUG
}

# Libs optionnelles si nécessaire
LIBS += -lm

# Copy de l'éxécutable dans ~/bin
#target.path = /Users/laurentdemaret/bin
#INSTALLS += xat
