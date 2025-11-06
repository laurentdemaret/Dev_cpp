TEMPLATE = lib
CONFIG += staticlib  # ou shared pour une librairie dynamique

# Chemin absolu
#DESTDIR = /Users/laurentdemaret/Dev_Cpp/MyLibs/Libs
DESTDIR = ../../../Libs
TARGET = algotools

SOURCES += \
	 utils.cpp \
	 matrix.cpp \
         symbolic.cpp \
         numeric.cpp \
         nl_opt.cpp \

HEADERS += \
    utils.h \
    matrix.h \
    symbolic.h \
    numeric.h \
    nl_opt.h \


# GiNaC via pkg-config
unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += ginac
}

unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += fftw3
}

