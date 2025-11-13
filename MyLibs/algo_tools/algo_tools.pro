TEMPLATE = lib
CONFIG += staticlib  # ou shared pour une librairie dynamique

# Chemin absolu
#DESTDIR = /Users/laurentdemaret/Dev_Cpp/MyLibs/Libs

message(Le chemin du projet algo_tools est: $$_PRO_FILE_PWD_)


# Choisir un sous-dossier selon le système
macx {
DESTDIR = $$_PRO_FILE_PWD_/../Libs/mac
}
else:unix {
    DESTDIR= $$_PRO_FILE_PWD_/../Libs/linux
#    DESTDIR     = $$BUILDDIR/linux/bin
}

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

