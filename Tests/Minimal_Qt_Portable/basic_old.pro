TEMPLATE = app
CONFIG += console c++23
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += sdk_no_version_check

message(Le chemin du projet basic est: $$_PRO_FILE_PWD_)
# Nom de l'éxécutable (binaire)
# TARGET = basic
# TARGET=~/bin/test2024_minimal

# Définir le chemin de destination
# DESTDIR = /Users/laurentdemaret/bin
# à des fins de portabilité
# pour une utilisation sous un nouvel environnement
DESTDIR = $$(MY_BIN_DIR)

# Dossier racine de build relatif au projet
BUILDDIR = $$_PRO_FILE_PWD_/build

# Choisir un sous-dossier selon le système
macx {
    OBJECTS_DIR = $$BUILDDIR/mac/obj
#    DESTDIR     = $$BUILDDIR/mac/bin
} else:unix {
    OBJECTS_DIR = $$BUILDDIR/linux/obj
#    DESTDIR     = $$BUILDDIR/linux/bin
}


SOURCES += \
      Edge2D.cpp \
      Point2D.cpp \
      main.cpp \
      test_functions.cpp

HEADERS += \
    Edge2D.h \
    Point2D.h \
    test_function.h

#INCLUDEPATH += /Users/laurentdemaret/Dev_Cpp/MyLibs/algo_tools
#LIBS += -L/Users/laurentdemaret/Dev_Cpp/MyLibs/Libs -lalgotools

INCLUDEPATH +=../../MyLibs/algo_tools
LIBS += -L../../MyLibs/Libs -lalgotools


# Copy de l'éxécutable dans ~/bin
#target.path = /Users/laurentdemaret/bin
#INSTALLS += target

# Copier le binaire après compilation
#QMAKE_POST_LINK += cp build/Desktop_Qt_6_9_1_clang_64_bit-Debug/basic $$DESTDIR/basic
