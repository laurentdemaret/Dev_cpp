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
#MY_BIN_DIR = $$system(echo $$ENV{MY_BIN_DIR})
#message(Le chemin de lexecutable basic est: $$MY_BIN_DIR)
#unix:!macx {
#    DESTDIR = /home/demaret/bin
#}
#macx {
#    DESTDIR = /Users/laurentdemaret/bin
#}
#DESTDIR = $$(MY_BIN_DIR)
DESTDIR = ~/bin

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

macx {
LIBS += -L../MyLibs/Libs/mac -lalgotools
} else:unix
{
LIBS += -L../MyLibs/Libs/linux -lalgotools
}


message(Le chemin absolu de include path est : $$absolute_path($$INCLUDEPATH, $$_PRO_FILE_PWD_))
message(Le chemin absolu de Libs est : $$absolute_path($$LIBS))

# Copy de l'éxécutable dans ~/bin
#target.path = /Users/laurentdemaret/bin
#INSTALLS += target

# Copier le binaire après compilation
#QMAKE_POST_LINK += cp build/Desktop_Qt_6_9_1_clang_64_bit-Debug/basic $$DESTDIR/basic
