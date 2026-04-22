# Project: xat
# Makefile created by Dev-C++ 4.9.9.2, and worked out to compile the xat project
# Author: Laurent Demaret

# ----------------------------machine dependent file
MODE = optim

# extern library and program variables
# compiler and linker variable 

#LIBS = -L/usr/X11R6/lib -lX11
LIBS = -lm    
CPP  = g++

# -Wall -Wno-sign-compare 
CPPFLAGS = -g -Wall -W -Wno-sign-compare -Wno-unused-label -MMD ${INCLUDES}

ifeq ($(MODE),debug)
else
	CPPFLAGS += -O3 -DNDEBUG
endif

#LDFLAGS = ${LIBS}

# rules
%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $*.cpp


OBJ  =  main.o IO.o  PGM.o coding.o Thinning.o fib.o bitstr.o arithcoder.o qsmodel.o rangecod.o\
matrix.o Edge.o Point2D.o PointGrid.o Triangle.o Triangulation.o Test.o $(RES) 
LINKOBJ  =  main.o IO.o PGM.o coding.o Thinning.o fib.o bitstr.o  arithcoder.o qsmodel.o rangecod.o\
matrix.o Edge.o Point2D.o PointGrid.o  Triangle.o Triangulation.o Test.o $(RES)
#BIN  = ~/bin/xat
RM = rm -f

ifndef MAIN
MAIN = main.cpp
endif

ifndef BIN
BIN = xat
endif


.PHONY: all all-before all-after clean clean-custom

all: all-before xat all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)
	${RM} *.d

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "xat" $(LIBS)

main.o: main.cpp
	$(CPP) -c main.cpp -o main.o $(CPPFLAGS)

IO.o: IO.cpp
	$(CPP) -c IO.cpp -o IO.o $(CPPFLAGS)

PGM.o: PGM.cpp
	$(CPP) -c PGM.cpp -o PGM.o $(CPPFLAGS)

coding.o: coding.cpp
	$(CPP) -c coding.cpp -o coding.o $(CPPFLAGS)

Thinning.o: Thinning.cpp
	$(CPP) -c Thinning.cpp -o Thinning.o $(CPPFLAGS)

fib.o: fib.cpp
	$(CPP) -c fib.cpp -o fib.o $(CPPFLAGS)

Logger.o: Logger.cpp
	$(CPP) -c Logger.cpp -o Logger.o $(CPPFLAGS)

bitstr.o: bitstr.cpp
	$(CPP) -c bitstr.cpp -o bitstr.o $(CPPFLAGS)

arithcoder.o: arithcoder.cpp
	$(CPP) -c arithcoder.cpp -o arithcoder.o $(CPPFLAGS)

qsmodel.o: qsmodel.cpp
	$(CPP) -c qsmodel.cpp -o qsmodel.o $(CPPFLAGS)

rangecod.o: rangecod.cpp
	$(CPP) -c rangecod.cpp -o rangecod.o $(CPPFLAGS)

matrix.o: matrix.cpp
	$(CPP) -c matrix.cpp -o matrix.o $(CPPFLAGS)

Point2D.o: Point2D.cpp
	$(CPP) -c Point2D.cpp -o Point2D.o $(CPPFLAGS)

Edge.o: Edge.cpp
	$(CPP) -c Edge.cpp -o Edge.o $(CPPFLAGS)

PointGrid.o: PointGrid.cpp
	$(CPP) -c PointGrid.cpp -o PointGrid.o $(CPPFLAGS)

Triangle.o: Triangle.cpp
	$(CPP) -c Triangle.cpp -o Triangle.o $(CPPFLAGS)

Triangulation.o: Triangulation.cpp
	$(CPP) -c Triangulation.cpp -o Triangulation.o $(CPPFLAGS)

Test.o: Test.cpp
	$(CPP) -c Test.cpp -o Test.o $(CPPFLAGS)







