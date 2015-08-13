CXX = mpicxx
TARGET = timetabling.ga.uk.2
CXXFLAGS = -Wall -ansi -O3 -fopenmp # -mpe=mpilog

OBJS = Problem.o Solution.o Random.o Timer.o Control.o util.o jsoncpp.o

all: ${TARGET}

timetabling.ga.uk.2: ga.cpp $(OBJS)
	${CXX} ${CXXFLAGS} -o $@ $^

clean:
	@rm -f *~ *.o ${TARGET} core DEADJOE
