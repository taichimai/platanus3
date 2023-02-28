CXX =g++
FLAGS=-O3  -fopenmp -pthread
LIBS= -I./src  -I/usr/local/opt/libomp/include
MAIN=main.cpp
PROGRAM = platanus3

all: ${PROGRAM}

${PROGRAM}: ${MAIN}
		$(CXX) $(FLAGS) ${LIBS}  -o $@ $^ 

clean:
	rm -f ${PROGRAM} 
 
