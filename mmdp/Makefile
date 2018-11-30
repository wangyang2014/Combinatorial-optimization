FLAGS=-O3
all:    init.o basicRoutines.o algo.o bestSolAndTexOutput.o
	g++ -Wall -o mmdp ./src/main.cpp init.o basicRoutines.o algo.o bestSolAndTexOutput.o $(FLAGS)  
%.o: src/%.cpp src/%.h
	g++ -c $< $(FLAGS)  
###Above target = same as list of below targets
#init.o: src/init.cpp src/init.h
#	g++ -c src/init.cpp $(FLAGS)
#basicRoutines.o: src/basicRoutines.cpp src/basicRoutines.h
#	g++ -c src/basicRoutines.cpp $(FLAGS)
clean: 
	rm *.o mmdp

