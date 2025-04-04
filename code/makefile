Compile=g++ -g -std=c++2a -O3 # improve run time e.g. 10s to 4s
r2_Ising: main.o r2_Ising.o
	$(Compile) -o r2_Ising main.o r2_Ising.o
main.o: main.cpp r2_Ising.h
	$(Compile) -o main.o -c main.cpp
r2_Ising.o: r2_Ising.cpp r2_Ising.h
	$(Compile) -o r2_Ising.o -c r2_Ising.cpp
clean:
	rm -f main.o r2_Ising.o