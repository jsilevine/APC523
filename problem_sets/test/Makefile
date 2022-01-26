main: main.o fact.o
	g++ -o main main.o fact.o

main.o: main.cpp fact.h
	g++ -c main.cpp

fact.o: fact.cpp fact.h
	g++ -c fact.cpp

clean:
	rm -f *.o
