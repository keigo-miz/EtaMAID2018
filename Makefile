CC = c++
ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --libs`
DEBUG = -Wall

bin/ana: src/main.cc lib/resonance.o
	$(CC) $(ROOTLIBS) $(ROOTFLAGS) $(DEBUG) -lMathMore src/main.cc lib/resonance.o -o bin/ana

lib/resonance.o: src/resonance.cc src/resonance.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/resonance.cc -o lib/resonance.o

clean:
	rm bin/ana
	rm lib/*.o
