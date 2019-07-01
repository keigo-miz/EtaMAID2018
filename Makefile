CC = c++
ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --libs`
DEBUG = -Wall

bin/ana: src/main.cc lib/hmaid.o
	$(CC) $(ROOTLIBS) $(ROOTFLAGS) $(DEBUG) -lMathMore src/main.cc lib/hmaid.o -o bin/ana

lib/hmaid.o: src/hmaid.cc src/hmaid.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/hmaid.cc -o lib/hmaid.o

clean:
	rm bin/ana
	rm lib/*.o
