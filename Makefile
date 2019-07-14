CC = c++
ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --libs`
DEBUG = -Wall

bin/ana: src/main.cc lib/resonance.o lib/born.o lib/regge.o lib/manager.o
	$(CC) $(ROOTLIBS) $(ROOTFLAGS) $(DEBUG) -lMathMore src/main.cc lib/resonance.o lib/born.o lib/regge.o lib/manager.o -o bin/ana

lib/resonance.o: src/resonance.cc src/resonance.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/resonance.cc -o lib/resonance.o

lib/born.o: src/born.cc src/born.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/born.cc -o lib/born.o

lib/regge.o: src/regge.cc src/regge.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/regge.cc -o lib/regge.o

lib/manager.o: src/manager.cc src/manager.h
	$(CC) $(ROOTFLAGS) $(DEBUG) -c src/manager.cc -o lib/manager.o

clean:
	rm bin/ana
	rm lib/*.o
