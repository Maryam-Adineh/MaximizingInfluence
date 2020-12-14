CC=g++
CFLAGS=-I.
DEPS = limit.h Graph.h IndependCascade.h AAPC.h EAPC.h

OBJ = Graph.o main.o IndependCascade.o AAPC.o EAPC.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

InfluenceMaximization: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
