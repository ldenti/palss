CC=gcc
CFLAGS=-g -O3 -Wall -Wno-unused-function -I..
LDFLAGS=-lz

.PHONY: all clean

all: extract_subgraph

extract_subgraph: extract_subgraph.o ../misc.o ../graph.o ../path.o ../segment.o
	@echo "* Linking $<"
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f test test.o ../misc.o ../graph.o ../path.o ../segment.o
