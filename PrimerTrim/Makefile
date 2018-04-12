# Makefile for CellDivide
CC = gcc
CFLAGS = -std=c99 -O2 -g
LIBS = -lz -lpthread
OBJECT = fastq.o hash.o index.o parse.o query.o dynamic.o main.o

primtrim-1.3.0: $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@ $(LIBS)

fastq.o: fastq.h utils.h
hash.o: hash.h
index.o: hash.h utils.h
parse.o: fastq.h utils.h
query.o: query.h
dynamic.o: dynamic.h
main.o: query.h

.PHONY : clean

clean:
	rm -f $(OBJECT)

