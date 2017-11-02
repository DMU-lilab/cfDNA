# Makefile for CellDivide
CC = gcc
CFLAGS = -std=c99 -O2
LIBS = -lz -lpthread
OBJECT = fastq.o hash.o index.o parse.o query.o main.o

primtrim-1.2: $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@ $(LIBS)

fastq.o: fastq.h
hash.o: hash.h
index.o: hash.h
parse.o: fastq.h
query.o: query.h
main.o: query.h

.PHONY : clean

clean:
	rm -f $(OBJECT)

