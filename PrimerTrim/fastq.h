#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <zlib.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <linux/limits.h>

#define FQLINE 256    // maxmum length for the fastq line
#define BUFNUM 2048   // default buf number for the fastq cache
#define BUFLINE 1024  // maxmum bufer line for file reading


/* fastq type [ read1 or read2 ]*/
enum TYPE { READ1 = 0, READ2 = 1 };

/*! @typedef arg_t
 @abstract structure for the comand line args
 @field help         [0 | 1] if 1: print the help infomation
 @field kmer         the length for primer index
 @field mismatch     mismatch allowed between the primer and sequence
 @field ampfile      the path of amplicon file [.csv]
 @field read1        the path of fastq file of R1 [_R1.fq]
 @field read2        the path of fastq file of R2 [_R2.fq]
 @field outdir       the path of the output directory [~/]
*/
typedef struct __arg_t {
    int help;
    int kmer;
    int mismatch;
    char ampfile[PATH_MAX];
    char read1[PATH_MAX];
    char read2[PATH_MAX];
    char outdir[PATH_MAX];
} arg_t;


/*! @typedef read_t
 @abstract structure for one fastq read group
 @field name         the seqname for the read
 @field seq          the sequence for the read
 @field mark         the fixed marker['+']
 @field qual         the quality for the read
*/
typedef struct __read_t {
    char name[FQLINE];
    char seq[FQLINE];
    char mark[FQLINE];
    char qual[FQLINE];
} read_t;


/*! @typedef fastq_t
 @abstract structure for the single fastq file
 @field in, out      the input and output pointer for the fastq [FILE *]
 @field read         a fastq read group
 @field bufnum       the bufnum for the cache
 @field cache        the buffer for the fastq cache
 @field inname       the input file name
 @field outname      the output file name
*/
typedef struct __fastq_t {
    gzFile in;
    FILE *out;
    read_t read;
    int bufnum;
    read_t cache[BUFNUM];
    char inname[PATH_MAX];
    char outname[PATH_MAX];
} fastq_t;


/* prototype function */
void FastqInit(fastq_t *, arg_t *, int);
int FastqRead(fastq_t *);
void FastqWrite(fastq_t *);

/* parse the comand line parmerters*/
void Usage(void);
arg_t *ParseOpt(int, char **);

#endif
