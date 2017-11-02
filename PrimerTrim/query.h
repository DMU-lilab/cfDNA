#ifndef QUERY_H
#define QUERY_H

#include "index.h"
#include "fastq.h"

/* Default search lenth for the primer seq
 * NOTE: Because primer always range from 20bp to 30bp, which means after 40bp,
 *   there is no need to search anymore. In this way, it coluld seed up the program
*/
#define SEARCHLEN 40


/* Default Illumina read length
 * NOTE: The value is used to decide whther your fragment length is larger 
 *   than your read length, if True, only 5' sequence will be consider, else,
 *   the program should take the reverse complement primer into consideration,
 *   and trim off both forward primer and reverse complement primer
*/
#define READLEN 150


/*! @typedef query_t
 @abstract structure for the return value of query function
 @field isfind      whether find the match primer
 @field ploc        the primer index in the amplicon list
 @field pstart      primer start index in the query sequecing
 @field pend        primer end index in the query sequecing
*/
typedef struct __query_t {
    int isfind;
    int ploc;
    int pstart;
    int pend;
} query_t;


/*! @typedef arginfo_t
 @abstract structure for the thread parmeters
 @field args            come from comand line parmeters
 @field fwdindex        foroward primer index
 @field revindex_list   reverse primer index table list
 @field fwdprim         foroward primer list
 @field revprim_list    reverse primer list
 @field fastq           read1 or read2 fastq container
 @field isread2         1=> read2 otherwise 0=> read1
*/
typedef struct __arginfo_t {
    arg_t *args;
    hash_t *fwdindex;
    hash_t **revindex_list;
    prim_t *fwdprim;
    prim_t **revprim_list;
    fastq_t fastq;
    int isread2;
} arginfo_t;


/* prototype function */
query_t PrimQuery(char *, arginfo_t *);
void PrimTrim(fastq_t *, query_t *);

#endif
