/*
*  PROGRAM: PrimerTrim
*  FUNCTION: Used to trim off the primer seq from the fastq file
*  INPUT:
*       (1) amplicon file [.csv]
*       (2) fastq1 and fastq2 [.fq or .fq.gz]
*       (3) output directory [path]
*  OUTPUT:
*       (1) Fastq file after trim off the primer seq [.fq]
*       (2) Amplicon count file [.ampcount]
*
*  AUTHOR: xiaolong Zhang
*  EMAIL:  xiaolongzhang2015@163.com
*  DATE:   2017-09-21 
*/

#include <pthread.h>
#include <time.h>

#include "query.h"

uint64_t Total = 0;
uint64_t Bad = 0;
pthread_mutex_t Lock = PTHREAD_MUTEX_INITIALIZER;


void *process(void *arginfo)
{
    query_t Q;
    int total=0, bad=0;
    int *ampcount, ampnum;
    arginfo_t *Arg = (arginfo_t *)arginfo;
    fastq_t *fq = &Arg->fastq;

    FastqInit(fq, Arg->args, Arg->isread2);

    ampnum = Arg->fwdprim->ampnum;
    ampcount = (int *)calloc(ampnum, sizeof(int));
    if (!ampcount)
        fprintf(stderr, "[Err:%s:%d]\
                Failed to alloc memory\n", __func__, __LINE__);

    while (FastqRead(fq)) {
        total++;
        if (fq->bufnum % BUFNUM == 0)
            FastqWrite(fq);
        
        Q = PrimQuery(fq->read.seq, Arg);
        if (Q.isfind)
            ampcount[Q.ploc]++;
        else bad++;

        PrimTrim(fq, &Q);
    } FastqWrite(fq);

    pthread_mutex_lock(&Lock); 
    amp_t *amp = Arg->fwdprim->amp;
    for (int i=0; i < ampnum; i++)
        amp[i].readnum += ampcount[i];
    
    Total += total; Bad += bad;
    pthread_mutex_unlock(&Lock);

    free(ampcount);
}


int main(int argc, char **argv)
{
    pthread_t tid[2];
    arginfo_t *arginfo1, *arginfo2;
    time_t start, end;
    arg_t *args;
   

    args = ParseOpt(argc, argv);
    if (args->help) 
        Usage();

    prim_t *fwdprim = GetPrim(args->ampfile);
    prim_t **revprim_list = RevPrim(fwdprim);
    hash_t *fwdindex = InitHash(fwdprim->ampnum << 5);
    PrimIndex(fwdprim, fwdindex, args->kmer);
    hash_t **revindex_list = RevIndex(revprim_list, fwdprim->ampnum, args->kmer);

    arginfo1 = (arginfo_t *)calloc(1, sizeof(arginfo_t));
    arginfo2 = (arginfo_t *)calloc(1, sizeof(arginfo_t));
    if (!arginfo1 || !arginfo2)
        goto _memerror;

    arginfo1->args = args;
    arginfo1->fwdprim = fwdprim;
    arginfo1->revprim_list = revprim_list;
    arginfo1->fwdindex = fwdindex;
    arginfo1->revindex_list = revindex_list;

    memcpy(arginfo2, arginfo1, sizeof(arginfo_t));
    arginfo1->isread2 = READ1; arginfo2->isread2 = READ2;

    time(&start);
    for (int i=0; i < 2; i++) {
        printf("[*] Processing the [thread: %d] ...\n", i+1);
        if (i == 0)
            pthread_create(&tid[i], NULL, process, (void*)arginfo1);
        else
            pthread_create(&tid[i], NULL, process, (void*)arginfo2);
    }

    for (int i=0; i < 2; i++)
        pthread_join(tid[i], NULL);

    AmpWrite(fwdprim, arginfo1->fastq.outname);

    time(&end);
    printf("Total time consume: %.2f\n", difftime(end,start));

    /* --------------------  summary ------------------- */

    fprintf(stderr, "\n----------------- Summary ------------------------\n");
    fprintf(stderr, "Total reads processed: %lu\n", Total);
    fprintf(stderr, "Reads can't find primer: %lu\n", Bad);
    fprintf(stderr, "Reads mismatch with primer: %.2f %%\n", (float)Bad/Total*100);

    return 0;

  _memerror:
      fprintf(stderr, \
          "[Err:%s:%d] Failed to alloc memory\n", __func__, __LINE__);
}


