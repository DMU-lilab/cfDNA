#include "fastq.h"

/*! @funciton: get the basename of path
 *   @parmeters:
 *   @    path       the pointer to the path [char *]
 *   @return:
 *   @               the pointer to the basename start
*/
static char *BaseName(char *path)
{
    int plen, i;

    plen = strlen(path);
    for (i=plen; i > 0; i--) {
        if (path[i] == '/') 
            break;
    }
    if (i)
        return (path+i+1);
    else
        return path;
}


/*! @funciton: fastq initialization
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @    args       the pointer to the arg_t structure [fastq_t *]
 *   @    type       the read1 or read2 for the fastq file
 *   @return:
 *   @    void
*/
void FastqInit(fastq_t *fq, arg_t *args, int type)
{
    char oname[128];

    type == READ1 
        ? strcpy(fq->inname, args->read1) 
        : strcpy(fq->inname, args->read2);

    fq->bufnum = 0;
    strcpy(oname, BaseName(args->read1));
    char *target = strstr(oname, "_R");
    if (!target)
        goto _invalidfq;

    type == READ1 
        ? strcpy(target, "_trim_R1.fq") 
        : strcpy(target, "_trim_R2.fq");
    strcpy(fq->outname, args->outdir);
    strcat(fq->outname, oname);

    fq->in = gzopen(fq->inname, "r"); // open the file [gz or fq]
    if (!fq->in) goto _inerror;
    gzbuffer(fq->in, 524288);
    fq->out = fopen(fq->outname, "w");

    return ;

  _inerror:
      fprintf(stderr, "[Err:%s] \
              Failed to open %s\n", __func__, fq->inname); exit(-1);
  _invalidfq:
      fprintf(stderr, "[Err:%s] \
              Invalid fastq name! [stander: *_R1.fq.gz]\n", __func__); exit(-2);
}


/*! @funciton: read a fastq group each time [seqname, seq, +, qual]
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @return:
 *   @    void
*/
int FastqRead(fastq_t *fq)
{
    char buf[BUFLINE];

    for (int i=0; i < 4; i++) {
        if (gzgets(fq->in, buf, BUFLINE)) {
            switch (i) {
                case 0: strcpy(fq->read.name, buf); break;
                case 1: strcpy(fq->read.seq,  buf); break;
                case 2: strcpy(fq->read.mark, "+\n"); break;
                case 3: strcpy(fq->read.qual, buf);
            }
        }
        else return 0;
    }
    return 1;
}


/*! @funciton: write the fastq group from cache
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @return:
 *   @    void
*/
void FastqWrite(fastq_t *fq)
{
    read_t *read;
   
    read = fq->cache;
    for (int i=0; i < fq->bufnum; i++) {
        for (int j=0; j < 4; j++) {
            switch (j) {
                case 0: fputs(read[i].name, fq->out); break;
                case 1: fputs(read[i].seq,  fq->out); break;
                case 2: fputs(read[i].mark, fq->out); break;
                case 3: fputs(read[i].qual, fq->out);
            }
        }
    } fq->bufnum = 0;
    return ;
}

