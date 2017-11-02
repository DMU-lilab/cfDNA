/* building index with hash table */

#include "index.h"
#include <linux/limits.h>


/*! @funciton: split the primer line
 *   @parmeters:
 *   @    bufline   the pointer to the primer line [char *]
 *   @    amp       the pointer to the amplicon line [amp_t *]
 *   @return:
 *   @    void         [void]
*/
static void LineSplit(char *bufline, amp_t *amp)
{
    char *pword;

    pword = strtok(bufline, ",");
    strcpy(amp->fwdprim, pword);

    for (int i=0; i < 2; i++) {
        pword = strtok(NULL, ",");

        if (!i)
            strcpy(amp->revprim, pword);
        else
            amp->insertlen = atoi(pword);
    } amp->readnum = 0;

    pword = strtok(NULL, "\n");
    strcpy(amp->auxinfo, pword);
}

/*! @funciton: split the primer line
 *   @parmeters:
 *   @    filename   the pointer to the filename of the primer [char *]
 *   @return:
 *   @    p          the pointer to the primer structure [prim_t *]
*/
prim_t *GetPrim(char *filename)
{
    char buf[BUFLINE];

    FILE *file = fopen(filename, "r");
    if (!file)
        goto _ferror;

    prim_t *p = (prim_t *)calloc(1, sizeof(prim_t));
    p->amp = (amp_t *)calloc(1, sizeof(amp_t));

    if (!p || !p->amp)
        goto _memerror;

    while (fgets(buf, BUFLINE, file)) {
        if (buf[0] != '#') {
            /* realloc necessary memory for amplicon list */
            if (p->ampnum % PRIMNUM == 0) {
                int newsize = p->ampnum + PRIMNUM;
                amp_t *tem = (amp_t *)realloc(p->amp, newsize*sizeof(amp_t));
                if (!tem) goto _memerror;
                p->amp = tem;
            }
            LineSplit(buf, &(p->amp[p->ampnum]));
            p->ampnum++;
        }
    } fclose(file);

    return p;

  _ferror:
      fprintf(stderr, "[Err:%s] Failed to open %s\n", __func__, filename);
      exit(-1);
  _memerror:
      fprintf(stderr, "[Err:%s] Failed to alloc memory\n", __func__);
}

/*! @funciton: get the suffix of path
 *   @parmeters:
 *   @    path       the pointer to the path [char *]
 *   @return:
 *   @               the pointer to the suffix start
*/
static char *GetSuffix(char *path)
{
    int plen, i;

    plen = strlen(path);
    for (i=plen; i > 0; i--) {
        /* get base path for out amplicon count */
        /*　_t => xlzh_trim_fq */
        if (path[i] == '_' && path[i+1] == 't')
            break;
    }
    return (path+i);
}


/*! @funciton: get the reverse and complement sequecne
 *   @parmeters:
 *   @    string     the pointer to target string [char *]
 *   @return:
 *   @    void       the pointer to the suffix start
*/
void RevComp(char *string)
{
    uint8_t tmp;
    char *last;

    for (last=string; *last; ++last) ;
    string--;
    while (++string <= --last) {
        tmp = *string;
        *string = _BASE[*last];
        *last = _BASE[tmp];
    }
    return ;
}


/*! @funciton: get the reverse primer structue
 *   @parmeters:
 *   @    P             the pointer to forward primer [prime_t *]
 *   @return:
 *   @   primlist       the pointer to the reverse primer list [prim_t **]
*/
prim_t **RevPrim(prim_t *P)
{
    prim_t **primlist;
    
    primlist = (prim_t **)calloc(P->ampnum, sizeof(prim_t*));
    if (!primlist)
        goto _memerror;

    for (int i=0; i < P->ampnum; i++) {
        primlist[i] = (prim_t *)malloc(sizeof(prim_t));
        primlist[i]->amp = (amp_t *)malloc(sizeof(amp_t));
        if (!primlist[i] || !primlist[i]->amp)
            goto _memerror;

        primlist[i]->ampnum = 1;
        primlist[i]->amp[0] = P->amp[i];

        RevComp(primlist[i]->amp[0].fwdprim);
        RevComp(primlist[i]->amp[0].revprim);
    }

    return primlist;

  _memerror:
      fprintf(stderr, "[Err:%s] Failed to alloc memory\n", __func__);
}


/*! @funciton: build index for the primer structure
 *   @parmeters:
 *   @    p          the pointer to primer [prime_t *]
 *   @    T          the hast talbe contain the index [hash_t *]
 *   @    kmer       the index length [hash_t *]
 *   @return:
 *   @    void
*/
void PrimIndex(prim_t *p, hash_t *T, int kmer)
{
    loc_t loc;
    char slic[KEYLEN];

    for (int i=0; i < p->ampnum; i++) {
        amp_t *amp = &p->amp[i];

        for (int j=0; j < 2; j++) {
            char *primseq = j ? amp->revprim : amp->fwdprim;
            int plen = strlen(primseq) -kmer + 1;

            for (int k=0; k < plen; k++) {
                strncpy(slic, primseq+k, kmer);
                slic[kmer] = '\0';
                loc.ploc = i; loc.frloc = j; loc.sloc = k;
                Insert(slic, &loc, T);
            }
        }   
    }
    return ;
}


/*! @funciton: get the index list of amplicon reverse primer
 *   @parmeters:
 *   @    P             the pointer to the primer list [prim_t **]
 *   @    ampnum        the amplicon count [int]
 *   @    kmer          the kmer for indexing [int]
 *   @return:
 *   @    revlist       [hash_t **]
*/
hash_t **RevIndex(prim_t **P, int ampnum, int kmer)
{
    hash_t **revlist;

    revlist = (hash_t **)calloc(ampnum, sizeof(hash_t*));
    if (!revlist)
        goto _memerror;

    for (int i=0; i < ampnum; i++) {
        revlist[i] = InitHash(1<<5);
        PrimIndex(P[i], revlist[i], kmer);
    }

    return revlist;

  _memerror:
      fprintf(stderr, "[Err:%s] Failed to alloc memory\n", __func__);
}


/*! @funciton: write out the amplicon infomation
 *   @parmeters:
 *   @    prim          the pointer to the forward primer [prim_t *]
 *   @    path          the pointer to the outpath [char *]
 *   @return:
 *   @    void        
*/
void AmpWrite(prim_t *prim, char *path)
{
    char fname[PATH_MAX];

    strcpy(fname, path);
    strcpy(GetSuffix(fname), ".ampcount");

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "Gene,Chrom,AmpStart,InsertStart,InsertEnd,AmpEnd,FwdPrim,RevPrim,AmpCount\n");
    for (int i=0; i < prim->ampnum; i++) {
        amp_t *amp = &prim->amp[i];
        fprintf(fp, "%s,%s,%s,%d\n", \
                amp->auxinfo, amp->fwdprim, amp->revprim, amp->readnum);
    } fclose(fp);

    return ;
}
