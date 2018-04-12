/* seq query function */

#include "query.h"
#include "utils.h"

/* buffer length for kmer search */
#define BLEN 6

static int _compscore(const void *a, const void *b) {
    /* sorted by descending order */
    return ((compo_t*)b)->score - ((compo_t*)a)->score;
}

/*! @funciton: check the mismatch base
 *   @parmeters:
 *   @    str1      the pointer to the string 1 [char *]
 *   @    str2      the pointer to the string 2 [char *]
 *   @    mismatch  the maxmum mismatch base allowed [int]
 *   @return:
 *   @    0     False
 *   @    1     True
*/
static int MisCheck(char *str1, char *str2, int mismatch)
{
    int misnum = 0;
    
    if (!*str1 || !*str2) return 0;

    while (*str1 && *str2) {
        if (*str1++ != *str2++)
            misnum++;

        if (misnum > mismatch)
            return 0;
    }
    return 1;
}

/*! @funciton: calculate the maximux primer length
 *   @parmeters:
 *   @    P        pointer to the primer struct [prim_t *]
 *   @return:
 *   @    max      the maximum primer length
*/
int MaxPrimLen(prim_t *P)
{
    int max=0, flen, rlen;
    amp_t *a;

    for (int i=0; i < P->ampnum; ++i) {
        a = &P->amp[i];
        
        /* calculate the length of forward and reverse primer */
        flen = strlen(a->fwdprim);
        rlen = strlen(a->revprim);

        /* obtain the maximux length of the both primer */
        max = flen > max ? flen : max;
        max = rlen > max ? rlen : max;
    }
    return max;
}

static void CalHits(hit_t *hit, node_t *node, int seqi)
{
    loc_t *loc;
    compo_t *c;

    for (int i=0; i < node->num; ++i) {
        loc = &node->loc[i];

        if (hit->hitnum % 0x8 == 0) {
            /* when n==0, realloc is similar to malloc */
            int n = hit->hitnum + 0x8;
            err_realloc(hit->compo, n, compo_t);
        }

        c = hit->compo;
        int flag = 0;
        for (int j=0; j < hit->hitnum; ++j) {
            if (loc->ploc == c[j].ploc && loc->frloc == c[j].frloc) {
                /* the loc has existed in the hit */
                c[j].score++; 
                flag = 1; break;
            }
        }

        if (!flag) { /* the loc not existed in the hit */
            c = &hit->compo[hit->hitnum++];
            c->ploc = loc->ploc; c->frloc = loc->frloc;
            c->sloc = loc->sloc; c->seqi = seqi;
            c->score = 1; /* the value may not be zero after relloc */
        }
    }
    return ;
}


static hit_t *GetHits(char *seq, hash_t *H, int kmer)
{
    char key[KEYLEN];
    status S;
    hit_t *hit;

    err_calloc(hit, 1, hit_t);
    int cutlen = strlen(seq);

    /* 'N' is imposible occured in primer index */
    for (int i=0; i < cutlen; ++i) {
        if (seq[i] == 'N') continue;

        strncpy(key, seq+i, kmer); key[kmer] = '\0';
        Search(key, H, &S);
        if (S.find == 0) continue;

        CalHits(hit, S.node, i);
    }

    if (hit->hitnum)
        /* sort the score by descending order */
        qsort(hit->compo, hit->hitnum, sizeof(compo_t), _compscore);

    return hit;
}

/* CORE FUNCTION 
 *      seq: AGAAATTTGCGGAGTAAGTTGCGCTGGGGCTTTCGGCGGCGGCGATTT
 *   primer:      TTTGCGGA
 *                |      |
 *               pstart  pend
 * */
static query_t SeqQuery(char *seq, hash_t *H, prim_t *P, int mis, int kmer)
{
    hit_t *hit;
    query_t Q = {0};

    hit = GetHits(seq, H, kmer);
    if (!hit->hitnum) {/* can't find even one kmer hit */
        err_free(hit); return Q;
    }

    /* only choose the first 2 high-score hit */
    int n = hit->hitnum > 1 ? 2 : 1;

    for (int i=0; i < 1; ++i) {
        compo_t *c = &hit->compo[i];
        int shift = c->seqi - c->sloc;
        
        char *string = seq;
        char *pstr = c->frloc ?
                P->amp[c->ploc].revprim : P->amp[c->ploc].fwdprim;

        if (shift >= 0)
            string = seq + shift;
        else if (shift < 0 && shift >= -1 * mis)
            pstr -= shift;
        else
            continue;

        if (MisCheck(string, pstr, mis)) {
            {
                Q.isfind = 1;
                Q.ploc = c->ploc;
                Q.pstart = shift < 0 ? 0 : shift;
                Q.pend = Q.pstart + strlen(pstr) -1;
            }
            err_free(hit->compo); err_free(hit);
            return Q;
        }
    } 
    err_free(hit->compo); err_free(hit);
    return Q;
}


/*! @funciton: get the insert fragment start and end (0-base)
 *   @parmeters:
 *   @    seq       the pointer to fastq sequence [char *]
 *   @    arg       all the necessary parmeters [arginfo_t *]
 *   @return:
 *   @    Q         the start and end of the insert fragment [query_t]
*/
query_t PrimQuery(char *seq, arginfo_t *arg)
{
    query_t f, r, Q ;
    char temseq[FQLINE];

    strncpy(temseq, seq, arg->maxpl+BLEN);
    temseq[arg->maxpl+BLEN] = '\0';


    f = SeqQuery(temseq, arg->fwdindex,
                arg->fwdprim, arg->args->mismatch, arg->args->kmer);
    if (!f.isfind) {
        Q.isfind = 0; return Q; 
    }

    amp_t *amp = &(arg->fwdprim->amp[f.ploc]);
    if (f.pend+amp->insertlen > arg->maxrl) {
        Q.isfind = 1;
        Q.ploc = f.ploc;
        Q.pstart = f.pend + 1;
        Q.pend = 0; // the revprimer is not sequenced!
        return Q;
    }

    int rstart = f.pend + amp->insertlen - BLEN;
    r = SeqQuery(seq+rstart, arg->revindex_list[f.ploc],
                arg->revprim_list[f.ploc], arg->args->mismatch, arg->args->kmer);
    if (!r.isfind) { 
        Q.isfind = 0; return Q; 
    }

    Q.isfind = 1;
    Q.ploc = f.ploc;
    Q.pstart = f.pend + 1;
    Q.pend = rstart + r.pstart -1;

    return Q;
}


/*! @funciton: trim off the primer seq
 *   @parmeters:
 *   @    fq        the pointer to fastq structure [fastq_t *]
 *   @    Q         the query result from PrimQuery [query_t *]
 *   @    arg       the args required [arginfo_t *]
 *   @return:
 *   @    void
*/
void PrimTrim(fastq_t *fq, query_t *Q, arginfo_t *arg)
{
    read_t *read;
    char *seq, *qual;
    int discard=0;

    read = &fq->cache[fq->bufnum++];
    seq = fq->read.seq; qual = fq->read.qual;

    strcpy(read->name, fq->read.name); 
    strcpy(read->mark, fq->read.mark);

    if (Q->isfind) { // find the primer sequence 
        strcpy(read->seq, seq+Q->pstart);
        strcpy(read->qual, qual+Q->pstart);

        if (Q->pend) {
            /* read through */
            int inlen = Q->pend - Q->pstart +1;
            strcpy(read->seq+inlen, "\n");
            strcpy(read->qual+inlen, "\n");
        }

        float q = MeanQuality(read->qual, arg->phred);
        if (q < arg->args->minqual) {
            Q->badqual = 1; discard = 1;
        }
    }
    else { // can't find the primer sequence
        if (arg->args->keep) {
            strcpy(read->seq, fq->read.seq);
            strcpy(read->qual, fq->read.qual);
        }
        else discard = 1;
    }
    if (discard) {
        strcpy(read->seq, "NNNNNNNNNNNNNNNNNNNN\n");
        strcpy(read->qual, "!!!!!!!!!!!!!!!!!!!!\n");
    }
    
    return ;
}


