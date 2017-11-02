/* seq query function */

#include "query.h"


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


/* CORE FUNCTION 
 *      seq: AGAAATTTGCGGAGTAAGTTGCGCTGGGGCTTTCGGCGGCGGCGATTT
 *   primer:      TTTGCGGA
 *                |      |
 *               pstart  pend
 * */
static query_t SeqQuery(char *seq, hash_t *H, prim_t *P, int mis, int kmer)
{
    status S;
    query_t Q = {1,0,0,0};
    char key[KEYLEN];

    int cutlen = strlen(seq) -kmer;

    /* 'N' is imposible occured in primer index */
    for (int i=0; i < cutlen; i++) {
        if (*(seq+i) == 'N')
            continue;

        strncpy(key, seq+i, kmer); key[kmer] = '\0';
        Search(key, H, &S);
        if (S.find == 0) // the key is not in the hash table
            continue;

        for (int j=0; j < S.node->num; j++) {
            loc_t *loc = &S.node->loc[j];
            int shift = i - loc->sloc;
            
            char *string = seq;
            char *pstr = loc->frloc ?
                    P->amp[loc->ploc].revprim : P->amp[loc->ploc].fwdprim;

            if (shift >= 0)
                string = seq + shift;
            else if (shift < 0 && shift >= -5) 
                pstr -= shift;
            else
                continue;

            if (MisCheck(string, pstr, mis)) {
                Q.ploc = loc->ploc;
                Q.pstart = shift < 0 ? 0 : shift;
                Q.pend = Q.pstart + strlen(pstr) - 1;
                return Q;
            }
        }
    } Q.isfind = 0;
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

    strncpy(temseq, seq, SEARCHLEN); 
    temseq[SEARCHLEN] = '\0';
    f = SeqQuery(temseq, arg->fwdindex, 
                arg->fwdprim, arg->args->mismatch, arg->args->kmer);
    if (!f.isfind) { 
        Q.isfind = 0; return Q; 
    }

    amp_t *amp = &(arg->fwdprim->amp[f.ploc]);
    if (f.pend+amp->insertlen > READLEN) {
        Q.isfind = 1;
        Q.ploc = f.ploc; 
        Q.pstart = f.pend + 1;
        Q.pend = 0; // the revprimer is not sequenced!
        return Q;
    }

    int rstart = f.pend + amp->insertlen - 4;
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
 *   @return:
 *   @    void
*/
void PrimTrim(fastq_t *fq, query_t *Q)
{
    read_t *read;
    char *seq, *qual;
    int inlen;

    read = &fq->read;
    seq = read->seq; qual = read->qual;

    if (Q->isfind && Q->pend) {
        inlen = Q->pend - Q->pstart +1;

        strncpy(seq, seq+Q->pstart, inlen); 
        seq[inlen] = '\n'; seq[inlen+1] = '\0';
        strncpy(qual, qual+Q->pstart, inlen);
        qual[inlen] = '\n'; qual[inlen+1] = '\0';
    }
    else if (Q->isfind && !Q->pend) {
        strcpy(seq, seq+Q->pstart);
        strcpy(qual, qual+Q->pstart);
    }
    else { // can't locate the primer seq
        strcpy(seq, "NNNNNNNNNNNNNNNNNNNN\n");
        strcpy(qual, "!!!!!!!!!!!!!!!!!!!!\n");
    }

    fq->cache[fq->bufnum++] = *read;

    return ;
}


