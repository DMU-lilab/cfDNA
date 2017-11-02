/* function: parse the comand line parmeters */

#include "fastq.h"

void Usage(void)
{
    char *usage =
        "\nUsage: primtrim [options]\n"
        "Version: 1.2\n"
        "\n"
        "Options:\n"
        "       -h|--help        print help infomation\n"
        "       -a|--ampfile     [required] input amplicon file[.csv]\n"
        "       -r|--read1       [required] read1 for fastq file[.fq|.gz]\n"
        "       -l|--read2       [required] read2 for fastq file[.fq|.gz]\n"
        "       -o|--outdir      [required] output directory for trimed fastq file\n"
        "       -k|--kmer        [optional] The kmer lenght for indexing [8]\n"
        "       -m|--mismatch    [optional] the maxmum mismatch for primer seq [3]\n\n";

    fprintf(stderr, "%s", usage); 
    exit(-1);
}

static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "ampfile", required_argument, NULL, 'a' },
    { "read1", required_argument, NULL, 'r' },
    { "read2", required_argument, NULL, 'l' },
    { "outdir", required_argument, NULL, 'o' },
    { "kmer", optional_argument, NULL, 'k' },
    { "mismatch", optional_argument, NULL, 'm' },
    { NULL, 0, NULL, 0 }
};


arg_t *ParseOpt( int argc, char **argv )
{
    int opt =0, opterr =0;

    arg_t *Arg = (arg_t *)calloc(1, sizeof(arg_t));
    if ( !Arg ) {
        fprintf(stderr, \
            "[Err::%s::%d]  Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    while ( (opt = getopt_long(argc, argv, "a:r:l:o:k:m:h", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'a': strcpy(Arg->ampfile, optarg); break;
            case 'r': strcpy(Arg->read1, optarg); break;
            case 'l': strcpy(Arg->read2, optarg); break;
            case 'o': strcpy(Arg->outdir, optarg); break;
            case 'k': Arg->kmer = atoi(optarg); break;
            case 'm': Arg->mismatch = atoi(optarg); break;
            case '?': fprintf(stderr, \
                            "[Err::%s::%d]  Option error occour!.\n", __func__, __LINE__);
                      Arg->help = 1;
        }
    } 
    if (!Arg->ampfile[0] || !Arg->read1[0] \
            || !Arg->read2[0] || !Arg->outdir[0]) {
        fprintf(stderr, \
            "[Err::%s::%d]  Please give the [requied] parmeters!\n", __func__, __LINE__);
        Arg->help = 1;
    }
    if (!Arg->kmer) Arg->kmer = 8;
    if (!Arg->mismatch) Arg->mismatch = 3;

    return Arg;
}

