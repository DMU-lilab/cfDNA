PrimerTrim
=========================
Used to trim off the primer sequence from amplicon fastq file


__PROGRAM: PrimerTrim__<br>
__PLATFORM: Linux__<br>
__COMPILER: gcc-4.8.5__<br>
__AUTHOR: xlzh__<br>
__EMAIL: xiaolongzhang2015@163.com__<br>
__DATE: 2017-09-21__<br>
__DEPENDENCE__<br>
* zlib-1.2.7<br>
#### NOTE
* The first you need to do are confirming the libraries above have been installed.<br />
* And the gcc compiler should be available on your server.<br /><br />


Description
=========================
The program is used to trim off the primer sequence of the target sequencing. <br>
The performing of k-mer indexing alogrithm makes it possible to deal with __thousands of amplicon primer pairs__ at the same time.<br>
Compared with other kinds of tools, this program could trim the primer sequence off directly from the fastq file, which could save you a lot of time.<br>


Building
=========================

The simple way to compile the program, is done as follows:
* cd ./cfDNA/PrimerTrim/
* make
* make clean

Then there will be generate an executable file ---> primtrim-1.2


Usage
========================
      Options:
         -h|--help        print help infomation
         -a|--ampfile     [required] input amplicon file[.csv]
         -r|--read1       [required] read1 for fastq file[.fq|.gz]
         -l|--read2       [required] read2 for fastq file[.fq|.gz]
         -o|--outdir      [required] output directory for trimed fastq file
         -k|--kmer        [optional] The kmer lenght for indexing [8]
         -m|--mismatch    [optional] the maxmum mismatch for primer seq [3]

Option
========================
#### \[-a|--ampfile]
      The path to amplicon file with csv format.
      There is an example in the "Example" directory, which will give you a guidance
      about how to specify each field of the amplicon file.
      [FIELDS]
            1. forwardprimer: the forward primer sequecne [required]
            2. reverseprimer: the reverse primer sequence [required]
            3. insertlength: the insert length between the primer pair [required]
            4. (gene, chorm, ...) : the auxiliary field  that used to describe the
               primer pair [ at least one field should be specified!]
      eg. chen_hg19_amplicon.csv

#### \[-r|--read1]
      The read1 of the fastq file or gziped fastq file [must have an marker of "_R"]
      eg. chen_WBC_R1.fq.gz

#### \[-r|--read2]
      The read2 of the fastq file or gziped fastq file [must have an marker of "_R"]
      eg. chen_WBC_R2.fq.gz

#### \[-o|--outdir]
      The output directory of the trimed fastq file

#### \[-k|--kmer]
      The kmer length for primer indexing, a larger kmer will perform an accurate locating
      to the primer sequence, which also means a lower posibility to find the primer 
      sequence. [8 is recommended]

#### \[-m|--mismatch]
      the maxmum mismatch allowed in the  primer matching [smaller than 4 is recommended]


OutPut
=========================
The program will generate 3 files:
1. \*_trim_R1.fq
2. \*_trim_R2.fq
3. \*.ampcount

 __\*_trim_R1.fq and \*_trim_R2.fq are the fastq file after primer trimed!__<br>
 __The last field of the \*.ampcount file is the reads belong to the amplicon.__<br>
