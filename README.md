Himapper: a Hi-C data analysis tool
===================

# Introduction


# Download

## Conda

## Github
  ```
  $ git clone https://github.com/hsinnan75/himapper.git
  ```
to download the package of Himapper.

# Dependencies

To compile Himapper, it requires libboost-all-dev, libbz2-dev, and liblzma-dev installed in your system.

Please download juicer_tools.jar (https://github.com/aidenlab/juicer/wiki/Download) for creating .hic files.

# Compiling

Please change to Himapper's folder and just type 'make' to compile himapper. If the compilation or the program fails, please contact me (arith@iis.sinica.edu.tw), Thanks.

# Usage

To index a reference genome, Himapper requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ bin/himapper index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
The above command is to index the reference (Ecoli.fa) with prefix of Ecoli.

To analyze Hi-C data, Himapper requires the the index files of the reference genome and read files of paired-end reads. Users should use -i to specify the prefix of the index files (including the directory path).

  ```
 $ bin/himapper -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -o alignment.txt
  ```
The above command is to align the paired-end reads in ReadFile1.fa and ReadFile2.fa with index files of ecoli and output the mapping result to alignment.txt.

Then use juicer_tools to create .hic file

  ```
 $java -jar juicer_tools.jar pre alignment.txt ecoli.hic ecoli
  ```

Then go to http://www.aidenlab.org/juicebox/ and upload the generated .hic file to produce the contact map.

# File formats

- Reference genome files

    All reference genome files should be in FASTA format.

- Read files

    All reads files should be in FASTA or FASTQ format. FASTQ files can be compressed with gzip format. We do not support FASTA files with gzip compression.
    Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

- Output file

    Output is in Medium format (https://github.com/aidenlab/juicer/wiki/Pre)
    A whitespace separated file that contains, on each line
    <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
    - str = strand (0 for forward, anything else for reverse)
    - chr = chromosome (must be a chromosome in the genome)
    - pos = position
    - frag = restriction site fragment
    - mapq = mapping quality score

# Parameter setting

 ```
-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq or fq.gz]

-f2 STR read filename2 [optional, fasta or fastq or fq.gz], f and f2 are files with paired reads

-o STR read alignment filename in txt format [alignment.txt]

  ```
