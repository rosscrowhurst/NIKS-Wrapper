NIKS-Wrapper
============

Wrapper script for NIKS pipeline

Modified version of the original NIKS pipeline kmerpipeline.sh script. The latter is part 
of the original NIKS pipeline code available from:

 	http://sourceforge.net/projects/niks/

DESCRIPTION
===========

The NIKS-kmerPipeline.sh NIKS wrapper script is an alternative way to 
run NIKS. Although you can run NIKS according to the README 
distributed with the NIKS code our wrapper does offer a few advantages
including:

1. It checks that the requisite programs are in your path and will provide
   an opportunity to run them from locations on your file system of your choice
   if you have multiple versions of the programs as we do in our software tree.

2. It incorporates the annotation steps that were in a separate script in the 
   original NIKS code distribution.

3. It includes code to be able to restart the analysis just by re-issuing the 
   original command line. It will check for the presences of expected output
   files and if found will not repeat the step required to produce those outputs.

4. It renames NIKS “contigs” to reflect the pairwise comparison undertaken
   enabling subsequent merging of results from multiple pairwise comparisons in 
   to a single results file without issues of replacement of similarly named 
   contigs from the different comparisons.

5. If you wish to leverage an A. thaliana accession genome reference then our 
   wrapper script will make a blast+ database from it if you have not done this 
   manually before running the pipeline.

6. It corrects a few minor code errors in the original annotation script from NIKS
   and incorporates these steps in to the pipeline directly.

7. Our wrapper writes all command lines to be used to a output script that you
   subsequently run. THis allows you to check the whole pipeline before running
   and to store the runnable in your repository for that analysis.
8. Our wrapper script captures the versions of all binaries and writes them to the
   runnable script - important for reproducibility. The original NIKS code 
   ommitted setting some variables in the annotation script. In this modified 
   version of the original script instead echo out the commands to a separate
   BASH shell script so the whole pipeline commands can be checked before running 
   which allows checking

IMPORTANT 
=========

1. When passing names for mutants to our wrapper only use an alphanumeric alphabet
   e.g ddi1, ddi2, ddi3 and do not include spaces or control characters

2. Unlike the original kmerPipeline.sh our alternative wrapper expects the
   FASTQ sequence files to be uncompressed. The original kmerPipeline.sh handles 
   both compressed and uncompressed. However, if you pass compressed files 
   that have been comptressed with gzip they will be automatically uncompressed
   for use in the pipeline

3. Although our code includes ability to re-start the pipeline this does not
   apply to the last (annotation) step as this is quick in any case.

4. Some of the more complex process substitutions in the original
   BASH kmerpipeline.sh script from (http://sourceforge.net/projects/niks/)
   were converted to sequential commands as the original code was failing
   when run under Open Lava where as the non process substiution code ran.
   The cause of failure under Open Lava may have been specific to our 
   installation setup
5. To permanently modify variables edit the parameters in the "USER SET VARIABLE SECTION"
    

This modified version can also be restarted and will automatically
skip steps that have already successfully been run.

Important note: original script allowed gz file input whereas this script uses uncompressed
fastq files. If gz compressed files are given on command line they will be decompressed
before use




USAGE
=====

${SCRIPTNAME} <options>

Where options include:

REQUIRED SWITCHES
-----------------
-DATADIR path
-WORKINGBASEDIR path (toplevel dir for pair analysis)
-MUTANTA_NAME name (e.g Dis9)
-MUTANTB_NAME name (e.g Dis15)
-MUTANTA_READ1_FILE filename (e.g. Dis9_ec_read_1.fastq)
-MUTANTA_READ2_FILE filename (e.g. Dis9_ec_read_2.fastq)
-MUTANTB_READ1_FILE filename (e.g. Dis15_ec_read_1.fastq)
-MUTANTB_READ2_FILE filename (e.g. Dis15_ec_read_2.fastq)

 If compressed then the read files will be decompressed before use
 
OPTIONAL SWITCHES
--------------------
-NIKS_SCRIPTS_DIR path (path to scripts dir for NIKS distribution if not default location)
-FASTQINSERT integer
-CORES integer
-KMERSIZE integer 
-MAXCUTOFF integer
-LMERSIZE integer
-MINSIZE integer
-HIGHBOUNDARY integer
-TOLERANCE integer
-JAVA_BIN path
-JELLYFISH_BIN path
-GENOME_REFERENCE_DIR path
-REFERENCE_GENOME_FASTA path
-NAMEBLASTDBGENOME string
-SHUFFLESEQUENCES_FASTQ path
-VELVETH path
-VELVETG path
-MAKEBLASTDB path
-BLASTN path

-DEBUG

