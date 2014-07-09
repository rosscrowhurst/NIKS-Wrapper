#!/bin/bash
# Modified version (by Ross Crowhurst ross.crowhurst@plantandfood.co.nz)
# of the original NIKS pipeline kmerpipeline.sh script. The latter is part 
# of the original NIKS pipeline code available from:
#
# 	(http://sourceforge.net/projects/niks/) 
#
# DESCRIPTION
#
# Our NIKS-kmerPipeline.sh NIKS wrapper script is an alternative way to 
# run NIKS. Although you can run NIKS according to the pipeline’s README 
# distributed with the NIKS code our wrapper does offer a few advantages
# including:
# 1. It checks that the requisite programs are in your path and will provide
#    an opportunity to run them from locations on your file system of your choice
#    if you have multiple versions of the programs as we do in our software tree.
# 2. It incorporates the annotation steps that were in a separate script in the 
#    original NIKS code distribution.
# 3. It includes code to be able to restart the analysis just by re-issuing the 
#    original command line. It will check for the presences of expected output
#    files and if found will not repeat the step required to produce those outputs.
# 4. It renames NIKS “contigs” to reflect the pairwise comparison undertaken
#    enabling subsequent merging of results from multiple pairwise comparisons in 
#    to a single results file without issues of replacement of similarly named 
#    contigs from the different comparisons.
# 5. If you wish to leverage an A. thaliana accession genome reference then our 
#    wrapper script will make a blast+ database from it if you have not done this 
#    manually before running the pipeline.
# 6. It corrects a few minor code errors in the original annotation script from NIKS
#    and incorporates these steps in to the pipeline directly.
# 7. Our wrapper writes all command lines to be used to a output script that you
#    subsequently run. THis allows you to check the whole pipeline before running
#    and to store the runnable in your repository for that analysis.
# 8. Our wrapper script captures the versions of all binaries and writes them to the
#    runnable script - important for reproducibility. The original NIKS code 
#    ommitted setting some variables in the annotation script. In this modified 
#    version of the original script instead echo out the commands to a separate
#    BASH shell script so the whole pipeline commands can be checked before running 
#    which allows checking
#
# IMPORTANT 
#
# 1. When passing names for mutants to our wrapper only use an alphanumeric alphabet
#    e.g ddi1, ddi2, ddi3 and do not include spaces or control characters
# 2. Unlike the original kmerPipeline.sh our alternative wrapper expects the
#    FASTQ sequence files to be uncompressed. The original kmerPipeline.sh handles 
#    both compressed and uncompressed. However, if you pass compressed files 
#    that have been comptressed with gzip they will be automatically uncompressed
#    for use in the pipeline
# 3. Although our code includes ability to re-start the pipeline this does not
#    apply to the last (annotation) step as this is quick in any case.
# 4. Some of the more complex process substitutions in the original
#    BASH kmerpipeline.sh script from (http://sourceforge.net/projects/niks/)
#    were converted to sequential commands as the original code was failing
#    when run under Open Lava where as the non process substiution code ran.
#    The cause of failure under Open Lava may have been specific to our 
#    installation setup
# 5. To permanently modify variables edit the parameters in the "USER SET VARIABLE SECTION"
#    
#  
#
# This modified version can also be restarted and will automatically
# skip steps that have already successfully been run.
#
# 
#
# Important note: original script allowed gz file input whereas this script uses uncompressed
# fastq files. If gz compressed files are given on command line they will be decompressed
# before use
#---------------------------------------------------------------------------------------------
# DATE OF LAST CODE MODIFICATION
LAST_UPDATED="28-05-2014 08:14"
# MODIFIED CODE VERSIONN
VERSION="1.1e"
# TURN ON DEBUGGING
DEBUG=yes

#---------------------------------------------------------------------------------------------
# START OF "USER SET VARIABLE SECTION"
#---------------------------------------------------------------------------------------------
# All variables in this section can be set via the command line and it is recommended that
# the command line options are used as opposed to changing the script itself
# 
# SCRIPTS:
#     Set to location of the "scripts" sub directory of the NIKS pipeline
#     code on your system or use the -NIKS_SCRIPTS_DIR command line option to set
SCRIPTS=/software/x86_64/NIKS/kmerPipeline/scripts

# The following are set to the defaults in the original kmerPipeline.sh
FASTQINSERT=200
#number of cores to use for blast searches
CORES=10
#KMERSIZE RD:61
KMERSIZE=31
# maximum value for a cutoff from the histogram
MAXCUTOFF=100
# l-mer for seed generation RD:55
LMERSIZE=25
# minimum number of kmers to make up a seed (minlength= MINSIZE+KMERSIZE-1) RD:51
MINSIZE=21
# seeds with a total coverage of this many kmers are discarded
HIGHBOUNDARY=10000
#tolerance in seed pairing
TOLERANCE=10

# Required binaries
JAVA=`which java`
#JELLYFISH=`which jellyfish-1.1.10`
JELLYFISH=`which jellyfish`
MAKEBLASTDB=`which makeblastdb`
BLASTN=`which blastn`
# Needed by "assemble.sh"
SHUFFLESEQUENCES_FASTQ=`which shuffleSequences_fastq.pl`
VELVETH=`which velveth`
VELVETG=`which velvetg` 

SHUFFLESEQUENCES_FASTA=`which shuffleSequences_fasta.pl`
#---------------------------------------------------------------------------------------------
# END OF "USER SET VARIABLE SECTION"
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
# DO NOT MODIFY CODE BELOW HERE
#---------------------------------------------------------------------------------------------

# Name of script as its running
SCRIPTNAME="$(basename "$(test -L "$0" && readlink "$0" || echo "$0")")"

# Setting debug output - on by default
set -x
if [ $DEBUG = yes ] ; then
        set -xv
else
        echo "[ok] DEBUG is $DEBUG"
fi


#---------------------------------------------------------------------------------------------
# USAGE
#---------------------------------------------------------------------------------------------
Usage()
{
        #[ "$DEBUG" = yes] && set -x
        echo "
Usage:  ${SCRIPTNAME} <options>

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
  (default $SCRIPTS)
-FASTQINSERT integer (default ${FASTQINSERT})
-CORES integer (default ${CORES})
-KMERSIZE integer (default ${KMERSIZE})
-MAXCUTOFF integer (default ${MAXCUTOFF})
-LMERSIZE integer (default ${LMERSIZE})
-MINSIZE integer (default ${MINSIZE})
-HIGHBOUNDARY integer (default ${HIGHBOUNDARY})
-TOLERANCE integer (default ${TOLERANCE})
-JAVA_BIN path (default ${JAVA})
-JELLYFISH_BIN path (default ${JELLYFISH})
-GENOME_REFERENCE_DIR path
-REFERENCE_GENOME_FASTA path
-NAMEBLASTDBGENOME string

-SHUFFLESEQUENCES_FASTQ path (default `which shuffleSequences_fastq.pl`)
-VELVETH path (default `which velveth`)
-VELVETG path (default `which velvetg`)

-MAKEBLASTDB path (default `which makeblastdb`)
-BLASTN path (default `which blastn`)

-DEBUG

LAST_UPDATED: ${LAST_UPDATED}
VERSION: ${VERSION}
"
        exit 1
}

# Default base working directory is the present working directory
WORKINGBASEDIR=`pwd`
WORKINGDIR=
DATADIR=
MUTANTA=
MUTANTB=

# If using reference genome and annotation. 
# Setting GENOME_REFERENCE_DIR etc will toogle DO_GENOME_ANNOTATION
DO_GENOME_ANNOTATION=no
GENOME_REFERENCE_DIR=
REFERENCE_GENOME_FASTA=
NAMEBLASTDBGENOME=
BLASTDBGENOME=temporary

#---------------------------------------------------------------------------------------------
# Do not change below
#
#---------------------------------------------------------------------------------------------
#==========================================================================#
# RETRIEVE COMMANDLINE VARIABLES
#--------------------------------------------------------------------------#
if [ $# -gt 0 ]; then
        while [ $# -gt 0 ]
        do
                case $1 in
                -DATADIR)       # Base directory for new sequences
                        DATADIR=$2
                        shift
                        ;;
                -MUTANTA_NAME)
                        MUTANTA=$2
                        shift;
                        ;;
                -MUTANTA_READ1_FILE)
                        MUTANTA_READ1_FILE=$2
                        shift;
                        ;;
                -MUTANTA_READ2_FILE)
                        MUTANTA_READ2_FILE=$2
                        shift;
                        ;;
                -MUTANTB_NAME)
                        MUTANTB=$2
                        shift;
                        ;;
                -MUTANTB_READ1_FILE)
                        MUTANTB_READ1_FILE=$2
                        shift;
                        ;;
                -MUTANTB_READ2_FILE)
                        MUTANTB_READ2_FILE=$2
                        shift;
                        ;;
                -WORKINGBASEDIR)
                        WORKINGBASEDIR=$2
                        shift;
                        ;;
                -NIKS_SCRIPTS_DIR)
                        SCRIPTS=$2
                        shift;
                        ;;
                -FASTQINSERT)
                        FASTQINSERT=$2
                        shift;
                        ;;
                -DEBUG)
                        DEBUG=$2
                        shift;
                        ;;
                -CORES)
                        CORES=$2
                        shift;
                        ;;
                -KMERSIZE)
                        KMERSIZE=$2
                        shift;
                        ;;
                -MAXCUTOFF)
                        MAXCUTOFF=$2
                        shift;
                        ;;
                -LMERSIZE)
                        LMERSIZE=$2
                        shift;
                        ;;
                -MINSIZE)
                        MINSIZE=$2
                        shift;
                        ;;
                -HIGHBOUNDARY)
                        HIGHBOUNDARY=$2
                        shift;
                        ;;
                -TOLERANCE)
                        TOLERANCE=$2
                        shift;
                        ;;
                -JAVA_BIN)
                        JAVA=$2
                        shift;
                        ;;
                -JELLYFISH_BIN)
                        JELLYFISH=$2
                        shift;
                        ;;
                -GENOME_REFERENCE_DIR)
                        GENOME_REFERENCE_DIR=$2
                        shift;
                        ;;
                -REFERENCE_GENOME_FASTA)
                        REFERENCE_GENOME_FASTA=$2
                        shift;
                        ;;
                -NAMEBLASTDBGENOME)
                        NAMEBLASTDBGENOME=$2
                        shift;
                        ;;
                -DO_GENOME_ANNOTATION)
                        DO_GENOME_ANNOTATION=$2
                        shift;
                        ;;
                -SHUFFLESEQUENCES_FASTQ)
                        SHUFFLESEQUENCES_FASTQ=$2
                        shift;
                        ;;
                -VELVETH)
                        VELVETH=$2
                        shift;
                        ;;
                -VELVETG)
                        VELVETG=$2
                        shift;
                        ;;
                -MAKEBLASTDB)
                        MAKEBLASTDB=$2
                        shift;
                        ;;
                -BLASTN)
                        BLASTN=$2
                        shift;
                        ;;
                esac
                shift
        done
else
        Usage
        exit 1
fi

# Reset our working directory to that passed by user
WORKINGDIR=${WORKINGBASEDIR}/NIKS/${MUTANTA}-${MUTANTB}

# Make the working directory
mkdir -p ${WORKINGDIR}
if [ -d ${WORKINGDIR} ] ; then
	echo "Using working directory: ${WORKINGDIR}"
else
	echo "Can not make directory ${WORKINGDIR}"
	exit 1
fi

# Set our current directory  
CURRENTDIR=`pwd`

# Set base name for output (scripts)
SHSCRIPTBASENAME=kmerPipeline.${MUTANTA}-${MUTANTB}.$$

# Set name of run script 
SHFILE=${WORKINGDIR}/${SHSCRIPTBASENAME}.sh

# Change to working directing
cd ${WORKINGDIR}

# Check and remove run script if it exists - not really requied as process id is used in script name - just a precaution
if [ -f ${WORKINGDIR}/${SHFILE} ] ; then
        rm -f ${WORKINGDIR}/${SHFILE}
fi

# Begin writting to run script - need to use bash as has process substitution
# which will fail if "sh" is not an alias for "bash" shell
echo "#!/bin/bash" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# PARAMETER SETTINGS" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}

# Check data directory exists
if [ -d ${DATADIR} ] ; then
        echo "DATADIR=${DATADIR}" >> ${SHFILE}
else
        echo "# A valid DATADIR name was not supplied" >> ${SHFILE}
        exit 1
fi

# Check NIKS scriptsdirectory exists
if [ -d ${SCRIPTS} ] ; then
        echo "SCRIPTS=${SCRIPTS}" >> ${SHFILE}
else
        echo "# A valid SCRIPTS name was not supplied" >> ${SHFILE}
        exit 1
fi

# Check JAVA var has been set
if [ -z "${JAVA}" ] ; then
	echo "ERROR: java path not set or in path" >> ${SHFILE}
        exit 1
fi
# Check java binary exists
if [ -e ${JAVA} ] ; then
        ${JAVA} -version 2>${WORKINGDIR}/jver.txt
        JAVA_VERSION=`grep version ${WORKINGDIR}/jver.txt`
        rm -f ${WORKINGDIR}/jver.txt
        echo "# using ${JAVA_VERSION}" >> ${SHFILE}
        echo "JAVA=${JAVA}" >> ${SHFILE}
else
        echo "ERROR: JAVA=${JAVA} does not exist" >> ${SHFILE}
        exit 1
fi
# Check jellyfish var has been set
if [ -z "${JELLYFISH}" ] ; then
	echo "ERROR: jellyfish path not set or in path" >> ${SHFILE}
        exit 1
fi
# Check jellyfish binary exists
if [ -e ${JELLYFISH} ] ; then
	JELLYFISH_VERSION=`${JELLYFISH} --version`
	echo "# Using ${JELLYFISH_VERSION}" >> ${SHFILE}
        echo "JELLYFISH=${JELLYFISH}" >> ${SHFILE}
        MAJOR_VERSION=`echo ${JELLYFISH_VERSION} | awk '{ split($2,version,"."); major_version=version[1]; print major_version}'`
        if [ ${MAJOR_VERSION} != 1 ] ; then
        	echo "# FATAL ERROR: JELLYFISH VERSION ERROR "
        	echo "# This ppiple need jellyfish version 1 currently "
        	echo "# Please supply the path to jellyfish version 1 via command line, e.g "
        	echo "# -JELLYFISH_BIN /software/x86_64/jellyfish-1.1.10/bin/jellyfish "
        	exit 1
        fi
else
        echo "ERROR: JELLYFISH=${JELLYFISH} does not exist" >> ${SHFILE}
        exit 1
fi
# Check shuffleseq var has been set
if [ -z "${SHUFFLESEQUENCES_FASTQ}" ] ; then
	echo "ERROR: shuffleSequences_fastq.pl not set or in path" >> ${SHFILE}
        exit 1
fi
# Check shuffleseq bin exists
if [ -e ${SHUFFLESEQUENCES_FASTQ} ] ; then
        echo "SHUFFLESEQUENCES_FASTQ=${SHUFFLESEQUENCES_FASTQ}" >> ${SHFILE}
else
        echo "ERROR: SHUFFLESEQUENCES_FASTQ=${SHUFFLESEQUENCES_FASTQ} does not exist" >> ${SHFILE}
        exit 1
fi
# Check velveth var has been set
if [ -z "${VELVETH}" ] ; then
	echo "ERROR: velveth not set or in path" >> ${SHFILE}
        exit 1
fi
# Check velveth binary exists
if [ -e ${VELVETH} ] ; then
        VELVETH_VERSION=`${VELVETH} | grep Version`
        echo "# using ${VELVETH_VERSION}" >> ${SHFILE}
        echo "VELVETH=${VELVETH}" >> ${SHFILE}
else
        echo "ERROR: VELVETH=${VELVETH} does not exist" >> ${SHFILE}
        exit 1
fi
# Check velvetg var has been set
if [ -z "${VELVETG}" ] ; then
	echo "ERROR: velvetg not set or in path" >> ${SHFILE}
        exit 1
fi
# Check velvetg binary exists
if [ -e ${VELVETG} ] ; then
        VELVETG_VERSION=`${VELVETG} | grep Version`
        echo "# using ${VELVETG_VERSION}" >> ${SHFILE}
        echo "VELVETG=${VELVETG}" >> ${SHFILE}
else
        echo "ERROR: VELVETG=${VELVETG} does not exist" >> ${SHFILE}
        exit 1
fi
# Check makeblastdb var is set
if [ -z "${MAKEBLASTDB}" ] ; then
	echo "ERROR: makeblastdb not set or in path" >> ${SHFILE}
        exit 1
fi
# Check makeblastdb bin exists
if [ -e ${MAKEBLASTDB} ] ; then
        ${MAKEBLASTDB} -version >${WORKINGDIR}/mkdb.txt
        MAKEBLASTDB_VERSION=`grep makeblastdb ${WORKINGDIR}/mkdb.txt`
        echo "# using ${MAKEBLASTDB_VERSION}" >> ${SHFILE}
        echo "MAKEBLASTDB=${MAKEBLASTDB}" >> ${SHFILE}
else
        echo "ERROR: MAKEBLASTDB=${MAKEBLASTDB} does not exist" >> ${SHFILE}
        exit 1
fi
# Check blastn var is set
if [ -z "${BLASTN}" ] ; then
	echo "ERROR: blastn not set or in path" >> ${SHFILE}
        exit 1
fi
# Check blastn bin exists
if [ -e ${BLASTN} ] ; then
        ${BLASTN} -version >${WORKINGDIR}/bltn.txt
        BLASTN_VERSION=`grep blastn ${WORKINGDIR}/bltn.txt`
        echo "# using ${BLASTN_VERSION}" >> ${SHFILE}
        echo "BLASTN=${BLASTN}" >> ${SHFILE}
else
        echo "ERROR: BLASTN=${BLASTN} does not exist" >> ${SHFILE}
        exit 1
fi

# Test KMERSIZE as max for jellyfish version 1 is 31
if [ ${KMERSIZE} -gt 31 ] ; then
	# Check version of jellyfish
	IS_V1=`echo $JELLYFISH_VERSION | grep -c "jellyfish 1"`
	if [ $IS_V1 = 1 ] ; then 
	        echo "# Using jellyfish version 1"
	        echo "# The maximum kmer supportedd by jellyfish version 1 is 31"
	        echo "# KMERSIZE must be 31 or less"
		exit 1
	fi
fi

echo "CURRENTDIR=${CURRENTDIR}" >> ${SHFILE}
echo "WORKINGDIR=${WORKINGDIR}" >> ${SHFILE}
echo "CORES=${CORES}" >> ${SHFILE}
echo "KMERSIZE=${KMERSIZE}" >> ${SHFILE}
echo "MAXCUTOFF=${MAXCUTOFF}" >> ${SHFILE}
echo "LMERSIZE=${LMERSIZE}" >> ${SHFILE}
echo "MINSIZE=${MINSIZE}" >> ${SHFILE}
echo "HIGHBOUNDARY=${HIGHBOUNDARY}" >> ${SHFILE}
echo "TOLERANCE=${TOLERANCE}" >> ${SHFILE}
echo "FASTQINSERT=${FASTQINSERT}" >> ${SHFILE}

KMERUTILS=${SCRIPTS}/kmerUtils.jar
if [ -e ${KMERUTILS} ] ; then
        echo "KMERUTILS=${KMERUTILS}" >> ${SHFILE}
else
        echo "# ERROR: KMERUTILS=${KMERUTILS} does not exist" >> ${SHFILE}
        exit 1
fi

# Stuff required if annotation to be done as part of pipeline
if [ -d ${GENOME_REFERENCE_DIR} ] ; then
        echo "GENOME_REFERENCE_DIR=${GENOME_REFERENCE_DIR}" >> ${SHFILE}
        REF=${GENOME_REFERENCE_DIR}/${REFERENCE_GENOME_FASTA}
        if [ -f ${REF} ] ; then
                echo "REFERENCE_GENOME_FASTA=${REFERENCE_GENOME_FASTA}" >> ${SHFILE}
        else
                echo "# ERROR: REFERENCE_GENOME_FASTA=${REFERENCE_GENOME_FASTA} does not exist and is required for genome annotation" >> ${SHFILE}
                exit 1
        fi
        # Check if there is a name for the blast db
        if [ -z ${NAMEBLASTDBGENOME} ] ; then
        	NAMEBLASTDBGENOME=`basename ${REFERENCE_GENOME_FASTA}`
        fi
        echo "NAMEBLASTDBGENOME=${NAMEBLASTDBGENOME}" >> ${SHFILE}
        BLASTDBGENOME=${GENOME_REFERENCE_DIR}/${NAMEBLASTDBGENOME}
        # Prepare reference blast+ database
        if [ ! -e ${BLASTDBGENOME}.nsq ] ; then
                echo "cd ${GENOME_REFERENCE_DIR}" >> ${SHFILE}
                echo "${MAKEBLASTDB} -in ${REFERENCE_GENOME_FASTA} -dbtype nucl -out ${NAMEBLASTDBGENOME}" >> ${SHFILE}
                echo "cd ${WORKINGDIR}" >> ${SHFILE}
        else
                echo "# Blast database for ${NAMEBLASTDBGENOME} already exists in ${GENOME_REFERENCE_DIR}" >> ${SHFILE}
                echo "# Skipping \"cd ${GENOME_REFERENCE_DIR}; ${MAKEBLASTDB} -in ${REFERENCE_GENOME_FASTA} -dbtype nucl -out ${NAMEBLASTDBGENOME}; cd ${WORKINGDIR}\" " >> ${SHFILE}
        fi
	DO_GENOME_ANNOTATION=yes        
fi
echo "DO_GENOME_ANNOTATION=${DO_GENOME_ANNOTATION}" >> ${SHFILE}

#Sample A fastq files
FASTQA1=${DATADIR}/${MUTANTA_READ1_FILE}
FASTQA2=${DATADIR}/${MUTANTA_READ2_FILE}
FILETYPE1=`file ${FASTQA1}`
FILETYPE2=`file ${FASTQA2}`
IS_GZIP1=`echo $FILETYPE1 | grep -c "gzip compressed"`
IS_GZIP2=`echo $FILETYPE2 | grep -c "gzip compressed"`
if [ $IS_GZIP1 = 1 ] ; then
	gzip -d ${FASTQA1}
	MUTANTA_READ1_FILE=`basename ${FASTQA1} .gz`
	FASTQA1=${DATADIR}/${MUTANTA_READ1_FILE}
fi
if [ $IS_GZIP2 = 1 ] ; then
	gzip -d ${FASTQA2}
	MUTANTA_READ2_FILE=`basename ${FASTQA2} .gz`
	FASTQA2=${DATADIR}/${MUTANTA_READ2_FILE}
fi
if [ -r ${FASTQA1} ] ; then
        echo "FASTQA1=${FASTQA1}" >> ${SHFILE}
else
        echo "# A valid MUTANTA_READ1_FILE name was not supplied" >> ${SHFILE}
        exit 1
fi
if [ -r ${FASTQA2} ] ; then
        echo "FASTQA2=${FASTQA2}" >> ${SHFILE}
else
        echo "# A valid MUTANTA_READ2_FILE name was not supplied" >> ${SHFILE}
        exit 1
fi
#Sample B fastq files
FASTQB1=${DATADIR}/${MUTANTB_READ1_FILE}
FASTQB2=${DATADIR}/${MUTANTB_READ2_FILE}
FILETYPE1=`file ${FASTQB1}`
FILETYPE2=`file ${FASTQB2}`
IS_GZIP1=`echo $FILETYPE1 | grep -c "gzip compressed"`
IS_GZIP2=`echo $FILETYPE2 | grep -c "gzip compressed"`
if [ $IS_GZIP1 = 1 ] ; then
	gzip -d ${FASTQB1}
	MUTANTB_READ1_FILE=`basename ${FASTQB1} .gz`
	FASTQB1=${DATADIR}/${MUTANTB_READ1_FILE}
fi
if [ $IS_GZIP2 = 1 ] ; then
	gzip -d ${FASTQB2}
	MUTANTB_READ2_FILE=`basename ${FASTQB2} .gz`
	FASTQB2=${DATADIR}/${MUTANTB_READ2_FILE}
fi
if [ -r ${FASTQB1} ] ; then
        echo "FASTQB1=${FASTQB1}" >> ${SHFILE}
else
        echo "# A valid MUTANTB_READ1_FILE name was not supplied" >> ${SHFILE}
        exit 1
fi
if [ -r ${FASTQB2} ] ; then
        echo "FASTQB2=${FASTQB2}" >> ${SHFILE}
else
        echo "# A valid MUTANTB_READ2_FILE name was not supplied" >> ${SHFILE}
        exit 1
fi

OUTPREFIXA=`echo $MUTANTA | tr '[A-Z]' '[a-z]'`
if [ -z "${OUTPREFIXA}" ] ; then
	echo "# OUTPREFIXA is empty - set variable MUTANTA" >> ${SHFILE}
	exit 1
fi
OUTPREFIXB=`echo $MUTANTB | tr '[A-Z]' '[a-z]'`
if [ -z "${OUTPREFIXB}" ] ; then
	echo "# OUTPREFIXB is empty - set variable MUTANTB" >> ${SHFILE}
	exit 1
fi

echo "OUTPREFIXA=${OUTPREFIXA}" >> ${SHFILE}
echo "OUTPREFIXB=${OUTPREFIXB}" >> ${SHFILE}

echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 1" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}

echo "# STEP 1A: JELLYFISH on Sample A" >> ${SHFILE}
if [ ! -f ${OUTPREFIXA}_kmers_jellyfish ] ; then
	echo "# Running STEP 1A" >> ${SHFILE}
        echo "cd ${WORKINGDIR}" >> ${SHFILE}
        echo "mkdir -p ${WORKINGDIR}/${OUTPREFIXA}_kmers" >> ${SHFILE}
        echo "${JELLYFISH} count -C -o ${OUTPREFIXA}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 10G < ${FASTQA1} ${FASTQA2}" >> ${SHFILE}
        echo "COUNT=\$(ls ${OUTPREFIXA}_kmers/tmp* |wc -l)" >> ${SHFILE}
        echo "if [ \$COUNT = 0 ] ; then" >> ${SHFILE}
        echo "	# JELLYFISH Sample A COUNT is Zero" >> ${SHFILE}
        echo "	exit 1" >> ${SHFILE}
        echo "fi" >> ${SHFILE}
        echo "if [ \$COUNT -eq 1 ] ; then"  >> ${SHFILE}
        echo "	mv ${OUTPREFIXA}_kmers/tmp_0 ${OUTPREFIXA}_kmers_jellyfish" >> ${SHFILE}
        echo "else"  >> ${SHFILE}
        echo "	${JELLYFISH} merge -o ${OUTPREFIXA}_kmers_jellyfish ${OUTPREFIXA}_kmers/tmp*"  >> ${SHFILE}
        echo "fi"  >> ${SHFILE}
        echo "rm -rf ${OUTPREFIXA}_kmers" >> ${SHFILE}
else
	echo "# ${OUTPREFIXA}_kmers_jellyfish already exists, skipping jellyfish count for ${OUTPREFIXA}" >> ${SHFILE}
fi

if [ ! -f ${OUTPREFIXA}.kmers.hist.csv ] ; then
        echo "${JELLYFISH} histo -f -o ${OUTPREFIXA}.kmers.hist.csv -t ${CORES} ${OUTPREFIXA}_kmers_jellyfish"  >> ${SHFILE}
        echo "awk '{print \$2\"\\t\"\$1}' ${OUTPREFIXA}.kmers.hist.csv > ${OUTPREFIXA}_tmp"  >> ${SHFILE}
        echo "mv ${OUTPREFIXA}_tmp ${OUTPREFIXA}.kmers.hist.csv" >> ${SHFILE}
        echo "if [ -e \"${OUTPREFIXA}.kmers.hist.csv\" ] ; then" >> ${SHFILE}
        echo "	echo \"# ${OUTPREFIXA}.kmers.hist.csv exists\" " >> ${SHFILE}
        echo "else" >> ${SHFILE}
        echo "	echo \" ${OUTPREFIXA}.kmers.hist.csv does not exist\"" >> ${SHFILE}
        echo "	exit 1" >> ${SHFILE}
        echo "fi" >> ${SHFILE}
else
	echo "# ${OUTPREFIXA}_kmers_jellyfish already exists, skipping jellyfish hist for ${OUTPREFIXA}" >> ${SHFILE}
fi
echo "# STEP 1A: End" >> ${SHFILE}

echo "# STEP 1B: JELLYFISH on Sample B" >> ${SHFILE}
if [ ! -f ${OUTPREFIXB}_kmers_jellyfish ] ; then
	echo "# Running STEP 1B" >> ${SHFILE}
        echo "cd ${WORKINGDIR}" >> ${SHFILE}
        echo "mkdir -p ${WORKINGDIR}/${OUTPREFIXB}_kmers"  >> ${SHFILE}
        echo "${JELLYFISH} count -C -o ${OUTPREFIXB}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 10G < ${FASTQB1} ${FASTQB2}" >> ${SHFILE}
        echo "COUNT=\$(ls ${OUTPREFIXB}_kmers/tmp* |wc -l)" >> ${SHFILE}
        echo "if [ \$COUNT = 0 ] ; then" >> ${SHFILE}
        echo "	echo \"#JELLYFISH B COUNT is Zero\" " >> ${SHFILE}
        echo "	exit 1" >> ${SHFILE}
        echo "fi" >> ${SHFILE}
        echo "if [ \$COUNT -eq 1 ] ; then" >> ${SHFILE}
        echo "	mv ${OUTPREFIXB}_kmers/tmp_0 ${OUTPREFIXB}_kmers_jellyfish" >>  ${SHFILE}
        echo "else" >> ${SHFILE}
        echo "	${JELLYFISH} merge -o ${OUTPREFIXB}_kmers_jellyfish ${OUTPREFIXB}_kmers/tmp*" >> ${SHFILE}
        echo "fi" >> ${SHFILE}
        echo "rm -rf ${OUTPREFIXB}_kmers" >> ${SHFILE}
else
	echo "# ${OUTPREFIXB}_kmers_jellyfish already exists, skipping jellyfish count for ${OUTPREFIXB}" >> ${SHFILE}
fi

if [ ! -f ${OUTPREFIXB}.kmers.hist.csv ] ; then
        echo "${JELLYFISH} histo -f -o ${OUTPREFIXB}.kmers.hist.csv -t ${CORES} ${OUTPREFIXB}_kmers_jellyfish" >> ${SHFILE}
        echo "awk '{print \$2\"\\t\"\$1}' ${OUTPREFIXB}.kmers.hist.csv > ${OUTPREFIXB}_tmp" >> ${SHFILE}
        echo "mv ${OUTPREFIXB}_tmp ${OUTPREFIXB}.kmers.hist.csv" >> ${SHFILE}
        echo "if [ -f \"${OUTPREFIXB}.kmers.hist.csv\" ] ; then" >> ${SHFILE}
        echo "	echo \"# ${OUTPREFIXB}.kmers.hist.csv exists\" " >> ${SHFILE}
        echo "else" >> ${SHFILE}
        echo "	echo \"${OUTPREFIXB}.kmers.hist.csv does not exist\" " >> ${SHFILE}
        echo "	exit 1" >> ${SHFILE}
        echo "fi" >> ${SHFILE}
else
	echo "# ${OUTPREFIXB}_kmers_jellyfish already exists, skipping jellyfish hist for ${OUTPREFIXB}" >> ${SHFILE}
fi
echo "# STEP 1B: End" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# End of Part 1" >> ${SHFILE}
echo "# " >> ${SHFILE}

# 2. find cutoff values
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 2" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 2A: FIND CUTOFF VALUES FOR Sample A" >> ${SHFILE}
if [ -f ${OUTPREFIXA}_cutoff.csv ] ; then
	CUTOFFA=`cat ${OUTPREFIXA}_cutoff.csv`
	echo "CUTOFFA=$CUTOFFA" >> ${SHFILE}
else
	echo "CUTOFFA=\$(awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk ${OUTPREFIXA}.kmers.hist.csv)" >> ${SHFILE}
	echo "echo \$CUTOFFA > ${OUTPREFIXA}_cutoff.csv" >> ${SHFILE}
fi	
echo "if [ \${CUTOFFA} -lt 1 ] ; then" >> ${SHFILE}
echo "	echo \"#CUTOFFA=\${CUTOFFA} not greater than 0\"" >> ${SHFILE}
echo "	exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# STEP 2B: FIND CUTOFF VALUES FOR Sample B" >> ${SHFILE}
if [ -f ${OUTPREFIXB}_cutoff.csv ] ; then
	CUTOFFB=`cat ${OUTPREFIXB}_cutoff.csv`
	echo "CUTOFFB=$CUTOFFB" >> ${SHFILE}
else
	echo "CUTOFFB=\$(awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk ${OUTPREFIXB}.kmers.hist.csv)" >> ${SHFILE}
	echo "echo \$CUTOFFB > ${OUTPREFIXB}_cutoff.csv" >> ${SHFILE}
fi	
echo "if [ \${CUTOFFB} -lt 1 ] ; then" >> ${SHFILE}
echo "	echo \"#CUTOFFB=\${CUTOFFB} not greater than 0\"" >> ${SHFILE}
echo "	exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 2" >> ${SHFILE}
echo "# " >> ${SHFILE}

# 3. contrast kmer-distributions
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 3" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 3: CONTRAST KMER-DISTRIBUTIONS" >> ${SHFILE}
WATCH_FILE=${WORKINGDIR}/niks_step3_completed
echo "WATCH_FILE=${WATCH_FILE}" >> ${SHFILE}
echo "if [ ! -e \"${WATCH_FILE}\" ] ; then" >> ${SHFILE}
echo "	bash ${SCRIPTS}/mergeCountFilesMarkFilterJellyfish.sh ${SCRIPTS} ${OUTPREFIXA}_kmers_jellyfish ${OUTPREFIXB}_kmers_jellyfish \${CUTOFFA} \${CUTOFFA} \${CUTOFFB} \${CUTOFFB} ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff" >> ${SHFILE}
echo "	if [ ! -e \"${OUTPREFIXA}_unique.kmerDiff\" ] ; then" >> ${SHFILE}
echo "		echo \"#${OUTPREFIXA}_unique.kmerDiff does not exists\" " >> ${SHFILE}
echo "		exit 1" >> ${SHFILE}
echo "	fi" >> ${SHFILE}
echo "	if [ ! -e \"${OUTPREFIXB}_unique.kmerDiff\" ] ; then" >> ${SHFILE}
echo "		echo \"# ${OUTPREFIXB}_unique.kmerDiff does not exists\"" >> ${SHFILE}
echo "		exit 1" >> ${SHFILE}
echo "	fi" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "touch ${WATCH_FILE}" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 3" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 4" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 4: GENERATE SEEDS AND REMOVE HIGHLY REPETITIVE SEQUENCES" >> ${SHFILE}
echo "# STEP 4A: " >> ${SHFILE}
WATCH_FILE=${WORKINGDIR}/niks_step4A_completed
echo "WATCH_FILE=${WATCH_FILE}" >> ${SHFILE}
echo "if [ ! -e \"${WATCH_FILE}\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -Xmx64G -Xms64G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXA}_seedsRep > /dev/null  #|gzip -c > ${OUTPREFIXA}_seedsRep_endCondition.csv.gz" >> ${SHFILE}
echo "	if [ ! -f ${OUTPREFIXA}_seedsRep_long.txt ] ; then" >> ${SHFILE}
echo "		echo \"# Error 4A - file does not exist ${OUTPREFIXA}_seedsRep_long.txt\"" >> ${SHFILE}
echo "		exit 1" >> ${SHFILE}
echo "	fi" >> ${SHFILE}
echo "	awk -f ${SCRIPTS}/longToInfo.awk ${OUTPREFIXA}_seedsRep_long.txt | gzip -c> ${OUTPREFIXA}_seedsRep_longTable.csv.gz" >> ${SHFILE}
echo "	gunzip -c ${OUTPREFIXA}_seedsRep_longTable.csv.gz| awk -v hb=${HIGHBOUNDARY} '\$3<hb && \$4<hb {print \$1}' > ${OUTPREFIXA}_seeds.gi" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} sub ${OUTPREFIXA}_seedsRep.fa ${OUTPREFIXA}_seeds.gi > ${OUTPREFIXA}_seeds.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "touch ${WATCH_FILE}" >> ${SHFILE}

echo "# STEP 4B: " >> ${SHFILE}
WATCH_FILE=${WORKINGDIR}/niks_step4B_completed
echo "WATCH_FILE=${WATCH_FILE}" >> ${SHFILE}
echo "if [ ! -e \"${WATCH_FILE}\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -Xmx64G -Xms64G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXB}_unique.kmerDiff ${OUTPREFIXA}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXB}_seedsRep > /dev/null  #|gzip -c > ${OUTPREFIXB}_seedsRep_endCondition.csv.gz" >> ${SHFILE}
echo "	awk -f ${SCRIPTS}/longToInfo.awk ${OUTPREFIXB}_seedsRep_long.txt |gzip -c> ${OUTPREFIXB}_seedsRep_longTable.csv.gz" >> ${SHFILE}
echo "	gunzip -c ${OUTPREFIXB}_seedsRep_longTable.csv.gz| awk -v hb=${HIGHBOUNDARY} '\$3<hb && \$4<hb {print $1}' > ${OUTPREFIXB}_seeds.gi" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} sub ${OUTPREFIXB}_seedsRep.fa ${OUTPREFIXB}_seeds.gi > ${OUTPREFIXB}_seeds.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "touch ${WATCH_FILE}" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 4.0" >> ${SHFILE}
echo "# " >> ${SHFILE}


echo "# STEP 4.1: CLEAN REVERSE COMPLEMENT" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_seeds.fa\" ] ; then" >> ${SHFILE}
echo "	echo \"# ${OUTPREFIXA}_seeds.fa does not exists\"" >> ${SHFILE}
echo "	exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_seeds.fa\" ] ; then" >> ${SHFILE}
echo "	echo \"# ${OUTPREFIXB}_seeds.fa does not exists\"" >> ${SHFILE}
echo "	exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# Find the smaller set and filter this. Identify the seeds in the other set which could match to the first filtered set. Filter those" >> ${SHFILE}
echo "if [ \$( grep -c '>' ${OUTPREFIXA}_seeds.fa ) -lt \$( grep -c '>' ${OUTPREFIXB}_seeds.fa ) ] ; then" >> ${SHFILE}
        #start with A
        echo "	PREFIX1=$OUTPREFIXA" >> ${SHFILE}
        echo "	PREFIX2=$OUTPREFIXB" >> ${SHFILE}
echo "else" >> ${SHFILE}
        #start with B
        echo "	PREFIX1=$OUTPREFIXB" >> ${SHFILE}
        echo "	PREFIX2=$OUTPREFIXA" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# Cleaning small set" >> ${SHFILE}
echo "if [ ! -e \"\${PREFIX1}_seeds_filtered.fa\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in \${PREFIX1}_seeds.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -query \${PREFIX1}_seeds.fa -db \${PREFIX1}_seeds.fa -num_threads ${CORES} -dust no -outfmt '6 std qlen slen' |${JAVA} -XX:-UseGCOverheadLimit -jar ${KMERUTILS} blastRevCompFilterPipe \${PREFIX1}_seeds.fa > \${PREFIX1}_seeds_filtered.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# Reducing large set" >> ${SHFILE}
echo "if [ ! -e \"\${PREFIX2}_seeds_toUse.fa\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in \${PREFIX2}_seeds.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query \${PREFIX1}_seeds_filtered.fa -db \${PREFIX2}_seeds.fa -outfmt 6 -task blastn | cut -f 2 | sort -u | ${JAVA} -jar ${KMERUTILS} subPipe \${PREFIX2}_seeds.fa > \${PREFIX2}_seeds_toUse.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# Cleaning reduced large set" >> ${SHFILE}
echo "if [ ! -e \"\${PREFIX2}_seeds_filtered.fa\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in \${PREFIX2}_seeds_toUse.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query \${PREFIX2}_seeds_toUse.fa -db \${PREFIX2}_seeds_toUse.fa -dust no -outfmt '6 std qlen slen' |${JAVA} -XX:-UseGCOverheadLimit -jar ${KMERUTILS} blastRevCompFilterPipe \${PREFIX2}_seeds_toUse.fa > \${PREFIX2}_seeds_filtered.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 4.1" >> ${SHFILE}
echo "# " >> ${SHFILE}

echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 5" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 5. CONTRAST SEEDS (BLASTN -TASK BLASTN) REMOVE SEQUENCES WITH MULTIPLE HITS" >> ${SHFILE}
echo "# STEP 5A. " >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in ${OUTPREFIXB}_seeds_filtered.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query ${OUTPREFIXA}_seeds_filtered.fa -db ${OUTPREFIXB}_seeds_filtered.fa -task blastn -dust no -outfmt '6 std qseq sseq'|gzip -c > ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# STEP 5B. " >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in ${OUTPREFIXA}_seeds_filtered.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query ${OUTPREFIXB}_seeds_filtered.fa -db ${OUTPREFIXA}_seeds_filtered.fa -task blastn -dust no -outfmt '6 std qseq sseq'| gzip -c > ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 5.0" >> ${SHFILE}
echo "# " >> ${SHFILE}

echo "# STEP 5.1. SNPs " >> ${SHFILE}
echo "# STEP 5.1A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_SNP_blastRep.csv\" ] ; then" >> ${SHFILE}
echo "	gunzip -c ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXA}_SNP_blastRep.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# STEP 5.1A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_SNP_blastRep.csv\" ] ; then" >> ${SHFILE}
echo "	gunzip -c ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXB}_SNP_blastRep.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_SNP_blast.csv\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} blastRemoveRepetitive ${OUTPREFIXA}_SNP_blastRep.csv ${OUTPREFIXB}_SNP_blastRep.csv ${OUTPREFIXA}_SNP_blast.csv ${OUTPREFIXB}_SNP_blast.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 5.1" >> ${SHFILE}
echo "# " >> ${SHFILE}


echo "# STEP 5.3 MERGE " >> ${SHFILE}
echo "# STEP 5.3A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_annotated_blast.csv\" ] ; then" >> ${SHFILE}
echo "	cat ${OUTPREFIXA}_SNP_blast.csv > ${OUTPREFIXA}_annotated_blast.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# STEP 5.3B" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_annotated_blast.csv\" ] ; then" >> ${SHFILE}
echo "	cat ${OUTPREFIXB}_SNP_blast.csv > ${OUTPREFIXB}_annotated_blast.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# End of Part 5.3" >> ${SHFILE}
echo "# " >> ${SHFILE}

# 5.4. Extract fasta file with all seeds
echo "# STEP 5.4. EXTRACT FASTA FILE WITH ALL SEEDS" >> ${SHFILE}
echo "# STEP 5.4A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_candidates.fa\" ] ; then" >> ${SHFILE}
echo "	cut -f 1 ${OUTPREFIXA}_SNP_blast.csv |sort -gu > ${OUTPREFIXA}_candidates.gi" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} sub ${OUTPREFIXA}_seeds_filtered.fa ${OUTPREFIXA}_candidates.gi > ${OUTPREFIXA}_candidates.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# STEP 5.4B" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_candidates.fa\" ] ; then" >> ${SHFILE}
echo "	cut -f 1 ${OUTPREFIXB}_SNP_blast.csv |sort -gu > ${OUTPREFIXB}_candidates.gi" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} sub ${OUTPREFIXB}_seeds_filtered.fa ${OUTPREFIXB}_candidates.gi > ${OUTPREFIXB}_candidates.fa" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 5.4" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 6" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 6. EXTEND CANDIDATE SEEDS" >> ${SHFILE}
echo "# STEP 6.1. KMERIZE CANDIDATE SEEDS" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_candidates.kmers\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} kmerize ${OUTPREFIXA}_candidates.fa ${KMERSIZE} |sort -u >  ${OUTPREFIXA}_candidates.kmers" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_candidates.kmers\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} kmerize ${OUTPREFIXB}_candidates.fa ${KMERSIZE} |sort -u >  ${OUTPREFIXB}_candidates.kmers" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# End of Part 6.1" >> ${SHFILE}
echo "# " >> ${SHFILE}

echo "# STEP 6.2. EXTRACT READS" >> ${SHFILE}
#A
##${JAVA} -jar ${KMERUTILS} kmerCoverage_extractCoveredZ ${FASTQA1} ${FASTQA2} ${OUTPREFIXA}_candidates.kmers 0 ${OUTPREFIXA}_candidates_reads

#when FASTQ files are not gzipped:
echo "if [ ! -e \"${OUTPREFIXA}_candidates_reads\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} kmerCoverage_extractCovered ${FASTQA1} ${FASTQA2} ${OUTPREFIXA}_candidates.kmers 0 ${OUTPREFIXA}_candidates_reads" >> ${SHFILE}
echo "fi" >> ${SHFILE}

#B
##${JAVA} -jar ${KMERUTILS} kmerCoverage_extractCoveredZ ${FASTQB1} ${FASTQB2} ${OUTPREFIXB}_candidates.kmers 0 ${OUTPREFIXB}_candidates_reads
#when FASTQ files are not gzipped:
echo "if [ ! -e \"${OUTPREFIXB}_candidates_reads\" ] ; then" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} kmerCoverage_extractCovered ${FASTQB1} ${FASTQB2} ${OUTPREFIXB}_candidates.kmers 0 ${OUTPREFIXB}_candidates_reads" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# End of Part 6.2" >> ${SHFILE}
echo "# " >> ${SHFILE}

echo "# STEP 6.3. EXTEND EACH SEED" >> ${SHFILE}
echo "# STEP 6.3A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_candidates_extData.tar.gz\" ] ; then" >> ${SHFILE}
echo "	mkdir ${OUTPREFIXA}_candidates" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} split ${OUTPREFIXA}_candidates.fa ${OUTPREFIXA}_candidates/seq 1" >> ${SHFILE}
echo "	for s in \`ls ${OUTPREFIXA}_candidates/seq*.fa |grep -v extended\`" >> ${SHFILE}
echo "	do" >> ${SHFILE}
echo "		bash ${SCRIPTS}/assemble.sh \${s} ${OUTPREFIXA}_candidates_reads_1.fastq ${OUTPREFIXA}_candidates_reads_2.fastq ${FASTQINSERT} ${KMERSIZE} ${SCRIPTS}" >> ${SHFILE}
echo "	done" >> ${SHFILE}
echo "	cat ${OUTPREFIXA}_candidates/seq*extended.fa > ${OUTPREFIXA}_candidates_ext.fa" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} lengths ${OUTPREFIXA}_candidates_ext.fa > ${OUTPREFIXA}_candidates_ext.lengths" >> ${SHFILE}
echo "	tar czf ${OUTPREFIXA}_candidates_extData.tar.gz ${OUTPREFIXA}_candidates" >> ${SHFILE}
echo "else" >> ${SHFILE}
echo "	echo \"# ${OUTPREFIXA}_candidates_extData.tar.gz already exists, skipping processing ${OUTPREFIXA}_candidates\"" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# STEP 6.3B" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_candidates_extData.tar.gz\" ] ; then" >> ${SHFILE}
echo "	mkdir ${OUTPREFIXB}_candidates" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} split ${OUTPREFIXB}_candidates.fa ${OUTPREFIXB}_candidates/seq 1" >> ${SHFILE}
echo "	for s in \`ls ${OUTPREFIXB}_candidates/seq*.fa |grep -v extended\`" >> ${SHFILE}
echo "	do" >> ${SHFILE}
echo "		bash ${SCRIPTS}/assemble.sh \${s} ${OUTPREFIXB}_candidates_reads_1.fastq ${OUTPREFIXB}_candidates_reads_2.fastq ${FASTQINSERT} ${KMERSIZE} ${SCRIPTS}" >> ${SHFILE}
echo "	done" >> ${SHFILE}
echo "	cat ${OUTPREFIXB}_candidates/seq*extended.fa > ${OUTPREFIXB}_candidates_ext.fa" >> ${SHFILE}
echo "	${JAVA} -jar ${KMERUTILS} lengths ${OUTPREFIXB}_candidates_ext.fa > ${OUTPREFIXB}_candidates_ext.lengths" >> ${SHFILE}
echo "	tar czf ${OUTPREFIXB}_candidates_extData.tar.gz ${OUTPREFIXB}_candidates" >> ${SHFILE}
echo "else" >> ${SHFILE}
echo "	echo \"# ${OUTPREFIXB}_candidates_extData.tar.gz already exists, skipping processing ${OUTPREFIXB}_candidates\"" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "# " >> ${SHFILE}
echo "# End of Part 6.3" >> ${SHFILE}
echo "# " >> ${SHFILE}

# 7. annotate position in extended seed
#  < seeds, extended seed, seed annotation file
#  > extended seed annotation file
#
echo "#---------------------------------------------------" >> ${SHFILE}
echo "# STEP 7" >> ${SHFILE}
echo "#---------------------------------------------------" >> ${SHFILE}

echo "# STEP 7. ANNOTATE POSITION IN EXTENDED SEED" >> ${SHFILE}
echo "# STEP 7A" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_candidates_ext_annotation.csv\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in ${OUTPREFIXA}_candidates_ext.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query ${OUTPREFIXA}_candidates.fa -db ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXA}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXA}_candidates_ext_annotation.csv" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# STEP 7B" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_candidates_ext_annotation.csv\" ] ; then" >> ${SHFILE}
echo "	${MAKEBLASTDB} -in ${OUTPREFIXB}_candidates_ext.fa -dbtype nucl" >> ${SHFILE}
echo "	${BLASTN} -num_threads ${CORES} -query ${OUTPREFIXB}_candidates.fa -db ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXB}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXB}_candidates_ext_annotation.csv" >>  ${SHFILE}
echo "else " >> ${SHFILE}
echo "	echo \"# Skipping production of ${OUTPREFIXB}_candidates_ext_annotation.csv - already exisits\" " >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_candidates_counts.csv\" ] ; then join -1 1 -2 1 <(gunzip -c ${OUTPREFIXA}_seedsRep_longTable.csv.gz |grep -Ff <(cut -f 1 ${OUTPREFIXA}_candidates.gi ) |sort -k1,1) <(cut -f 1 ${OUTPREFIXA}_candidates.gi |sort -k1,1) |awk '{print \$1\"\\t\"\$3\"\\t\"\$4}' |sort -k1,1g > ${OUTPREFIXA}_candidates_counts.csv; fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_candidates_counts.csv\" ] ; then join -1 1 -2 1 <(gunzip -c ${OUTPREFIXB}_seedsRep_longTable.csv.gz |grep -Ff <(cut -f 1 ${OUTPREFIXB}_candidates.gi ) |sort -k1,1) <(cut -f 1 ${OUTPREFIXB}_candidates.gi |sort -k1,1) |awk '{print \$1\"\\t\"\$3\"\\t\"\$4}' |sort -k1,1g > ${OUTPREFIXB}_candidates_counts.csv; fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_summarySelf.csv\" ] ; then ${JAVA} -jar ${KMERUTILS} leftjoin ${OUTPREFIXA}_SNP_blast.csv ${OUTPREFIXA}_candidates.gi ${OUTPREFIXA}_candidates_ext.lengths ${OUTPREFIXA}_candidates_ext_annotation.csv ${OUTPREFIXA}_candidates_counts.csv 0 |awk 'NR>1' > ${OUTPREFIXA}_summarySelf.csv; fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXB}_summarySelf.csv\" ] ; then ${JAVA} -jar ${KMERUTILS} leftjoin ${OUTPREFIXB}_candidates.gi ${OUTPREFIXB}_candidates_ext.lengths ${OUTPREFIXB}_candidates_ext_annotation.csv ${OUTPREFIXB}_candidates_counts.csv 0 |awk 'NR>1' > ${OUTPREFIXB}_summarySelf.csv; fi" >> ${SHFILE}
echo "if [ ! -e \"${OUTPREFIXA}_summaryBOTH.csv\" ] ; then" >> ${SHFILE}
echo "	join -1 2 -2 1 <(sort -k 2,2 ${OUTPREFIXA}_summarySelf.csv) <(sort -k 1,1 ${OUTPREFIXB}_summarySelf.csv) |awk -vKMER=${KMERSIZE} -v OFS=\"\\t\" ' {print \$2,\$1,\$4,\$5,\$6,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,\$23,\$18+\$22,\$17+\$23,\$4-length(\$16)-2*(KMER-1)}'|awk -vOFS=\"\\t\" -vKMER=${KMERSIZE} 'BEGIN {print \"A.id\\tB.id\\talignment length\\tmismatches\\tgaps\\tA.length\\tA.annotation\\tA.EMS\\tA.Acount\\tA.Bcount\\tB.length\\tB.annotation\\tB.EMS\\tB.Acount\\tB.Bcount\\tmirror count\\tsupport count\\ttolerance\\tquality\"} {qual=\"\"} \$16<=KMER && \$18>-4 {qual=\"high quality\"} {print \$0,qual}' > ${OUTPREFIXA}_summaryBOTH.csv" >>  ${SHFILE}
echo "	cp ${OUTPREFIXA}_summaryBOTH.csv renamed_${OUTPREFIXA}_summaryBOTH.csv" >> ${SHFILE}
echo "	perl -p -i -e 's/^(\d+)\t(\d+)/${OUTPREFIXA}_\$1\t${OUTPREFIXB}_\$2/' renamed_${OUTPREFIXA}_summaryBOTH.csv" >> ${SHFILE}
echo "else " >> ${SHFILE}
echo "	echo \"# Skipping production of ${OUTPREFIXA}_summaryBOTH.csv - already exisits\" " >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "if [ ! -f ${OUTPREFIXA}_candidates_ext.fa ] ; then" >> ${SHFILE}
echo "	echo \"# Error: File does not exist: ${OUTPREFIXA}_candidates_ext.fa\"" >> ${SHFILE}
echo "	exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "if [ ! -f ${OUTPREFIXB}_candidates_ext.fa ] ; then" >> ${SHFILE}
echo "  # Error: File does not exist: ${OUTPREFIXB}_candidates_ext.fa" >> ${SHFILE}
echo "  exit 1" >> ${SHFILE}
echo "fi" >> ${SHFILE}

echo "if [ $DO_GENOME_ANNOTATION = yes ]; then " >> ${SHFILE}
echo "  ${BLASTN} -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv" >> ${SHFILE}
echo "  ${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXA}_candidates_ext_annotation.csv 2 > ${OUTPREFIXA}_blastMutatedHits" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_blastMutatedHits | ${JAVA} -jar ${KMERUTILS} topResultPipe 1 > ${OUTPREFIXA}_topResultPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_topResultPipe | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXA}_candidates_ext_annotation.csv 1 > ${OUTPREFIXA}_blastTransferAnnotationToChrPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_blastTransferAnnotationToChrPipe | sort -k 3,3 -g | sort -k 2,2 -s > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv | awk '{print \$2\"_\"\$3\"\\t\"\$1} ' > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.trans" >> ${SHFILE}
echo "  rm -f ${OUTPREFIXA}_candidates_ext_mutated_raw.csv" >> ${SHFILE}
echo "  ${BLASTN} -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv" >> ${SHFILE}
echo "  ${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXB}_candidates_ext_annotation.csv 2 > ${OUTPREFIXB}_blastMutatedHits" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_blastMutatedHits | ${JAVA} -jar ${KMERUTILS} topResultPipe 1 > ${OUTPREFIXB}_topResultPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_topResultPipe | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXB}_candidates_ext_annotation.csv 1 > ${OUTPREFIXB}_blastTransferAnnotationToChrPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_blastTransferAnnotationToChrPipe | sort -k 3,3 -g | sort -k 2,2 -s > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv | awk '{print \$2\"_\"\$3\"\\t\"\$1} ' > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.trans" >> ${SHFILE}
echo "  rm -f ${OUTPREFIXB}_candidates_ext_mutated_raw.csv" >> ${SHFILE}
echo "  perl -p -i.orig -e 's/^/${OUTPREFIXA}_/' ${OUTPREFIXA}_blastTransferAnnotationToChrPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_blastTransferAnnotationToChrPipe | sort -k 3,3 -g | sort -k 2,2 -s > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv | awk '{print \$2\"_\"\$3\"\\t\"\$1} ' > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.trans" >> ${SHFILE}
echo "  perl -p -i.orig -e 's/^/${OUTPREFIXB}_/' ${OUTPREFIXB}_blastTransferAnnotationToChrPipe" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_blastTransferAnnotationToChrPipe | sort -k 3,3 -g | sort -k 2,2 -s > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv" >> ${SHFILE}
echo "  cat ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv | awk '{print \$2\"_\"\$3\"\\t\"\$1} ' > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.trans" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv >> ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv" >> ${SHFILE}
echo "  grep \"mitochondria\" ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.mitochondria.csv" >> ${SHFILE}
echo "  sort -nk 3,3 ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.mitochondria.csv > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.mitochondria.sorted.csv" >> ${SHFILE}
echo "  grep \"chloroplast\" ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chloroplast.csv" >> ${SHFILE}
echo "  sort -nk 3,3 ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chloroplast.csv > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chloroplast.sorted.csv" >> ${SHFILE}
echo "  grep -v \"chloroplast\" ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv | grep -v \"mito\" > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chromosome.csv" >> ${SHFILE}
echo "  sort -nk 2,2 -nk 3,3 ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chromosome.csv > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chromosome.sorted.csv" >> ${SHFILE}
echo "  cat ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.chromosome.sorted.csv | awk -F\" \" '{ print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5}' > ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.sorted.chromosome.sum.cvs" >> ${SHFILE}
echo "  #cp ${OUTPREFIXA}-${OUTPREFIXB}_candidates_ext_genomeAnnotation.sorted.chromosome.sum.cvs ../" >> ${SHFILE}
echo "  perl -p -i.orig -e 's/^>(\d+)/>${OUTPREFIXA}_\$1_${OUTPREFIXB}/' ${OUTPREFIXA}_candidates_ext.fa" >> ${SHFILE}
echo "  perl -p -i.orig -e 's/^>(\d+)/>${OUTPREFIXB}_\$1_${OUTPREFIXA}/' ${OUTPREFIXB}_candidates_ext.fa" >> ${SHFILE}
echo "	if [ ! -f ${OUTPREFIXA}_candidates_ext.fa ] ; then" >> ${SHFILE}
echo "  	echo \"# Error: File does not exist: ${OUTPREFIXA}_candidates_ext.fa" >> ${SHFILE}
echo "  	echo \"# please check pipeline command manually" >> ${SHFILE}
echo "  	exit 1" >> ${SHFILE}
echo "	fi" >> ${SHFILE}
echo "	if [ ! -f ${OUTPREFIXB}_candidates_ext.fa ] ; then" >> ${SHFILE}
echo "  	echo \"# Error: File does not exist: ${OUTPREFIXB}_candidates_ext.fa" >> ${SHFILE}
echo "  	echo \"# please check pipeline command manually" >> ${SHFILE}
echo "  	exit 1" >> ${SHFILE}
echo "	fi" >> ${SHFILE}
echo "fi" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# End of Part 7" >> ${SHFILE}
echo "# " >> ${SHFILE}
echo "# " 
echo "# " 
echo "####================================================================####"
echo "#### To run on command line do:                                     ####"
echo "cd ${WORKINGDIR}; sh ${SHFILE} 1>${SHSCRIPTBASENAME}.log 2>${SHSCRIPTBASENAME}.err &"
echo "####----------------------------------------------------------------####"
echo "####" 
echo "####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++####"
echo "#### To run on Open Lava specifying larger servers explicitly do:   ####"
echo "cd ${WORKINGDIR}; NAME=${MUTANTA}-${MUTANTB}; bsub -n ${CORES} -m \"genome7.pfr.co.nz genome10.pfr.co.nz\" -J \"\${NAME}\" -o ${WORKINGDIR}/ol.\${NAME}.niks.stdout -e ${WORKINGDIR}/ol.\${NAME}.niks.stderr < ${SHFILE}"
echo "####----------------------------------------------------------------####"
echo "####" 
echo "####================================================================####"
echo "#### To run on Open Lava specifying servers by resource limits (memory > 300 Gb) do: ####"
echo "cd ${WORKINGDIR}; NAME=${MUTANTA}-${MUTANTB}; bsub -n ${CORES} -R \"mem > 300000 order[ut]\" -J \"\${NAME}\" -o ${WORKINGDIR}/ol.\${NAME}.niks.stdout -e ${WORKINGDIR}/ol.\${NAME}.niks.stderr < ${SHFILE}"
echo "####----------------------------------------------------------------####"

exit 0
