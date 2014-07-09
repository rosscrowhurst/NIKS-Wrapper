#!/usr/bin/perl -w
use strict;

=head1 DESCRIPTION

Warning this script is designed to work with TAIR10 and annotation + fasta
from 19 genomes. It has only been tested with ler-0.

This script takes output file from

1. Our NIKS- kmerPipeline.sh wrapper script which is an alternative way to run the 
   orignal NIKS pipeline (available from http://sourceforge.net/projects/niks/) from
   the following paper:
   
   Nordström KJ1, Albani MC, James GV, Gutjahr C, Hartwig B, Turck F, Paszkowski U, 
   Coupland G, and Schneeberger K 
   
   Mutation identification by direct comparison of whole-genome sequencing data from 
   mutant and wild-type individuals using k-mers.
   
   Nat Biotechnol. 2013 Apr;31(4):325-30. doi: 10.1038/nbt.2515.
   
   
2. The descriptions from the TAIR10 protein fasta file

3. The locations of locus from 19 genomes consolidated 
   (e.g. Ler-0) annotations


This script outputs a file with the NIKS identified locations, the Arabidopsis genome 
AGI code and the AGI locus description.

The current implementation works only on a single chromosome at a time.

=head1 SYNOPSIS

	niksMergeAndGetDescriptionFromTairDesc.pl -c=chrNum -tair_desc=tairDecFile \
		-refChr=ChrFasta -niks_locations=file -ler0Anno=ler0_consolidated_annotations

=head2 Example Bash Shell Workflow

MUTANTA=dis9
MUTANTB=dis58
TAIR_DESC_FILE=TAIR10.desc
REFERENCE_DIR=/input/ComparativeDataSources/Arabidopsis/thaliana/Landsberg/erecta/Ler-0

for chr in 1 2 3 4 5
do
	niksMergeAndGetDescriptionFromTairDesc.pl \
	-c=${chr} \
	-tair_desc=${TAIR_DESC_FILE} \
	-refChr=${REFERENCE_DIR}/Chr${CHR}.fa \
	-niks_locations=${MUTANTA}-${MUTANTB}_candidates_ext_genomeAnnotation.sorted.chromosome.sum.cvs 
	-ler0Anno=../../consolidated_annotation.Ler_0.gff3 \
	-niks_summaryBoth=renamed_${MUTANTA}_summaryBOTH.csv \
	-o=NIKS_CHR${CHR}_${MUTANTA}-${MUTANTB}_SUMMARY.txt \
	-refMutantName=${MUTANTA} \
	-subjMutantName=${MUTANTB}

=head1 NIKS files

=head2 candidates_ext_genome 

This file will have been produced after running the NIKS pipeline with our alternate wrapper. 

It is named according to the following format:

	${MUTANTA}-${MUTANTB}_candidates_ext_genome

The format of the file is 5 columns as follows where
        column 1 is the NIKS contig name
        column 2 is chromosome number
        column 3 is location of variation on chromosome
        column 4 and 5 are reference and alternative alleles relative to the named contig

=head3 Example NIKS location file format

        dis9_1845       2       1890    G       T
        dis9_348798     2       1890    G       T
        dis9_9492       2       1930    T       A
        dis9_13458      2       2653    A       C
        dis15_25744     2       2654    C       A
        dis9_13458      2       2654    C       T
        dis15_25744     2       2655    T       C
        dis9_2082168    2       3410    G       T



The input data will be mangled such that the order of the alleles will be reference followed by alternate
where the reference allele will be the genome reference allele.

=head2 ${MUTANTA}_topResultPipe

Part of original NIKS distribution outputs.

This file will contain data like the following format:

11      chloroplast     99.07   107     1       0       1       107     106438  106544  1.0E-47 188.0
14      chloroplast     98.25   57      1       0       1       57      100658  100714  7.0E-21 98.7


=head1 renamed_${MUTANTA}_summaryBOTH.csv

An output file specific to our alternate NIKS wrapper. It really is just renaming of the
original NIKS pipeline output file to reflect the names of mutants in the contig names.

Note: for this MUTANTA=dis9 and MUTANTB=dis58

This file should contain the following format:

        dis9_86404      dis15_10013     55      2       0       64      27:A>G,28:C>T   RM      238     0       278     117:A>G,118:C>T RM      47      64      47      302     -7
        dis9_19649      dis15_10016     56      2       0       66      26:T>G,31:A>T   UU      406     0       58      26:G>T,31:T>A   UU      98      136     98      542     -6

=head1 SINGLE CHROMOSOME REFERENCE FASTA

Only use one chromosome or scaffold as a reference. Do not use
a multi sequence fasta file - it will not work!!

=head TAIR 10 DESCRIPTION FILE

Obtain the appropriate data set from 

	ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/

=head2 Example download and extract descriptions

	#!/bin/sh
	cd /input/ComparativeDataSources/Arabidopsis/thaliana/TAIR10
	for file in TAIR10_cdna_20101214_updated TAIR10_cdna_20110103_representative_gene_model_updated TAIR10_cds_20101214_updated TAIR10_cds_20110103_representative_gene_model_updated TAIR10_exon_20101028 TAIR10_intergenic_20101028 TAIR10_intron_20101028 TAIR10_pep_20101214_updated TAIR10_pep_20110103_representative_gene_model_updated TAIR10_seq_20101214_updated TAIR10_seq_20110103_representative_gene_model_updated TAIR10_3_utr_20101028 TAIR10_5_utr_20101028 TAIR10_bac_con_20101028
	do
		wget -c ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/$file
		grep ">" $file > $file.desc
	        perl -p -i -e 's/^>//' $file.desc
		mv $file ${file}.fasta
	done

=head2 Extract descriptions manually

To product this file simple do the following:

        grep ">" TAIR10.mRNA.fasta > TAIR10.desc
        perl -p -i -e 's/^>//' TAIR10.desc

This should produce the following type of output:

        AT1G51370.2 | Symbols: |F-box/RNI-like/FBD-like domains-containing protein| chr1:19045615-19046748 FORWARD LENGTH=346
        AT1G50920.1 | Symbols: |Nucleolar GTP-binding protein| chr1:18870555-18872570 FORWARD LENGTH=671

This must have the chr1:19045615-19046748 style location data for parsing.

=head1 19 Genome Project files

=head2 Genome sequence files

Obtain genome sequence for ecotype from:

	http://mus.well.ox.ac.uk/19genomes/fasta/UNMASKED/

e.g

	http://mus.well.ox.ac.uk/19genomes/fasta/UNMASKED/ler_0.v7.fas

=head3 Obtain associated annotation file for ecotype

Get from:

	http://mus.well.ox.ac.uk/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/

e.g

	http://mus.well.ox.ac.uk/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/consolidated_annotation.Ler_0.gff3.bz2

Uncompress before use

=head1 refMutantExtFastaFile and subjMutantExtFastaFile

These are the two candidates_ext.fa files produced by the NIKS pipeline
for the mutant used as a reference and the mutant used as the subject in
each pairwise comparison.

=head1 Use at your own risk



=cut

#------------------------------------
# Defaults are for our local system - change or pass on command line 
#------------------------------------
my $tairDescriptionFile = 
	"/input/ComparativeDataSources/Arabidopsis/thaliana/TAIR10/TAIR10_cdna_20101214_updated.desc";
my $ler0ConsolidatedAnnotationsGFF3 = 
	"/input/ComparativeDataSources/Arabidopsis/thaliana/Landsberg/erecta/Ler-0/19genomes/annotations/consolidated_annotation_9.4.2011/consolidated_annotation.Ler_0.gff3";


#------------------------------------
# Vars
#------------------------------------
my $niksLocationsFile = "";
my $renamedSummaryFile = "";
my $chr = "";
my %summaryBothCVSData = ();
my %tair10 = ();
my %niksData = ();
my $outfile = "";
my $refMutantName = "";
my $subjMutantName = "";
my %chrLocationToContigsLookup = ();
my $referenceGenomeFastaFile = "";
my $refMutantExtFastaFile = "";
my $subjMutantExtFastaFile = "";
my %refGenomeFastaData = ();
my %refMutantExtFastaData = ();
my %subjMutantExtFastaData = ();
my $mutantAtopResultPipeFile = "";
my %reftopResultPipeData = ();

#------------------------------------
# Command line args
#------------------------------------
foreach my $arg (@ARGV)
{
        ($arg =~ m/^-c=(\d+)$/) and $chr = $1;
        ($arg =~ m/^-tair_desc=(.+)$/) and $tairDescriptionFile = $1;
        ($arg =~ m/^-ler0Anno=(.+)$/) and $ler0ConsolidatedAnnotationsGFF3 = $1;
        ($arg =~ m/^-refChr=(.+)$/) and $referenceGenomeFastaFile = $1;
        ($arg =~ m/^-niks_locations=(.+)$/) and $niksLocationsFile = $1;
        ($arg =~ m/^-niks_summaryBoth=(.+)$/) and $renamedSummaryFile = $1;
        ($arg =~ m/^-o=(.+)$/) and $outfile = $1;
        ($arg =~ m/^-refMutantName=(.+)$/) and $refMutantName = $1;
        ($arg =~ m/^-subjMutantName=(.+)$/) and $subjMutantName = $1;
        ($arg =~ m/^-refMutantExtFastaFile=(.+)$/) and $refMutantExtFastaFile = $1;
        ($arg =~ m/^-subjMutantExtFastaFile=(.+)$/) and $subjMutantExtFastaFile = $1;
        ($arg =~ m/^-mutantAtopResultPipeFile=(.+)$/) and $mutantAtopResultPipeFile = $1;
}

#------------------------------------
# Check opts
#------------------------------------
if (($#ARGV < 1) or (!$refMutantName) or (!$subjMutantName))
{
        print "USAGE: $0 -c=chrNum -tair_desc=tairDecFile -refChr=ChrFasta -niks_locations=file \\ \n";
        print "  -ler0Anno=ler0_consolidated_annotations -niks_summaryBoth=renamed_disX_summaryBOTH.csv \\ \n";
        print "  -o=outfile -refMutantName=name -subjMutantName=name -mutantAtopResultPipeFile=name \\ \n";
        print "  -refMutantExtFastaFile=refMutantExtFastaFile -subjMutantExtFastaFile=subjMutantExtFastaFile\n";
        exit(0);
}
($outfile) or $outfile = $niksLocationsFile . "parsed.out";
$refMutantName = lc($refMutantName);
$subjMutantName = lc($subjMutantName);

my $oldSep = $/;
$/=">";
#--------------------------------------------------------
# Load refMutantExtFastaFile File (x_candidates_ext.fa)
#--------------------------------------------------------
select STDERR;
print STDERR "Loaded $refMutantExtFastaFile\n";
open(IN, "<$refMutantExtFastaFile") or die "Can not open refMutantExtFastaFile: $refMutantExtFastaFile $!\n";
while( my $record = <IN>)
{
        chomp $record;
        next if $record eq "";
        next if $record eq ">";

        my ($nameLine, @sequence) = split/\n/, $record;
        my $sequence = join("", @sequence);
        ($sequence) or die "no sequence for $record\n";
        $sequence =~ s/>//gs;
        my ($name, @rest) = split/\s+/, $nameLine;
        my $description = join(" ", @rest);
        # Name will be of the format dis9_906_dis58 (only want dis9_906 part)
        my @name = split/_/, $name;
        $name = $refMutantName . "_" . $name[1];
        $refMutantExtFastaData{$name}{'desc'} = $description;
        $refMutantExtFastaData{$name}{'seq'} = $sequence;
}
close(IN);
$/ = $oldSep;
#--------------------------------------------------------
# Load subjMutantExtFastaFile File
#--------------------------------------------------------
$/=">";
select STDERR;
print STDERR "Loaded $subjMutantExtFastaFile\n";
open(IN, "<$subjMutantExtFastaFile") or die "Can not open subjMutantExtFastaFile: $subjMutantExtFastaFile $!\n";

while( my $record = <IN>)
{
        chomp $record;
        next if $record eq "";
        next if $record eq ">";

        my ($nameLine, @sequence) = split/\n/, $record;
        my $sequence = join("", @sequence);
        ($sequence) or die "no sequence for $record\n";
        $sequence =~ s/>//gs;
        my ($name, @rest) = split/\s+/, $nameLine;
        my $description = join(" ", @rest);
        # Name will be of the format dis1_906_dis9 (only want dis1_906 part)
        my @name = split/_/, $name;
        $name = $subjMutantName . "_" . $name[1];
        $subjMutantExtFastaData{$name}{'desc'} = $description;
        $subjMutantExtFastaData{$name}{'seq'} = $sequence;
}
close(IN);

$/ = $oldSep;
#--------------------------------------------------------
# Load reference genome File
#--------------------------------------------------------
$/=">";
select STDERR;
print STDERR "Loaded $referenceGenomeFastaFile\n";
open(IN, "<$referenceGenomeFastaFile") or die "Can not open referenceGenomeFastaFile: $referenceGenomeFastaFile $!\n";

while( my $record = <IN>)
{
        chomp $record;
        next if $record eq "";
        next if $record eq ">";

        my ($nameLine, @sequence) = split/\n/, $record;
        my $sequence = join("", @sequence);
        ($sequence) or die "no sequence for $record\n";
        $sequence =~ s/>//gs;
        my ($name, @rest) = split/\s+/, $nameLine;
        my $description = join(" ", @rest);
        my @refBases = split('', $sequence);
        my $splitBases = scalar( @refBases);
        $refGenomeFastaData{$name}{'desc'} = $description;
        $refGenomeFastaData{$name}{'seq'} = $sequence;
        $refGenomeFastaData{$name}{'seqLen'} = length($sequence);
        $refGenomeFastaData{$name}{'refBases'} = \@refBases;
        $refGenomeFastaData{$name}{'splitBases'} = $splitBases;
	select STDERR;
	print STDERR "Splitting $name in to bases (expect $refGenomeFastaData{$name}{'seqLen'})\n";
	print STDERR "Split $name in to $splitBases bases (expected $refGenomeFastaData{$name}{'seqLen'})\n";
}
close(IN);

$/ = $oldSep;
#--------------------------------------------------------
# Load Summary File
#--------------------------------------------------------
select STDERR;
print STDERR "Loading $renamedSummaryFile file\n";
open(IN, "<$renamedSummaryFile") or die "Can no open renamedSummaryFile $renamedSummaryFile file $!\n";
my $lineCount = 0;
while (my $line = <IN>)
{
        chomp $line;
        $lineCount++;
        next if ($line =~ m/^A.id\t/);
        my @lineData = split /\t/, $line;
        my $key = $lineData[0] . $lineData[1];
        if (defined $summaryBothCVSData{$key})
        {
                $summaryBothCVSData{$key} .= "\t$line";
        }
        else
        {
                $summaryBothCVSData{$key} = $line;
        }
}
close(IN);

$/ = $oldSep;

#--------------------------------------------------------
# Load mutantAtopResultPipeFile File (e.g dis1_topResultPipe)
#--------------------------------------------------------
select STDERR;
print STDERR "Loading $mutantAtopResultPipeFile file\n";
open(IN, "<$mutantAtopResultPipeFile") or die "Can not open $mutantAtopResultPipeFile $!\n";
while (my $line = <IN>)
{
        chomp $line;
        my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore) = split/\t/, $line;
        my $queryContigName = $refMutantName . "_" . $queryId;
        $reftopResultPipeData{$queryContigName}{'queryId'} = $queryId;
        $reftopResultPipeData{$queryContigName}{'subjectId'} = $subjectId;
        $reftopResultPipeData{$queryContigName}{'percIdentity'} = $percIdentity;
        $reftopResultPipeData{$queryContigName}{'alnLength'} = $alnLength;
        $reftopResultPipeData{$queryContigName}{'mismatchCount'} = $mismatchCount;
        $reftopResultPipeData{$queryContigName}{'gapOpenCount'} = $gapOpenCount;
        $reftopResultPipeData{$queryContigName}{'queryStart'} = $queryStart;
        $reftopResultPipeData{$queryContigName}{'queryEnd'} = $queryEnd;
        $reftopResultPipeData{$queryContigName}{'subjectStart'} = $subjectStart;
        $reftopResultPipeData{$queryContigName}{'subjectEnd'} = $subjectEnd;
        $reftopResultPipeData{$queryContigName}{'eVal'} = $eVal;
        $reftopResultPipeData{$queryContigName}{'bitScore'} = $bitScore;

        my @slice = @refBases[($subjectStart -1) .. ($subjectEnd -1)];
        $reftopResultPipeData{$queryContigName}{'genome_sequence'} = join("", @slice);


}
close(IN);

#--------------------------------------------------------
# Load Single Sequence Reference Fasta File
#--------------------------------------------------------
select STDERR;
print STDERR "Loading $niksLocationsFile\n";
open(IN, "<$niksLocationsFile") or die "Can not open NIKS location file $niksLocationsFile $!\n";
while( my $line = <IN>)
{
        chomp $line;
        #dis15_16138     1       11270679        G       A
        #dis9_663        1       11270679        A       G

        #my ($contigName, $atChr, $chrLocation, $niksSubjectBase, $niksReferenceBase) = split /\t/, $line;
        my ($contigName, $atChr, $chrLocation, $niksReferenceBase, $niksSubjectBase) = split /\t/, $line;
        next if ($atChr != $chr);
        $niksReferenceBase = uc( $niksReferenceBase);
        $niksSubjectBase = uc($niksSubjectBase);

        if ($contigName =~ m/^$refMutantName/)
        {
                $chrLocationToContigsLookup{$chrLocation}{'reference_contig'} = $contigName;
        }
        elsif ($contigName =~ m/^$subjMutantName/)
        {
                $chrLocationToContigsLookup{$chrLocation}{'subject_contig'} = $contigName;
        }
        else
        {
                die "$contigName does not begin with $refMutantName or $subjMutantName\n";
        }

        #next if ($contigName =~ m/^Dis15D/);
        #$contigName =~ s/Dis15//;
        $niksData{$contigName}{'location'} = $chrLocation;
        my $genomeBase = uc($refBases[$chrLocation - 1]);
       # my $genomeBaseReverseComp = $genomeBase;
       # $genomeBaseReverseComp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

        my $niksReferenceBaseRevComp = uc($niksReferenceBase);
        #$niksReferenceBaseRevComp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

        my $niksSubjectBaseRevComp = uc($niksSubjectBase);
        #$niksSubjectBaseRevComp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

        $niksData{$contigName}{'niksSubjectBase'} = $niksReferenceBase;
        $niksData{$contigName}{'niksReferenceBase'} = $niksSubjectBase;
        $niksData{$contigName}{'genomeBase'} = $genomeBase;
        $niksData{$contigName}{'orientation'} = ".";

}
close(IN);
my %lerAnnotation = ();

select STDERR;
print STDERR "Loading $ler0ConsolidatedAnnotationsGFF3\n";

open(IN, "<$ler0ConsolidatedAnnotationsGFF3") or die "Can not open ler0ConsolidatedAnnotationsGFF3 location file $ler0ConsolidatedAnnotationsGFF3 $!\n";
while( my $line = <IN>)
{
        chomp $line;
        next if ($line =~ m/^#/);
        my @gff3 = split/\t/, $line;
        my $localChr = "Chr$chr";
        next if ($gff3[0] ne $localChr);
        if ($gff3[2] eq "mRNA")
        {
                my $agi = "";

                if ($gff3[8] =~ m/^ID=Transcript:(AT\dG\d+\.\d+);Parent=/)
                {
                        $agi = $1;
                }
                elsif ($gff3[8] =~ m/;Parent=Gene:(AT\dG\d+)/)
                {
                        $agi = $1 . ".1";
                }

                if (! $agi)
                {
                        die "No AGI from $line\n";
                }
                foreach my $i ($gff3[3] .. $gff3[4])
                {
                        $lerAnnotation{$i} = $agi
                }
        }
}
##gff-version 3
##Seqid Source  Type    Start   End     Score   Phase   Attributes
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   gene    3644    5912    .       +       .       ID=Gene:AT1G01010;explain="take_identical "
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   mRNA    3644    5912    1       +       .       ID=Transcript:AT1G01010.1;Parent=Gene:AT1G01010;transcript_origin=1;explain="take_identical "
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   five_prime_UTR  3644    3772    .       +       .       ID=five_prime_UTR:AT1G01010.1.1;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     3773    3926    .       +       .       ID=CDS:AT1G01010.1.1;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     4009    4289    .       +       .       ID=CDS:AT1G01010.1.2;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     4499    4618    .       +       .       ID=CDS:AT1G01010.1.3;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     4719    5108    .       +       .       ID=CDS:AT1G01010.1.4;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     5187    5339    .       +       .       ID=CDS:AT1G01010.1.5;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   CDS     5452    5643    .       +       .       ID=CDS:AT1G01010.1.6;Parent=Transcript:AT1G01010.1
#Chr1   mgene/tair10 consolidated Ler_0 (nss)   three_prime_UTR 5644    5912    .       +       .       ID=three_prime_UTR:AT1G01010.1.1;Parent=Transcript:AT1G01010.1

select STDERR;
print STDERR "Loading $tairDescriptionFile\n";
open(IN, "<$tairDescriptionFile") or die "Can not open $tairDescriptionFile\n";
$lineCount = 0;
while (my $line = <IN>)
{
        #print $line;

        chomp $line;
        $lineCount++;
        #AT4G05475.1|Symbols:|RNI-like superfamily protein| chr4:2765962-2767957 REVERSE LENGTH=309
        #AT4G15165.1|Symbols:|N-terminal nucleophile aminohydrolases (Ntn hydrolases) superfamily protein| chr4:8649460-8650212 FORWARD LENGTH=208
        #AT4G15160.2|Symbols:|Bifunctional inhibitor/lipid-transfer protein/seed storage 2S albumin superfamily protein| chr4:8646192-8647019 FORWARD LENGTH=193
        #AT4G22305.1|Symbols:|alpha/beta-Hydrolases superfamily protein| chr4:11789546-11791055 REVERSE LENGTH=228
        my ($agi, $symbols, $description, $chromDesc) = split/\|/, $line;
        $agi =~ s/^\s+//;
        $agi =~ s/\s+$//;
        $description  =~ s/^\s+//;
        $description  =~ s/\s+$//;
        $chromDesc =~ s/^\s+//;
        $chromDesc  =~ s/\s+$//;
select STDERR;
print STDERR "AGI:$agi\n";
print STDERR "DESCRIPTION:$description\n";
print STDERR "LINECOUNT: $lineCount\n";

        my ($chrPosData, @junk) = split / /, $chromDesc;
        $tair10{$agi} = $description;
        my $start = 0;
        my $stop = 0;

}
close(IN);
select STDERR;
print STDERR "Writing ouptut to $outfile\n";

open(OUT, ">$outfile") or die "Can not open outfile\n";
select STDOUT;
foreach my $contigName (sort {$niksData{$a}{'location'} <=> $niksData{$a}{'location'}} keys %niksData)
{
        my $outLine = "";
        my $chrLocation = $niksData{$contigName}{'location'};
        my $reference_contig = "";
        my $reference_contig_sequence = "";
        my $subject_contig = "";
        my $subject_contig_sequence = "";
        my $summaryData = "";
        if (exists $chrLocationToContigsLookup{$chrLocation})
        {
                if ((exists $chrLocationToContigsLookup{$chrLocation}{'reference_contig'}) and (exists $chrLocationToContigsLookup{$chrLocation}{'subject_contig'}))
                {
                        $reference_contig = $chrLocationToContigsLookup{$chrLocation}{'reference_contig'};
                        $subject_contig = $chrLocationToContigsLookup{$chrLocation}{'subject_contig'};
                        $reference_contig_sequence = $refMutantExtFastaData{$reference_contig}{'seq'};
                        $subject_contig_sequence = $subjMutantExtFastaData{$subject_contig}{'seq'};
                        my $key = $chrLocationToContigsLookup{$chrLocation}{'reference_contig'} . $chrLocationToContigsLookup{$chrLocation}{'subject_contig'};
                        if (exists $summaryBothCVSData{$key})
                        {
                                $summaryData = $summaryBothCVSData{$key};
                        }
                }

        }
        if (exists $lerAnnotation{$chrLocation})
        {
                my $agi = $lerAnnotation{$chrLocation};
                $outLine = join("\t",
                        $chr,
                        $chrLocation,
                        $agi,
                        $contigName,
                        $niksData{$contigName}{'niksReferenceBase'},
                        $niksData{$contigName}{'niksSubjectBase'},
                        $niksData{$contigName}{'genomeBase'},
                        $niksData{$contigName}{'orientation'},
                        $tair10{$agi},
                        $reference_contig_sequence,
                        $subject_contig_sequence,
                        $reftopResultPipeData{$contigName}{'genome_sequence'}

                );

        }
        else
        {
                $outLine = join("\t",
                        $chr,
                        $chrLocation,
                        "not in gene",
                        $contigName,
                        $niksData{$contigName}{'niksReferenceBase'},
                        $niksData{$contigName}{'niksSubjectBase'},
                        $niksData{$contigName}{'genomeBase'},
                        $niksData{$contigName}{'orientation'},
                        "-",
                        $reference_contig_sequence,
                        $subject_contig_sequence,
                        $reftopResultPipeData{$contigName}{'genome_sequence'}
                );
        }
        ($outLine) or die "No outline\n";
        ($summaryData ne "") and $outLine .= "\t$summaryData";
        select OUT;
        print OUT "#\n";
        print OUT "$outLine\n";
}
close(OUT);
select STDERR;
print STDERR "Completed\n";

exit(0);
