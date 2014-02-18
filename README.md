#Introduction
_pecnv_ is a collection of C++ programs used for the detection of copy-number variation (CNV) by means of paired-end mapping of short reads.  

This pipeline corresponds to the clustering algorithm described and used in

>Rogers, R. L. _et al._ Landscape of standing variation for tandem duplications in Drosophila yakuba and Drosophila simulans

A preprint of the manuscript is available on [arXiv](http://arxiv.org/abs/1401.7371).

This software also contains the C++ portions of the transposable element (TE) detection pipeline descrbed in:
> Cridland, J.M., S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two \textit{Drosophila} QTL Mapping Resources  Molecular Biology and Evolution 30: 2311-2327, which is available [here](http://mbe.oxfordjournals.org/content/30/10/2311.full)

The full version of the TE detection pipeline, along with test data, is available from the Thornton lab [website](http://www.molpopgen.org/tepipeline/line99_example.tar.gz).  Users interested primarily in TE detection are encouraged to download that archive and study how it works.

This code has been used by the lab for detecting tandem duplications from short-read data, and the accompanying perl/shell scripts are intended to walk a user through doing such an analysis.

#Installation

##Dependencies

This software requires the following libraries to compile:

1. [libsequence](http://www.github.com/molpopgen/libsequence) - version 1.8.0 or greater
2. [boost](http://www.boost.org) - version 0.5.3 or greater

##Compilation

On most systems, simply type "make" to compile the programs.  On systems where the dependencies may be installed in funny locations, edit the Makefile accordingly.

##Installation

Either copy the binaries to somewhere in your user's path, or add the directory where you compiled them to your user's path.  To copy them into your user's ~/bin directory on a Linux machine, this will do the trick:

find . -perm -111 -type f -maxdepth 1 -exec cp "{}" ~/bin \;

The above command is executed via the install2home.sh script included with this package

##Using the programs to detect CNVs (not TEs, though -- see above for our pipeline to do that).

The easiest way to use the programs will be to run the master script, _pecnv.pl_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.
2. [samtools](http://samtools.sourceforge.net/)

__Note:__ the _pecnv.pl_ script uses ForkManager to attempt to make more effective use of CPU resources during several of the steps.  Make sure that this perl module is installed on your system!

###Processing the reference genome

I assume that you have a reference genome in fasta format.  Let's call it reference.fasta.  Different systems have different conventions for how chromosomes are named, so we'll just rename them 0 through n - 1 using the perl script _processRef.pl_, which works as follows:

processRef.pl reference.fasta new_reference_name.fasta

The script also runs "bwa index" on new_reference_name.fasta.

###Running the master script:

_pecnv.pl_ takes the following options:

> pecnv.pl -outdir pecnv_output -minqual 30 -mismatches 3 -gaps 0 -infile (infilename) -sample 0 -cpu 32 -ref (reference_fasta_filename)

In the above line, default values are shown for each option where the exist and values in parentheses must be provided by the user.  The definition of each argument is:

1. outdir = the name of the directory to write the pipeline output.  This is ALL files, including bam files, etc.
2. minqual = minimum mapping quality to consider a read in the CNV clustering
3. mismatches = max number of mismatches to allow in a read alignment to consider it for CNV clustering
4. gaps = max number of gaps to allow in a read alignment to consider it for CNV clustering
5. infile = an input file containing names of FASTQ files.  They must be in order of left read, right read, left read, right read, for all read files from a single sample.
6. sample = an integer that will be used to identify this sample.  Must be non-negative.
7. cpu = number of CPU to use
8. ref = name of the fasta file containing the renamed reference.  This should be "new_reference_name.fasta" in the above example (or the full path to that reference if it is not present in pwd/cwd).

###Running the pipeline on an Open Grid Engine (OGE) system (formerly known as Sun Grid Engine, or SGE)

The _pecnv.pl_ master script described above attempts to make the best possible use of CPU resources on a multi-core desktop machine.  However, the pipeline can be sped up dramatically with the use of a proper compute cluster with a good scheduling system.  At UCI, we use Open Grid Engine on our [cluster](http://hpc.oit.uci.edu).  I have written a basic [tutorial](http://hpc.oit.uci.edu/~krthornt/BioClusterGE.pdf) on how to use such a system.  The tutorial is somewhat specific to the UCI cluster in the details, but the main concepts are generic with respect to OGE systems.

The pecnv packages comes with a perl script called _gridify.pl_.  This script creates a set of OGE scripts that will execute the pipeline.  The scripts are designed such that, when submitted to the queue, later steps are held in the queue until previous steps are finished, and each step requests the resources necessary to complete it.  This allows for sample with multiple fastq files to use multiple nodes (or parts of nodes) for alignment, or you can submit the scripts for multiple samples to the queue.

This script is known to produce a workflow that runs successfuly on the UCI cluster, which is maintained by a group of excellent IT guys who have configured the OGE setup nicely.  The script assumes the following about your OGE system (it probably also assumes things that I've forgotten to list here...):

1. If an OGE scripts assumes that cwd and pwd are the same as the directory from which the script was submitted, unless the script specifically executes a _cd_ command
2. The system allows for name-based job dependencies via the -N and -hold_jid options to the OGE _qsub_ command

###The output

The (interesting) output from _pecnv.pl_ is the following:

1. div.gz
2. par.gz
3. ul.gz

The above correspond to clusters of reads mapping in divergent orientation, parallel orientation, and read pairs mapping to different chromosomes, respectively.

The format of the output files is as follows:

id = Event identification number (arb. integer)<br>
chrom1 = Chromosome number in reference where the first read cluster is<br>
coverage = Number of read pairs supporting the event<br>
strand1 = Strand of first read cluster.  0 = plus, 1 = minus<br>
start1 = Start position of first read cluster.  
stop1 = Stop position of second read cluster. 
chrom2 = Chromosome number in reference where the second read cluster is<br>
start2 = Start position of second read cluster.  
stop2 = Stop position of second read cluster. 
reads = Pipe-separated (the pipe is the | character) list of the read pairs supporting the event.  Format is readPairName;start,stop,strand,start,stop,strand, where the last two values are for the two reads in the pair.


