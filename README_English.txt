Procedure to perform mapping and following analysis of rearrangements in organelle

This procedure includes all the programs used to carry out the analysis in the papre.

"Ultra-deep sequencing reveals dramatic alteration of organellar genomes in
Physcomitrella patens due to biased asymmetric recombination"

by Odahara, Nakamura, Oshima and Sekine.

The file suffix appears in the paper .mhr and .mhmr are replaced with .hr and .hmr, 
respectively in this document and with the actual programs.
Likewise, program names in the paper such as comp_mhr and comp_mhmr 
are replaced with comp_hr and comp_hmr, respectively.

It includes many analysis.
The main part of the rearrangement analysis including the drawing of links, can 
be carried out through the section 0.0-3.2.

Lines with '%' in the first column is the command line to input,
and /DBdir/ is the directory you put the reference files.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

0. Preparing the mapping: Indexing and Reference list

---------------------------------------------------------------------------------------------

0.0 Indexing reference sequence (mkindex)

Our mapping software immap uses simple fixed length index.
To make the index, use 'mkindex' with specified index length.
For instance, to make an index file for the reference file 'foo.fasta' with index length 14, 
type

% mkindex -14 foo.fasta

An index file 'foo.fasta.index14' will be created.
The index file should stay in the same directory as the reference file (foo.fasta).
The most appropreate index length depends on the length of reference (10~15).

---------------------------------------------------------------------------------------------

0.1 Prepareing reference list file (ref.fasta)

Our mapping software immap only takes single fasta files as input. 
When the reference sequence consists of multiple sequence,
such as in the case of multiple chromosomes of eukaryotic genomes, 
we need to prepare a list file of single fasta files as follows.

========================== ref.list =============================
#Title line
0       /DBdir/CH_MOD.fasta                   /DBdir/CH_MOD.gb
0       /DBdir/MT_MOD.fasta                   /DBdir/MT_MOD.gb
0       /DBdir/Ppatens.fa.out.fasta
=================================================================

The first line starts with # and is a title line.
Following lines starts with 0, and are followed by one or two space separated colums.
The first column is the reference single fasta file name, and the second column is the 
Genbank format file that includes the description/annotation of the reference file.
The genbank file (.gb) may be omitted.
In our analysis, the first 2 files are for the genome sequences of organelles 
(chroloplast and mitochondorion),
and the third file contains nuclear chromosome sequences.
The multiple nuclear choromosome sequences are connected with 500 N's as linkers to form a 
single fasta format.
The index files for each fasta files should be prepared according to the procedure 
0.0 prior to the mapping.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

1 Mapping (immap, chopfastq)

---------------------------------------------------------------------------------------------

1.1 Mapping original read (100~300bp) allowing many mismatches to detect junction reads

To carry out the mapping, use immap as follows.

% immap -a 16 -i 14 -h 75 -r 10 -job i14h75r10 ref.list PpWT.fastq

Here the last argument (PpWT.fastq) is the read file, 
and the second last argument ( ref.list) is the list of reference fasta files 
described in the previous section(0.1).
Paired end reads should be merged into a single file in the order of R1-R2.
Following is the description of options with '-'.

-a   # :  Number of threads for multi thread calculations.
-i   # :  Index length. An appropreate index file should be prepared according to 0.0
-h   # :  Maximum number of mismatch per read. The value should be as many as half the 
          length of reads for junction analysis.
-r   # :  Number of search per read shifting index position. h+1 if possible. 
          10~20 if h is large.
-job ~ :  Job name. We usually put the mapping parameters as shown in the example.

The result of mapping is output in the .hit file. .hit file consists of 2 or 3 columns.
The first column is the read id. (sequentical number), and the second column is the 
mapped position of the read.
The negative value for the mapped position indicates the read is mapped as reverse 
complement, and the absolute value indicate the mapped position of the 5' end of 
the reverse complement read on the reference.
The third column (if exist) is the mapping depth which can be used to display 
the read alignment.

immap generates .hit output file for each reference file in the reference list file 
with index number starting from 0.
For instance, there are 3 reference files (CH_MOD.fasta, MT_MOD.fasta and 
Ppatens.fa.out.fasta) in the ref.list described in 0.1.
In this example immap outputs, PpWT.fastq_i14h75r10_0.hit for the hits on the first 
reference CH_MOD.fasta, and PpWT.fastq_i14h75r10_1.hit for the hits on the second 
reference MT_MOD.fasta.

By adding -map option, immap generate a text file (.map) that shows the mapping status 
of each read.
The width of .map file is 300 by default and can be altered with -mw option. 
Add option '-mw 100' to make the width of the map file 100.

---------------------------------------------------------------------------------------------

1.2 Mapping of truncated reads (chopfastq)

In our study, 150bp reads were used for the junction read analysis.
Besides, reads truncated to 50bp were used for the analysis of mapping depth, 
and rearrangement analysis using paired end reads. 
To prepare the truncated reads, we used 'chopfastq'.

% chopfastq -l 50 PpWT.fastq

This generates chp50_PpWT.fastq that contains truncated reads.
The generated reads are mapped using immap allowing upto 2 bp mismatch.

% immap -a 16 -i 14 -h 2 -r 3 -job i14h2r3  ref.list chp50_PpWT.fastq

Options and arguments are as same as the example in 1.0.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

2   Depth plot, and preparation of the .peb file used in circos-like plot 
    described later(2,3) (ispmap) Depth plot with calibration will be described in 4.0.

---------------------------------------------------------------------------------------------

2.1 Generage depth plot (histogram) and generate '.peb' file for circmaps

To generate a postscript file describing mapping depth, and to generate .peb files 
that contains mapped read population along the genome sequence, used for the 
circos-like plot (3) and calibrated depth plot (4), we use 'ispmap' as follows.

% ispmap -peb 1 -hd 50 -ds 0.02 -job i14h2r3_0 /DBdir/CH_MOD.fasta chp50_PpWT.fastq
% ispmap -peb 1 -hd 50 -ds 0.02 -job i14h2r3_1 /DBdir/MT_MOD.fasta chp50_PpWT.fastq

The last two arguments are the reference file and read file.
Unlike 'immap', the reference file is specified as a single fasta file, 
not the list of reference files.  
In our example ispmap is carried out for both chroloplast (CH_MOD.fasta) and 
mitochondorion (MT_MOD.fasta).
Please note that the suffix of the job names (string after -job) are '_0' 
for chroloplast and '_1' for mitochondorion as these suffix are applied by immap 
for each '.hit' file.

As of other options, '-peb 1' specifies the generation of read population for every 
'1' base of reference, that is used for circos-like plot.

Options '-hd 50' and '-ds 0.02' are used for the postscript depth display from 
ispmap and does not affect .peb file for circos-like plot. 
-hd   #   : Window size of depth  depth count along reference sequence
-ds   #.# : Scale factor of the Y axis (depth) of the postscript plot

---------------------------------------------------------------------------------------------

2.2 Generate base resolution mapping plot and text .map file.

The former example generates depth histogram output.
Another option can display each reads mapped underneath the reference sequence, 
showing the mismatch base as red dots. To generate this type of output, 
use '-bw 1' and '-cc' instead of '-hd 50' and '-ds 0.02'.

% ispmap -bw 1 -cc -map -job i14h75r10_0 /DBdir/CH_MOD.fasta PpWT.fastq

.map file can be generated with ispmap, if you forgot to make it with immap.
'-cc' is a coloring option (no other choice available at the moment).

Please note the output postscript file by -bw 1 option can be quite large, 
if the reference sequence is very long.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

3   Rearrangement analysis and visualization
    Paired-end read analysis, Junction read analysis, and Generation of Circos-like map. 
    (midhr, circmaps, idss, exclude_id_pairs, sort_pdist)

---------------------------------------------------------------------------------------------

3.1.0 Paired-end read analysis  (midhr, idss, exclude_id_pairs, sort_pdist, circmaps)
    
To perform rearrangement analysis using 'midhr' we use the mapping results 
of truncated reads (2.1).

% midhr -ptd 1000 -job i14h3r4_0 /DBdir/CH_MOD.fasta chp50_PpWT.fastq
% midhr -ptd 1000 -job i14h3r4_1 /DBdir/MT_MOD.fasta chp50_PpWT.fastq

Here the last two arguments are reference sequence and read data as same as ispmap (2.0).
Please note that the job name has the same suffix as ispmap 
(_0 for chloroplast and _1 for mitochondorion).
The value specified with the ptd option is the distance threshold to distinguish 
a pair of reads facing each other is normal (closer) or abnormal (distant).
There are four output files, .pdist, .pdist2, .pdist3, .pdist4.

pdist:  Mapped positions of pair reads are more than -ptd value appart. 
        (Abnormal pair, regardless of their directions)
pdist2: Mapped positions of pair reads are within -ptd value appart, 
        and faced to each other. (Normal pair)
pdist3: Mapped positions of pair reads are within -ptd value appart, 
        and pointing the same direction. (Abnormal pair)
pdist4: All abnormal pairs, pdist + pdist3.

---------------------------------------------------------------------------------------------

3.1.1 Detection of long (more than 50bp) repeats in the reference

Some of the abnormal pairs in .pdist (.pdist4) file can be artifacts caused by 
the presence of long identical sequence pairs in reference genome. 
To decrease this kind of artifact, first we prepare the list of identical sequences 
longer than threshfold length (50 bp in this example).

% idss  -i 8 -m 50 CH_MOD.fasta > CHi8m50.ids
% idss  -i 8 -m 50 MT_MOD.fasta > CHi8m50.ids

The threshold value is specified with '-m' option, and '-i' option specifies 
the index value used for the search.
Index value should be smaller than 15, and also, smaller than the threshold value. 
The idss outputs the results into standard output, so it should be redirected 
into a file (CHi8m50.ids).

---------------------------------------------------------------------------------------------

3.1.2 Removal of artifacts caused by repeats in reference

Based on the ids file obtained, we can remove some artifact using 'exclude_id_pairs'.

% exclude_id_pairs chp50_PpWT_0.pdist CHi8m50.ids > chp50_PpWT_0.excl.pdist
% exclude_id_pairs chp50_PpWT_1.pdist MTi8m50.ids > chp50_PpWT_1.excl.pdist

Here again, the output is redirected. 
The positions of read pairs suggesting the same rearrangement may have similar values, 
but usually not exactly identical.
Therefore we need to cluster the read pairs suggesting the same rearrangement 
based on the mapped positions of pair reads.

% sort_pdist chp50_PpWT_0.excl.pdist
% sort_pdist chp50_PpWT_1.excl.pdist

'sort_pdist' generates a .pe file (here we get chp50_PpWT_0.excl.pdist.pe and 
chp50_PpWT_1.excl.pdist.pe), that contains the list of clusters of the paired reads.

---------------------------------------------------------------------------------------------

3.1.3 Circos-type drawing of links by paired end information

Using this .pe file and .peb file obtained by ispmap (2.0), 
we can generate a Circos-like map as follows.

% circmaps -rpb chp50_PpWT.fastq_i14h2r3_0.peb \
           -rgb /DBdir/CH_MOD.gb \
           -rfa /DBdir/CH_MOD.fasta \
           -pe \
           -bls 1.0 \
           -psc 2 \ 
           -drm 5 \ 
           chp50_PpWT.fastq_i14h2r3_0.excl.pdist.pe

% circmaps -rpb chp50_PpWT.fastq_i14h2r3_1.peb \
           -rgb /DBdir/MT_MOD.gb \
           -rfa /DBdir/MT_MOD.fasta \
           -pe \
           -bls 1.0 \
           -psc 2 \ 
           -drm 5 \ 
           chp50_PpWT.fastq_i14h2r3_1.excl.pdist.pe

These are 2 long single line commands separated by backslash (\).
The argument in the last line is the .pe file just generated(3.1.2).
Following is the description of the options.

-rpb ~ :  Specifies the .peb file name created in 2.0 
-rgb ~ :  Specifies the gen bank(.gb) file name. Used to plot the gene position 
-rfa ~ :  Specifies the reference genome fasta file.
-pe    :  Paired end analysis reading the .pe file
-bls   :  Scale factor for the plot of predicted depth increase/decrease by links. 
          (default 1.0)
-psc   :  Scale factor for the depth plot by .peb. (default 1.0)
-drm   :  Draw method. 5 is the only option for .pe, at the moment.

This will generate a postscript file with the Circos-like diagram.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

3.2 Junction read analysis (midhr, circmaps)

---------------------------------------------------------------------------------------------

3.2.1 Detection of Junction reads, Formation of Junction read clusters and Rearrangement Prediction

Based on the mapping result of untruncated reads (1.0), 
'midhr' detects the junction reads/cluster and rearrangements suggested by them.

% midhr -cl 13 -cd 3 -job i14h75r10_0 /DBdir/CH_MOD.fasta PpWT.fastq

The detection of junction and clustering junction reads is somewhat complicated.
Please refer to the description in our paper, or contact us for further information.
You can control the sensitivity by 2 criteria.

-cl  # : Minimum length of the mismatch consensus sequence 
         (the longer the more sensitive)
-cd  # : Minimum number of reads to form mismatch consensus 
         (the more the more sensitive)
-job ~ : jobname with suffix (_0 for the first, _1 for the second reference)

The last two arguments are the file names of the reference sequence and the 
read sequence.

Output files with different suffix contains the different kind of mutation.

.del : Simple deletion less than 50 bps.
.ins : Simple insertion less than 50 bps.
.hr  : Homologous rearrangement with more than 2bp microhomology
.hmr : Homologous rearrangement (in both direction).
.pal : Homologous rearrangement caused by (quasi) palindrome less than 70bp. 
       (Hairpin turn rearrangement)
.unk : Rearrangement without microhomology
.arr : .hr + .unk + .del + .ins

Output files with following suffix contains the information about the status 
of junction clusters.

.jcr : Junction read information. Shows alignment with reference and mismatch 
       start point call.
.jcc : Junction cluster information. Shows alignment of mismatch consensus 
       sequences of all junction reads in cluster.
.jcs : List of junction clusters with mismatch starting position and search 
       results of the mismatch consensus sequence on reference.

---------------------------------------------------------------------------------------------

3.2.1.1 Rebuilding .hmr file

The .hmr file generated in the previous section 3.2.1 include the junction read that may include 
several mismatches. To improve the reliability of .hmr output, we can use program h2h to generate
.hmr file from .hr file which is more reliable. With the .hr file PpWT.fastq_i14h75r10_0.hr, type

% h2h PpWT.fastq_i14h75r10_0.hr

This will generate a .hmr file PpWT.fastq_i14h75r10_0.hr.hmr
Sorting the second column of this file in numerical order will provide an outpu which is easier to
recognize.

% sort -k 2 PpWT.fastq_i14h75r10_0.hr.hmr > PpWT_CH.hmr

---------------------------------------------------------------------------------------------

3.2.2 Drawing Circos-like map by junction read analysis

Using circmaps we can generate Circos-like diagram as similar fashion as 
paired-end analysis (3.1.3).
This time run it without -pe option, and specify 6 for -drm option, 
and put .hr, .arr, or .unk file as the last argument.

circmaps -rpb chp50_PpWT.fastq_i14h2r3_0.peb \
           -rgb /DBdir/CH_MOD.gb \
           -rfa /DBdir/CH_MOD.fasta \
           -bls 1.0 \
           -psc 2 \
           -drm 6 \
           PpWT.fastq_i14h75r10_0.arr

This will generate a postscript file of the diagram. PpWT.fastq_i14h75r10_0.arr.ps

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

4. Comparison of rearrangements/mutations (links) among experiments
   (comp_del, comp_ins, comp_hr, comp_hmr, comp_pal, comp_unk, comp_pe)

Some of the rearrangements are shared among experiments, 
and some are found only in some experiments.
In order to compare .del, .ins, .hr, .hmr, .pal, .unk, .pe files among experiments 
and organize the identical rearrangements/mutations, we prepared small programs:comp_del, 
comp_ins, comp_hr, comp_hmr, comp_pal, comp_unk, comp_pe, respectively.

Their usages are similar, and here we show how to use comp_del.
First we need to prepare a list of .del files as follows.

========================== del.list =============================
PpWT.fq_i11h35r36immap_0.del   40223458
Mut1.fq_i11h35r36immap_0.del   37033756
Mut2.fq_i11h35r36immap_0.del   42008745
=================================================================

The list consists of two column. The left column is the .del file name 
and the right column is the integer number
to calibrate the intensity among experiments (such as the number of total/mapped reads).
The first line is the standard and typically is the wild type experiment. 

% comp_del del.list

Program 'comp_del' takes the prepared list file name (del.list) as the only argument.
It generates two outputs with suffixes '.out' and '.out2'.
As for .ins, .del, .pal, .hmr, the output files are cvs (comma separated), 
and for .hr, .unk, .pe, the output files are fixed length space separated files.
Columns in right hand side are the number of junction reads indicating 
the rearrangement/mutaion (.out) for each experiment, 
and the numbers calibrated by the second column of the list file (.out2).
For instance, the example del.list file contains three lines, PpWT, Mut1 and Mut2, 
the 3 column from the right corresponds to the PpWT, Mut1, and Mut2, respectively.

Definitions of columns on the left side are as follows.

del : 1. Genome position, 2. Deletion length
ins : 1. Genome position, 2. Insertion length, 3. Inserted sequence
pal : 1. Genome position(left) 2. Genome position (right) 
      3. (quasi) Palindrome sequence
hmr : 1. Genome position(homologous sequence 1 left) 
      2. Genome position (homologous sequence 1 right) 
      3. Genome position(homologous sequence 2 left) 
      4. Genome position (homologous sequence 2 right) 
      5. Homologous sequence

unk : 1. Genome position (rearrange from) 
      2. Genome position (rearrange to)
hr  : 1. Genome position (rearrange from) 
      2. Genome position (rearrange to) 
      3. Homologous sequence
pe  : 1. Genome position (rearrange from) 
      2. Genome position (rearrange to) 
(For pe, the positions are the mapped positions of the seed pair of the cluster, 
indicating the same rearrangement.)

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

5. Depth plot (calibrated in comparison with WT experiment)
   (pop_comp)

Mapping depth plot of each experiment was obtained with 'ispmap' (2.1).
'ispmap' also creates .peb file that contains the number of read mapped 
at each position of reference sequence as a text file.
Here we describe the procedure to create a plot calibrated by standare 
(wild type) experiment.
In our case, we used the number of reads mapped on nuclear chromosome 
for each experiment to calculate the scale factor.
We assumed, gene knowckout mutant may have decreased/increased read count 
for each organella, while the number of reads mapped on nuclear chromosome 
should stay constant. Therefore we divide the number of reads mapped 
on nuclear chromosome of wild type, by the number of reads mapped on nuclear 
chromosome of knock out mutant, to obtain the scale factor (0.588 in our example).

% pop_comp -sf 0.588 WT.fastq_i14h2r3_0.peb MUT1.fastq_i14h2r3_0.peb
% pop_comp -sf 0.588 WT.fastq_i14h2r3_1.peb MUT1.fastq_i14h2r3_1.peb

The value after '-sf' option is the scale factor, and the last two arguments are 
.peb files of wild type and mutant.
This will create a postscript file 'MUT1.fastq_i14h2r3_0.peb.ps' that contains, 
depth plot of wild type, mutant, and the normalized plot.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

6. Histogram of number of repeats longer/shorter than a threshold
   (reptile)

Count the number of repeats longer than -min value and shorter than -max value
and draw a histgram along the genome sequence.

% reptile -ph -ps -idx 10 -min 11 -max 100 MT.fasta > MT_11_100.out

The program generates output file with suffix .pos, .dst, .hmp.
Here we use MT.fasta.pos to draw a histgram using R.

Launch R and load science library.
> library(MASS)
Load .pos file.
> MT_11_100 = scan("/Users/.../MT.fasta.pos")
Specify output ps file name.
> postscript("MT_11_100.ps")
Draw histogram
> truehist(MT_11_100,prob=FSLSE,nbins = 50)
Create .ps file
> dev.off()

'MT.fasta.ps' is generated.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

7. Plot histogram of rearrangement position along genome, based on the midhr output (.hr)
   (hr_count)

Based on the junction analysis of 'midhr' (.hr output), 'hr_count' creates a histgram of 
the number of homologous rearrangement along the genome sequence.

% hr_count Mut1.fastq_i14h75r10_0.hr

This will generate three files Mut1.fastq_i14h75r10_0.hr.3_5.pos, 
Mut1.fastq_i14h75r10_0.hr.6_10.pos and Mut1.fastq_i14h75r10_0.hr.11_100.pos
The output file Mut1.fastq_i14h75r10_0.hr.3_5.pos contains the positions of 
rearrangement with which the length of homologous sequnce is from 3 to 5 bps.
The other output files contains positions of rearrangements with the length 
of homologous sequences, from 6 to 10, and 11 to 100.

To draw the histogram from one of these outputs 'Mut1.fastq_i14h75r10_0.hr.3_5.pos'
Launch R and load science library.
> library(MASS)
Load .pos file.
> Mut1_3_5 = scan("/Users/.../Mut1.fastq_i14h75r10_0.hr.3_5.pos")
Specify output ps file name.
> postscript("Mut1_3_5.ps")
Draw histogram
> truehist(Mut1_3_5,prob=FALSE,nbins = 50)
Create .ps file
> dev.off()

'Mut1_3_5.ps' is generated.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

8. Reorganize .hmr from .hr output of midhr and obtain link position file .psd and number of 
   reads consiting the links .psd2
   (h2h)
   
.hmr file generated from midhr may include junction reads excluded in the formation of regular
junction cluster in midhr. Here we reconstruct .hmr file only from junction clusters in .hr file
using program h2h.

% h2h PpWT.fastq_i14h75r10_1.hr

This generates a file with extra .hmr suffix, as PpWT.fastq_i14h75r10_1.hr.hmr
This also generates .psd file consisting of three columns position of link, the number of 
links toward upstream, and toward downstream at the position (PpWT.fastq_i14h75r10_1.hr.psd).
h2h also generates .psd2 file with which the total number of junction reads consisting each 
each at the position (PpWT.fastq_i14h75r10_1.hr.psd2).
 
---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

9. Using .hmr output from the previous section 8, generate a histgram of the number of links.
   (h2h_histo)

First sort the .psd file obtained in the previous section 8, in the order of genome position.

% sort -k 1 -n PpWT.fastq_i14h75r10_1.hr.psd > PpWT_MT.psd

Then run h2h_histo with two additional arguments, the window size (1000 in the following 
example) and the length of the reference genome (105,340 in the example).

% h2h_histo 1000 105340 PpWT_MT.psd  >  PpWT_MT.hst

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

10. Using .hmr output from the previous section 8, calculate Morisita index regarding the
    distribution of the links.
    (h2h_morisita)

Calculate the Morisita index by using h2h_morisita with same arguments as the previous command
h2h_histo.

% h2h_morisita 1000 105340 PpWT_MT.psd

This outputs two numbers id1 and id2, to standard output.
id1 is the Morisita index for the links toward upstream, and id2 is the index for the links
toward downstream.

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
