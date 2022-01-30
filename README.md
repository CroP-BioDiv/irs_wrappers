# IRs wrappers
Wrapper scripts for chloroplast inverted repeats annotation tools.

Scripts have same interface:
- input is filename with chloroplast sequence in fasta or GenBank format,
- outputs are four numbers of IRs locations (IRa start, IRa end, IRb start, IRb end).

Each wrapper script is described with tool specifics and requirenmnets.

Biopython


## Chloë

Web: https://chloe.plantenergy.edu.au/ 
Source code: https://githubSource date.com/ian-small/chloe


## Chloroplot

The R package Chloroplot wrapped the functions for visualizing the organelle genomes.

**Source code:** https://github.com/shuyuzheng/Chloroplot

**Web page:** https://irscope.shinyapps.io/chloroplot/

**Wrapper script:** [chloroplot.py](chloroplot.py)

**Requirements:**
- R

**Notes:**


## GeSeq

GeSeq is a Web application used for annotation of organelle genomes, in particular chloroplast genomes.

**Source code:** -

**Web page:** https://chlorobox.mpimp-golm.mpg.de/geseq.html

**Wrapper script:** [ge_seq.py](ge_seq.py)

**Requirements:** -

**Notes:**
- sequence has to be processed with Web application and wrapper extracts IRs location from result GenBank file.



## ORG.Annotate

A pipeline for annotating Chloroplast genomes.

**Source code:** https://git.metabarcoding.org/org-asm/org-annotate 

**Web page:** -

**Wrapper script:** [org_annotate.py](org_annotate.py)

**Requirements:**
- repseek executable has to be on the PATH or location has to be set with environment variable REPSEEK_EXE.
- blastn has to be on the PATH.

**Notes:** 


## PGA

Plastid Genome Annotator

**Source code:** https://github.com/quxiaojian/PGA

**Web page:** -

**Wrapper script:** [pga.py](pga.py)

**Requirements:**
- PGA.pl script has to be on the PATH.
- Check PGA's [requirenments](https://github.com/quxiaojian/PGA) (Perl, blastn)

**Notes:** 


## Plann

Plann is a command-line application for annotating chloroplast sequences.

**Source code:** https://github.com/daisieh/plann

Our fork: https://github.com/CroP-BioDiv/plann

**Web page:** -

**Wrapper script:** [plann.py](plann.py)

**Requirements:**
- plann.pl script has to be on the PATH or location has to be set with environment variable PLANN_SCRIPT.
= Check plann's [requirenments](https://github.com/daisieh/plann) (Perl, blastn)

**Notes:**
- Plann leaks blastn result files in temp folder. If run on lot of files it can fill partition.
Files have name of length 10. Delete them with:
```
cd /tmp
find . -name '??????????' -delete
```


## Airpg

Automatically accessing the inverted repeats of archived plastid genomes

**Source code:** https://github.com/michaelgruenstaeudl/airpg

**Web page:** -

**Wrapper script:** [airpg.py](airpg.py)

**Requirements:** -
