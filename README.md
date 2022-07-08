# IRs wrappers

Wrapper scripts for chloroplast inverted repeats annotation tools. Scripts are implemented in python3 and use [Biopython](https://biopython.org/) package.

List of annotation tools:

| Tool | Script | Source code | Web application |
| ----------- | ----------- | ----------- | ----------- |
| Chloë | [chloe.py](chloe.py) | https://github.com/ian-small/chloe | https://chloe.plantenergy.edu.au/ |
| Chloroplot | [chloroplot.py](chloroplot.py) | https://github.com/shuyuzheng/Chloroplot | https://irscope.shinyapps.io/chloroplot/ |
| GeSeq | [ge_seq.py](ge_seq.py) | - | https://chlorobox.mpimp-golm.mpg.de/geseq.html |
| ORG.Annotate | [org_annotate.py](org_annotate.py) | https://git.metabarcoding.org/org-asm/org-annotate | - |
| PGA | [pga.py](pga.py) | https://github.com/quxiaojian/PGA | - |
| Plann | [plann.py](plann.py) | https://github.com/daisieh/plann | - |
| Airpg | [airpg.py](airpg.py) | https://github.com/michaelgruenstaeudl/airpg | - |


## Usage

Scripts have same interface:
- input is a filename with chloroplast sequence in fasta or GenBank format,
- outputs are four numbers of IRs locations (IRa start, IRa end, IRb start, IRb end) or none if IRs were not identified.

Helper script [run_more.py](run_more.py) is used to run wrapper script on more sequences. Usage:
```
python3 run_more.py <script_name> filename [filename]+
```


## Implementation specifics

### Chloë

Organelle Annotator.

**Notes:**
- Script does not run Julia program, instead it fetches annotation through Chloë web API which is of type
https://chloe.plantenergy.edu.au/annotate-ncbi?ncid={accession_number}&force_circular=true


### Chloroplot

The R package Chloroplot wrapped the functions for visualizing the organelle genomes.

**Requirements:**
- R programming language (Rscript)


### GeSeq

GeSeq is a Web application used for annotation of organelle genomes, in particular chloroplast genomes.

**Notes:**
- sequence has to be processed with Web application and wrapper extracts IRs location from result GenBank file.


### ORG.Annotate

A pipeline for annotating Chloroplast genomes.

**Requirements:**
- repseek executable has to be on the PATH or location has to be set with environment variable REPSEEK_EXE.
- blastn has to be on the PATH.


### PGA

Plastid Genome Annotator

**Requirements:**
- PGA.pl script has to be on the PATH.
- Check PGA's [requirenments](https://github.com/quxiaojian/PGA) (Perl, blastn)


### Plann

Plann is a command-line application for annotating chloroplast sequences.

Our fork: https://github.com/CroP-BioDiv/plann

**Requirements:**
- plann.pl script has to be on the PATH or location has to be set with environment variable PLANN_SCRIPT.
= Check plann's [requirenments](https://github.com/daisieh/plann) (Perl, blastn)

**Notes:**
- Plann leaks blastn result files in temp folder. If run on lot of files it can fill partition.
Files have name of length 10. Delete them with:
```
find /tmp -name '??????????' -delete
```

### Airpg

Automatically accessing the inverted repeats of archived plastid genomes



## Research

Scripts are result of research paper "Chloroplast genome annotation tools: Challenges and recommendations".
Similar scripts, implemented in [ZCItools environment](https://github.com/CroP-BioDiv/zcitools), are used to obtain results for the research.
Reproduction is described on this [page](https://github.com/CroP-BioDiv/zcitools/blob/master/docs/irs_annotation_tools.md).
