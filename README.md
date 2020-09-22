# Evaluation of extended fastidious clustering implemented in GeFaST

Scripts for the evaluation of GeFaST's extended fastidious-clustering capabilities as presented in [*GeFaST: An improved method for OTU assignment by generalising Swarm's fastidious clustering approach*](https://doi.org/10.1186/s12859-018-2349-1).


## How to use it

To set up the necessary programs (except one, see the comment in the makefile for more information), run `make` in the top-level directory of the repository.


The analyses are launched collectively by `./run.sh`.
This will download and prepare the data sets, perform and evaluate all the clusterings and visualise the results.

To recreate the visualisations from existing analysis results, run `./revisualise.sh`.


## Required software
 * [GeFaST](https://github.com/romueller/gefast) (version 1.0.0)
 * [Swarm](https://github.com/torognes/swarm) (version 1.2.3 and version 2.1.13)
 * [USEARCH](http://www.drive5.com/usearch/download.html) (version 10)
 * [VSEARCH](https://github.com/torognes/vsearch) (version 2.7.1)
 * [CD-HIT](https://github.com/weizhongli/cdhit) (version 4.6.8)
 * [DNACLUST](http://dnaclust.sourceforge.net/) (release 3)
 * [Sumaclust](https://git.metabarcoding.org/obitools/sumaclust/wikis/home) (version 1.0.31)
 * [seqtk](https://github.com/lh3/seqtk) (version 1.2-r95-dirty)
 * GCC (version 4.9.2 or higher)
 * python (version 2.7 or higher; with Biopython)
 * perl (version 5.20.2 or higher)
 * make (version 4.0 or higher)
 * R / Rscript (version 3.3.2 or higher; with packages dplyr, ggplot2, grid, gridExtra, RColorBrewer and reshape2)
 * /usr/bin/time (external command)
 * common Unix tools (awk, bash, bzip2, cat, cut, wget, sed, shuf)

_Notes:_
The scripts in this repository are not designed for more recent versions of GeFaST (2.0.0 or higher).   
Older versions of GCC, python, perl etc. might also work but have not been tested. 
Makefile takes care of GeFaST, Swarm, VSEARCH, CD-HIT, DNACLUST and Sumaclust, 
while the other software prerequisites have to be satisfied by the user.  
The location of the binaries of software such as seqtk has to be in `PATH`.
