# GeFaST: An improved method for OTU assignment by generalising Swarm's fastidious clustering approach (submitted) 

Scripts for the evaluation of GeFaST's extended fastidious-clustering capabilities.


## How to use it

To set up the necessary programs (except one, see the comment in the makefile for more information), run `make` in the top-level directory of the repository.


The analyses are launched collectively by `./run.sh`.
This will download and prepare the data sets, perform and evaluate all the clusterings and visualise the results.

To recreate the visualisations from existing analysis results, run `./revisualise.sh`.


## Required software
 * [GeFaST](https://github.com/romueller/gefast) (version 1.0.0)
 * [Swarm](https://github.com/torognes/swarm) (version 1.2.3 and version 2.1.13)
 * [USEARCH](http://www.drive5.com/usearch/download.html) (version 7)
 * [seqtk](https://github.com/lh3/seqtk) (version 1.2-r95-dirty)
 * GCC (version 4.9.2 or higher)
 * python (version 2.7 or higher; with Biopython)
 * perl (version 5.20.2 or higher)
 * make (version 4.0 or higher)
 * R / Rscript (version 3.3.2 or higher; with packages dplyr, ggplot2, grid, gridExtra and reshape2)
 * /usr/bin/time (external command)
 * common Unix tools (bash, bzip2, cat, wget)

_Notes:_
Older versions of GCC, python, perl etc. might also work, but have not been tested. 
Makefile takes care of GeFaST and Swarm, while the other software prerequisites have to be satisfied by the user.  
The location of the binaries of software such as seqtk has to be in `PATH`.
