#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Replaces identifiers (obtained from FASTA headers) in a taxonomic assignment
# with the taxonomy corresponding to the identifier.



import sys
from optparse import OptionParser

# Parse arguments from command line.
def option_parser():
    des = """This program reads (1) a table of taxonomic assignments 
    (e.g. produced with USEARCH or VSEARCH) and (2) a taxonomy
    file (e.g. from GreenGenes or Silva). It replaces the identifiers of 
    the representative sequences in (1) with the taxonomy information from (2)."""

    parser = OptionParser(usage = "Usage: %prog -t <file> -r <file>",
                          description = des)

    parser.add_option("-t", "--taxonomic_assignments",
                      action = "store", dest = "taxonomic_assignments")

    parser.add_option("-r", "--rep_taxonomy",
                      action = "store", dest = "rep_taxonomy")

    (options, args) = parser.parse_args()
    return options.taxonomic_assignments, options.rep_taxonomy


# Parse taxonomy file (e.g. gg_13_8_otu/taxonomy/97_otu_taxonomy.txt)
def parse_otu_taxonomy(otu_taxonomy):
    id2taxonomy = dict()
    with open(otu_taxonomy, "rU") as otu_taxonomy:
        for line in otu_taxonomy:
            id, taxon = line.split("\t")
            taxon = "".join(taxon.split(" "))
            id2taxonomy[id] = taxon.rstrip()
    return id2taxonomy


# Parse taxonomic assignments
def parse_assignments(taxonomic_assignments):
    assignments = []
    with open(taxonomic_assignments, "rU") as taxonomic_assignments:
        for line in taxonomic_assignments:
            assignments.append(line.split("\t"))
    return assignments


# Look up the taxonomy and replace the identifier in the assignment with the taxonomy.
def replace_identifiers(assignments, id2taxonomy):
    for a in assignments:
        a[1] = id2taxonomy[a[1]]
    return assignments


if __name__ == '__main__':

    ## Parse command-line arguments
    taxonomic_assignments, rep_taxonomy = option_parser()

    ## Parse taxonomy file and assignments
    id2taxonomy = parse_otu_taxonomy(rep_taxonomy)
    assignments = parse_assignments(taxonomic_assignments)

    ## Replace the identifiers with the taxonomy
    assignments = replace_identifiers(assignments, id2taxonomy)

    ## Print relabelled assignments
    for a in assignments:
        sys.stdout.write("\t".join(a))

sys.exit(0)
