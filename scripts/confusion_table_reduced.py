#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Reproduced from:
# Mah� F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8



"""
Parse taxonomic assignments and swarm results to produce a confusion
table.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2013/09/21"
__version__ = "$Revision: 0.1"

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
from optparse import OptionParser


def option_parser():
    """
    Parse arguments from command line.
    """
    desc = """This program reads a table of taxonomic assignments and
    a swarm clustering file and outputs a confusion table (OTUs vs
    taxonomic assignments)."""

    parser = OptionParser(usage="usage: %prog -t FILENAME -s FILENAME",
                          description=desc,
                          version="%prog version 0.1")

    parser.add_option("-t", "--taxonomic_assignments",
                      metavar="<FILENAME>",
                      action="store",
                      dest="taxonomic_assignments",
                      help="set <FILENAME> as input.")

    parser.add_option("-s", "--swarm_OTUs",
                      metavar="<FILENAME>",
                      action="store",
                      dest="swarm_OTUs",
                      help="set <FILENAME> as input.")

    (options, args) = parser.parse_args()
    return options.taxonomic_assignments, options.swarm_OTUs


def parse_taxonomy(taxonomic_assignments):
    """
    Parse taxonomy assignments
    """
    amplicon2taxonomy = dict()
    taxa = []
    with open(taxonomic_assignments, "rU") as taxonomic_assignments:
        for line in taxonomic_assignments:
            amplicon, taxon = line.split()[0:2]
            amplicon, abundance = amplicon.split("_")
            amplicon2taxonomy[amplicon] = (abundance, taxon)
            taxa.append(taxon)
    return amplicon2taxonomy, taxa


def dereplicate_taxa(taxa):
    """
    Dereplicate taxa (list of unique taxon names, and initialized dictionary)
    """
    taxa_list = list(set(taxa))
    taxa_dict = dict(zip(taxa_list, [0] * len(taxa_list)))
    return taxa_list, taxa_dict


def OTU_parser(taxa_dict, taxa_list, amplicon2taxonomy, amplicons):
    """
    Parse each OTU and output one line of the confusion table.
    """
    for amplicon in amplicons:
        try:
            amplicon, abundance = amplicon.split("_")
        except ValueError:
            amplicon, abundance = amplicon, 1
        try:
            abundance, taxon = amplicon2taxonomy[amplicon]
        except KeyError:
            taxon = "Unassigned"
        taxa_dict[taxon] += int(abundance)
    OTU_abundance_per_taxa = [str(taxa_dict[taxon]) for taxon in taxa_list]
    return OTU_abundance_per_taxa


if __name__ == '__main__':

    ## Parse command line arguments
    taxonomic_assignments, swarm_OTUs = option_parser()

    ## Parse taxonomy assignments
    amplicon2taxonomy, taxa = parse_taxonomy(taxonomic_assignments)

    ## Dereplicate taxa
    taxa_list, taxa_dict = dereplicate_taxa(taxa)

    ## Output the table header
    print("OTUs_vs_Taxa", "\t".join(taxa_list), sep="\t", file=sys.stdout)

    ## Parse swarm OTUs
    with open(swarm_OTUs, "rU") as swarm_OTUs:
        for i, line in enumerate(swarm_OTUs):
            amplicons = line.strip().split()
            OTU_abundance_per_taxa = OTU_parser(taxa_dict.copy(), taxa_list,
                                                amplicon2taxonomy, amplicons)
            print(str(i+1), "\t".join(OTU_abundance_per_taxa), sep="\t",
                  file=sys.stdout)

sys.exit(0)
