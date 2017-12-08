#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Reproduced from:
# Mah� F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8



"""
    Detect and break chains of amplicons in a swarm.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/01/30"
__version__ = "$Revision: 1.1"

import os
import sys
import tempfile
import itertools
import subprocess
from operator import itemgetter
from optparse import OptionParser

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Detect and break chains of amplicons in a swarm. That
    script will search for the swarm binary in /usr/bin/. If swarm is
    installed at a different location, please modify the corresponding
    line in the function run_swarm."""

    parser = OptionParser(usage="usage: %prog -f filename -s filename",
                          description=desc,
                          version="%prog version 1.1")

    parser.add_option("-b", "--binary",
                      metavar="<BINARY>",
                      action="store",
                      default="/usr/bin/swarm",
                      dest="binary",
                      help="swarm binary location. Default is /usr/bin/swarm")

    parser.add_option("-f", "--fasta_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="fasta_file",
                      help="set <FILENAME> as fasta file.")

    parser.add_option("-s", "--swarm_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="swarm_file",
                      help="set <FILENAME> as swarm file.")

    parser.add_option("-d", "--differences",
                      metavar="<THRESHOLD>",
                      action="store",
                      type="int",
                      default=1,
                      dest="threshold",
                      help="set local clustering <THRESHOLD>. Default is 1")

    (options, args) = parser.parse_args()
    return options.binary, options.fasta_file, options.swarm_file, options.threshold


def fasta_parse(fasta_file):
    """
    List amplicon ids, abundances and sequences, make a list and a dictionary
    """
    with open(fasta_file, "rU") as fasta_file:
        all_amplicons = dict()
        for line in fasta_file:
            if line.startswith(">"):
                amplicon, abundance = line.strip(">\n").split("_")
            else:
                sequence = line.strip()
                all_amplicons[amplicon] = (int(abundance), sequence)
        return all_amplicons


def swarm_parse(swarm_file):
    """
    List amplicons contained in each swarms, sort by decreasing
    abundance. Sort the list of swarms by decreasing mass and
    decreasing size.
    """
    with open(swarm_file, "rU") as swarm_file:
        swarms = list()
        for line in swarm_file:
            amplicons = [(amplicon.split("_")[0], int(amplicon.split("_")[1]))
                         for amplicon in line.strip().split(" ")]
            # Sort amplicons by decreasing abundance and alphabetical order
            amplicons.sort(key=itemgetter(1, 0), reverse=True)
            top_amplicon, top_abundance = amplicons[0]
            swarm_size = len(amplicons)
            swarm_mass = sum([amplicon[1] for amplicon in amplicons])
            swarms.append([top_amplicon, swarm_mass, swarm_size,
                           top_abundance, amplicons])
        # Sort swarms on mass, size and seed name
        swarms.sort(key=itemgetter(1, 2, 0), reverse=True)
        return swarms


def run_swarm(binary, all_amplicons, swarm, threshold):
    """
    Write temporary fasta files, run swarm and collect the graph data
    """
    swarm_command = [binary, "-b", "-d", str(threshold)]
    with open(os.devnull, "w") as devnull:
        with tempfile.SpooledTemporaryFile() as tmp_fasta_file:
            with tempfile.SpooledTemporaryFile() as tmp_swarm_results:
                for amplicon, abundance in swarm:
                    sequence = all_amplicons[amplicon][1]
                    print(">", amplicon, "_", str(abundance), "\n", sequence,
                          sep="", file=tmp_fasta_file)
                tmp_fasta_file.seek(0)  # rewind to the begining of the file
                proc = subprocess.Popen(swarm_command,
                                        stderr=tmp_swarm_results,
                                        stdout=devnull,
                                        stdin=tmp_fasta_file,
                                        close_fds=True)
                proc.wait()  # usefull or not?
                tmp_swarm_results.seek(0)  # rewind to the begining of the file
                graph_data = [line.strip().split("\t")[1:4]
                              for line in tmp_swarm_results
                              if line.startswith("@")]
                return graph_data


def build_graph(graph_data):
    """
    List pairwise relations in a swarm. Note that not all pairwise
    relations are stored. That's why the graph exploration must always
    start from the most abundant amplicon, and must be reiterated for
    sub-swarms after a breaking.
    """
    graph = dict()
    for line in graph_data:
        ampliconA, ampliconB, differences = line
        if ampliconA in graph:
            graph[ampliconA] += [ampliconB]
        else:
            graph[ampliconA] = [ampliconB]
    return graph


def find_path(graph, start, end, path=[]):
    """
    Recursively explore the graph and find all paths connecting two
    amplicons (http://www.python.org/doc/essays/graphs.html). As the
    graph is not complete, some pairs of amplicon cannot be linked.
    """
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath:
                return newpath
    return None


def graph_breaker(amplicons, graph, all_amplicons, ABUNDANT):
    """
    Find deep valleys and cut the graph
    """
    # High peaks to test (starting and ending points)
    top_amplicons = [amplicon[0] for amplicon in amplicons
                     if amplicon[1] >= ABUNDANT]
    # Ending peak is RATIO times higher than the valley
    RATIO = 50
    # Debugging
    print("## OTU ", top_amplicons[0], "\n", "# List potential bridges",
          sep="", file=sys.stderr)
    # Initialize the list of new seeds
    new_swarm_seeds = [top_amplicons[0]]
    # Break if there is no second peak
    if len(top_amplicons) < 2:
        return new_swarm_seeds, graph
    # Loop over the list of top amplicons
    pairs_of_peaks = itertools.combinations(top_amplicons, 2)
    for pair_of_peaks in pairs_of_peaks:
        start_amplicon, end_amplicon = pair_of_peaks
        path = find_path(graph, start_amplicon, end_amplicon)
        # Path can be empty if the relation have been deleted
        if path and len(path) > 1:
            abundances = [int(all_amplicons[node][0]) for node in path]
            # Find the weakest spot
            lowest = min(abundances)
            if lowest != abundances[-1]:
                # LOW VALLEY MODEL (CHANGE HERE)
                if (abundances[-1] / lowest > RATIO / 2 and abundances[0] / abundances[-1] < 10) or abundances[-1] / lowest >= RATIO:
                    # Debugging
                    print(abundances, "\tBREAK!", file=sys.stderr)
                    # Find the rightmost occurence of the lowest point
                    index = len(abundances) - (abundances[::-1].index(lowest) + 1)
                    left_amplicon = path[index-1]
                    right_amplicon = path[index]
                    # Delete the relation from the graph
                    graph[left_amplicon].remove(right_amplicon)
                    # Remove the graph entry if the relation is now empty
                    if not graph[left_amplicon]:
                        del graph[left_amplicon]
                    # Lowest point will be a new swarm seed
                    new_swarm_seeds.append(right_amplicon)
                else:
                    print(abundances, file=sys.stderr)
    return new_swarm_seeds, graph


def swarmer(graph, seed, path=[]):
    """
    Recursively explore the graph and find all amplicons linked to the
    seed
    """
    path = path + [seed]
    if seed in graph:
        for node in graph[seed]:
            path = swarmer(graph, node, path)
    return path


def swarm_breaker(binary, all_amplicons, swarms, threshold):
    """
    Recursively inspect and break the newly produced swarms
    """
    # ARBITRARY PARAMETERS
    ABUNDANT = 100
    # Deal with each swarm
    for swarm in swarms:
        top_amplicon, swarm_mass, swarm_size, top_abundance, amplicons = swarm
        if swarm_size > 2 and top_abundance > ABUNDANT:
            # Run swarm to get the pairwise relationships
            graph_raw_data = run_swarm(binary, all_amplicons, amplicons, threshold)
            # Build the graph of pairwise relationships
            graph = build_graph(graph_raw_data)
            new_swarm_seeds, graph = graph_breaker(amplicons, graph,
                                                   all_amplicons, ABUNDANT)
            # Explore the graph and find all amplicons linked to the seeds
            observed = 0
            new_swarms = list()
            for seed in new_swarm_seeds:
                new_swarm = swarmer(graph, seed)
                observed += len(new_swarm)
                # Give to the new swarms the same structure and
                # re-order them by decreasing abundance
                amplicons = [(amplicon, all_amplicons[amplicon][0])
                             for amplicon in new_swarm]
                amplicons.sort(key=itemgetter(1), reverse=True)
                top_amplicon, top_abundance = amplicons[0]
                swarm_size = len(amplicons)
                swarm_mass = sum([amplicon[1] for amplicon in amplicons])
                new_swarms.append([top_amplicon, swarm_mass,
                                   swarm_size, top_abundance, amplicons])
            # Deal with the new swarms (no need to treat again the
            # first swarm). There will always be at least one swarm in
            # new_swarms.
            print(" ".join(["_".join([amplicon[0], str(amplicon[1])])
                            for amplicon in new_swarms[0][4]]), file=sys.stdout)
            new_swarms.pop(0)
            if new_swarms:
                # Sort the rest of the new swarms by decreasing mass
                # and size. Inject them into swarm_breaker.
                new_swarms.sort(key=itemgetter(1, 2), reverse=True)
                swarm_breaker(binary, all_amplicons, new_swarms, threshold)
        else:
            # Output the swarm
            print(" ".join(["_".join([amplicon[0], str(amplicon[1])])
                            for amplicon in amplicons]), file=sys.stdout)
    return None


def main():
    """
    Hypothesis: chain of amplicons happen among the most abundant
    amplicons of the swarm. The number of chains in a swarm is
    small. The abundances of each sub-swarm centroids are
    comparable. The "valleys" are deep compared to the "peaks". Swarm
    graphs are acyclical, so there is only one path joining two
    amplicons.

    Synopsis: Break bridges as you discover them. Find the weakest
    point in the chain, break on the left of that point and mark it as
    the seed of a new swarm. Repeat the process with the nth most
    abundant amplicon, until all amplicons in the arbitrary range have
    been treated.
    """
    # Parse command line options.
    binary, fasta_file, swarm_file, threshold = option_parse()
    # Load all amplicon ids, abundances and sequences
    all_amplicons = fasta_parse(fasta_file)
    # Load the swarming data
    swarms = swarm_parse(swarm_file)
    # Deal with each swarm
    swarm_breaker(binary, all_amplicons, swarms, threshold)


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
