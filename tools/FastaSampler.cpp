/*
 * FastaSampler
 *
 * Copyright (C) 2017 Robert Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <vector>


typedef std::pair<std::string, std::string> Entry; // used as (header, sequence)

char convert[128] = { // lower-case for A-Z, everything else remains the same
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 91,  92,  93,  94,  95,
    96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
};


void lowerCase(std::string& s) {

    for (auto i = 0; i < s.size(); i++) {
        s[i] = convert[s[i]];
    }

}

// reads all entries from a FASTA file as header-sequence pairs
void readInput(const std::string fileName, std::vector<Entry>& entries) {

    std::ifstream iStream(fileName);
    if (!iStream.good()) {

        std::cerr << "ERROR: File '" << fileName << "' not opened correctly. No sequences are read from it." << std::endl;
        return;

    }

    std::string line, dl, seq;
    bool first = true;

    while (std::getline(iStream, line).good()) {

        if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

        if (line[0] == '>') { // header line

            if (first) { // first entry found (no previous entry to finish), simply get header

                dl = line;
                first = false;

            } else { // finish and store previous entry, then get header of new entry

                lowerCase(seq);
                entries.push_back(std::make_pair(dl, seq));


                seq.clear();
                dl = line;

            }

        } else { // still the same entry, continue to collect sequence
            seq += line;
        }
    }

    if (!first) { // ensures that last entry (if any) is written to file.

        lowerCase(seq);
        entries.push_back(std::make_pair(dl, seq));

    }

}

// randomly chooses a given number of samples
std::vector<Entry> sample(std::vector<Entry>& entries, size_t n) {

    // create indices 0, ..., |entries| - 1
    std::vector<size_t> indices(entries.size());
    std::iota(indices.begin(), indices.end(), 0);

    // shuffle indices randomly
    std::random_device rng;
    std::mt19937 urng(rng());
    std::shuffle(indices.begin(), indices.end(), urng);

    // choose those entries corresponding to the first n indices
    std::vector<Entry> s;
    for (auto i = 0; i < n; i++) {
        s.push_back(entries[indices[i]]);
    }

    return s;

};

// writes the given entries as a FASTA file to disk
void writeSample(const std::string fileName, std::vector<Entry>& sample) {

    std::ofstream oStream(fileName);

    for (auto& s : sample) {
        oStream << s.first << std::endl << s.second << std::endl;
    }

    oStream.close();

}

void printHelp() {

    std::cout << "Expected inputs: <FASTA file> <percentages-csv> <repetitions> <output-stem>" << std::endl;
    std::cout << "\t<FASTA file>: input file to sample from" << std::endl;
    std::cout << "\t<percentages>: comma-separated list of percentages (integer values, e.g. 50)" << std::endl;
    std::cout << "\t<repetitions>: number of samples to obtain per percentage" << std::endl;
    std::cout << "\t<output-stem>: path and prefix of output files (e.g. /home/user/sample), " << std::endl;
    std::cout << "\t\t completed by the percentage and repetition number (e.g. /home/user/sample_50_0.fasta)" << std::endl;

}


/**
 * Reads the entries from a FASTA file and samples it according to
 * the given percentages and repetitions.
 * Comments in the input FASTA file are not transferred to the output FASTA files.
 * Sequences in the output FASTA files are one-liners (regardless of the input FASTA file).
 *
 * Expected inputs: <FASTA file> <percentages-csv> <repetitions> <output-stem>
 */
int main(int argc, const char* argv[]) {

    if (argc < 5) {

        std::cout << "Error: Not enough arguments!" << std::endl;
        printHelp();
        return 1;

    }

    std::string inputFile = argv[1];
    std::string percentagesStr = argv[2]; // comma-separated percentages
    int reps = std::stoi(argv[3]); // number of samples per percentage
    std::string outputStem = std::string(argv[4]);

    std::vector<int> percentages;
    size_t start = 0;
    size_t end = percentagesStr.find(',');
    while (end != std::string::npos) {

        percentages.push_back(std::stoi(percentagesStr.substr(start, end - start)));
        start = end + 1;
        end = percentagesStr.find(',', start);

    }

    percentages.push_back(std::stoi(percentagesStr.substr(start, end)));


    std::cout << "FASTA file: " << inputFile << std::endl;
    std::cout << "Percentages: ";
    for (auto p : percentages) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
    std::cout << "Repetitions: " << reps << std::endl;
    std::cout << "Stem of output files: " << outputStem << std::endl << std::endl;

    std::cout << "Reading FASTA file..." << std::flush;
    std::vector<Entry> entries;
    readInput(inputFile, entries);
    std::cout << "DONE" << std::endl;

    std::cout << "Sampling..." << std::endl;
    for (auto p : percentages) {

        std::cout << p << " %: 0 / " << reps << " completed\r" << std::flush;
        for (auto i = 0; i < reps; i++) {

            std::vector<Entry> s = sample(entries, size_t(ceil((p / 100.0) * entries.size())));
            writeSample(outputStem + "_" + std::to_string(p) + "_" + std::to_string(i) + ".fasta", s);

            std::cout << p << " %: "<< (i + 1) << " / " << reps << " completed\r" << std::flush;

        }

        std::cout << std::endl;

    }

    std::cout << "\nAll samples obtained!" << std::endl;


}