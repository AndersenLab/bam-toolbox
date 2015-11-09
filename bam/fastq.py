#! /usr/bin/env python
"""
usage:
  bam fastq [options] <fq>

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --header                    Print header
  --tsv                       Print output in tsv format
"""
from docopt import docopt
from clint.textui import colored, indent, puts_err
import os
from collections import OrderedDict
from subprocess import Popen, PIPE
import gzip
import re
from pprint import pprint as pp
from eav import *
from itertools import groupby as g

old_illumina_header = ["instrument",
                       "flowcell_lane",
                       "flowcell_number",
                       None,  # x-tile
                       None,  # y-tile
                       "barcode",  # barcode - fetched later
                       "pair"]

illumina_header = ["instrument",
                   "run_id",
                   "flowcell_id",
                   "flowcell_lane",
                   None,  # tile number
                   None,  # x-tile
                   None,  # y-tile
                   "pair",
                   "filtered",
                   "control_bits",
                   "barcode"]  # barcode/index sequence; fetched later.

SRR_header = ['SRR']

stat_header = ["Total_Reads",
         "Unique_Reads",
         "Percent_Unique", 
         "Most_Abundant_Sequence",
         "Most_Abundant_Frequency",
         "Percentage_Unique_fq"]

def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def most_common(L):
    # Fetch most common item from a list.
  try:
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
  except:
    return ""

class fastq:
    # Simple class for reading fastq files.
    def __init__(self, filename):
        self.filename = filename
        # Get fastq information
        header_lines = [x["info"] for x in self.read(1)]
        header = re.split(r'(\:|#|/| )',header_lines[0])[::2]
        fetch_barcode = True

        if len(header) == 11:
            # Use new header format.
            use_header = illumina_header
        elif len(header) == 6:
            # Use old header format.
            use_header = old_illumina_header
        elif len(header) == 7:
            # Use old header and add pair if available.
            use_header = old_illumina_header + ["pair"]
        elif header[0].startswith("@SRR"):
            # Setup SRR Header
            header = header[0].split(".")[0:1]
            use_header = SRR_header
            fetch_barcode = False
        else:
            # If unknown header, enumerate
            use_header = ["h" + str(x) for x in range(0,len(header))]
            fetch_barcode = False

        if fetch_barcode == True:
            # Fetch index
            index_loc = use_header.index("barcode")
            fetch_index = [re.split(r'(\:|#|/| )',x["info"])[::2][index_loc] for x in self.read(1000)]
            self.barcode = most_common(fetch_index)


        self.header = {}
        # Set remaining attributes.
        for attr, val in zip(use_header, header):
            if attr is not None:
                val = autoconvert(val) # Set variable type
                self.header[attr] = val
                setattr(self, attr, val)


    def read(self, n=-1):
        """
            Iterate through gzipped fastq file and put
            yield sequence+info in dictionary.
        """
        if self.filename.endswith(".gz"):
            open_file = gzip.open(self.filename, 'rb')
        else:
            open_file = open(self.filename, 'r') 
        with open_file as f:
            dna = {}
            for linenum, line in enumerate(f):
                dna["info"] = line.strip()
                dna["seq"] = f.next().strip()
                f.next()
                dna["qual"] = f.next().strip()
                if linenum < n or n == -1:
                    yield dna
                else:
                    break

    def calculate_fastq_stats(self):
        if self.filename.endswith(".fq"):
            # Read if not zipped.
            awk_read = "cat"
        else:
            awk_read = "gunzip -c"
        awk_one_liner =  """ awk '((NR-2)%4==0){ 
                                read=$1;total++;count[read]++;
                                print $0;
                             }
                             END{
                             for(read in count){
                             if(!max||count[read]>max) {
                                 max=count[read];
                                 maxRead=read};
                                 if(count[read]==1){
                                    unique++
                                 }
                             };

                             print "#AWK",
                                   total,
                                   unique,
                                   unique*100/total,
                                   maxRead,
                                   count[maxRead],
                                   count[maxRead]*100/total}'"""
        awk = awk_read + " " + self.filename + " | " + awk_one_liner
        out = Popen([awk], shell = True, stdout = PIPE, stderr = PIPE)
        min_length = ""
        max_length = None
        cum_length = 0
        A, T, C, G, N = [0]*5
        for n, line in enumerate(out.stdout):
            if line.startswith("#AWK"):
                stats = OrderedDict(zip(stat_header, line.strip().split(" ")[1:]))
            else:
                line = line.strip()
                length = len(line)
                cum_length += length
                A += line.count("A")
                T += line.count("T")
                C += line.count("C")
                G += line.count("G")
                N += line.count("N")
                if length < min_length:
                    min_length = length
                if length > max_length:
                    max_length = length
        stats["cum_length"] = cum_length # Stat line at end + 1
        stats["A_count"] = A
        stats["T_count"] = T
        stats["C_count"] = C
        stats["G_count"] = G
        stats["N_count"] = N 
        stats["bases"] = A + T + C + G
        stats["GC_count"] = (G + C) / float(cum_length)
        stats["min_length"] = min_length
        stats["avg_length"] = cum_length / float(stats["Total_Reads"])
        stats["max_length"] = max_length
        return stats


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='BAM-Toolbox v0.1',
                  options_first=False)
    if not os.path.exists(args["<fq>"]):
        with indent(4):
            puts_err("\n" + colored.red(args["<fq>"] + " does not exist") + "\n")

    fq = fastq(args["<fq>"])
    fq_name = os.path.basename(args["<fq>"])
    stats = fq.calculate_fastq_stats()
    dict_to_eav(fq_name, stats)

    







