#! /usr/bin/env python
"""
usage:
  bam coverage header [--eav]
  bam coverage <bam> [--mtchr=<mtchr>] [--eav]
  bam coverage <bam> <chrom:start-end> [--eav]

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --eav                       Print eav representation of data.

"""
from docopt import docopt
import sys
from collections import defaultdict, OrderedDict
from clint.textui import colored, indent, puts_err
import os
from eav import *
import re
from subprocess import Popen, PIPE
from pysam import AlignmentFile



def get_contigs(bam):
    header, err = Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    # Extract contigs from header and convert contigs to integers
    contigs = OrderedDict()
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
        contigs[x[0]] = int(x[1])
    return contigs

def coverage(bam, mtchr = None):
    # Check to see if file exists
    if os.path.isfile(bam) == False and bam != "-":
        raise Exception("Bam file does not exist")
    contigs = get_contigs(bam)


    # Guess mitochondrial chromosome
    if mtchr == None:
        mtchr = [x for x in contigs if x.lower().find("m") == 0]
        if len(mtchr) != 1:
            mtchr = None
        else:
            mtchr = mtchr[0]
            with indent(4):
                puts_err(colored.blue("\nGuessing Mitochondrial Chromosome: " + mtchr + "\n"))

    coverage_dict = OrderedDict()
    for c in contigs.keys():
        command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
        coverage_dict[c] = {}
        coverage_dict[c]["Bases_Mapped"], coverage_dict[c]["Sum_of_Depths"] = map(int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))
        coverage_dict[c]["Breadth_of_Coverage"] = coverage_dict[c]["Bases_Mapped"] / float(contigs[c])
        coverage_dict[c]["Depth_of_Coverage"] = coverage_dict[c]["Sum_of_Depths"] / float(contigs[c])
        coverage_dict[c]["Length"] = int(contigs[c])

    # Calculate Genome Wide Breadth of Coverage and Depth of Coverage
    genome_length = float(sum(contigs.values()))
    coverage_dict["genome"] = {}
    coverage_dict["genome"]["Length"] = int(genome_length)
    coverage_dict["genome"]["Bases_Mapped"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Sum_of_Depths"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Breadth_of_Coverage"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    coverage_dict["genome"]["Depth_of_Coverage"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)

    if mtchr != None:
        # Calculate nuclear breadth of coverage and depth of coverage
        ignore_contigs = [mtchr, "genome", "nuclear"]
        coverage_dict["nuclear"] = {}
        coverage_dict["nuclear"]["Length"] = sum([x["Length"] for k,x in coverage_dict.iteritems() if k not in ignore_contigs ])
        coverage_dict["nuclear"]["Bases_Mapped"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Sum_of_Depths"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Breadth_of_Coverage"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])
        coverage_dict["nuclear"]["Depth_of_Coverage"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])

        # Calculate the ratio of mtDNA depth to nuclear depth
        coverage_dict["genome"]["mt_ratio"] = coverage_dict[mtchr]["Depth_of_Coverage"] / float(coverage_dict["nuclear"]["Depth_of_Coverage"])

    # Flatten Dictionary 
    coverage = []
    for k,v in coverage_dict.items():
        for x in v.items():
            coverage += [(k,x[0], x[1])]
    return coverage

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  options_first=False)

    if args["<chrom:start-end>"]:
        chrom, start, end = re.split("[:-]", args["<chrom:start-end>"]) 
        start, end = int(start), int(end)
        bamfile = AlignmentFile(args["<bam>"])

        # If end extends to far, adjust for chrom
        chrom_len = bamfile.lengths[bamfile.gettid(chrom)]
        if end > chrom_len:
            with indent(4):
                puts_err(colored.yellow("\nSpecified chromosome end extends beyond chromosome length. Set to max of: " + str(chrom_len) + "\n"))
                end = chrom_len


        region = bamfile.pileup(chrom, start, end+1, truncate = True, max_depth = 1e8)
        cum_depth = 0
        pos_covered = 0
        for n,i in enumerate(region):
            pos_covered += 1
            cum_depth += i.nsegments
            print n, pos_covered, i.nsegments
        coverage = cum_depth / float(n+1)
        breadth  = pos_covered / float(end - start + 1)
        if args["--eav"]:
            print eav(args["<bam>"], OrderedDict([("ATTR", "depth_of_coverage"), ("chrom", chrom), ("start", start), ("end", end)]), coverage)
            print eav(args["<bam>"], OrderedDict([("ATTR", "breadth_of_coverage"), ("chrom", chrom), ("start", start), ("end", end)]), breadth)
        else:
            print args["<bam>"] + "\tdepth_of_coverage\t" + str(coverage)
            print args["<bam>"] + "\tbreadth_of_tcoverage\t" + str(breadth)

    elif args["header"]:
        if args["--eav"]:
            print(eav.header)
        else:
            print("CONTIG\tATTR\tVALUE")
    elif args["<bam>"]:
        bam = args["<bam>"]
        cov = coverage(bam, args["--mtchr"])
        if args["--eav"]:
            for i in cov:
                out = eav(bam, {"CONTIG": i[0], "ATTR": i[1]}, i[2])
                print(out)
        else:
            for i in cov:
                print '\t'.join(map(str,i))



