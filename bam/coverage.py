#! /usr/bin/env python
"""
usage:
  bam coverage <bam> [options] [--mtchr=<mtchr>]
  bam coverage <bam> [options] <chrom:start-end>...
  bam coverage <bam> [options] --window=<size>
  bam coverage <bam> [options] --regions=<gff/bed>

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --header                    print(header

"""
from docopt import docopt
from collections import OrderedDict
from clint.textui import colored, indent, puts_err
import os
from output import *
import re
from pysam import AlignmentFile
from subprocess import Popen, PIPE


def get_contigs(bam):
    header, err = Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    # Extract contigs from header and convert contigs to integers
    contigs = {}
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
        contigs[x[0]] = int(x[1])
    return contigs


def iterate_window(contigs, size):
    for chrom, size in contigs:
        for i in xrange(1, size, window):
            if i + window > size:
                end = size
            else:
                end = i + window - 1
            yield "{chrom}:{i}-{end}".format(**locals())


def calc_coverage(bamfile, regions=None, mtchr=None):
    from pybedtools.cbedtools import Interval
    depths = []
    for region in regions:
        output_dir = OrderedDict()
        if type(region) == Interval:
            # Add one to start as starts are 0 based; ends are 1 based.
            rchrom = str(region.chrom)
            chrom, start, end = rchrom, region.start + 1, region.stop
            output_dir["name"] = region.name
        else:
            chrom, start, end = re.split("[:-]", region)
            start, end = int(start), int(end)
        output_dir["chrom"] = chrom
        output_dir["start"] = start
        output_dir["end"] = end

        # If end extends to far, adjust for chrom
        chrom_len = bamfile.lengths[bamfile.gettid(chrom)]
        if end > chrom_len:
            m = "\nSpecified chromosome end extends beyond chromosome length. Set to max of: "
            with indent(4):
                puts_err(colored.yellow(m + str(chrom_len) + "\n"))
                end = chrom_len

        region = bamfile.pileup(chrom,
                                start,
                                end + 1,
                                truncate=True,
                                max_depth=1e9)
        cum_depth = 0
        pos_covered = 0
        for n, i in enumerate(region):
            pos_covered += 1
            cum_depth += i.nsegments
        length = end - start + 1
        coverage = cum_depth / float(length)
        breadth = pos_covered / float(length)
        output_dir["ATTR"] = "bases_mapped"
        print(output_line(bam_name, output_dir, cum_depth, args["--header"]))
        output_dir["ATTR"] = "depth_of_coverage"
        print(output_line(bam_name, output_dir, coverage))
        output_dir["ATTR"] = "breadth_of_coverage"
        print(output_line(bam_name, output_dir, breadth))
        output_dir["ATTR"] = "length"
        print(output_line(bam_name, output_dir, length))
        output_dir["ATTR"] = "pos_mapped"
        print(output_line(bam_name, output_dir, pos_covered))
        depths.append({"chrom": chrom,
                       "bases_mapped": cum_depth,
                       "pos_covered": pos_covered,
                       "depth_of_coverage": coverage})
    return depths


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='VCF-Toolbox v0.1',
                  options_first=False)
    if args["<bam>"]:
        # Add check for file here
        bam_name = os.path.basename(args["<bam>"]).replace(".bam", "")
        bamfile = AlignmentFile(args["<bam>"])

    if args["<chrom:start-end>"]:
        """
            Calculate coverage in a given region or regions
        """
        calc_coverage(bamfile, args["<chrom:start-end>"])
    elif args["--window"]:
        """
            Calculate coverage across a window of given size.
        """
        contigs = [(x["SN"], x["LN"]) for x in bamfile.header["SQ"]]
        window = int(args["--window"])
        regions = iterate_window(contigs, window)
        calc_coverage(bamfile, regions)

    elif args["--regions"]:
        """
            Calculate coverage in specified regions
        """
        from pybedtools import BedTool
        bed = BedTool(args["--regions"])
        calc_coverage(bamfile, bed[:])
    elif args["<bam>"]:
        """
            Calculate coverage genome wide
        """
        bam = args["<bam>"]
        ncontig = len(bamfile.references)
        chroms = [
            "{0}:{1}-{2}".format(*x) for x in zip(bamfile.references, ["1"] * ncontig, map(str, bamfile.lengths))]
        if not args["--mtchr"]:
            mtchr = [x for x in bamfile.references if x.lower().find("m") == 0]
        if len(mtchr) != 1:
            mtchr = None
        else:
            mtchr = mtchr[0]
            with indent(4):
                puts_err(colored.blue("\nGuessing Mitochondrial Chromosome: " + mtchr + "\n"))
        depths = []
        cov = calc_coverage(bamfile, chroms, mtchr)

        # Genomewide depth
        output_dir = {}
        genome_length = sum([x for x in bamfile.lengths])
        output_dir["start"] = 1
        output_dir["end"] = genome_length
        output_dir["chrom"] = "genome"

        bases_mapped = sum([x["bases_mapped"] for x in cov])
        output_dir["ATTR"] = "bases_mapped"
        print(output_line(bam_name, output_dir, bases_mapped))
        output_dir["ATTR"] = "depth_of_coverage"
        coverage = bases_mapped / float(genome_length)
        print(output_line(bam_name, output_dir, coverage))

        output_dir["ATTR"] = "breadth_of_coverage"
        breadth = sum([x["pos_covered"] for x in cov]) / float(genome_length)
        print(output_line(bam_name, output_dir, breadth))

        output_dir["ATTR"] = "positions_mapped"
        pos_mapped = sum([x["pos_covered"] for x in cov])
        print(output_line(bam_name, output_dir, pos_mapped))

        if mtchr:
            # Nuclear
            nuclear_length = sum([y for x, y in zip(bamfile.references, bamfile.lengths) if x != mtchr])
            output_dir["end"] = nuclear_length
            output_dir["chrom"] = "nuclear"
            bases_mapped = sum([x["bases_mapped"] for x in cov if x["chrom"] != mtchr])
            output_dir["ATTR"] = "bases_mapped"
            print(output_line(bam_name, output_dir, bases_mapped))

            output_dir["ATTR"] = "depth_of_coverage"
            coverage = bases_mapped / float(genome_length)
            print(output_line(bam_name, output_dir, coverage))

            output_dir["ATTR"] = "breadth_of_coverage"
            breadth = sum([x["pos_covered"] for x in cov if x["chrom"] != mtchr]) / float(genome_length)
            print(output_line(bam_name, output_dir, breadth))

            output_dir["ATTR"] = "positions_mapped"
            pos_mapped = sum([x["pos_covered"] for x in cov if x["chrom"] != mtchr])
            print(output_line(bam_name, output_dir, pos_mapped))

            # mt:nuclear ratio
            output_dir = {"start": 1, "end": genome_length, "chrom": "genome", "ATTR": "mt_nuclear_ratio"}
            mt_nuc = [x for x in cov if x["chrom"] == mtchr][0]["depth_of_coverage"] / coverage
            print(output_line(bam_name, output_dir, mt_nuc))
