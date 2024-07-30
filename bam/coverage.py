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
  --header                    print header

"""
import sys
from collections import OrderedDict
import os
import re
from subprocess import Popen, PIPE
from datetime import datetime
from collections import OrderedDict

from docopt import docopt
from colorama import Fore, just_fix_windows_console, init
from termcolor import colored

just_fix_windows_console()
init(autoreset=True)

class output_line:

    """
        Entity-Attributes-Value Model
    """
    header_out = False

    def __init__(self, entity, attributes, value, header=False):
        self.entity = entity
        if type(attributes) in [dict, OrderedDict]:
            attributes = [k + "=" + str(v) for k, v in attributes.items()]
        elif type(attributes) != list:
            attributes = [attributes]
        self.attributes = attributes
        self.value = value
        if not output_line.header_out and header:
            print("bam\tcontig\tstart\tend\tproperty\tvalue")
            output_line.header_out = True

    def __setattr__(self, name, value):
        # Value is attribute
        if name == "add_attr" or name == "set_attr":
            if type(value) in [dict, OrderedDict]:
                value = [k + "=" + v for k, v in value.items()]
            elif type(value) != list:
                value = [value]
            if name == "add_attr":
                self.__dict__["attributes"].extend(value)
            else:
                self.__dict__["attributes"] = value
        else:
            self.__dict__[name] = value

    def __repr__(self):
        attributes = '\t'.join(map(str, [x.split("=")[1] for x in self.attributes]))
        out = [self.entity, attributes, self.value]
        output = map(str, out)
        return '\t'.join(output)


class bam_file:

    def __init__(self, fname, mtchr = None):
        self.fname = fname
        self.mtchr = mtchr
        self.parse_header()

    def parse_header(self):
        header, err = Popen(["samtools", "view", "-H", self.fname], stdout=PIPE, stderr=PIPE).communicate()
        if err != b"":
            raise Exception(err)
        self.header = header
        contigs = OrderedDict()
        contig_regions = []
        for x in re.findall(b"@SQ\t[A-Za-z0-9]*SN:(?P<chrom>[A-Za-z0-9_]*)[A-Za-z0-9]*\tLN:(?P<length>[0-9]+)", header):
            contigs[x[0].decode('utf-8')] = int(x[1])
            region = "%s:%s-%s" % (x[0].decode('utf-8'), "1", x[1].decode('utf-8'))
            contig_regions.append(region)
        self.contigs = contigs
        self.contig_regions = contig_regions

        mtchr = [x for x in self.contigs.keys() if x.lower().find("m") == 0]
        if len(mtchr) == 1:
            self.mtchr = mtchr[0]
            print(Fore.BLUE + "\n    Guessing Mitochondrial Chromosome: " + self.mtchr + "\n",
                  file=sys.stderr)

        self.genome_length = sum(contigs.values())
        if mtchr:
            self.nuclear_length = sum([x for x in self.contigs.values() if x != self.contigs[self.mtchr]])


    def sum_coverage(self, region=None):
        for n, i in enumerate(region):
            comm = Popen(["samtools", "depth", "-r", region, self.fname], stdout=PIPE, stderr=PIPE)
            pos_covered = 0
            cum_depth = 0
            for row in comm.stdout:
                chrom, pos, depth = row.strip().split(b"\t")
                pos_covered += 1
                cum_depth += int(depth)
            return pos_covered, cum_depth


def iterate_window(bamfile, size):
    for chrom, size in bamfile.contigs.items():
        for i in range(1, size, window):
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
        chrom_len = bamfile.contigs[chrom]
        if end > chrom_len:
            m = "\n    Specified chromosome end extends beyond chromosome length. Set to max of: "
            print(Fore.YELLOW + m + str(chrom_len) + "\n", file=sys.stderr)
            end = chrom_len

        region = "{c}:{s}-{e}".format(c=chrom, s=start, e=end + 1)
        pos_covered, cum_depth = bamfile.sum_coverage(region)

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
                  version='BAM-Toolbox v0.1',
                  options_first=False)
    if args["<bam>"]:
        # Add check for file here
        bam_name = os.path.basename(args["<bam>"]).replace(".bam", "")
        b = bam_file(args["<bam>"], args["--mtchr"])

    if args["<chrom:start-end>"]:
        """
            Calculate coverage in a given region or regions
        """
        calc_coverage(b, args["<chrom:start-end>"])
    elif args["--window"]:
        """
            Calculate coverage across a window of given size.
        """
        window = int(args["--window"])
        regions = iterate_window(b, window)
        calc_coverage(b, regions)

    elif args["--regions"]:
        """
            Calculate coverage in specified regions
        """
        from pybedtools import BedTool
        bed = BedTool(args["--regions"])
        calc_coverage(b, bed[:])
    elif args["<bam>"]:
        """
            Calculate coverage genome wide
        """
        bam = args["<bam>"]
        cov = calc_coverage(b, b.contig_regions)

        # Genomewide depth
        output_dir = {}
        output_dir["chrom"] = "genome"
        output_dir["start"] = 1
        output_dir["end"] = b.genome_length

        bases_mapped = sum([x["bases_mapped"] for x in cov])
        output_dir["ATTR"] = "bases_mapped"
        print(output_line(bam_name, output_dir, bases_mapped))
        output_dir["ATTR"] = "depth_of_coverage"
        coverage = bases_mapped / float(b.genome_length)
        print(output_line(bam_name, output_dir, coverage))

        output_dir["ATTR"] = "breadth_of_coverage"
        breadth = sum([x["pos_covered"] for x in cov]) / float(b.genome_length)
        print(output_line(bam_name, output_dir, breadth))

        output_dir["ATTR"] = "positions_mapped"
        pos_mapped = sum([x["pos_covered"] for x in cov])
        print(output_line(bam_name, output_dir, pos_mapped))

        if b.mtchr:
            # Nuclear
            output_dir["chrom"] = "nuclear"
            output_dir["end"] = b.nuclear_length
            bases_mapped = sum([x["bases_mapped"] for x in cov if x["chrom"] != b.mtchr])
            output_dir["ATTR"] = "bases_mapped"
            print(output_line(bam_name, output_dir, bases_mapped))

            output_dir["ATTR"] = "depth_of_coverage"
            coverage = bases_mapped / float(b.nuclear_length)
            print(output_line(bam_name, output_dir, coverage))

            output_dir["ATTR"] = "breadth_of_coverage"
            breadth = sum([x["pos_covered"] for x in cov if x["chrom"] != b.mtchr]) / float(b.nuclear_length)
            print(output_line(bam_name, output_dir, breadth))

            output_dir["ATTR"] = "positions_mapped"
            pos_mapped = sum([x["pos_covered"] for x in cov if x["chrom"] != b.mtchr])
            print(output_line(bam_name, output_dir, pos_mapped))

            # mt:nuclear ratio
            output_dir = {"chrom": "genome", "start": 1, "end": b.nuclear_length, "ATTR": "mt_nuclear_ratio"}
            mt_nuc = [x for x in cov if x["chrom"] == b.mtchr][0]["depth_of_coverage"] / coverage
            print(output_line(bam_name, output_dir, mt_nuc))
