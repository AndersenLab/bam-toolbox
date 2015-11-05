#! /usr/bin/env python
"""
usage:
  bam readgroups <bam>
  bam readgroups <bam> (--include=<args>...|--exclude=<args>...) 

options:
  -h --help                   Show this screen.
  --version                   Show version.
  --eav                       Print eav representation of data.

"""
from docopt import docopt
import sys
from clint.textui import colored, indent, puts_err
from pysam import AlignmentFile
import os
from subprocess import Popen
import tempfile

from pprint import pprint as pp


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='BAM-Toolbox v0.1',
                  options_first=False)
    if not os.path.exists(args["<bam>"]):
        with indent(4):
            puts_err("\n" + colored.red(args["<bam>"] + " does not exist") + "\n")
    bamfile = AlignmentFile(args["<bam>"])
    readgroups = []
    for x in [x for x in bamfile.header["RG"]]:
        readgroups.append(x["ID"])
    if len(args["--include"] + args["--exclude"]) == 0:
        with indent(4):
            puts_err("Readgroups:")
            with indent(4):
                exit(puts_err(colored.green("\n" + '\n'.join(readgroups)) + "\n"))
    elif args["--include"]:
        include = args["--include"][0].split(",")
        for x in include:
            if x not in readgroups:
                with indent(4):
                    puts_err(colored.red("\n" + x + " not a Readgroup"))

        readgroups = [x for x in readgroups if x not in include]
    else:
        exclude = args["--exclude"][0].split(",")
        for x in exclude:
            if x not in readgroups:
                with indent(4):
                    puts_err(colored.red("\n" + x + " not a Readgroup"))
        
        readgroups = [x for x in readgroups if x in exclude]

    for i in readgroups:
        print i






