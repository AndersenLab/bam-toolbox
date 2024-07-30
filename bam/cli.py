#! /usr/bin/env python
"""
usage:
  bam <command> [<args>...]
  bam -h | --help
  bam --version

commands:
  coverage

"""
from subprocess import call, check_output, CalledProcessError
import sys
import os

from docopt import docopt
from colorama import Fore, just_fix_windows_console, init
from termcolor import colored

import bam


just_fix_windows_console()
init(autoreset=True)

debug = None
if len(sys.argv) == 1:
    debug = [""]

def getScriptPath():
    return os.path.dirname(bam.__file__)


# Stack Overflow: 377017
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def main():
    args = docopt(__doc__, 
                  version='bam-toolbox v1.0',
                  argv = debug,
                  options_first=True)
    argv = [args['<command>']] + args['<args>']
    program_list = {"samtools":"samtools"}
    if args["<command>"] == "setup":
        """
            Use Homebrew to install programs!
        """
        program_installed = program_list.keys()
        for install_name, program in program_list.items():
            check_output(["brew", "tap", "homebrew/science"])
            try:
                print(Fore.BLUE + "    Installing " + install_name)
                check_output(["brew", "install", install_name])
                program_installed.remove(install_name)
            except CalledProcessError:
                try:
                    check_output(["which", program])
                    print(Fore.BLUE + "    " + program + " previously installed")
                    program_installed.remove(install_name)
                except CalledProcessError:
                    print(Fore.RED + "    Error installing " + install_name)
        if len(program_installed) == 0:
            print(Fore.BLUE + "    Programs successfully installed!")
        else:
            print(Fore.RED + "    Error: Not all programs successfully installed: "  + ", ".join(program_installed))
    elif args["<command>"] == "":
        print(__doc__)
        for prog in program_list.values():
            try:
                check_output(["which", prog])
            except CalledProcessError:
                print(Fore.RED + "    " + prog + " not installed. Use a package manager to install or try using 'tb.py setup'\n")
    elif args['<command>'] in ['coverage', 'readgroups', 'fastq']:
        comm = ['python', getScriptPath() + '/' + args["<command>"] + ".py"] + argv
        exit(call(comm))

if __name__ == '__main__':
    main()

