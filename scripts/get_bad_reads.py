#!/usr/bin/env python3

import argparse
import sys
"""
James M. Ferguson (j.ferguson@garvan.org.au)
Genomic Technologies
Garvan Institute
Copyright 2025

ut2.py


----------------------------------------------------------------------------
MIT License

Copyright (c) 2025 James M. Ferguson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def get_readset(readlist):
    '''
    get read set from read list
    '''
    with open(readlist, 'r') as f:
        read_set = set()
        for line in f:
            read_set.add(line.strip())
    return read_set


def read_fastq(infile):
    '''
    read fastq file
    '''
    with open(infile, 'r') as f:
        c = 0
        read = []
        for line in f:
            c += 1
            read.append(line.strip("\n"))
            if c >= 4: 
                yield read
                c = 0
                read = []



def main():
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:

    
    """

    VERSION = "0.0.1"

    parser = MyParser(description="Get parent ID for bad read-id",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args
    parser.add_argument("input",
                        help="fastq to get parent ID using ")
    parser.add_argument("-l", "--list",
                        help="read-id list")


    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    

    read_set = get_readset(args.list)
    sys.stderr.write(read_set[:5])
    reads = read_fastq(args.input)
    # process reads
    for read in reads:
        line1 = read[0].strip().split(" ")
        readID = line1[0][1:]
        parentID = line1[1][15:]
        sys.stderr.write("{}\t{}\n".format(readID, parentID))
        if readID in read_set:
            print(parentID)


    sys.stderr.write("Done!\n")

if __name__ == '__main__':
    main()