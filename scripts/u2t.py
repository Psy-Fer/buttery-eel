#!/usr/bin/env python3

import argparse
import sys
import re
"""
James M. Ferguson (j.ferguson@garvan.org.au)
Genomic Technologies
Garvan Institute
Copyright 2024

ut2.py


----------------------------------------------------------------------------
MIT License

Copyright (c) 2024 James M. Ferguson

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

def get_file_type(infile):
    '''
    find if file is fastq or sam
    '''
    if infile == "-":
        return "si"
    
    if infile.split(".")[-1] in ["fastq", "fq"]:
        return "fastq"

    if infile.split(".")[-1] in ["sam",]:
        return "sam"
    
    return None


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
                

def stream_fastq():
    '''
    stream a fastq file from stdin
    '''
    c = 0
    read = []
    for line in sys.stdin:
        c += 1
        read.append(line.strip("\n"))
        if c >= 4: 
            yield read
            c = 0
            read = []

def get_header(infile):
    '''
    get header of sam file
    '''
    header = []
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == "@":
                header.append(line)
    return header

def read_sam(infile):
    '''
    read sam file
    '''
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == "@":
                continue
            read = line.strip("\n")
            yield read

def stream_sam():
    '''
    stream a sam file from stdin
    '''
    head = False
    header = []
    for line in sys.stdin:
        if not head:
            if line[0] == "@":
                header.append(line)
                continue
            else:
                head = True
                yield header
        read = line.strip("\n")
        yield read



def main():
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:

    python3 u2t.py output.sam -> new_output.sam
    python3 u2t.py output.fastq -> new_output.fastq

    """

    VERSION = "0.0.1"

    parser = MyParser(description="Replace U (uracil) with T (thymine) in a sam or fastq file",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args
    parser.add_argument("input",
                        help="fastq, sam file or - for STDIN/pipe")
    # parser.add_argument("output",
    #                     help="output path to write pass/fail files")
    parser.add_argument("-f", "--format", choices=["fastq", "sam"],
                        help="if using - for stdin, give format as fastq or sam")


    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # detect input
    ftype = get_file_type(args.input)

    if ftype is None:
        sys.stderr.write("ERROR: cannot derive file type of input. must be - for stdin or one of the following: {}\n".format(["fastq", "fq", "sam"]))
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if ftype == "fastq":
        # if fastq, do fastq things
        reads = read_fastq(args.input)
        # process reads
        for read in reads:
            read[1] = re.sub("U", "T", read[1])
            for l in read:
                print(l)
    
    elif ftype == "sam":
        # if sam, do sam things
        header = get_header(args.input)
        reads = read_sam(args.input)
        #write header
        for l in header:
            print(l)
        # process sam
        for read in reads:
            r = read.split("\t")
            r[9] = re.sub("U", "T", r[9])
            # write record
            print("\t".join(r))
    
    elif ftype == "si":
        # do detection on first line
        # if stdin, detect type, then do the things
        if not args.format:
            sys.stderr.write("ERROR: please provide format of stdin with -f/--format fastq or sam\n")
            parser.print_help(sys.stderr)
            sys.exit(1)
        if args.format == "fastq":
            reads = stream_fastq()
            # process reads
            for read in reads:
                read[1] = re.sub("U", "T", read[1])
                for l in read:
                    print(l)
        elif args.format == "sam":
            reads = stream_sam()
            #write header
            for l in header:
                print(l)
            # process sam
            for read in reads:
                r = read.split("\t")
                r[9] = re.sub("U", "T", r[9])
                # write record
                print("\t".join(r))
        else:
            sys.stderr.write("ERROR: please provide format of stdin with -f/--format fastq or sam\n")
            parser.print_help(sys.stderr)
            sys.exit(1)
    


    sys.stderr.write("Done!\n")

if __name__ == '__main__':
    main()
