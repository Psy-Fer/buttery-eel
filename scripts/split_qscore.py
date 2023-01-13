#!/usr/bin/env python3

import argparse
import sys
import numpy as np

"""
James M. Ferguson (j.ferguson@garvan.org.au)
Genomic Technologies
Garvan Institute
Copyright 2022

split_qscore.py


----------------------------------------------------------------------------
MIT License

Copyright (c) 2022 James M. Ferguson

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

def find_score(read, ftype):
    '''
    attempt to find score in fastq/sam record
    '''
    if ftype == "fastq":
        header = read[0].split(" ")
        for i in header:
            if i[:4] == "qs:i":
                score = int(i[5:])
            elif i.split("=")[0] in ["qscore", "mean_qscore"]:
                score = int(i.split("=")[1])
            else:
                score = None
    
    elif ftype == "sam":
        sread = read.split("\t")
        for i in sread:
            if i[:4] == "qs:i":
                score = int(i[5:])
            else:
                score = None
    return score

def calculate_qscore(qstring):
    '''
    calculate a qscore from a qstring
    '''
    qs = (np.array(qstring, 'c').view(np.uint8) - 33)
    mean_err = np.exp(qs * (-np.log(10) / 10.)).mean()
    score = -10 * np.log10(max(mean_err, 1e-4))
    return score

def test_read(qscore, score):
    '''
    Compare qscore and score
    '''
    if score >= qscore:
        return True
    return False

def write_fastq(f, read):
    '''
    write fastq file
    '''
    for line in read:
        f.write(line)
        f.write("\n")


def write_sam_header(f, header):
    '''
    write sam file
    '''
    for line in header:
        f.write(line)

def write_sam(f, read):
    '''
    write sam file
    '''
    f.write(read)
    f.write("\n")



def main():
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:

    python3 split_qscore.py -q 9 in.fastq output/path/

    # outputs 2 files, labelled in.pass_9.fastq and out.fail_9.fastq
    """

    VERSION = "0.0.1"

    parser = MyParser(description="split a fastq or sam file by qscore value into pass/fail files",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args
    parser.add_argument("input",
                        help="fastq, sam file or - for STDIN/pipe")
    parser.add_argument("output",
                        help="output path to write pass/fail files")
    parser.add_argument("-p", "--prefix", default="reads",
                        help="filename prefix to give pass/fail output files, eg -p example would give example.pass.fastq & example.fail.fastq")
    parser.add_argument("-q", "--qscore", type=int, default=9,
                        help="qscore to split on: pass=>qscore")
    parser.add_argument("-f", "--format", choices=["fastq", "sam"],
                        help="if using - for stdin, give format as fastq or sam")


    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    pass_reads = 0
    fail_reads = 0
    total_reads = 0
    # hist_list = [0,]*60

    # detect input
    ftype = get_file_type(args.input)
    if ftype is None:
        sys.stderr.write("ERROR: cannot derive file type of input. must be - for stdin or one of the following: {}\n".format(["fastq", "fq", "sam"]))
        parser.print_help(sys.stderr)
        sys.exit(1)
    if ftype == "fastq":
        # if fastq, do fastq things
        reads = read_fastq(args.input)
        # create new pass/fail files
        out_pass = open("{}/{}.pass.fastq".format(args.output, args.prefix), 'w')
        out_fail = open("{}/{}.fail.fastq".format(args.output, args.prefix), 'w')
        # process reads
        for read in reads:
            score = find_score(read, ftype)
            if score is None:
                score = calculate_qscore(read[3])
            if test_read(args.qscore, score):
                pass_reads += 1
                write_fastq(out_pass, read)
            else:
                fail_reads += 1
                write_fastq(out_fail, read)
    elif ftype == "sam":
        # if sam, do sam things
        header = get_header(args.input)
        reads = read_sam(args.input)
        # write new sam files with headers
        out_pass = open("{}/{}.pass.sam".format(args.output, args.prefix), 'w')
        out_fail = open("{}/{}.fail.sam".format(args.output, args.prefix), 'w')
        # write the sam header for pass and fail
        write_sam_header(out_pass, header)
        write_sam_header(out_fail, header)
        # process sam
        for read in reads:
            score = find_score(read, ftype)
            if score is None:
                score = calculate_qscore(read.split("\t")[10])
            if test_read(args.qscore, score):
                pass_reads += 1
                write_sam(out_pass, read)
            else:
                fail_reads += 1
                write_sam(out_fail, read)
    elif ftype == "si":
        # do detection on first line
        # if stdin, detect type, then do the things
        if not args.format:
            sys.stderr.write("ERROR: please provide format of stdin with -f/--format fastq or sam\n")
            parser.print_help(sys.stderr)
            sys.exit(1)
        if args.format == "fastq":
            reads = stream_fastq()
            out_pass = open("{}/{}.pass.fastq".format(args.output, args.prefix), 'w')
            out_fail = open("{}/{}.fail.fastq".format(args.output, args.prefix), 'w')
            for read in reads:
                score = find_score(read, args.format)
                if score is None:
                    score = calculate_qscore(read[3])
                if test_read(args.qscore, score):
                    pass_reads += 1
                    write_fastq(out_pass, read)
                else:
                    write_fastq(out_fail, read)
                    fail_reads += 1
        elif args.format == "sam":
            reads = stream_sam()
            out_pass = open("{}/{}.pass.sam".format(args.output, args.prefix), 'w')
            out_fail = open("{}/{}.fail.sam".format(args.output, args.prefix), 'w')
            # get the header, as first output from generator
            header = next(reads)
            # write the sam header for pass and fail
            write_sam_header(out_pass, header)
            write_sam_header(out_fail, header)
            for read in reads:
                score = find_score(read, args.format)
                if score is None:
                    score = calculate_qscore(read.split("\t")[10])
                if test_read(args.qscore, score):
                    pass_reads += 1
                    write_sam(out_pass, read)
                else:
                    fail_reads += 1
                    write_sam(out_fail, read)
        else:
            sys.stderr.write("ERROR: please provide format of stdin with -f/--format fastq or sam\n")
            parser.print_help(sys.stderr)
            sys.exit(1)
    
    # close files
    out_pass.close()
    out_fail.close()

    # Do some stats
    
    total_reads = pass_reads + fail_reads
    pass_fraction = round((pass_reads/total_reads)*100, 2)

    sys.stderr.write("Done!\n")
    print("~~~ STATS ~~~:\n")
    print("pass_reads: {}".format(pass_reads))
    print("fail_reads: {}".format(fail_reads))
    print("total_reads: {}".format(total_reads))
    print("pass fraction: {}%".format(pass_fraction))

    # write out pass/fail files
    # write a summary showing basic statistics of number of reads in each
    # also give a basic histogram of read qscore distribution
    # make this a separate stat command to get that first to then split on approptiate qscore


if __name__ == '__main__':
    main()
