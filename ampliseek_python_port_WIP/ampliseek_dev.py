#!/usr/bin/env python3

from os import close, name
import sys
import argparse
import gzip
import subprocess
from Bio import SeqIO



# Arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="Ampliseek dev",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')

    # input forward fastq.gz file
    required.add_argument('-f', '--forward_reads', action='store',
                          required=True,
                          help='Forward fastq.gz file')

    # input reverse fastq.gz file
    required.add_argument('-r', '--reverse_reads', action='store',
                          required=True,
                          help='Reverse fastq.gz file')
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    return args


# validate forward and reverse fq.gz files
def validate_input():
    args = get_args(),
    with gzip.open(args.forward_reads, "rt") as handle:
        fasta_content = list(SeqIO.parse(handle, "fastq"))
        assert len(fasta_content) >= 1, 'No records found in forward fq.gz file'

    with gzip.open(args.reverse_reads, "rt") as handle:
        fasta_content = list(SeqIO.parse(handle, "fastq"))
        assert len(fasta_content) >= 1, 'No records found in reverse fq.gz file'


# get filename for writing files - currently relies on standard output from Illumina name_<read#>.fq.gz/.fastq.gz
def get_filename(): 
    args = get_args()
    name = args.forward_reads.split('_1')
    filename = (name[0].split('/'))[-1]
    return(filename)
    

def interleave_reads():
    args = get_args()
    filename = get_filename()
    reformat_command = [
        "reformat.sh",
        str("in1=" + args.forward_reads),
        str("in2=" + args.reverse_reads),
        str("out=" + (str(filename + '_interleaved.fq.gz')))]
    print(str(filename + '_interleaved.fq.gz')) # name works locally but not on analysis1 - python version?
    p = subprocess.Popen(
        reformat_command,
        stderr=subprocess.PIPE)
    err = p.communicate()  # read stderr data
    err = str(err)
    o = open(str(filename + '_stderr'), 'w')  # create output stderr file - format needs cleaning
    o.write(err)  
    o.close()  
    
interleave_reads()




