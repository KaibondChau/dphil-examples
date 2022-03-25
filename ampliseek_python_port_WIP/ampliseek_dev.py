#!/usr/bin/env python3


import sys
import argparse
import gzip
import subprocess
import time
import os
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
def validate_input(args):
    with gzip.open(args.forward_reads, "rt") as handle:
        fasta_content = list(SeqIO.parse(handle, "fastq"))
        assert len(fasta_content) >= 1, 'No records found in forward fq.gz file'

    with gzip.open(args.reverse_reads, "rt") as handle:
        fasta_content = list(SeqIO.parse(handle, "fastq"))
        assert len(fasta_content) >= 1, 'No records found in reverse fq.gz file'


# todo wrap below 2 chunks into if intermediate and logs wanted

# get filename for writing files - currently relies on standard output from Illumina name_<read#>.fq.gz/.fastq.gz # todo could add as argument instead
def get_filename(args): 
    name = args.forward_reads.split('_1')
    filename = (name[0].split('/'))[-1]
    return(filename)

# get directories to write optional intermediate outputs and logs    
def get_dirs(filename):
    log_dir = str(os.getcwd() + "/ampliseek_work/" + filename +"/logs/")
    output_dir = str(os.getcwd() + "/ampliseek_work/" + filename +"/outputs/")
    cov_dir = str(os.getcwd() + "/ampliseek_work/" + filename +"/")
    return(log_dir, output_dir, cov_dir)


# make directories
def prep_dirs(filename):
    os.makedirs(str(os.getcwd() + "/ampliseek_work/" + filename +"/logs"))
    os.makedirs(str(os.getcwd() + "/ampliseek_work/" + filename +"/outputs"))
    

# todo consider bash.sh commands as params to input into universal function

# Interleaves forward and reverse reads using reformat.sh
def interleave_reads(args, filename, log_dir, output_dir):
    reformat_command = [
        "reformat.sh",
        str("in1=" + args.forward_reads),
        str("in2=" + args.reverse_reads),
        str("out=" + output_dir + filename + '_interleaved.fq.gz')]
    p_interleave = subprocess.run(
        reformat_command,
        stdout=subprocess.DEVNULL, # supress stdout
        stderr=subprocess.PIPE, # capture stderr
        universal_newlines=True) # capture as str not bytes (supersedes decode)
    e = open(str(log_dir + filename + '_interleave_stderr'), 'w')  # create output stderr file  
    e.write(p_interleave.stderr)
    e.close()
    

# Trims interleaved reads using bbduk2.sh and standard Illumina adapters
def trim_reads(filename, log_dir, output_dir):
    bbduk_command= [
        "bbduk2.sh",
        str("in="+ str(output_dir + filename + '_interleaved.fq.gz')),
        str("out="+ output_dir + filename + '_trimmmed.fq.gz'),
        "mink=6",
        "ktrim=r",
        "k=19",
        "hdist=1",
        "edist=0",
        "ref=adapters.fa", # standard illumina adapters from BBTools #todo consider param
        "minlength=75",
        "qin=33"]
    p_bbduk = subprocess.run(
        bbduk_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    e = open(str(log_dir + filename + '_trim_stderr'), 'w')  # create output stderr file 
    e.write(p_bbduk.stderr)  
    e.close()  


# Merge trimmed reads using bbmerge-auto.sh 
def merge_reads(filename, log_dir, output_dir):
    bbmerge_command= [
        "bbmerge-auto.sh",
        str("in="+ str(output_dir + filename + '_trimmmed.fq.gz')),
        str("out="+ output_dir + filename + '_merged.fq.gz'),
        "k=62",
        "extend2=50",
        "ecct"]
    p_bbmerge = subprocess.run(
        bbmerge_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    e = open(str(log_dir + filename + '_merge_stderr'), 'w')  # create output stderr file 
    e.write(p_bbmerge.stderr)  
    e.close()  


# Map merged reads against Ampliseq AMR panel targets using bbmapskimmer.sh semiperfect mode
def map_reads(filename, log_dir, output_dir):
    bbmap_command= [
        "bbmapskimmer.sh",
        str("in="+ str(output_dir + filename + '_merged.fq.gz')),
        "out=stdout.sam",
        "ref=ampliseq_targets_only.fasta", # ampliseq AMR panel targets - argument for generalisation
        "ambig=all",
        "minid=0.80", # added options for relaxed mapping (minid needs to be set as 1 below idfilter)
        "idfilter=0.98", # arg here for 0.9 or 0.98
        "minscaf=73",
        "saa=f",
        "sam=1.3",
        #"semiperfectmode=t", # comment out for relaxed mapping
        "int=f"]
    p_bbmap = subprocess.run(
        bbmap_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    e = open(str(log_dir + filename + '_mapping_stderr'), 'w')  # create output stderr file 
    e.write(p_bbmap.stderr)  
    e.close()  
    o = open(str(output_dir + filename + '_mapped.sam'), 'w')  # create output stdout file 
    o.write(p_bbmap.stdout)  
    o.close() 

# Output coverage using pileup.sh
def pileup_reads(filename, log_dir, output_dir, cov_dir):
    pileup_command= [
        "pileup.sh",
        str("in="+ str(output_dir + filename + '_mapped.sam')),
        "out=stdout"]
    p_pileup = subprocess.run(
        pileup_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    e = open(str(log_dir + filename + '_pileup_stderr'), 'w')  # create output stderr file 
    e.write(p_pileup.stderr)  
    e.close()  
    o = open(str(cov_dir + filename + '_coverage.txt'), 'w')  # create output stdout file 
    o.write(p_pileup.stdout)  
    o.close()  

def main_function():
    start = time.time()
    args = get_args()
    filename = get_filename(args)
    log_dir, output_dir, cov_dir = get_dirs(filename)
    #validate_input(args) # todo works locally but not on analysis1
    prep_dirs(filename)
    print("Interleaving reads...")
    interleave_reads(args, filename, log_dir, output_dir)
    print("Completed")
    print("Trimming reads...")
    trim_reads(filename, log_dir, output_dir)
    print("Completed")
    print("Merging reads...")
    merge_reads(filename, log_dir, output_dir)
    print("Completed")
    print("Mapping reads...")
    map_reads(filename, log_dir, output_dir)
    print("Completed")
    print("Pileup")
    pileup_reads(filename, log_dir, output_dir, cov_dir)
    end = time.time()
    print(str("All done in " + str((round(end - start, 2))) + " seconds"))
    
main_function()


