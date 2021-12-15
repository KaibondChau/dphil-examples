#!/bin/bash

############ Define Help function
function Help()
{
   # Display Help
   echo "Copy number estimation as per" 
   echo "Stoesser N, Phan HTT, Seale AC, Aiken Z, Thomas S, Smith M, Wyllie D, George R, Sebra R," 
   echo "Mathers AJ, Vaughan A, Peto TEA, Ellington MJ, Hopkins KL, Crook DW, Orlek A, Welfare W," 
   echo "Cawthorne J, Lenney C, Dodgson A, Woodford N, Walker AS; TRACE Investigators' Group." 
   echo "Genomic Epidemiology of Complex, Multispecies, Plasmid-Borne blaKPC Carbapenemase in Enterobacterales in the United Kingdom from 2009 to 2014." 
   echo "Antimicrob Agents Chemother. 2020 Apr 21;64(5):e02244-19. doi: 10.1128/AAC.02244-19. PMID: 32094139; PMCID: PMC7179641."
   echo 
   echo "Maps short reads onto reference genome for estimating mean depth at positions of interest"
   echo "Currently only generates sorted indexed BAM for manual searching via e.g."
   echo "samtools depth -r <chromosomeID>:<position-position> <sorted indexed BAM> > <output.txt>"
   echo 
   echo "Syntax: scriptTemplate [-r|a|b|h]"
   echo "options:"
   echo "r     Reference.fa"
   echo "a     Forward read.fq.gz"
   echo "b     Reverse read.fq.gz"
   echo "h     Print this Help"
   echo
}



############ Main program

while getopts r:a:b:"h" option
do
    case "$option" in
        r) ref=${OPTARG};;
        a) read1=${OPTARG};;
        b) read2=${OPTARG};;
        h) # display Help
         Help
         exit;;
       \?) # incorrect option
         echo "Error: Invalid option or missing argument | Use -h for help"
         exit;;
    esac
done

# Exit if no options specified
if ((OPTIND == 1))
then
    echo "No options specified | Use -h for help"
    exit
fi

# Exit if partially missing options
if [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$ref" ]
then
    echo "Missing required options | Use -h for help"
    exit
fi


# Index reference
bwa index -p indexed_ref_bwamem $ref

# Map short reads to reference
bwa mem indexed_ref_bwamem $read1 $read2 > mapped_reads_bwamem.sam

# Convert to BAM
samtools view -S -b mapped_reads_bwamem.sam > mapped_reads_bwamem.bam

# Sort BAM
samtools sort mapped_reads_bwamem.bam -o mapped_reads_bwamem_sorted.bam

# Index BAM
samtools index mapped_reads_bwamem_sorted.bam

