#!/bin/bash

############ Define Help function
function Help()
{
   # Display Help
   echo "Ampliseq amplicon mapper v1"
   echo 
   echo "Syntax: scriptTemplate [-f|r|h]"
   echo "options:"
   echo "f     Forward read.fq.gz"
   echo "r     Reverse read.fq.gz"
   echo "o     Output for coverage"
   echo "h     Print this Help"
   echo
}


############ Main program
# Currently writes all intermediate files in addition to specified output.txt (-o)

# Define options
while getopts f:r:o:"h" option
do
    case "$option" in
        f) read1=${OPTARG};;
        r) read2=${OPTARG};;
        o) output=${OPTARG};;
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
if [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$output" ]
then
    echo "Missing required options | Use -h for help"
    exit
fi


# Interleave reads
reformat.sh in1=$read1 in2=$read2 out=interleaved.fq.gz  

# Trim
bbduk2.sh in=interleaved.fq.gz out=trimmed.fq.gz mink=6 ktrim=r k=19 hdist=1 edist=0 ref=adapters.fa minlength=75 qin=33

# Merge trimmed treads
bbmerge-auto.sh in=trimmed.fq.gz out=merged.fq k=62 extend2=50 ecct

# Map merged reads against reference
bbmapskimmer.sh in=merged.fq outm=mapped.sam ref=ampliseq_targets_only.fasta ambig=all saa=f sam=1.3 semiperfectmode=t int=f

# Pileup coverage
pileup.sh in=mapped.sam out=$output