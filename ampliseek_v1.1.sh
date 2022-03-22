#!/bin/bash

############ Define Help function
function Help()
{
   # Display Help
   echo "Ampliseq amplicon mapper v1"
   echo 
   echo "Syntax: scriptTemplate [-f|r|o|l|h]"
   echo "options:"
   echo "f     Forward read.fq.gz"
   echo "r     Reverse read.fq.gz"
   echo "o     Output for coverage"
   echo "l     Silent mode and write logs"
   echo "h     Print this Help"
   echo
}

############ Define and validate options
logs=FALSE
# Define options
while getopts f:r:o:"lh" option
do
    case "$option" in
        f) read1=${OPTARG};;
        r) read2=${OPTARG};;
        o) output=${OPTARG};;
        l) logs=TRUE;;
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

# Exit if either input fq.gz is empty
if [ ! -s "$read1" ] || [ ! -s "$read2" ]
then
    echo "Specified fastq.gz does not exist or is empty"
    exit
fi

############ Main program

if [ "$logs" = TRUE ];
then
    echo "Printing stderr to files"
    echo "Interleaving reads..."
    # Interleave reads
    reformat.sh in1=$read1 in2=$read2 out=interleaved.fq.gz 1>reformat.txt 2>&1 
    echo "Complete" # timer to add to each step i.e. completed in ...
    echo "Trimming..."
    # Trim
    bbduk2.sh in=interleaved.fq.gz out=trimmed.fq.gz mink=6 ktrim=r k=19 hdist=1 edist=0 ref=adapters.fa minlength=75 qin=33 1>trim.txt 2>&1
    echo "Complete"
    echo "Merging trimmed reads..."
    # Merge trimmed treads
    bbmerge-auto.sh in=trimmed.fq.gz out=merged.fq.gz k=62 extend2=50 ecct 1>merge.txt 2>&1
    echo "Complete"
    echo "Mapping merged reads..."
    # Map merged reads against reference
    bbmapskimmer.sh in=merged.fq.gz outm=mapped.sam ref=ampliseq_targets_only.fasta ambig=all minscaf=73 saa=f sam=1.3 semiperfectmode=t int=f 1>bbmapskimmer.txt 2>&1
    echo "Complete"
    echo "Writing output..."
    # Pileup coverage
    pileup.sh in=mapped.sam out=$output # need to silence here
    echo "Complete"
    echo "All finished" # timer to add total run time
else
    # Same as above but no stderr direction - # Todo abstract command into a separate function (duplication of code hinders maintainability) 
    reformat.sh in1=$read1 in2=$read2 out=interleaved.fq.gz 
    bbduk2.sh in=interleaved.fq.gz out=trimmed.fq.gz mink=6 ktrim=r k=19 hdist=1 edist=0 ref=adapters.fa minlength=75 qin=33 
    bbmerge-auto.sh in=trimmed.fq.gz out=merged.fq.gz k=62 extend2=50 ecct 
    bbmapskimmer.sh in=merged.fq.gz outm=mapped.sam ref=ampliseq_targets_only.fasta ambig=all minscaf=73 saa=f sam=1.3 semiperfectmode=t int=f
    pileup.sh in=mapped.sam out=$output
fi
