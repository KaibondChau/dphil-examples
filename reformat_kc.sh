#!/bin/bash

############ Define Help function
function Help()
{
   # Display Help
   echo "Simple reformat.sh script"
   echo 
   echo "Syntax: scriptTemplate [-i]"
   echo "options:"
   echo "i     Input .fq"
   echo
}
############ Define and validate options
# Define options
while getopts i:"h" option
do
    case "$option" in
        i) fq=${OPTARG};;
        h) # display Help
         Help
         exit;;
       \?) # incorrect option
         echo "Error: Invalid option or missing argument | Use -h for help"
         exit;;
    esac
done

filename="${fq%.*}"

# Main

reformat.sh in=$fq out1=$filename"_r1.fq" out2=$filename"_r2.fq"
mv $fq raw
gzip $filename"_r1.fq"
gzip $filename"_r2.fq"
mv $filename"_r1.fq.gz" input
mv $filename"_r2.fq.gz" input
