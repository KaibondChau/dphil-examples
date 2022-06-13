#!/bin/bash

############ Define Help function
function Help()
{
   # Display Help
   echo "Simple VSEARCH subsampling script"
   echo 
   echo "Syntax: scriptTemplate [-i]"
   echo "options:"
   echo "i     Input .bam"
   echo
}
############ Define and validate options
# Define options
while getopts i:"h" option
do
    case "$option" in
        i) bam=${OPTARG};;
        h) # display Help
         Help
         exit;;
       \?) # incorrect option
         echo "Error: Invalid option or missing argument | Use -h for help"
         exit;;
    esac
done

filename="${bam%.*}"

# Main

samtools bam2fq $bam > $filename".fq"
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_1000000.fq" --sample_size 2000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_2000000.fq" --sample_size 4000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_4000000.fq" --sample_size 8000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_6000000.fq" --sample_size 12000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_10000000.fq" --sample_size 20000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_20000000.fq" --sample_size 40000000
/home/ndm.local/kevinchau/vsearch/vsearch-2.21.1-linux-x86_64-static/bin/vsearch --fastx_subsample $filename".fq" --fastqout $filename"_40000000.fq" --sample_size 80000000

mkdir raw
mkdir input
mv $filename".fq" $bam raw

for f in *.fq; do bash ~/vscode/kc_repos/dphil-misc/reformat_kc.sh -i $f; done