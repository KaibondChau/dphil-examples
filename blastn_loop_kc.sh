#!/bin/bash

# simple script for looped blastn similarity search of specified dir using local sequences as db
# testing transformation of long options which can't be natively handled by getopts

############ Define Help function
function Help()
{
   # Display Help
  
   echo "Looped blastn local query - KC"
   echo
   echo "--ref      Fasta of reference sequences"
   echo "--db       Desired name of BLAST db"
   echo "--dir      Target directory of fasta files"
   echo "--help     Print this Help"
   echo
}

############ Define options
# Transform long options to short 
for arg in "$@"; do
  shift
  case "$arg" in
    "--help")   set -- "$@" "-h" ;;
    "--ref")    set -- "$@" "-r" ;;
    "--db")     set -- "$@" "-d" ;;
    "--dir")    set -- "$@" "-i" ;;
    *)          set -- "$@" "$arg"
  esac
done

# Parse short options
OPTIND=1
while getopts r:d:i:"h" opt
do
  case "$opt" in
  "h") Help; exit 0 ;;
    r) ref=${OPTARG} ;;
    d) db=${OPTARG} ;;
    i) dir=${OPTARG} ;;
   \?) echo "Error: Invalid option or missing argument | Use --help for help"
         exit;;
  esac
done

# Exit if no options specified
if ((OPTIND == 1))
then
    echo "No options specified | Use --help for help"
    exit
fi

# Exit if partially missing options
if [ -z "$ref" ] || [ -z "$db" ] || [ -z "$dir" ]
then
    echo "Missing required options | Use --help for help"
    exit
fi

shift $(expr $OPTIND - 1) # remove options from positional parameters


############ Main program

# Generate database using (local) sequences
makeblastdb -in $ref -out $db -parse_seqids -dbtype nucl

for f in $dir*; do
echo "${f##*/}"
blastn -query $f -db $db -outfmt 6 >> blastn_output.csv
done 
echo "Finished"
