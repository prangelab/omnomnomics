#! /bin/bash

## There's a bug in HOMER that can cause bedGraph2BigWig to fail. Result is a .bw file tha's in fact a .bg
## This scripts checks if the HOMER generated .bw file is binary, and therefore prolly a bigwig, or is text mode, and therefore prolly still a bedgraph.
## In case of the latter: remove some junk, sort the file, and make it a bigwig afterall!


#--------------------------------------------------------------------------------------------------------------
# Function definitions
# Usage function
function display_help {
	echo ""
	echo "Usage: $0 -i <.bw> -g GENOME -c <chrom_size_dir> -h"
	echo ""
	echo "	-i:			Input bigWig file"
	echo "				Required argument."
	echo "	-g:			Genome"
	echo "				Default: hg38"
	echo "	-c:			Directory containing <genome>_chrom_sizes.2_column"
	echo "				Required argument."
	echo ""
	echo "	-h:			Print this help message"
	echo ""
	exit 1
}

#--------------------------------------------------------------------------------------------------------------
# Display help and exit if not enough arguments

if [ $# -lt 1 ]
then
	display_help
fi

#--------------------------------------------------------------------------------------------------------------
# Define options

while getopts ":i:g:c:h" opt; do
  case $opt in
	i)
		MYBW=$OPTARG
		;;
	g)
		MYGENOME=$OPTARG
		;;
	c)
		CHROMDIR=$OPTARG
		;;
	h)
		display_help
		;;
    \?)
      	echo "Invalid option: -$OPTARG" >&2
      	exit 1
    	  ;;
    :)
      	echo "Option -$OPTARG requires an argument." >&2
      	exit 1
      	;;
  esac
done

# Set defaults
if [ -z "$MYGENOME" ]
then
	MYGENOME="hg38"
fi

# Check vars
if [ -z "$MYBW" ]
then
	echo "Option -i is required!"
	exit 2
fi

if [ -z "$CHROMDIR" ]
then
	echo "Option -c is required!"
	exit 2
fi

if [ ! -f  $CHROMDIR/$MYGENOME"_chrom_sizes.2_column" ]
then
	echo "Genome $MYGENOME is not supported! (Add a chrom size file)"
	exit 3
fi

FILE_PATH="$MYBW"
FILE_OUTPUT=$(file "$FILE_PATH")

# Check filetype
if [ $(file $MYBW | grep -c text) != 0 ]
then
	# Get a random number to name the temp file
	MYRND=$RANDOM
	let "MYRND %= 9999"

	# Clean up the unwanted contigs
	grep -v 'chrUn\|random\|alt\|GL\|KI' $MYBW > tmp.$MYRND.bg

	# Remove the trackline from the tmp file
	sed -i '1d' tmp.$MYRND.bg

	# Sort the tmp file
	sort -k1,1 -k2,2n -k3,3n tmp.$MYRND.bg > tmp.sorted.$MYRND.bg

	# Make bigwig
	bedGraphToBigWig tmp.sorted.$MYRND.bg $CHROMDIR/$MYGENOME"_chrom_sizes.2_column" $MYBW

	# Clean up
	rm tmp.$MYRND.bg
	rm tmp.sorted.$MYRND.bg
else
	echo "Already not a bedGraph file! (well, not a text file at least... Is it already a bigwig?)"
	echo "Assuming we have a working HOMER installation, the file $FILE_PATH is already of type $FILE_OUTPUT. No need to do anything."
fi
