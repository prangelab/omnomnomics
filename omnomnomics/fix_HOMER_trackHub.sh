#! /bin/bash



#--------------------------------------------------------------------------------------------------------------
# Function definitions
# Usage function
function display_help {
	echo ""
	echo "Usage: $0 -i <TRACKHUB_DIR> -g GENOME -h"
	echo ""
	echo "	-i <DIR>:		Input trackhub"
	echo "				Required argument."
	echo "	-g <DIR>:		Genome"
	echo "				Default: hg38"
	echo ""
	echo "	-h:			Print this help message"
	echo ""
	exit 1
}

#--------------------------------------------------------------------------------------------------------------
# Display help and exit if not enough arguments

if [ $#  -lt 1 ]
then
	display_help
fi

#--------------------------------------------------------------------------------------------------------------
# Define options

while getopts ":i:g:h" opt; do
  case $opt in
	i)
		MYHUB=$OPTARG
		;;
	g)
		MYGENOME=$OPTARG
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
if [ -z "$MYHUB" ]
then
	echo "Option -i is required!"
	exit 2
fi

if [ ! -f  $OMNOM_HOME/genomes/$MYGENOME"_chrom_sizes.2_column" ]
then
	echo "Genome $MYGENOME is not supported! (Add a chrom size file)"
	exit 3
fi

FILE_PATH="$MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirpos.ucsc.bigWig"
FILE_OUTPUT=$(file "$FILE_PATH")
# Check filetype
if [ $(file $MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirpos.ucsc.bigWig | grep -c text) != 0 ]
then
	# Get a random number to name the temp file
	MYRND=$RANDOM
	let "MYRND %= 9999"

	# Clean up the unwanted contigs
	grep -v 'chrUn\|random\|alt\|GL\|KI' $MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirpos.ucsc.bigWig > pos_tmp.$MYRND.bg &
	grep -v 'chrUn\|random\|alt\|GL\|KI' $MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirneg.ucsc.bigWig > neg_tmp.$MYRND.bg &
	wait

	# Remove the trackline from the tmp files
	sed -i '1d' pos_tmp.$MYRND.bg
	sed -i '1d' neg_tmp.$MYRND.bg

	# Sort the tmp files
	sort -k1,1 -k2,2n -k3,3n pos_tmp.$MYRND.bg > pos_tmp.sorted.$MYRND.bg
	sort -k1,1 -k2,2n -k3,3n neg_tmp.$MYRND.bg > neg_tmp.sorted.$MYRND.bg

	# Make bigwigs
	bedGraphToBigWig pos_tmp.sorted.$MYRND.bg $OMNOM_HOME/genomes/$MYGENOME"_chrom_sizes.2_column" $MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirpos.ucsc.bigWig &
	bedGraphToBigWig neg_tmp.sorted.$MYRND.bg $OMNOM_HOME/genomes/$MYGENOME"_chrom_sizes.2_column" $MYHUB/$MYGENOME/$(basename $MYHUB .hub).HOMER_tagDirneg.ucsc.bigWig &
	wait

	# Set vars as we like them
	sed -i 's/maxHeightPixels [^ ]\+/maxHeightPixels 128:64:8/' $MYHUB/$MYGENOME/trackDb.txt
	sed -i 's/color 255,150,150/color 255,128,0/' $MYHUB/$MYGENOME/trackDb.txt
	sed -i 's/color 255,180,180/color 0,128,255/' $MYHUB/$MYGENOME/trackDb.txt

	# Clean up
	rm pos_tmp.$MYRND.bg &
	rm neg_tmp.$MYRND.bg &
	rm pos_tmp.sorted.$MYRND.bg &
	rm neg_tmp.sorted.$MYRND.bg &
	wait

else
	echo "Already not a bedGraph file! (well, not a text file at least... Is it already a bigwig?)"
  	echo "Assuming we have a working HOMER installation, the file $FILE_PATH is already of type $FILE_OUTPUT. No need to do anything."
fi

