#! /bin/bash

##============================##
## run_quant_peaks.sh		  ##
## Copyright Koen Prange 2023 ##
## v0.2						  ##
##============================##

###---------------------------------------------------------------------------------------------------------------
### Wrapper script to dispatch a job to count reads under peaks using sbatch
###---------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#Function definitions
#Usage function
function display_help {
	echo ""
	echo "Usage: $0  -d <DIR> -p <BED> -g genome -r regex -h"
	echo "Takes a number of ChIP-seq (or ATAC-seq) experiments as input, in HOMER tag directory format. Also takes a peak file in BED format."
	echo "Counts reads under the peaks."
	echo ""
	echo "	-d <DIR>:	Directory containing HOMER tag directories"
	echo "				Required argument"
	echo ""
	echo "	-p <BED>:	BED format peak file of genomic locations to quantify"
	echo "				Required argument"
	echo ""
	echo "	-g	genome:	Genome to use. Check HOMER config for availability."
	echo "				Default: hg38"
	echo "	-r	regex:	SED command to clean up the sample names in the table header. If not defined, will attempt to clean up standard pipeline output names."
	echo "				Default: 's|_L00[0-9]_R1\.trimmed\.sorted\.dups_marked\.filtered\.HOMER_tagDir ([0-9]\+\.[0-9]\+ total)||g;s/_[ACTG]\+-//g'"
	echo ""
	echo "	-h:		Print this help message"
	echo ""
	exit 1
}
#--------------------------------------------------------------------------------------------------------------
#Display help and exit if not enough arguments

if [ $# -lt 1 ]
then
	display_help
fi

#--------------------------------------------------------------------------------------------------------------
#Parse options

while getopts ":d:p:g:r:h" opt; do
  case $opt in
  	d)
  		MYDIR=$OPTARG
  		;;
  	p)
  		MYPEAKS=$OPTARG
  		;;
  	g)
  		MYGENOME=$OPTARG
  		;;
  	r)
  		SUBRE=$OPTARG
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

#--------------------------------------------------------------------------------------------------------------
#Set defaults or die
if [ ! -d "$MYDIR" ]
then
	echo "No directory (-d) given! This is a required argument." >&2
	exit 1
fi

if [ -z "$MYPEAKS" ]
then
	echo "No peak file (-p) given! This is a required argument." >&2
	exit 1
fi

if [ ! -f "$MYPEAKS" ]
then
	echo "Peak file (-p) Does not exist!" >&2
	exit 1
fi

if [[ $MYPEAKS != *.bed ]]
then
	echo "Peak file (-p) has to be in BED format!" >&2
	exit 1
fi

if [ -z "$MYGENOME" ]
then
	MYGENOME="hg38"
fi

if [ -z "$SUBRE" ]
then
	SUBRE="Default"
fi

#---------------------------------------------------------------------------------------------------------------
#Set variables
#Remove slash from the end of the directory if necessary.
[[ $MYDIR == */ ]] && MYDIR=${MYDIR%?}

#Set the number of cores per node
MYCORES=16

#Set job
THEJOB="$OMNOM_HOME/jobs/quant_peaks.job"

#Set extension
MYEXT="tagDir"
if [[ $(ls $MYDIR | grep -c "$MYEXT$") == 0 ]]
then
	MYEXT="tagDir.tar.gz"
    if [[ $(ls $MYDIR | grep -c "$MYEXT$") == 0 ]]
	then
		echo "No HOMER tag dirs found in $MYDIR! Exiting..."
		exit 1
	else
		echo "Compressed tagDirs found..."
	fi
else
	echo "Uncompressed tagDirs found..."
fi

#---------------------------------------------------------------------------------------------------------------
#Distribute files over jobs and submit them

# Get number of files of type MYEXT in MYDIR, and put their names into an array
NUMFILES=$(ls -d $MYDIR/*$MYEXT | wc -l | sed 's/^\W\+//' | cut -f1)

# Report some stats and create the batch job
echo "Running quant peaks"
echo "FILES:	$NUMFILES"
echo "DIR:	$MYDIR"
echo "SED:	$SUBRE"
echo ""

sbatch --export=MYDIR=$MYDIR,MYPEAKS=$MYPEAKS,MYGENOME=$MYGENOME,SUBRE=$SUBRE $THEJOB