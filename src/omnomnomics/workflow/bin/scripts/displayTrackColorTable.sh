#! /bin/bash

# Visualize a 'color table' file as made by createTrackcolorTable.sh and used by make_ChIP_hubs.sh.

#--------------------------------------------------------------------------------------------------------------
# Function definitions
# Usage function
function display_help {
	echo ""
	echo "Usage: $0  -i <COLOR_TABLE> -h"
	echo ""
	echo "	-i:			Color table file to display"
	echo "				Mandatory argument"
	echo "	-t:			Display color test"
	echo "				Display a spectrum to see if your terminal supports true color"
	echo "	-h:			Print this help message"
	echo ""
	echo "Examples:	~/bin/color_data_for_hubs"
	exit 1
}

#--------------------------------------------------------------------------------------------------------------
# Parse options

while getopts ":i:th" opt; do
  case $opt in
	i)
		THEFILE=$OPTARG
		;;
	t)
		THETEST=1
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

# Run test and exit, if required
if [ ! -z "$THETEST" ]
then
	echo "Color Test: Following line should be a rainbow. If not your terminal does not support true color."
	COLUMNS=$(tput cols)
	awk -v ncols=$COLUMNS 'BEGIN{
		s="Full spectrum // "; s=s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s s;
		for (colnum = 0; colnum<ncols; colnum++) {
			r = 255-(colnum*255/(ncols-1));
			g = (colnum*510/(ncols-1));
			b = (colnum*255/(ncols-1));
			if (g>255) g = 510-g;
			printf "\033[48;2;%d;%d;%dm", r,g,b;
			printf "\033[38;2;%d;%d;%dm", 255-r,255-g,255-b;
			printf "%s\033[0m", substr(s,colnum+1,1);
		}
		printf "\n";
	}'
	exit 1
fi


# Set defaults
if [ ! -f "$THEFILE" ]
then
	echo "No color table file (-i) specified!"
	exit 1
fi


#--------------------------------------------------------------------------------------------------------------
# Display color table
echo "Displaying color table:"
echo "$THEFILE"
NUMCOLS=$(wc -l $THEFILE | cut -f 1 -d " ")
COLUMNS=$(tput cols)
echo "Number of colors: $NUMCOLS"


# Display color table
awk 'BEGIN{FS=","}{
	r = $1; g = $2; b = $3;
	printf "%d %d %d", r,g,b
	printf "\033[48;2;%d;%d;%dm", r,g,b;
    printf "\n";
}' $THEFILE

# Reset terminal colors
tput sgr0 