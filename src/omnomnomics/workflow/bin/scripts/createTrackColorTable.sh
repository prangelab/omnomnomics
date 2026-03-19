#! /bin/bash

#Given a base color, construct a collection of shades and tints for that color and store it in a 'color table' as used by make_ChIP_hubs.sh

#--------------------------------------------------------------------------------------------------------------
#Function definitions
#Usage function
function display_help {
	echo ""
	echo "Usage: $0  -c BASE_COLOR -2 SECOND_COLOR -3 THIRD_COLOR -f FRACTION -s SKEW -n NAME -P <DIR> -d DIRECTION -r -v -t -h"
	echo ""
	echo "	-c:			Base color to work shade (darken) or tint (lighten). {RRR,GGG,BBB} {0-255} or R color name."
	echo "				Default: random!"
	echo "	-2:			Second color to  shade (darken) or tint (lighten) to. {RRR,GGG,BBB} {0-255} or R color name."
	echo "				Used with type 'mix'. Default: random."
	echo "	-3:			Third color to  shade (darken) or tint (lighten) to. {RRR,GGG,BBB} {0-255} or R color name."
	echo "				Used with type 'mix3'. Default: random."
	echo "	-f:			Fraction to shade or tint by. {0-1}"
	echo "				Default: 0.25"
	echo "	-s:			Multiplier to skew the fraction with after each round. {>1}"
	echo "				Default: 1"
	echo "	-n:			Name given to the output table"
	echo "				Default: mycolor"
	echo "	-P <DIR>:			Path to a folder with color tables in which to save the output table."
	echo "				Default: ~/bin/color_data_for_hubs"
	echo "	-d:			Direction to manipulate in: 'shade', 'tint', 'both', 'mix', or 'mix3'"
	echo "				Default: tint."
	echo "	-r:			Reverse the order of the colors in the palette."
	echo ""
	echo "	-v:			Display the results"
	echo ""
	echo "	-t:			Trash the results!"
	echo "     				Only display the colors, do not keep the file!"
	echo ""
	echo "	-h:			Print this help message"
	echo ""
	echo "Direction"
	echo "	shade:	Shades color (-c) to black."
	echo "	tint:	Tints  color (-c) to white."
	echo "	both:	Make a gradient from black via color (-c) to white."
	echo "	mix:	Make a gradient between colors (-c) and (-2)."
	echo "	mix3:	Make a gradient from color (-c) via (-2) to (-3)."
	echo ""
	echo "Palette length"
	echo "	Fraction:	Determines the step size bewteen each output color."
	echo "			The smaller the fraction, the more steps will be needed to reach the end color."
	echo "	Skew:		Controls if and by how much the step size will be increased after each round of color picking."
	echo "			If skew is 1, step size will not increase."
	echo "			If skew is > 1, step size will increase and the end color will be reached faster." 
	echo "			Useful if the colors in the beginning of a palette are visually separated better than at the end (e.g. when reaching black the last few colors tend to be to dark to distinguish)."
	exit 1
}

#--------------------------------------------------------------------------------------------------------------
#Parse options

while getopts ":c:2:3:f:s:n:P:d:rvth" opt; do
  case $opt in
	c)
		THECOL=$OPTARG
		;;
	2)
		THECOL2=$OPTARG
		;;	
	3)
		THECOL3=$OPTARG
		;;	
	f)
		THEFRACTION=$OPTARG
		;;
	s)
		FRAC_SKEW=$OPTARG
		;;
	n)
		THENAME=$OPTARG	
		;;
	P)
		THEFOLDER=$OPTARG	
		;;
	d)
		THEDIRECTION=$OPTARG
		;;
	r)
		REVERSE=1
		;;	
	v)
		SHOW_TABLE=1
		;;
	t)
		TRASH_TABLE=1
		SHOW_TABLE=1
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

#Set defaults
if [ -z "$THECOL" ]
then
	THECOL="$(( RANDOM % 256 )),$(( RANDOM % 256 )),$(( RANDOM % 256 ))"
fi

if [ -z "$THEFRACTION" ]
then
	THEFRACTION=0.25
fi

if [ -z "$THENAME" ]
then
	THENAME="mycolor"
fi

if [ -z "$THEFOLDER" ]
then
	THEFOLDER="$OMNOM_HOME/bin/color_data_for_hubs"
	echo $THEFOLDER
fi

if [ ! -d "$THEFOLDER" ]
then
	echo "Color table folder (-P) does not exist!"
	echo "Exiting..."
	exit 1
fi

if [ -z "$THEDIRECTION" ]
then
	THEDIRECTION="tint"
fi

if [ -z "$REVERSE" ]
then
	REVERSE=0
fi

if [ -z "$SHOW_TABLE" ]
then
	SHOW_TABLE=0
fi

if [ -z "$TRASH_TABLE" ]
then
	TRASH_TABLE=0
fi

if [ -z "$FRAC_SKEW" ]
then
	FRAC_SKEW=1
fi

if [[ $FRAC_SKEW < 1 ]]
then
	FRAC_SKEW=1
fi

#Remove slash from the end of the color table path if necessary.
[[ $THEFOLDER == */ ]] && THEFOLDER=${THEFOLDER%?}	

#Check if color file exists already (do not overwrite)
if [ -f "$THEFOLDER/$THENAME.color.table" ]
then
	echo "File $THEFOLDER/$THENAME.color.table already exists! Please use a different name."
	exit 1
fi

#Check for valid color
#Set color table file path
R_COLOR_TABLE="$OMNOM_HOME/bin/R.color.table"

#If we have an RGB color, do nothing
if [[ !($THECOL =~ ^([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5])$) ]]
then
	#If not, check if we have a valid R color name
	COLEXISTS=$(grep -wc $THECOL $R_COLOR_TABLE)
	if [[ "$COLEXISTS" == 0 ]]
	then
		#If not, exit
		echo "$THECOL is not a valid RGB color or R color name!"
		exit 1
	else
		#Retrieve the RGB values for the R color name
		THECOL=$(grep -w $THECOL $R_COLOR_TABLE | cut -f2)
	fi
fi

#Separate the color components into vars
R=$(echo $THECOL | cut -f 1 -d ",")
G=$(echo $THECOL | cut -f 2 -d ",")
B=$(echo $THECOL | cut -f 3 -d ",")

#If we are mixing, make sure we have a (proper) second color
if [[ "$THEDIRECTION" =~ "mix" ]]
then
	if [ -z "$THECOL2" ]
	then
		THECOL2="$(( RANDOM % 256 )),$(( RANDOM % 256 )),$(( RANDOM % 256 ))"
	fi
	
	#If we have an RGB color, do nothing
	if [[ !($THECOL2 =~ ^([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5])$) ]]
	then
		#If not, check if we have a valid R color name
		COLEXISTS=$(grep -wc $THECOL2 $R_COLOR_TABLE)
		if [[ "$COLEXISTS" == 0 ]]
		then
			#If not, exit
			echo "$THECOL2 is not a valid RGB color or R color name!"
			exit 1
		else
			#Retrieve the RGB values for the R color name
			THECOL2=$(grep -w $THECOL2 $R_COLOR_TABLE | cut -f2)
		fi
	fi

	#Separate the color components into vars
	R2=$(echo $THECOL2 | cut -f 1 -d ",")
	G2=$(echo $THECOL2 | cut -f 2 -d ",")
	B2=$(echo $THECOL2 | cut -f 3 -d ",")
	
	#And if we are doing a 3 way mix, check the third color also
	if [ "$THEDIRECTION" == "mix3" ]
	then
		if [ -z "$THECOL3" ]
		then
			THECOL3="$(( RANDOM % 256 )),$(( RANDOM % 256 )),$(( RANDOM % 256 ))"
		fi
	
		#If we have an RGB color, do nothing
		if [[ !($THECOL3 =~ ^([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5]),([0-9]|[1-9][0-9]|1[0-9][0-9]|2[0-4][0-9]|25[0-5])$) ]]
		then
			#If not, check if we have a valid R color name
			COLEXISTS=$(grep -wc $THECOL3 $R_COLOR_TABLE)
			if [[ "$COLEXISTS" == 0 ]]
			then
				#If not, exit
				echo "$THECOL3 is not a valid RGB color or R color name!"
				exit 1
			else
				#Retrieve the RGB values for the R color name
				THECOL3=$(grep -w $THECOL3 $R_COLOR_TABLE | cut -f2)
			fi
		fi

		#Separate the color components into vars
		R3=$(echo $THECOL3 | cut -f 1 -d ",")
		G3=$(echo $THECOL3 | cut -f 2 -d ",")
		B3=$(echo $THECOL3 | cut -f 3 -d ",")
	
	fi
fi

#--------------------------------------------------------------------------------------------------------------
#Generate the color palette

#Tinting only
if [ "$THEDIRECTION" == "tint" ]
then	
	THECOL2="255,255,255"
	#Separate the color components into vars
	R=$(echo $THECOL | cut -f 1 -d ",")
	G=$(echo $THECOL | cut -f 2 -d ",")
	B=$(echo $THECOL | cut -f 3 -d ",")
	R2=$(echo $THECOL2 | cut -f 1 -d ",")
	G2=$(echo $THECOL2 | cut -f 2 -d ",")
	B2=$(echo $THECOL2 | cut -f 3 -d ",")
	THEDIRECTION="mix"
		
#Shading only	
elif [ "$THEDIRECTION" == "shade" ]
then
	THECOL2="0,0,0"
	#Separate the color components into vars
	R=$(echo $THECOL | cut -f 1 -d ",")
	G=$(echo $THECOL | cut -f 2 -d ",")
	B=$(echo $THECOL | cut -f 3 -d ",")
	R2=$(echo $THECOL2 | cut -f 1 -d ",")
	G2=$(echo $THECOL2 | cut -f 2 -d ",")
	B2=$(echo $THECOL2 | cut -f 3 -d ",")
	THEDIRECTION="mix"
	
#Tinting and shading	
elif [ "$THEDIRECTION" == "both" ]
then
	THECOL2=$THECOL
	THECOL="0,0,0"
	THECOL3="255,255,255"
	#Separate the color components into vars
	R=$(echo $THECOL | cut -f 1 -d ",")
	G=$(echo $THECOL | cut -f 2 -d ",")
	B=$(echo $THECOL | cut -f 3 -d ",")
	R2=$(echo $THECOL2 | cut -f 1 -d ",")
	G2=$(echo $THECOL2 | cut -f 2 -d ",")
	B2=$(echo $THECOL2 | cut -f 3 -d ",")
	R3=$(echo $THECOL3 | cut -f 1 -d ",")
	G3=$(echo $THECOL3 | cut -f 2 -d ",")
	B3=$(echo $THECOL3 | cut -f 3 -d ",")
	THEDIRECTION="mix3"	
fi


#Mix
if [ "$THEDIRECTION" == "mix" ]
then
	echo "Mixing $THECOL and $THECOL2"

	#As long as we are not fully color 2 yet
	while [[ ($R != $R2) || ($G != $G2) || ($B != $B2) ]]
	do
		#Write the color to the palette
		echo "$R,$G,$B" >> "$THEFOLDER/$THENAME.color.table"
		echo "Color: $R,$G,$B Color2: $R2,$G2,$B2 Fraction: $THEFRACTION"

		#And calculate the next step
		#First check if we are over, adjust, and then make sure we don't go under!
		if [[ "$R" -gt "$R2" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$R" -lt "$R2" && R=$R2
		fi
		if [[ "$G" -gt "$G2" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$G" -lt "$G2" && G=$G2
		fi
		if [[ "$B" -gt "$B2" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$B" -lt "$B2" && B=$B2
		fi
		
		#And vice versa
		if [[ "$R" -lt "$R2" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$R" -gt "$R2" && R=$R2
		fi
		if [[ "$G" -lt "$G2" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$G" -gt "$G2" && G=$G2
		fi
		if [[ "$B" -lt "$B2" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$B" -gt "$B2" && B=$B2
		fi
		
		#Add some skew to the fraction
		THEFRACTION=$(echo "scale=5;$THEFRACTION * $FRAC_SKEW" | bc)
	done
	#Add color 2 to the palette
	echo "$R2,$G2,$B2" >> "$THEFOLDER/$THENAME.color.table"
	echo "Color: $R2,$G2,$B2 Color2: $R2,$G2,$B2 Fraction: $THEFRACTION"
	
#Mix 3 colors
elif [ "$THEDIRECTION" == "mix3" ]
then
	echo "Mixing $THECOL, $THECOL2, and $THECOL3"
	#Save the frac
	OLDFRAC=$THEFRACTION
	
	#First tint to a temp file
	MYRND=$RANDOM
	let "MYRND %= 9999"
	
	#As long as we are not fully color 2 yet
	while [[ ($R != $R2) || ($G != $G2) || ($B != $B2) ]]
	do
		#Write the color to the palette
		echo "$R,$G,$B" >> "$THEFOLDER/$THENAME.color.table"
		echo "Color: $R,$G,$B Color2: $R2,$G2,$B2 Fraction: $THEFRACTION"

		#And calculate the next step
		#First check if we are over, adjust, and then make sure we don't go under!
		if [[ "$R" -gt "$R2" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$R" -lt "$R2" && R=$R2
		fi
		if [[ "$G" -gt "$G2" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$G" -lt "$G2" && G=$G2
		fi
		if [[ "$B" -gt "$B2" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$B" -lt "$B2" && B=$B2
		fi
		
		#And vice versa
		if [[ "$R" -lt "$R2" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$R" -gt "$R2" && R=$R2
		fi
		if [[ "$G" -lt "$G2" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$G" -gt "$G2" && G=$G2
		fi
		if [[ "$B" -lt "$B2" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$B" -gt "$B2" && B=$B2
		fi
		
		#Add some skew to the fraction
		THEFRACTION=$(echo "scale=5;$THEFRACTION * $FRAC_SKEW" | bc)
	done
	
	#Reset the color and fraction
	R=$R2
	G=$G2
	B=$B2
	THEFRACTION=$OLDFRAC
	
	#As long as we are not fully color 3 yet
	while [[ ($R != $R3) || ($G != $G3) || ($B != $B3) ]]
	do
		#Write the color to the palette
		echo "$R,$G,$B" >> "$THEFOLDER/$THENAME.color.table"
		echo "Color: $R,$G,$B Color3: $R3,$G3,$B3 Fraction: $THEFRACTION"

		#And calculate the next step
		#First check if we are over, adjust, and then make sure we don't go under!
		if [[ "$R" -gt "$R3" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$R" -lt "$R3" && R=$R3
		fi
		if [[ "$G" -gt "$G3" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$G" -lt "$G3" && G=$G3
		fi
		if [[ "$B" -gt "$B3" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 * (1 - frac)))}')
			test "$B" -lt "$B3" && B=$B3
		fi
		
		#And vice versa
		if [[ "$R" -lt "$R3" ]]
		then
			R=$(echo "$R" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$R" -gt "$R3" && R=$R3
		fi
		if [[ "$G" -lt "$G3" ]]
		then
			G=$(echo "$G" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$G" -gt "$G3" && G=$G3
		fi
		if [[ "$B" -lt "$B3" ]]
		then
			B=$(echo "$B" | awk -v frac=$THEFRACTION '{printf("%d\n", ($1 + ((255 - $1) * frac)) + 1)}')
			test "$B" -gt "$B3" && B=$B3
		fi
		
		#Add some skew to the fraction
		THEFRACTION=$(echo "scale=5;$THEFRACTION * $FRAC_SKEW" | bc)
	done
	
#Intercept invalid direction argument
else
	echo "Invalid direction $THEDIRECTION (-d): please use 'tint', 'shade', 'both', 'mix', or 'mix3'"
	exit 1
fi

#Reverse the output lines
if [[ $REVERSE == 1 ]]
then
	MYRND=$RANDOM
	let "MYRND %= 9999"
	tac "$THEFOLDER/$THENAME.color.table" > $MYRND.col.tmp
	mv $MYRND.col.tmp "$THEFOLDER/$THENAME.color.table"
fi

#Display color table
if [[ $SHOW_TABLE == 1 ]]
then
	echo "$(wc -l $THENAME.color.table | cut -f1 -d" ") colors in this table"
	awk 'BEGIN{FS=","}{
		r = $1; g = $2; b = $3;
		printf "\033[48;2;%d;%d;%dm", r,g,b;
		printf "\n";
	}' "$THEFOLDER/$THENAME.color.table"

	#Reset terminal colors
	tput sgr0 
fi

#Trash color table
if [[ $TRASH_TABLE == 1 ]]
then
	if [ -f "$THEFOLDER/$THENAME.color.table" ]
	then
		rm "$THEFOLDER/$THENAME.color.table"
	fi
fi