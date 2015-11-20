#!/bin/bash

# Contains the following commands:
#
# autoGNUPlot - plots contents of text files using gnuplot
# gevgenScript - runs gevgen and gntpc
# mkspl - makes a new spline
# mkspl_all - calls mkspl 8 time for all combos of R, C, and NS
# mkspl_rungevgen - makes a new spline, then runs gevgen
# set_RPA_Coul_Supp - sets config file for Nieves XSec method
# setVar - sets single parameter in config file for Nieves XSec method
# xmlGetData - puts xsec data from an xml file for CCQE into 2 columns

# ---------------------------------------------------------------------------

# Function to script creation of a png file using gnuplot
#
# Arguments
# 1 File to plot
# 2 output file name
# 3 first column
# 4 second column
# 5 xlabel
# 6 ylabel
# 7 plot title
# 8 xmin
# 9 xmax
# 10 ymin
# 11 ymax
# put nothing to use defaults for arguments 5 and up
autoGNUPlot() {
    out=$2
    if [ -z "$2" ]
	then
	out="out.png"
    fi
    echo "Storing graph of data in $1 to png file $out"

    # Assemble instructions to send to gnuplot
    instruc=""
    if [ -n "$5" ]
	then
	instruc="$instruc set xlabel \"$5\";"
    fi
    if [ -n "$6" ]
	then
	instruc="$instruc set ylabel \"$6\";"
    fi
    if [ -n "$7" ]
	then
	instruc="$instruc set title \"$7\";"
    fi
    
    # Commands to set graph axis bound
    if [[ -n "$8" && -n "$9" ]]
	then
	instruc="set xrange[$8:$9];"
	echo "xmin=$8, xmax=$9"
    else
	echo "Using default domain"
    fi
    
    if [[ -n "${10}" && -n "${11}" ]]
	then
	instruc="$instruc set yrange[${10}:${11}];"
	echo "ymin=${10}, ymax=${11}"
    else
	echo "Using default range"
    fi
    
    # Commands to plot
    instruc="$instruc set term png; set output \"$out\";"
    instruc="$instruc plot"
    for file in $1; do
	instruc="$instruc \"$file\" using $3:$4,"
    done
    instruc=${instruc%?}
    instruc="$instruc;"
#echo $instruc
    gnuplot <<EOF
$instruc
EOF
}

# ---------------------------------------------------------------------------

#Runs gevgen and converts the output to a gst.root file. 
#
#Arguments:
#$1 Target, eg C12 or Pb208
#$2 Target PDG code, eg 1000060120
#$3 Run number
#$4 Number of events to generate
#$5 Neutrino energy
#$6 spline
#$7 Messenger settings file (defaults to Messenger_inuke_verbose.xml)
gevgenScript() {
    if [ -n "$7" ]
	then
	messenger=$7
    else
	messenger="/usr/GENIEdevel/config/Messenger_inuke_verbose.xml"
    fi
    
    if [ -n $6 ]
	then
	gevgen -n $4 -p 14 -e $5 --event-generator-list CC \
	    --message-thresholds $messenger --cross-sections $6 \
	    -t $2 -r $3
	
	gntpc -f gst -i gntp.$3.ghep.root -o "numu_"$1"_run"$3"_noFSI.gst.root"
    else
	echo -e "spline "$6" does not exist.\nNot running gevgen or gntpc."
    fi

}

# ---------------------------------------------------------------------------

# Function to make a new spline
# $1 - name of directory to create files in
# $2, $3 - Change RPA from $2 to $3 (true or false)
# $4,$5 - Change Coul from $4 to $5
# $6,$7 - Change Suppression from $6 to $7
# $8 - Name of partial xml spline that does not include the QEL spline
# $9 - Messenger config file (Messenger_inuke_verbose.xml by default);
function mkspl() {
    newdir=$1
    partialSpline=$8

    CWD=$(pwd)
    mkdir $newdir
    cd $newdir

    echo -e "top directory: "$CWD >> "status_"$newdir".txt"
    echo -e "working in: "$CWD/$newdir >> "status_"$newdir".txt"
    echo -e "running gmkspl" >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    messenger="/usr/GENIEdevel/config/Messenger_inuke_verbose.xml"
    if [ -n "$9" ]
	then
	messenger=$9
    fi

    set_RPA_Coul_Supp $2 $3 $4 $5 $6 $7
    gmkspl -p 14 -t 1000060120 --event-generator-list CCQE \
	--message-thresholds $messenger
    mv xsec_splines.xml $newdir"_only.xml"

    echo -e 'N LFG gmkspl finished for numu_C12' >> "status_"$newdir".txt"
    echo -e 'Generating graph and full spline' >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    xmlGetData $newdir"_only.xml" "data_"$newdir".txt"
    autoGNUPlot "data_"$newdir".txt" "output_"$newdir".png" 1 2 \
        'Neutrino Energy GeV' 'Cross Section' 'C12 Spline'\
	0 40 0 1.7e-10 

    tail -n+6 $newdir"_only.xml" > temp.txt
    cat $partialSpline temp.txt > $newdir"_spline.xml"
    rm temp.txt

    date >> "status_"$newdir".txt"
    echo -e 'xmlGetData, plot, and full spline finished for numu_C12' \
	>> "status_"$newdir".txt"

    cd $CWD
}

# ---------------------------------------------------------------------------

# Generate all 8 splines possible from combinations of RPA, Coul, and Supp
mkspl_all() {
    mkspl noR_noC_S_spline true false true false false true \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl R_noC_S_spline false true true false false true \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl R_C_S_spline false true false true false true \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl noR_C_S_spline true false false true false true \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl noR_noC_noS_spline true false true false true false \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl R_noC_noS_spline false true true false true false \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl R_C_noS_spline false true false true true false \
	/data/jpj13/scripts/important_files/partial_spline.xml
    sleep 5

    mkspl noR_C_noS_spline true false false true true false \
	/data/jpj13/scripts/important_files/partial_spline.xml
}

# ---------------------------------------------------------------------------

# Function to make a new spline and the run gevgen
# $1 - name of directory to create files in
# $2, $3 - Change RPA from $2 to $3 (true or false)
# $4,$5 - Change Coul from $4 to $5
# $6,$7 - Change Suppression from $6 to $7
# $8 - Name of partial xml spline that does not include the QEL spline
# $9 - gevgen run number
# $10 - Number of events to generate
# $11 - Neutrino energy
# $12 - Messenger settings (Messenger_inuke_verbose.xml by default)
#
# Sample call:
# mkspl_rungevgen noR_noC_S true false true false false true \
#     /data/jpj13/scripts/partial_spline.xml \
#     1 1000000 1.0 &
function mkspl_rungevgen() {
    mkspl $1 $2 $3 $4 $5 $6 $7 $8 ${12}

    # Return to the new directory
    CWD=$(pwd)
    newdir=$1
    cd $newdir

    # Run gevgen
    echo -e 'running gevgen and gntpc' >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    messenger="/usr/GENIEdevel/config/Messenger_inuke_verbose.xml"
    if [ -n "${12}" ]
	then
	messenger=${12}
    fi
    set_RPA_Coul_Supp $2 $3 $4 $5 $6 $7
    gevgenScript C12 1000060120 $9 ${10} ${11} $newdir"_spline.xml"

    echo -e 'gevgen and gntpc finished' >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    echo -e 'mkspl_rungevgen.sh finished\n\n' >> "status_"$newdir".txt"
    cd $CWD
}

# ---------------------------------------------------------------------------

# Set var changes the variable $1 from $2 to $3 in Nieves XSec config file
setVar(){
    sed -i 's/\(.*\)'$1'\(.*\)'$2'\(.*\)/\1'$1'\2'$3'\3/' \
	/usr/GENIEdevel/config/NievesQELCCPXSec.xml
}

# ---------------------------------------------------------------------------

# Set RPA, then Coulomb, then Suppression
function set_RPA_Coul_Supp() {
    setVar RPA $1 $2
    setVar Coulomb $3 $4
    setVar Suppression $5 $6
}

# ---------------------------------------------------------------------------

# Gets data from an xml file created by gmkspl. Assumes only one interaction
# is included in the xml file, and thus gets all numbers, with E in the first
# column and xsec in the second. First line contains # followed by the number 
# of knots. # Second line is #E xsec to label columns. Remaining lines are 
# data
#
# Arguments
# $1 xml file to get data from
# $2 name of output file

xmlGetData() {
    out=$2
    if [ -z $2 ]
	then
	out="data.txt"
    fi
    sed -n 's/.*nknots="\([0-9]*\).*/#\1/p' $1 > $out
    echo -e "#E\txsec" >> $out
    sed -n 's/\t<knot> <E> *\([0-9.]*\) <\/E> <xsec> *\([0-9.e-]*\).*/\1\t\2/p'\
	$1 >> $out
}


# ---------------------------------------------------------------------------