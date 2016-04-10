#!/bin/bash

# Contains the following commands:
#
# mkspl_rungevgen - makes a new spline, then runs gevgen
# autoGNUPlot - plots contents of text files using gnuplot
# combine_spl - combines two splines
# gevgenScript - runs gevgen and gntpc
# mkspl - makes a new spline
# mkspl_all - calls mkspl 8 time for all combos of R, C, and RFG/LFG
# set_RPA_Coul_FG - sets config file for Nieves XSec method
# setVar - sets single parameter in config file for Nieves XSec method
# xmlGetData - puts xsec data from an xml file for CCQE into 2 columns
# ---------------------------------------------------------------------------
# Sample commands:
#
# Code to create a spline and run gevgen for CCQE, C12, RPA, Coul, LFG
# mkspl_rungevgen C12_R_C_LFG false true false true \
#     LFGNuclearModel FGMBodekRitchie CC 1000060120 C12 \
#     1 1000000 1.0
#
# Code to create and combine a CCQE spline with a CC spline:
# 
#
# ---------------------------------------------------------------------------
# Function to make a new spline and the run gevgen
# $1 - name of directory to create files in
# $2,$3 - Change RPA from $2 to $3 (true or false)
# $4,$5 - Change Coul from $4 to $5 (true or false)
# $6,$7 - Change Nuclear Model from $6 to $7 
#         (FGMBodekRitchie or LFGNuclearModel)
# $8 - type of spline to create (eg CC, CCQE)
# $9 - target PDG code, eg 1000060120
# $10 - target name, eg C12
# $11 - gevgen run number
# $12 - Number of events to generate
# $13 - Neutrino energy
# $14,$15 - change xsec model from $15 to $16
# $16 - Messenger settings (same default as gevgenScript)

#
# Sample call:
# mkspl_rungevgen noR_noC_RFG true false true false \
#     LFGNuclearModel FGMBodekRitchie CC 1000060120 C12 \
#     1 1000000 1.0

function mkspl_rungevgen() {
    mkspl $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${14} ${15} ${16}

    # Return to the new directory
    CWD=$(pwd)
    newdir=$1
    cd $newdir

    # Run gevgen
    echo -e 'running gevgen and gntpc' >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    messenger="$GENIE/config/Messenger_inuke_verbose.xml"
    if [ -n "${16}" ]
	then
	messenger=${16}
    fi
    set_RPA_Coul_FG $2 $3 $4 $5 $6 $7
    setVar 'QEL-CC' ${14} ${15} UserPhysicsOptions.xml
    gevgenScript ${10} $9 ${11} ${12} ${13} $8 $1'.xml' ${16}

    echo -e 'gevgen and gntpc finished' >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    echo -e 'mkspl_rungevgen.sh finished\n\n' >> "status_"$newdir".txt"
    cd $CWD
}

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
# Take a partial spline and combine it with a newly created spline
# eg, take a spline with all CC events except CCQE, and combine it
# with a complete CCQE spline
#
# $1 - first spline, without the interactions in the second spline
# $2 - second spline
# $3 - name of the combined spline
function combine_spl() {
    tail -n+6 $2 > temp.txt
    cat $1 temp.txt > $3
    rm temp.txt
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
#$6 event generator list (eg CC, CCQE)
#$7 spline
#$8 Messenger settings file (defaults to Messenger_inuke_verbose.xml)
gevgenScript() {
    if [ -n "$8" ]
	then
	messenger=$8
    else
	messenger="$GENIE/config/Messenger_inuke_verbose.xml"
    fi
    
    if [ -n $7 ]
	then
	gevgen -n $4 -p 14 -e $5 --event-generator-list $6 \
	    --message-thresholds $messenger --cross-sections $7 \
	    -t $2 -r $3
	
	gntpc -f gst -i gntp.$3.ghep.root \
	    -o "numu_"$1"_"$6"_run"$3".gst.root"
    else
	echo -e "spline "$7" does not exist.\nNot running gevgen or gntpc."
    fi

}
# ---------------------------------------------------------------------------
# Function to make a new spline
# $1 - name of directory to create files in
# $2,$3 - Change RPA from $2 to $3 (true or false)
# $4,$5 - Change Coul from $4 to $5
# $6,$7 - Change Suppression from $6 to $7
# $8 - type of spline to create (eg CC, CCQE)
# $9 - target PDG code, eg 1000060120
# $10 - target name, eg C12
# $11,$12 - change xsec model from $11 to $12
# $13 - Messenger config file (Messenger_inuke_verbose.xml by default);
#
# created spline is named $1.xml

function mkspl() {
    newdir=$1

    CWD=$(pwd)
    mkdir $newdir
    cd $newdir

    echo -e "top directory: "$CWD >> "status_"$newdir".txt"
    echo -e "working in: "$CWD/$newdir >> "status_"$newdir".txt"
    echo -e "running gmkspl" >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    messenger="$GENIE/config/Messenger.xml"
    if [ -n "$13" ]
	then
	messenger=${13}
    fi

    set_RPA_Coul_FG $2 $3 $4 $5 $6 $7
    setVar 'QEL-CC' ${11} ${12} UserPhysicsOptions.xml
    gmkspl -p 14 -t $9 --event-generator-list $8 \
	--message-thresholds $messenger
    mv xsec_splines.xml $1".xml"

    echo -e 'N LFG gmkspl finished for '$10'_'$8 >> "status_"$newdir".txt"
    echo -e 'Generating png graph' >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    xmlGetData $1".xml" "data_"$1".txt"
    autoGNUPlot "data_"$1".txt" $1".png" 1 2 \
        'Neutrino Energy GeV' 'Cross Section' $10'_'$8\
	0 40 0 1.7e-10 

    echo -e 'mkspl complete' >> "status_"$newdir".txt"
    date >> "status_"$newdir".txt"
    echo -e "\n" >> "status_"$newdir".txt"

    cd $CWD
}
# ---------------------------------------------------------------------------
# NOT UPDATED. CURRENTLY WILL NOT WORK.
# Generate all 8 splines possible from combinations of RPA, Coul, and Supp
function mkspl_all() {
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
# setVar changes the variable $1 from $2 to $3 in $GENIE/config/$4
setVar(){
    sed -i 's/\(.*\)'$1'\(.*\)'$2'\(.*\)/\1'$1'\2'$3'\3/' \
	$GENIE/config/$4
}
# ---------------------------------------------------------------------------
# Set RPA, then Coulomb, then Fermi Gas Model
function set_RPA_Coul_FG() {
    setVar RPA $1 $2 NievesQELCCPXSec.xml
    setVar Coulomb $3 $4 NievesQELCCPXSec.xml
    setVar NuclearModel $5 $6 UserPhysicsOptions.xml
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
