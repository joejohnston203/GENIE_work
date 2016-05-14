#!/bin/bash

# ---------------------------------------------------------------------------
# Function to take the output of the modified qe.f code and modified GENIE
# code, which both contain the values of Q2, form factors, tensors, and more
# from each iteration.
#
# Arguments
# 1 modified GENIE code output file Furmanski
# 2 modified GENIE code output file Nieves
# 3 prefix for saved files
function compare_xsec() {
    files="$1 $2"
    # Add a line for each stored variable to make a plot
    autoGNUPlot "$files" 'q2_xsec.png' 1 2 'Q2' 'xsec' 'Q2 vs xsec'
    #autoGNUPlot "$files" 'q2_dq.png' 1 3 'Q2' 'dq' 'Q2 vs dq'
    #autoGNUPlot "$files" 'q2_xsec.png' 1 4 'Q2' 'xsec' 'Q2 vs xsec'
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