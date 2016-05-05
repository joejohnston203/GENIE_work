#!/bin/bash

# ---------------------------------------------------------------------------
# Function to take the output of the modified qe.f code and modified GENIE
# code, which both contain the values of Q2, form factors, tensors, and more
# from each iteration.
#
# Arguments
# 1 modified qe.f output file (out of qe-gen-intOmega-tensors.o)
# 2 modified GENIE code output file
# 3 prefix for saved files
function comp_tensors_choose_kine() {
    files="$1 $2"
    # Add a line for each stored variable to make a plot
    autoGNUPlot "$files" 'kine_'$3'_q2_q0.png' 1 2 'Q2' 'q0' 'Q2 vs q0'
    autoGNUPlot "$files" 'kine_'$3'_q2_dq.png' 1 3 'Q2' 'dq' 'Q2 vs dq'
    autoGNUPlot "$files" 'Amunu_'$3'_q2_axx.png' 1 4 'Q2' 'axx' 'Q2 vs axx'
    autoGNUPlot "$files" 'Amunu_'$3'_q2_azz.png' 1 5 'Q2' 'azz' 'Q2 vs azz'
    autoGNUPlot "$files" 'Amunu_'$3'_q2_a0z.png' 1 6 'Q2' 'a0z' 'Q2 vs a0z'
    autoGNUPlot "$files" 'Amunu_'$3'_q2_a00.png' 1 7 'Q2' 'a00' 'Q2 vs a00'
    autoGNUPlot "$files" 'Amunu_'$3'_q2_axy.png' 1 8 'Q2' 'axy' 'Q2 vs axy'
    autoGNUPlot "$files" 'pc_'$3'_q2_CT.png' 1 9 'Q2' 'CT' 'Q2 vs CT' # fact
    autoGNUPlot "$files" 'pc_'$3'_q2_CL.png' 1 10 'Q2' 'CL' 'Q2 vs CL' # facl
    autoGNUPlot "$files" 'pc_'$3'_q2_CN.png' 1 11 'Q2' 'CN' 'Q2 vs CN' # f00
    autoGNUPlot "$files" 'kine_'$3'_q2_Tl.png' 1 12 'Q2' 'Tl' 'Q2 vs Tl'
    autoGNUPlot "$files" 'kine_'$3'_q2_xlind.png' 1 13 'Q2' 'xlind' 'Q2 vs xlind'
    autoGNUPlot "$files" 'ff_'$3'_q2_f1v.png' 1 14 'Q2' 'f1v' 'Q2 vs f1v'   
    autoGNUPlot "$files" 'ff_'$3'_xmuf2v.png' 1 15 'Q2' 'xmuf2v' 'Q2 vs xmuf2v'
    autoGNUPlot "$files" 'ff_'$3'_gaq.png' 1 16 'Q2' 'gaq' 'Q2 vs gaq'
    autoGNUPlot "$files" 'ff_'$3'_gpq.png' 1 17 'Q2' 'gpq' 'Q2 vs gpq'

    # Make plots vs tmugev
    #autoGNUPlot "$files" 'kine_'$3'_Tl_q0.png' 12 2 'Tmu' 'q0' 'Tmu vs q0'
    #autoGNUPlot "$files" 'kine_'$3'_Tl_dq.png' 12 3 'Tmu' 'dq' 'Tmu vs dq'
    #autoGNUPlot "$files" 'Amunu_'$3'_Tl_axx.png' 12 4 'Tmu' 'axx' 'Tmu vs axx'
    #autoGNUPlot "$files" 'Amunu_'$3'_Tl_azz.png' 12 5 'Tmu' 'azz' 'Tmu vs azz'
    #autoGNUPlot "$files" 'Amunu_'$3'_Tl_a0z.png' 12 6 'Tmu' 'a0z' 'Tmu vs a0z'
    #autoGNUPlot "$files" 'Amunu_'$3'_Tl_a00.png' 12 7 'Tmu' 'a00' 'Tmu vs a00'
    #autoGNUPlot "$files" 'Amunu_'$3'_Tl_axy.png' 12 8 'Tmu' 'axy' 'Tmu vs axy'
    #autoGNUPlot "$files" 'pc_'$3'_Tl_CT.png' 12 9 'Tmu' 'CT' 'Tmu vs CT' # fact
    #autoGNUPlot "$files" 'pc_'$3'_Tl_CL.png' 12 10 'Tmu' 'CL' 'Tmu vs CL' # facl
    #autoGNUPlot "$files" 'pc_'$3'_Tl_CN.png' 12 11 'Tmu' 'CN' 'Tmu vs CN' # f00
    #autoGNUPlot "$files" 'pc_'$3'_Tl_q2.png' 12 1 'Tl' 'Q2' 'Tl vs Q2'
    #autoGNUPlot "$files" 'pc_'$3'_Tl_xlind.png' 12 13 'Tl' 'xlind' 'Tl vs xlind'
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
