#!/bin/bash

Enu='1'
cosl='0.5'
RPA='1'
rfrac='0'

gfortran qe_iterate_kine.f -o qe_iterate_kine.o
echo -e $Enu "\n" $cosl "\n" $RPA "\n" $rfrac "\n" | ./qe_iterate_kine.o
mv fort.61 'fort.RPA_E'$RPA'_ctl'$cosl'_r'$rfrac
