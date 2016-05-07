#!/bin/bash

gfortran qe_iterate_kine.f -o qe_iterate_kine.o
echo -e 1 "\n" 0 "\n" 1 "\n" 0 "\n" | ./qe_iterate_kine.o
mv fort.61 fort.RPA_E1_ctl0_r0
