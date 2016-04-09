#!/bin/bash

# Script to run qe-gen-integral-dOmega.f and generate data files 
CWD=$(pwd)

target='Pb208'
for RPA in 1 0; do
for Coul in 1 0; do
for Enu in 0.2 1.0 5.0; do
    mkdir $target"_"$Enu"_R"$RPA"_C"$Coul
    cp fort.4 $target"_"$Enu"_R"$RPA"_C"$Coul/fort.4
    cd $target"_"$Enu"_R"$RPA"_C"$Coul
    # Integrate over all angles
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" -1.0 1.0 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_all_angles.txt"
    # Integrate over small subsets
    
    cd $CWD
done
done
done
