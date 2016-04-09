#!/bin/bash

# Script to run qe-gen-integral-dOmega.f and generate data files 
CWD=$(pwd)

target='C12'
for RPA in 1 0; do
for Coul in 1 0; do
for Enu in 0.2 1 5; do
    mkdir $target"_"$Enu"_R"$RPA"_C"$Coul
    cp fort.4 $target"_"$Enu"_R"$RPA"_C"$Coul/fort.4
    cd $target"_"$Enu"_R"$RPA"_C"$Coul
    # Integrate from -1.0 to 1.0
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" -1.0 1.0 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_-1_1.txt"
    # Integrate from -1.0 to -0.9
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" -1.0 -0.9 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_-1_-0.9.txt"
    # Integrate from -0.05 to 0.05
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" -0.05 0.05 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_-0.05_0.05.txt"
    # Integrate from 0.4 to 0.5
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" 0.4 0.5 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_0.4_0.5.txt"
    # Integrate from 0.7 to 0.8
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" 0.7 0.8 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_0.7_0.8.txt"
    # Integrate from 0.9 to 0.95
    echo -e $Enu "\n" 1.0, $RPA, $Coul "\n" -1.0 1.0 "\n" | ../qe-gen-integral-dOmega.o
    mv fort.60 $target"_"$Enu"_R"$RPA"_C"$Coul"_0.9_0.95.txt"
    
    # return to the working directory
    cd $CWD
done
done
done
