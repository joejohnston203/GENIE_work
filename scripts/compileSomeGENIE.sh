#!bin/bash

source /usr/GENIEdevel/setup_genie

CWD=$(pwd)

echo -e "LlewelynSmith:\n"
cd /usr/GENIEdevel/src/LlewellynSmith/
make

echo -e "\n\nQEL:\n"
cd /usr/GENIEdevel/src/QEL/
make

echo -e "\n\nCrossSections:\n"
cd /usr/GENIEdevel/src/CrossSections/
make

echo -e "\n\nNuclear:\n"
cd /usr/GENIEdevel/src/Nuclear/
make

cd $CWD


#echo -e "\n\nEVGModules:\n"
#cd /usr/GENIEdevel/src/EVGModules/
#make


#echo -e "\n\nUtils:\n"
#cd /usr/GENIEdevel/src/Utils/
#make
