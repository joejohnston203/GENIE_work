#!/bin/bash

. /local/xrootd/a/users/genie/GENIEtrunk290/setup_genie
. ~/scripts/core.sh

prefix='2016_05_20'
target_pdg='1000060120'
target='C12'
numev='100000'
neutrinoE='0.2'

# LwlynSmith, RFG
#mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_LS_RFG' true false true false LocalFGM FGMBodekRitchie CCQE $target_pdg $target 0 $numev $neutrinoE Nieves LwlynSmith 

# LwlynSmith, LFG
mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_LS_LFG' true false true false FGMBodekRitchie LocalFGM CCQE $target_pdg $target 0 $numev $neutrinoE Nieves LwlynSmith 

# Nieves, RFG, No RPA, No Coul
#mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_N_noRPA_noCoul_RFG' true false true false LocalFGM FGMBodekRitchie CCQE $target_pdg $target 0 $numev $neutrinoE LwlynSmith Nieves

# Nieves, LFG, No RPA, No Coul
mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_N_noRPA_noCoul_LFG' true false true false FGMBodekRitchie LocalFGM CCQE $target_pdg $target 0 $numev $neutrinoE LwlynSmith Nieves

# Nieves, RFG, RPA, Coul
#mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_N_RPA_Coul_RFG' false true false true LocalFGM FGMBodekRitchie CCQE $target_pdg $target 0 $numev $neutrinoE LwlynSmith Nieves

# Nieves, LFG, RPA, Coul
mkspl_rungevgen $prefix'_'$target'_'$neutrinoE'GeV_N_RPA_Coul_LFG' false true false true FGMBodekRitchie LocalFGM CCQE $target_pdg $target 0 $numev $neutrinoE LwlynSmith Nieves

