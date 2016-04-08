#!/bin/bash

. ~/path_setup.sh

function copy_LFG_gentr_to_home {
    cp $gentr/src/CrossSections/QELXSec.h ~/LFG_Nieves_files/LFG/src/CrossSections/
    cp $gentr/src/CrossSections/QELXSec.cxx ~/LFG_Nieves_files/LFG/src/CrossSections/
    cp $gentr/src/EVGModules/FermiMover.h ~/LFG_Nieves_files/LFG/src/EVGModules/
    cp $gentr/src/EVGModules/FermiMover.cxx ~/LFG_Nieves_files/LFG/src/EVGModules/
    cp $gentr/src/EVGModules/PauliBlocker.h ~/LFG_Nieves_files/LFG/src/EVGModules/
    cp $gentr/src/EVGModules/PauliBlocker.cxx ~/LFG_Nieves_files/LFG/src/EVGModules/
    cp $gentr/src/Interaction/Target.h ~/LFG_Nieves_files/LFG/src/Interaction/
    cp $gentr/src/Interaction/Target.cxx ~/LFG_Nieves_files/LFG/src/Interaction/
    cp $gentr/src/LlewellynSmith/LwlynSmithQELCCPXSec.cxx ~/LFG_Nieves_files/LFG/src/LlewellynSmith/
    cp $gentr/src/LlewellynSmith/LwlynSmithQELCCPXSec.h ~/LFG_Nieves_files/LFG/src/LlewellynSmith/
    cp $gentr/src/Nuclear/LinkDef.h ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/LocalFGM.cxx ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/LocalFGM.h ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/NuclearModel.h ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/NuclearModelI.h ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/NuclearModelMap.h ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/Nuclear/NuclearModelMap.cxx ~/LFG_Nieves_files/LFG/src/Nuclear/
    cp $gentr/src/QEL/QELKinematicsGenerator.cxx ~/LFG_Nieves_files/LFG/src/QEL/
    cp $gentr/src/Utils/NuclearUtils.cxx ~/LFG_Nieves_files/LFG/src/Utils/
    cp $gentr/config/master_config.xml ~/LFG_Nieves_files/LFG/config/
    cp $gentr/config/LFGNuclearModel.xml ~/LFG_Nieves_files/LFG/config/
    cp $gentr/config/Messenger* ~/LFG_Nieves_files/LFG/config/
}

function copy_LFG_home_to_gentr {
    cp ~/LFG_Nieves_files/LFG/src/CrossSections/QELXSec.h $gentr/src/CrossSections/
    cp ~/LFG_Nieves_files/LFG/src/CrossSections/QELXSec.cxx $gentr/src/CrossSections/
    cp ~/LFG_Nieves_files/LFG/src/EVGModules/FermiMover.h $gentr/src/EVGModules/
    cp ~/LFG_Nieves_files/LFG/src/EVGModules/FermiMover.cxx $gentr/src/EVGModules/
    cp ~/LFG_Nieves_files/LFG/src/EVGModules/PauliBlocker.h $gentr/src/EVGModules/
    cp ~/LFG_Nieves_files/LFG/src/EVGModules/PauliBlocker.cxx $gentr/src/EVGModules/
    cp ~/LFG_Nieves_files/LFG/src/Interaction/Target.h $gentr/src/Interaction/
    cp ~/LFG_Nieves_files/LFG/src/Interaction/Target.cxx $gentr/src/Interaction/
    cp ~/LFG_Nieves_files/LFG/src/LlewellynSmith/LwlynSmithQELCCPXSec.cxx $gentr/src/LlewellynSmith/
    cp ~/LFG_Nieves_files/LFG/src/LlewellynSmith/LwlynSmithQELCCPXSec.h $gentr/src/LlewellynSmith/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/LinkDef.h $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/LocalFGM.cxx $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/LocalFGM.h $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/NuclearModel.h $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/NuclearModelI.h $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/NuclearModelMap.h $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/Nuclear/NuclearModelMap.cxx $gentr/src/Nuclear/
    cp ~/LFG_Nieves_files/LFG/src/QEL/QELKinematicsGenerator.cxx $gentr/src/QEL/
    cp ~/LFG_Nieves_files/LFG/src/Utils/NuclearUtils.cxx $gentr/src/Utils/
    cp ~/LFG_Nieves_files/LFG/config/master_config.xml $gentr/config/
    cp ~/LFG_Nieves_files/LFG/config/LFGNuclearModel.xml $gentr/config/
    cp ~/LFG_Nieves_files/LFG/config/Messenger* $gentr/config/
}
