The files in this folder are used to validate the tensor elements
calculated in NievesQELCCPXSec.cxx against the tensor elements
calculated in Nieves' fortran code. In order to make validation
plots follow these steps:

1. Use Nieves' fortran code to output the necessary values. Either
(temporarily) copy the versions of the NievesQEL code in this folder into
$GENIE/src/LwlynSmith/, or uncomment the testing code in the files
already there. Edit the input values in the CompareNievesTensors, then
compile the code and run it. This will generate a file with data in columns.

2. Edit compile_run.sh and fort.4 to have the same input values that you
selected in GENIE. Then run compile_run.sh to create a file with the same
data in columns.

3. Use compare_tensors_choose_kine.sh to create plots comparing the tensor
elements, as well as the form factors, polarization constants, coulomb
factor, and various kinematics.
