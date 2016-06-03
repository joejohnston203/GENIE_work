In order to create the plots in this directory showing that without
RPA, Nieves fortran code and my genie code give the same results,
I had to change several things:

1. Several constants in Nieves' code were changed slightly to match
   genie

2. genie must use LFG as the nuclear model (or else the densities will
   not have the proper r-dependence).

3. Genie must use DipoleELFormFactors. For some reason this class in
   genie sets the GEn form factor to 0, so I had to change GEn to
   0 inside Nieves' fortran code.
 
