# orca_led
                  ************************************************************** 
                  ***                                                        *** 
                  ***                        ORCA_LED                        *** 
                  ***                      iharden 07/21                     *** 
                  ***                                                        *** 
                  ***                                                        *** 
                  ************************************************************** 

The purpose of this little program is to facilitate the LED analysis as implemented in ORCA 
ORCA_LED does not help you with the calculations. It works only with ORCA output files 
Its main purpose is to print the essential results of LED + it performs the decomposition of the interaction energy 
However, it might be useful when you want to compare the outcomes of many LED analyzes 

Some general remarks:
Without using the compare-mode (--compare) I will only calculate the interaction energies 
To include the electronic preparation I need DLPNO-CCSD(T) calculations of the fragments in the dimer geometry 
To include the geometric preparation I need DLPNO-CCSD(T) calculations of the fragments in their equilibrium geometries 
Type 'orca_led --help' to get a list of options and some command line examples 

Installation: In /binaries there are precompiled binaries for Ubuntu and Windows 10. 

(The source files can easily be compiled by yourself. However, orca_led makes use of third party libraries (boost and fmt), so both libraries
need to be installed properly and the include paths need to be set in the appropriate manner.)
