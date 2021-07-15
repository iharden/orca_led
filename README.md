# orca_led
Prints main results from Local Energy Decomposition analysis and decomposes the binding energy between the monomers.

## Description
The purpose of this little program is to facilitate the LED analysis as implemented in ORCA (https://orcaforum.kofo.mpg.de/app.php/portal).
`orca_led` does not help you with the calculations. It works only with ORCA output files.
Its main purpose is to print the essential results of LED + it performs the decomposition of the binding energy.
However, it might be useful when you want to compare the outcomes of many LED analyzes. 

Some general remarks:
To include the electronic preparation it needs DLPNO-CCSD(T) calculations of the fragments in the dimer geometry. 
To include the geometric preparation it needs DLPNO-CCSD(T) calculations of the fragments in their equilibrium geometries. 

The program produces an output file `basename.led` and a `basename.csv` file that can be opened with e.g. Microsoft Excel.

Example output files can be found at [examples](examples/). 

Snippet:
```
******* DECOMPOSITION OF HARTREE-FOCK REFERENCE ENERGY ******* 

E(HF) for Fragment 1:                                                  -192.35466 
E(HF) for Fragment 2:                                                  -377.06253 

Electrostatic Interaction between fragments 2 and 1:                   -0.55525 

Exchange Interaction between fragments 2 and 1:                        -0.20502 

Consistency check:
Total Reference energy:                                                -570.17746 
Sum of the fragments HF reference and interaction energies:            -570.17746 
```
or
```
************************************************************** 
*******                Approach ONE                    ******* 
************************************************************** 



******* GEOMETRIC PREPARATION ENERGY ******* 

Geometric preparation for Fragment 1:                                  13.40991 
Geometric preparation for Fragment 2:                                  8.01998 
Sum:                                                                   21.42989 



******* HF-Interaction energy ******* 

Electronic preparation for Fragment 1:                                 296.86135 
Electronic preparation for Fragment 2:                                 183.61209 
Sum:                                                                   480.47344 

Electrostatic Interaction between fragments 2 and 1:                   -348.42518 

Exchange Interaction between fragments 2 and 1:                        -128.65306 

Consistency check:
Total HF interaction energy:                                           3.39520 
Sum of elec.preparation/electrostatic- and exchange Interaction:       3.39520 
```

## Usage
`orca_led --help` prints a help message.

`orca_led --dimer example.out --monomers <list of monomers out files>` standard use for printing main LED results decomposing binding energies.

`orca_led --compare file1.led file2.led` compares two .led files and calculates the differences of the various contributions to the binding energy.

## Requirements
Tested with `ORCA 5.0`. Precompiled binaries for Ubuntu and Windows 10 are statically linked and do not have further dependencies. When compiling the source code from scratch, Third-Party libraries (Boost and fmt) are needed. Boost.program_options is not header-only and therefore requires a proper installation of Boost. Note: fmt was used in header_only mode so the FMT_HEADER_ONLY macro needs to be set.

## Contributor
contributed by Ingolf Harden
