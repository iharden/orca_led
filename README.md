# orca_led
Prints main results from Local Energy decomposition analysis and decomposes the interaction energy.

## Description
The purpose of this little program is to facilitate the LED analysis as implemented in ORCA (https://orcaforum.kofo.mpg.de/app.php/portal).
ORCA_LED does not help you with the calculations. It works only with ORCA output files.
Its main purpose is to print the essential results of LED + it performs the decomposition of the interaction energy.
However, it might be useful when you want to compare the outcomes of many LED analyzes. 

Some general remarks:
Without using the compare-mode (--compare) it will only calculate the interaction energies. 
To include the electronic preparation it needs DLPNO-CCSD(T) calculations of the fragments in the dimer geometry. 
To include the geometric preparation it needs DLPNO-CCSD(T) calculations of the fragments in their equilibrium geometries. 

The program produces an output file `basename.led` and a `basename.csv` file that can be opend with e.g. Microsoft Excel.

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
******* I will start now with the Interaction Energies ******* 
************************************************************** 

You gave me 4 monomer files and 2 fragments are present in the Dimer-file! 
Name of Dimer-file 
   endo.mpi8.out 

Name of Monomer-files 
   frag1.mpi8.out 
   frag2.mpi8.out 
   frag1_opt.mpi8.out 
   frag2_opt.mpi8.out 
                              Dimer         Monomer  1    Monomer  2    

E(HF):                        -570.17746    -192.82774    -377.35513    
E(CCSD):                      -572.32448    -193.65822    -378.64950    
E(CCSD(T)):                   -572.43223    -193.69929    -378.71098    
E(CORR):                      -2.14702      -0.83048      -1.29437      
E(Triples):                   -0.10776      -0.04107      -0.06147      

I will calculate the geometric and the electronic preparation 

If not stated otherwise, energies from now on are in kcal/mol !!! 



******* Geometric preparation energy ******* 

Geometric preparation for Fragment 1:                                  13.40991 
Geometric preparation for Fragment 2:                                  8.01998 
Sum:                                                                   21.42989 
```

## Usage
`orca_led --help` prints a help message.

`orca_led --dimer example.out --monomers <list of monomers out files>` standard use for decomposting interaction energies.

`orca_led --compare file1.led file2.led` compares two .led files and calculates the differences of the interaction energies.

## Requirements
Tested with `ORCA 5.0`. Precompiled binaries for Ubuntu and Windows 10 are statically linked and do not have further dependencies. When compiling the source code from scratch, Third-Party libraries (Boost and fmt) are needed. Note: FMT was used in header_only mode so the FMT_HEADER_ONLY macro needs to be set.

## Contributor
contributed by Ingolf Harden
