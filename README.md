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

Installation: In /binaries there are precompiled and statically linked binaries for Ubuntu and Windows 10. 
You can just download the executables and use them right away.


## Usage
`orca_led --help` prints a help message
`orca_led --dimer example.out --monomers <list of monomers out files>` standard use for decomposting interaction energies.
`orca_led --compare file1.led file2.led` compares two .led files and calculates the differences of the interaction energies.

## Requirements
Tested with `ORCA 5.0`. Precompiled binaries for Ubuntu and Windows 10 are statically linked and do not have further dependencies. When compiling the source code from scratch, Third-Party libraries (Boost and fmt) are needed. Note: FMT was used in header_only mode so the FMT_HEADER_ONLY macro needs to be set.

## Contributor
contributed by Ingolf Harden
