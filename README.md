# AMMCR
Ab initio model for mobility and conductivity calculation by using Rode Algorithm (AMMCR).
It is written in C++. Currently it is interfaced with Vienna Ab initio Simulation Package (VASP). It is very simple to use. Once compiled, it produces an executable called AMMCR. The executable can be run in a directory where four VASP output files: OUTCAR, EIGENVAL, PROCAR and DOSCAR are present. The specification of the simulation conditions such as temperature, doping concentration, external electric field etc., as well as the materials constants are provided in a file called “input.dat”. 

There is also an youtube video about how to compile and run the AMMCR code:
https://youtu.be/II-7I68iPnI
