# AMMCR
Ab initio model for mobility and conductivity calculation by using Rode Algorithm (AMMCR), solves Boltzmann transport equation within Rode's iterative scheme.
At present it works for only electrons.

It is written in C++. Currently it is interfaced with Vienna Ab initio Simulation Package (VASP). It is very simple to use. Once compiled, it produces an executable called AMMCR. The executable can be run in a directory where four VASP output files: **OUTCAR, EIGENVAL, PROCAR and DOSCAR** are present. 

If the band structure calculation is done in spin-polarized mode, the code automatically detetcs it, and produces two sets of results: For majority and 
minority spin electrons.


The specification of the simulation conditions such as temperature, doping concentration, external electric field etc., as well as the materials constants are provided in a file called **input.dat**. 

There is also an youtube video about how to compile and run the AMMCR code:
https://youtu.be/II-7I68iPnI

If you use AMMCR for your research, please cite the following article:

Anup Kumar Mandia, Bhaskaran Muralidharan, Jung-Hae Choi, Seung-Cheol Lee and Satadeep Bhattacharjee,
**AMMCR: Ab initio model for mobility and conductivity calculation by using Rode Algorithm,**
Computer Physics Communications,
Volume 259,
2021,
107697,
ISSN 0010-4655,
https://doi.org/10.1016/j.cpc.2020.107697.
