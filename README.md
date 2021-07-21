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

**AMMCR flowchart:**

We have implemented the Rode algorithm in the AMMCR code, which is valid for the low electric field.
For low electric field, the total distribution can be written as
                                         
                                                fk=f0(k)+xg(k)

Where f0(k) is the equilibrium distribution function and g(k) is the perturbation (non-equilibrium part of the distribution function).  The factor x defines the angle between the electric field and the direction of the crystal momentum k.

The flow-chart of the code is shown here. First, one has to calculate all the required inputs using the first -principles method. The typical inputs are energy band dispersion curve, density of states, the phonon frequencies and the elastic constants.

We then perform the analytical fitting of the band structure to obtain smooth curves to calculate the group velocity. For a given doping concentration, the Fermi energy is calculated. Then we calculate various scattering rates (discussed below). Finally a loop is used to obtain the non-equilibrium part, g(k) of the distribution function.

![image](https://user-images.githubusercontent.com/68414451/126433576-e6e68c38-5f22-4769-b3a1-37280aef6dfd.png)

The various scattering mechanisms which are implemented in the code are:
Defect induced scatterings from:
Ionized ion impurities
Neutral impurities
Alloy 
Dislocations
Scattering mechanisms involving lattice:
Piezoelectric
Intravalley acoustic deformation potential
Polar optical phonons.
 
Once the loop converges for g(k), all the transport co-efficients are calculated.
