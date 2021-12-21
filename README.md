# AMMCR
Ab initio model for mobility and conductivity calculation by using Rode Algorithm (AMMCR), solves Boltzmann transport equation within Rode's iterative scheme.
At present it works for only electrons.

It is written in C++. Currently it is interfaced with Vienna Ab initio Simulation Package (VASP). It is very simple to use. Once compiled, it produces an executable called AMMCR. The executable can be run in a directory where four VASP output files: **OUTCAR, EIGENVAL, PROCAR and DOSCAR** are present. 

If the band structure calculation is done in spin-polarized mode, the code automatically detects it, and produces two sets of results: For majority and 
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
For low electric field, the total distribution function can be written as
                                         
                                                fk=f0(k)+xg(k)

Where f0(k) is the _equilibrium distribution function_ and g(k) is the _perturbation_ (out of equilibrium part of the distribution function).  The factor x defines the angle between the electric field and the direction of the crystal momentum k.

The flow-chart of the code is shown here. First, one has to calculate all the required inputs using the first -principles method. The typical inputs are energy band dispersion curve, density of states, the phonon frequencies and the elastic constants.

We then perform the analytical fitting of the band structure to obtain smooth curves to calculate the group velocity. For a given doping concentration, the Fermi energy is calculated. Then we calculate various scattering rates (discussed below). Finally a loop is used to obtain the out of equilibrium part, g(k) of the distribution function.

![image](https://user-images.githubusercontent.com/68414451/126433576-e6e68c38-5f22-4769-b3a1-37280aef6dfd.png)

_The various scattering mechanisms which are implemented in the code are:_

**Defect induced scatterings from:**

Ionized ion impurities

Neutral impurities

Alloy 

Dislocations

**Scattering mechanisms involving lattice:**

Piezoelectric

Intravalley acoustic deformation potential

Polar optical phonons.

 
Once the loop converges for g(k), all the transport co-efficients are calculated.

**Studies where AMMCR was used:**
1. Mandia, A. K., Patnaik, R., Muralidharan, B., Lee, S. C., & Bhattacharjee, S. (2019). 
Ab initio semi-classical electronic transport in ZnSe: the role of inelastic scattering mechanisms. 
Journal of Physics: Condensed Matter, 31(34), 345901, https://doi.org/10.1088/1361-648X/ab229b.

![11](https://user-images.githubusercontent.com/68414451/146901858-f4073509-1d2c-4540-81c5-5dae83c1edc0.jpg)
![12](https://user-images.githubusercontent.com/68414451/146901867-7862c11d-b97e-4623-be01-ada7f81f9d59.jpg)
![13](https://user-images.githubusercontent.com/68414451/146901870-9728a86f-35d2-4f12-af7e-1c626ca27383.jpg)
![14](https://user-images.githubusercontent.com/68414451/146901871-1175deed-527f-413e-b43c-2ee47d7567df.jpg)


2. Chakrabarty, S., Mandia, A. K., Muralidharan, B., Lee, S. C., & Bhattacharjee, S. (2019). 
Semi-classical electronic transport properties of ternary compound AlGaAs2: role of different scattering mechanisms. 
Journal of Physics: Condensed Matter, 32(13), 135704, https://doi.org/10.1088/1361-648X/ab5edf.  

![21](https://user-images.githubusercontent.com/68414451/146901988-1494180b-4d50-445f-b2f8-cd5c8ee5680a.jpg)
![22](https://user-images.githubusercontent.com/68414451/146901998-7288c334-dd96-4ac1-9d36-2f0df2852ddd.jpg)
![23](https://user-images.githubusercontent.com/68414451/146902003-a10742bb-4d18-47c7-9ba1-1ce04fd62722.jpg)
![24](https://user-images.githubusercontent.com/68414451/146902006-e4f43abe-0f5f-4f45-a619-e6b3fde29e5b.jpg)

3. Kumar, U., Nayak, S., Chakrabarty, S., Bhattacharjee, S., & Lee, S. C. (2020). 
Gallium–Boron–Phosphide ($$\hbox {GaBP} _ {2} $$ GaBP 2): a new III–V semiconductor for photovoltaics.
Journal of Materials Science, 55(22), 9448-9460, https://doi.org/10.1007/s10853-020-04631-5.

