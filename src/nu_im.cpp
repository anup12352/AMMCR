#include"main.h"

double nu_im(double k_dum,int counter, double epsilon_s,double N_im,double v)
// Neutral impurity scattering rate according to equation 4.25 from Garick Ng thesis this is write expression
{
    double im = (N_im*1e6)*(80*pi*epsilon_s*epsilon_0*(h_bar*e)*pow((v/100),2))/(e*e*pow((k_dum*1e9),2));

    return im;
// 1e6 is to convert cm-3 to m-3
// e with hbar is to convert eV-s to J-s
// v(counter) is divided with 100 to convert cm/s to m/s
// 1e9 with k_dum to convert 1/nm to 1/m

}
