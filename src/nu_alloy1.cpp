
#include"main.h"

double nu_alloy1(double k, double v)
{
    double alloy = (3*pi*pow((k*1e9),2)*pow(Uall,2)*(V0*1e-27)*xx*(1-xx))/(16*h_bar*h_bar*(v*1e-2) );
    // Equation No. 17 from Ramu paper 1
    // The coefficient of conversion is to get 1/s
    return alloy;
}
