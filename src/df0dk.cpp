#include"main.h"

double df0dk(double k_grid[],double k,double T,double E_f, double coefficients[5][7],double kindex[], int a[])
{
    double DkFermi;
    double En = conduction_dispersion(k,coefficients,kindex,a);
    if (f0(En,E_f,T) < 1e-300)
        DkFermi=0;
    else
        DkFermi=-1*exp((En-E_f)/(k_B*T))/(k_B*T*pow((exp((En-E_f)/(k_B*T))+1),2))*dedk(k,coefficients,kindex,a);
        // Based on chain rule: df0/dk=df0/de*de/dk
        // unit (nm)
    return DkFermi;
}



