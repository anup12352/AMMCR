#include"main.h"
                
double df0dk(double k,double T,double e_f, double coefficients[5][7],double kindex[], int aa[])
{
    double DkFermi;
    double En = conduction_dispersion(k,coefficients,kindex,aa);
    if (f0(En,e_f,T) < 1e-300)
        DkFermi=0;
    else
        DkFermi=-1*exp((En-e_f)/(k_B*T))/(k_B*T*pow((exp((En-e_f)/(k_B*T))+1),2))*dedk(k,coefficients,kindex,aa);
        // Based on chain rule: df0/dk=df0/de*de/dk
        // unit (nm)
    return DkFermi;
}



