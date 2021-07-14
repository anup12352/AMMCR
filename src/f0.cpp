#include"main.h"

double f0(double E1, double E_f, double T)
{
    double fermi = 1/(1+exp((E1-E_f)/(k_B*T)));
    return fermi;
}

