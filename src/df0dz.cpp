
#include"main.h"

double df0dz(double k, double e_f, double T, double df0dz_integral, double coefficients[5][7], double kindex[], int aa[])
{
    // It gives the driving force for transport under dTdz in unit of 1/cm,
    // required for thermopower calculations

    // According to equation (54) in Semiconductors and Semimetals, volume 10 (Rode's chapter):

    double E1 = conduction_dispersion(k,coefficients,kindex,aa);
    double df0dz_temp = f0(E1,e_f,T)*(1-f0(E1,e_f,T))*(E1/(k_B*T)-df0dz_integral)*1/T*dTdz;
    return df0dz_temp;
    // According to equation (54) in Rode's book. Unit is 1/length consistent with units of dTdz

}
