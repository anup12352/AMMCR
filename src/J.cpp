
#include"main.h"

double J(double T,double m,double g_th[], int points)
// It gives the current density in units of (???) from dT/dz as a driving force in
// thermopower calculation
// According to equation (A10) of PRB, doi: 10.1103/PhysRevB.84.075315
{


    double integral_dum = 0;
    double dk;

    for (int counter = 0;counter<=points-2;counter++)
    {
        dk = k_grid[counter+1]-k_grid[counter];
        integral_dum = integral_dum + dk*pow((k_grid[counter]/pi),2)*v_n[counter]*g_th[counter];
        // J=e/3*integral of [v(En)*g_th(En)*DOS(En)
    }

    double current_density = e/3*integral_dum*1e21;
    return current_density;

    // The conversion coefficient is to bring J to the units of A/cm^2
    // consistent with the rest of the units used here
}
