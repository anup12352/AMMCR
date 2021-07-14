#include"main.h"

double betaplus(int counter,double k_grid[],double omega,double epsilon_s, double epsilon_inf,
                double energy[], double v[], int points )
// gives \betaplus [1/s] in equations for inelastic optical phonon scattering; equation (118) of Rode's book
{
    double k_plus = kplus(counter,k_grid,omega,energy,points);

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_plus);
    int plus_index =FindMinInd(arr,points);

    double bp = (e*e*omega*k_plus)/(4*pi*h_bar*k_grid[counter]*v[plus_index])*
    (1/(epsilon_inf*epsilon_0)-1/(epsilon_s*epsilon_0))*3.895643846e28*1.60217657/1e8;

    return bp;
}

