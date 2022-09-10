#include"main.h"

double kplus(int counter, double omega, int points, double energy[])
{
    double E0 = energy[counter];
    double E1 = E0 + h_bar*omega;

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(energy[i] - E1);
    int index =FindMinInd(arr,points);
    double kk = k_grid[index];
    if (abs(k_grid[counter]-kk)==0)
        kk = kk+1e-4;

    return kk;

}


