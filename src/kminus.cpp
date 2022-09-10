
#include"main.h"

double kminus(int counter, double omega, int points, double energy[])
{
    double E0 = energy[counter];
    double E1 = E0 - h_bar*omega;

    double arr[points];
    double kk;
    if (E0>h_bar*omega)
    {
        for (int i=0;i<points;i++)
            arr[i] = abs(energy[i] - E1);
        int index =FindMinInd(arr,points);
        kk = k_grid[index];
    }
    else
        kk = k_grid[0]/2;

    if (abs(k_grid[counter]-kk)==0)
        kk = kk+1e-4;

    return kk;

}

