#include"main.h"

double Aminus(int counter,double k_grid[],double omega, double a[],double c[],double energy[], int points)
{
    double k_minus =kminus(counter,k_grid,omega,energy,points);
    double AA;
    if (energy[counter] < h_bar*omega)
        AA =0;
    else
    {
        double arr[points];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k_minus);
        int minus_index =FindMinInd(arr,points);
        double k = k_grid[counter];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k);
        int index =FindMinInd(arr,points);

        AA = a[index]*a[minus_index]+(k_minus*k_minus+k*k)/(2*k_minus*k)*c[index]*c[minus_index];
    }
    return AA;
}
