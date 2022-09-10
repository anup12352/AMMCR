#include "main.h"

double conduction_dispersion(double k,double coefficients[5][7],double kindex[4],int a[2])
{
    double l = a[1];
    int segment_of_the_curve = 0;
    for (int j=0;j<l;j++)    // l contains length of kindex
    {
         if (k <= kindex[j])
         {
            segment_of_the_curve = j;
            break;
         }
         else
            segment_of_the_curve = l;
    }
    double energy=0;
    double *powers;
    powers = linspace(a[0],0,a[0]+1);
    for(int i=0;i<a[0]+1;i++)
        energy = energy + coefficients[segment_of_the_curve][i]*(pow(k,powers[i]));

    return energy;
}

