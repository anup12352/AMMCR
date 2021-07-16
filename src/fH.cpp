
#include"main.h"

double fH(double k,double E_f,double T,double  coefficients[5][7],double kindex[],double g[], double h[],int points, int aa[])
{

// gives the electron distribution by calling g as a function (when iteration==0 f==f0)

    double E1=conduction_dispersion(k,coefficients,kindex,aa);
    int n=0;
    double min1=100000;
    //here we assumed that x, cosine of the angle between vector k and the
    // electric field, is 1 and they are in the same direction %%%

    for (int i=0;i<points;i++)
    {
            if (abs(k_grid[i]-k)<min1)
            {
                n=i;
                min1=abs(k_grid[i]-k);
            }
    }

	double distribution=1/(1+exp((E1-E_f)/(k_B*T)))+ g[n] + h[n];
    return distribution;
}
