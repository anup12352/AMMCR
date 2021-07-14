
#include"main.h"

double lambda_o_plus(int counter,double k_grid[],double omega,double A_plus,double a[],double c[],
        double epsilon_s,double epsilon_inf,double energy[],double v[],int points)
// gives Lambda+_o in equations for inelastic optical phonon scattering; equation (117) of Rode's book
{
    //cout<<endl<<"Inside lambda_o_plus"<<endl;
    double k_plus = kplus(counter,k_grid,omega,energy,points);

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_plus);
    int plus_index =FindMinInd(arr,points);

    double l;

    double k = k_grid[counter];
    if ((k_plus==k))
        l = 0;
    else
    {
        double aa = betaplus(counter,k_grid,omega,epsilon_s, epsilon_inf,energy,v, points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*(A_plus*A_plus*log(abs((k_plus+k)/(k_plus-k)))
            -A_plus*c[counter]*c[plus_index]-a[counter]*a[plus_index]*c[counter]*c[plus_index]);

    }
    //cout<<"l = "<<l<<endl;
    //cout<<"End of  lambda_o_plus"<<endl;
    return l;
}
