
#include"main.h"

double lambda_o_minus(int counter,double k_grid[],double omega,double A_minus,double a[],double c[],
        double epsilon_s,double epsilon_inf,double energy[],double v[],int points)
// gives Lambda-_o in equations for inelastic optical phonon scattering; equation (117) of Rode's book
{
    //cout<<endl<<"Inside lambda_o_minus"<<endl;

    double k_minus = kminus(counter,k_grid,omega,energy,points);
    //cout<<"k_minus = "<<k_minus<<endl;

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_minus);
    int minus_index =FindMinInd(arr,points);
    //cout<<"minus_index = "<<minus_index<<endl;

    double l;

    double k = k_grid[counter];
    //cout<<"k = "<<k<<endl;

    if ((energy[counter]<h_bar*omega)||(k_minus==k))
        l = 0;
        // If the energy of electron is lower than the phonon, there will be no emission
        // so lambda_minus terms will be zero (page 40 of Semiconductors and Semimetals, volume 10)
    else
    {   double aa = betaminus(counter,k_grid,omega,epsilon_s, epsilon_inf, energy,v,points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*(A_minus*A_minus*log(abs((k_minus+k)/(k_minus-k)))-
         A_minus*c[counter]*c[minus_index]-a[counter]*a[minus_index]*c[counter]*c[minus_index]);

    }
    //cout<<"l = "<<l<<endl;
    //cout<<"End of lambda_o_minus"<<endl;

    return l;

}
