
#include"main.h"

double lambda_i_minus(int counter,double k_grid[],double omega,double A_minus,double c[],double epsilon_s,
                     double epsilon_inf,double energy[],double v[],int points)
// gives Lambda-_i in equations for inelastic optical phonon scattering; equation (123) of Rode's book
{
    double l;
    double k_minus = kminus(counter,k_grid,omega,energy,points);

    //cout<<endl<<"Inside lambda_i_minus"<<endl;
    //cout<<"k_minus = "<<k_minus<<endl;

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_minus);
    int minus_index =FindMinInd(arr,points);

    //cout<<"minus_index = "<<minus_index<<endl;

    double k = k_grid[counter];
    //cout<<"k = "<<k<<endl;
    //cout<<"c[counter] = "<<c[counter]<<endl;
    //cout<<"c[minus_index] = "<<c[minus_index]<<endl;

    if ((energy[counter]<h_bar*omega)||(k_minus==k))
        l = 0;
    //If the energy of electron is lower than the phonon, there will be no elmission so
    // lambda_minus terms will be zero (page 40 of Semiconductors and Semimetals, volume 10)
    else
    {
        double aa;
        aa = betaminus(counter,k_grid,omega, epsilon_s, epsilon_inf, energy,v,points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*((k_minus*k_minus+k*k)/(2*k_minus*k)*A_minus*A_minus*log(abs((k_minus+k)/(k_minus-k)))-
         A_minus*A_minus-(pow(c[counter],2))*(pow(c[minus_index],2)/3));
    }
    //cout<<"l  =" <<l<<endl;
    //cout<<"End of  lambda_i_minus"<<endl;
    return l;
}


