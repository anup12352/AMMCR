
#include"main.h"

double lambda_i_minus(int counter,double omega,double A_minus, double epsilon_s, double epsilon_inf,int points)
// gives Lambda-_i in equations for inelastic optical phonon scattering; equation (123) of Rode's book
{
    double l;
    double k_minus = kminus_grid_pop[0][counter];

    //cout<<endl<<"Inside lambda_i_minus"<<endl;
    //cout<<"k_minus = "<<k_minus<<endl;

    int minus_index = minus_index_pop[0][counter];;

    //cout<<"minus_index = "<<minus_index<<endl;

    double k = k_grid[counter];
    //cout<<"k = "<<k<<endl;
    //cout<<"c[counter] = "<<c[counter]<<endl;
    //cout<<"c[minus_index] = "<<c[minus_index]<<endl;

    if ((energy_n[counter]<h_bar*omega)||(k_minus==k))
        l = 0;
    //If the energy of electron is lower than the phonon, there will be no elmission so
    // lambda_minus terms will be zero (page 40 of Semiconductors and Semimetals, volume 10)
    else
    {
        double aa;
        aa = betaminus(counter, omega, epsilon_s, epsilon_inf, points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*((k_minus*k_minus+k*k)/(2*k_minus*k)*A_minus*A_minus*log(abs((k_minus+k)/(k_minus-k)))-
         A_minus*A_minus-(pow(c_n[counter],2))*(pow(c_n[minus_index],2)/3));
    }
    //cout<<"l  =" <<l<<endl;
    //cout<<"End of  lambda_i_minus"<<endl;
    return l;
}


