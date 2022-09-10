
#include"main.h"

double lambda_o_plus(int counter,double omega,double A_plus, double epsilon_s,double epsilon_inf,int points)
// gives Lambda+_o in equations for inelastic optical phonon scattering; equation (117) of Rode's book
{
    //cout<<endl<<"Inside lambda_o_plus"<<endl;
    double k_plus = kplus_grid_pop[0][counter];

    int plus_index = plus_index_pop[0][counter];

    double l;

    double k = k_grid[counter];
    if ((k_plus==k))
        l = 0;
    else
    {
        double aa = betaplus(counter, omega,epsilon_s, epsilon_inf, points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*(A_plus*A_plus*log(abs((k_plus+k)/(k_plus-k)))
            -A_plus*c_n[counter]*c_n[plus_index]-a_n[counter]*a_n[plus_index]*c_n[counter]*c_n[plus_index]);

    }
    //cout<<"l = "<<l<<endl;
    //cout<<"End of  lambda_o_plus"<<endl;
    return l;
}
