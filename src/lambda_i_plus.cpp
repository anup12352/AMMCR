
#include"main.h"

double lambda_i_plus(int counter,double omega,double A_plus,double epsilon_s, double epsilon_inf,int points)
// gives Lambda+_i in equations for inelastic optical phonon scattering; equation (123) of Rode's book
{
    double l;

    double k_plus = kplus_grid_pop[0][counter];
    //cout<<endl<<"Inside lambda_i_plus"<<endl;
    //cout<<"A_plus = "<<A_plus<<endl;
    //cout<<"omega = "<<omega<<endl;
    //cout<<"k_plus = "<<k_plus<<endl;
    //cout<<"c_n[counter] = "<<c_n[counter]<<endl;

    //A_plus =  47.695069008345868;
    //omega =  3.694512960621596e+13;

    int plus_index = plus_index_pop[0][counter];

    //cout<<"plus_index = "<<plus_index<<endl;
    //cout<<"c[plus_index] = "<<c[plus_index]<<endl;

    double k = k_grid[counter];
    // if abs(k_plus-k_grid(length(k_grid)))<1e-5;
    //cout<<" k = "<<k<<endl;
    if (k_plus == k)
            l = 0;
    else
    {
        //double yy = ;
        //cout<<"aa = "<<aa<<endl;
        //aa =      6.115495758768291e+16;
        l = betaplus(counter,omega,epsilon_s,epsilon_inf,points)*
        ((pow(k_plus,2)+k*k)/(2*k_plus*k) * A_plus*A_plus*log(abs((k_plus+k)/(k_plus-k))) -
                A_plus*A_plus - (c_n[counter]*c_n[counter])*(pow(c_n[plus_index],2)/3));
    }
    //cout<<"l = "<<l<<endl;
    //cout<<"End of lambda_i_plus"<<endl;
    return l;
}
