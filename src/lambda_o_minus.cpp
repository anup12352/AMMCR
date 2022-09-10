
#include"main.h"


double lambda_o_minus(int counter,double omega,double A_minus, double epsilon_s,double epsilon_inf,int points)
// gives Lambda-_o in equations for inelastic optical phonon scattering; equation (117) of Rode's book
{
    //cout<<endl<<"Inside lambda_o_minus"<<endl;
	
    double k_minus = kminus_grid_pop[0][counter];

    //cout<<"counter = "<<counter<<endl;    
    //cout<<"k_minus = "<<k_minus<<endl;

    int minus_index = minus_index_pop[0][counter];;
    //cout<<"minus_index = "<<minus_index<<endl;

    double l;

    double k = k_grid[counter];
    //cout<<"k = "<<k<<endl;

    if ((energy_n[counter]<h_bar*omega)||(k_minus==k))
        l = 0;
        // If the energy of electron is lower than the phonon, there will be no emission
        // so lambda_minus terms will be zero (page 40 of Semiconductors and Semimetals, volume 10)
    else
    {   double aa = betaminus(counter,omega,epsilon_s, epsilon_inf, points);
        //cout<<"aa = "<<aa<<endl;
        l = aa*(A_minus*A_minus*log(abs((k_minus+k)/(k_minus-k)))-
         A_minus*c_n[counter]*c_n[minus_index]-a_n[counter]*a_n[minus_index]*c_n[counter]*c_n[minus_index]);

    }
    
    /*	
    cout<<"A_minus = "<<A_minus<<endl;
    cout<<"c_n = "<<c_n[counter]<<endl;
    cout<<"a_n = "<<a_n[counter]<<endl;
    cout<<"c_n minus_index = "<<c_n[minus_index]<<endl;
    cout<<"a_n  minus_index = "<<a_n[minus_index]<<endl;
    	
    cout<<"l = "<<l<<endl;
    cout<<"End of lambda_o_minus"<<endl;
    getchar();
    	*/
    	
    return l;

}
