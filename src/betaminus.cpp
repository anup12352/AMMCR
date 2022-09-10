#include"main.h"

double betaminus(int counter, double omega,double epsilon_s, double epsilon_inf, int points)
// gives \betaplus [1/s] in equations for inelastic optical phonon scattering; equation (118) of Rode's book
{
    double k_minus = kminus_grid_pop[0][counter];
	
    int minus_index;
    	
    if (energy_n[counter]<h_bar*omega)
    {
        k_minus = k_grid[counter];
	double arr[points];

	for (int i=0;i<points;i++)
	arr[i] = abs(k_grid[i] - k_minus);
	minus_index =FindMinInd(arr,points);
    }
    else		
    	minus_index = minus_index_pop[0][counter];

    double bp = (e*e*omega*k_minus)/(4*pi*h_bar*k_grid[counter]*v_n[minus_index])*
    (1/(epsilon_inf*epsilon_0)-1/(epsilon_s*epsilon_0))*3.895643846e28*1.60217657/1e8;

    //cout<<"k_minus  = "<<k_minus<<endl;
    //cout<<"minus_index = "<<minus_index<<endl;
	
    return bp;
}

