
#include"main.h"

double mu_overall_time(double e_f,double T,double coefficients[5][7],double kindex[], double g[],double nu_el[],int points,int aa[], int j)
// It gives the overall mobility in units of (cm^2/V.s)
{
    // According to Equation (46) in Semiconductors and Semimetals volume1 10
    // (Rode's chapter), but with calculated group velocity from band structure and DOS both calculated from DFT:
    double integral_numerator = 0;
    double integral_denominator = 0;
    double ddf,k_step,de,dv,ds,df;
    int factor = 1;

	double f11;
		
    if (T < T_trans)
    {
        for (int counter = 0;counter<points-1;counter++)
        {            
            g[counter] = (-1)*Efield_time[j]/(h_bar*nu_el[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)*1e-7);
            // The last number is the conversion from convensional units to cancel out to be unitless (as in g)
            //end
        }
    }


    if (free_e ==1)
    {
        for (int counter =0;counter<=points-2;counter++)
        {
            dv = (v_n[counter+1] - v_n[counter])/factor;
            k_step = (k_grid[counter+1] - k_grid[counter])/factor;
					
                df = (f0(energy_n[counter+1],e_f,T) - f0(energy_n[counter],e_f,T))/factor;

		if(energy_n[counter+1]==0)
			df = 0;


            for (int i = 0;i<=factor-1;i++)
            {
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v_n[counter]+i*dv)*g[counter]/Efield_time[j];
                        // =1/Efield*int[g(En)*DOS(En)*v(En)*dEn]

                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy_n[counter],e_f,T)+i*df);
                    // =int[f(En)*DOS(En)*$

            }
            /*
            if (counter==199||counter==399||counter==599||counter==points-2)
            {
                cout<<"dv = "<<dv<<endl;
                cout<<"k_step = "<<k_step<<endl;
                cout<<"df = "<<df<<endl;
                cout<<"integral_numerator = "<<integral_numerator<<endl;
                cout<<"integral_denominator = "<<integral_denominator<<endl;
                getchar();
            }
            */


        }
    }
    else
    {
        for (int counter = 0;counter<=points-2;counter++)
        {
            de = (energy_n[counter+1] - energy_n[counter]);
            if(de<0)
            	de = 0;
            	
            integral_numerator = integral_numerator+de*(Ds_n[counter]/volume1)*v_n[counter]*g[counter]/Efield_time[j];
            
            //cout<<"counter = "<<counter<<"   integral_numerator = "<<integral_numerator<<endl;
                    // =1/Efield*int[g(En)*DOS(En)*v(En)*dEn]

            integral_denominator = integral_denominator + de*(Ds_n[counter]/volume1)*f0(energy_n[counter],e_f,T);
                    // =int[f(En)*DOS(En)*dEn]
                    // =int[f(En)*DOS(En)*dEn]
        }
    }
	
	
    double mobility_overall = (1/3.0)*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm$
    return mobility_overall;
}
