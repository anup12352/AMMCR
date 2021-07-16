
#include"main.h"

double mu_overall(double E_f,double T,double coefficients[5][7],double kindex[], double g[],double nu_el[],int points,int aa[])
// It gives the overall mobility in units of (cm^2/V.s)
{
    // According to Equation (46) in Semiconductors and Semimetals volume1 10
    // (Rode's chapter), but with calculated group velocity from band structure and DOS both calculated from DFT:
    double integral_numerator = 0;
    double integral_denominator = 0;
    double ddf,k_step,de,dv,ds,df;
    int factor = 10;


    if (T < T_trans)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            ddf = (df0dk(k_grid[counter+1],T,E_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa))/factor;
            //for i = 0:factor-1
            g[counter] = (-1)*e*E/(h_bar*nu_el[counter])*(df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa)+ddf*9)*6.241509324e11;
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

            if (T < T_trans)
                df = (f0(energy_n[counter+1],E_f,T) - f0(energy_n[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],E_f,T,coefficients,kindex,g,points,aa)-f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa))/factor;


            for (int i = 0;i<=factor-1;i++)
            {
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v_n[counter]+i*dv)*g[counter]/E;
                        // =1/E*int[g(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy_n[counter],E_f,T)+i*df);
                    // =int[f(En)*DOS(En)*$
                else
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa)+i*df);  ////// =int[f(En)*DOS(En)*$

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
            integral_numerator = integral_numerator+de*(Ds_n[counter]/volume1)*v_n[counter]*g[counter]/E;
                    // =1/E*int[g(En)*DOS(En)*v(En)*dEn]
            if (T < T_trans)
                integral_denominator = integral_denominator + de*(Ds_n[counter]/volume1)*f0(energy_n[counter],E_f,T);
                    // =int[f(En)*DOS(En)*dEn]
            else
                integral_denominator = integral_denominator + de*(Ds_n[counter]/volume1)* f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa);
                    // =int[f(En)*DOS(En)*dEn]
        }
    }

    double mobility_overall = (1/3.0)*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm$
    return mobility_overall;
}
