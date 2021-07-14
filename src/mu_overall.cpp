
#include"main.h"

double mu_overall(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],double Ds[],
                  double v[],double g[],double nu_el[],int points,int a[])
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
            ddf = (df0dk(k_grid,k_grid[counter+1],T,E_f,coefficients,kindex,a)-df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a))/factor;
            //for i = 0:factor-1
            g[counter] = (-1)*e*E/(h_bar*nu_el[counter])*(df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a)+ddf*9)*6.241509324e11;
            // The last number is the conversion from convensional units to cancel out to be unitless (as in g)
            //end
        }
    }


    if (free_e ==1)
    {
        for (int counter =0;counter<=points-2;counter++)
        {
            dv = (v[counter+1] - v[counter])/factor;
            k_step = (k_grid[counter+1] - k_grid[counter])/factor;

            if (T < T_trans)
                df = (f0(energy[counter+1],E_f,T) - f0(energy[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],k_grid,E_f,T,coefficients,kindex,g,points,a)-f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a))/factor;


            for (int i = 0;i<=factor-1;i++)
            {
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v[counter]+i*dv)*g[counter]/E;
                        // =1/E*int[g(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy[counter],E_f,T)+i*df);
                    // =int[f(En)*DOS(En)*$
                else
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a)+i*df);  ////// =int[f(En)*DOS(En)*$

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
            de = (energy[counter+1] - energy[counter]);
            integral_numerator = integral_numerator+de*(Ds[counter]/volume1)*v[counter]*g[counter]/E;
                    // =1/E*int[g(En)*DOS(En)*v(En)*dEn]
            if (T < T_trans)
                integral_denominator = integral_denominator + de*(Ds[counter]/volume1)*f0(energy[counter],E_f,T);
                    // =int[f(En)*DOS(En)*dEn]
            else
                integral_denominator = integral_denominator + de*(Ds[counter]/volume1)* f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a);
                    // =int[f(En)*DOS(En)*dEn]
        }
    }

    double mobility_overall = (1/3.0)*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm$
    return mobility_overall;
}
