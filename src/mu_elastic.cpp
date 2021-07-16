
#include"main.h"

double mu_elastic(double E_f,double T,double coefficients[5][7],double kindex[], double nu_elastic[],double g[],int points, int aa[])
// It gives effect of an elastic scattering mechanism on mobility in units of (cm^2/V.s)
{
    // According to Equation
    // (A7) in	PRB(2011) doi: 10.1103/PhysRevB.84.075315  ,
    // but with calculated group velocity from band structure and DOS both calculated from DFT:

    double g_elastic[points]={0};
    double integral_numerator = 0;
    double integral_denominator = 0;
    int dos_intervals = points;
    int factor = 10;

    double k_step,dv,ddf,df,de,ds;
    if (free_e ==1)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            k_step = (k_grid[counter+1]-k_grid[counter])/factor;
            dv = (v_n[counter+1]-v_n[counter])/factor;


            ddf = (df0dk(k_grid[counter+1],T,E_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa))/factor;
            if (T < T_trans)
                df = (f0(energy_n[counter+1],E_f,T)-f0(energy_n[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],E_f,T,coefficients,kindex,g,points,aa)-f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa))/factor;

            for (int i=0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*e*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa)+i*ddf)*6.241509324e11;
                // The last number is the conversion from convensional units to cancel out to be unitless (as in g)
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v_n[counter]+i*dv)*g_elastic[counter]/E;
                // =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy_n[counter],E_f,T)+i*df);
                    ////// =int[f(En)*DOS(En)*dEn]
                else
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa)+i*df);
            }
        }
    }
    else
    {
        for (int counter = 0;counter<dos_intervals-1;counter++)
        {
            de = (energy_n[counter+1]-energy_n[counter])/factor;
            ds = (Ds_n[counter+1]-Ds_n[counter])/factor;
            dv = (v_n[counter+1]-v_n[counter])/factor;
            ddf = (df0dk(k_grid[counter+1],T,E_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa))/factor;
            if (T < T_trans)
                df = (f0(energy_n[counter+1],E_f,T)-f0(energy_n[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],E_f,T,coefficients,kindex,g,points,aa)-f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa))/factor;

            for (int i = 0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*e*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid[counter],T,E_f,coefficients,kindex,aa)+i*ddf)*6.241509324e11;
                //// The last number is the conversion from convensional units to cancel out to be unitless (as in g)
                integral_numerator = integral_numerator+de*(Ds_n[counter]+i*ds)/volume1*(v_n[counter]+i*dv)*g_elastic[counter]/E;
                //// =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+de*(Ds_n[counter]+i*ds)/volume1*(f0(energy_n[counter],E_f,T)+i*df);
                    // =int[f(En)*DOS(En)*dEn]
                else
                    integral_denominator = integral_denominator+de*(Ds_n[counter]+i*ds)/volume1*(f(k_grid[counter],E_f,T,coefficients,kindex,g,points,aa)+i*df);
                    // =int[f(En)*DOS(En)*dEn]
            }
        }
    }
    //cout<<"integral_numerator = "<<integral_numerator<<endl;
    //cout<<"integral_denominator = "<<integral_denominator<<endl;

    double mobility_elastic = 1/3.0*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm/s) and E (V/cm)
    return mobility_elastic;
}

