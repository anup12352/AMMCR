
#include"main.h"

double mu_elastic(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],
                  double Ds[],double v[],double nu_elastic[],double g[],int points, int a[])
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
            dv = (v[counter+1]-v[counter])/factor;


            ddf = (df0dk(k_grid,k_grid[counter+1],T,E_f,coefficients,kindex,a)-df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a))/factor;
            if (T < T_trans)
                df = (f0(energy[counter+1],E_f,T)-f0(energy[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],k_grid,E_f,T,coefficients,kindex,g,points,a)-f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a))/factor;

            for (int i=0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*e*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a)+i*ddf)*6.241509324e11;
                // The last number is the conversion from convensional units to cancel out to be unitless (as in g)
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v[counter]+i*dv)*g_elastic[counter]/E;
                // =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy[counter],E_f,T)+i*df);
                    ////// =int[f(En)*DOS(En)*dEn]
                else
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a)+i*df);
            }
        }
    }
    else
    {
        for (int counter = 0;counter<dos_intervals-1;counter++)
        {
            de = (energy[counter+1]-energy[counter])/factor;
            ds = (Ds[counter+1]-Ds[counter])/factor;
            dv = (v[counter+1]-v[counter])/factor;
            ddf = (df0dk(k_grid,k_grid[counter+1],T,E_f,coefficients,kindex,a)-df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a))/factor;
            if (T < T_trans)
                df = (f0(energy[counter+1],E_f,T)-f0(energy[counter],E_f,T))/factor;
            else
                df = (f(k_grid[counter+1],k_grid,E_f,T,coefficients,kindex,g,points,a)-f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a))/factor;

            for (int i = 0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*e*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid,k_grid[counter],T,E_f,coefficients,kindex,a)+i*ddf)*6.241509324e11;
                //// The last number is the conversion from convensional units to cancel out to be unitless (as in g)
                integral_numerator = integral_numerator+de*(Ds[counter]+i*ds)/volume1*(v[counter]+i*dv)*g_elastic[counter]/E;
                //// =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
                if (T < T_trans)
                    integral_denominator = integral_denominator+de*(Ds[counter]+i*ds)/volume1*(f0(energy[counter],E_f,T)+i*df);
                    // =int[f(En)*DOS(En)*dEn]
                else
                    integral_denominator = integral_denominator+de*(Ds[counter]+i*ds)/volume1*(f(k_grid[counter],k_grid,E_f,T,coefficients,kindex,g,points,a)+i*df);
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

