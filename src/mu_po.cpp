#include"main.h"

double mu_po(double e_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double g[],double nu_el[],int points, int aa[], double energy[], double v[], double Ds[])
// It gives the effct of polar optical phonon on mobility in units of (cm^2/V.s)
// According to Equation (46) in Semiconductors and Semimetals volume1 10 (Rode's chapter),
// but with calculated group velocity from band structure and DOS both calculated from DFT:
{
    double integral_numerator = 0;
    double integral_denominator = 0;
    int factor = 1;

    double ddf,df,dv,k_step,de;
    //{
    if (T < T_trans)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            ddf = (df0dk(k_grid[counter+1],T,e_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa))/factor;
            for (int i = 0;i<=factor-1;i++)
                g[counter] = (-1)*e*E/(h_bar*nu_el[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)+i*ddf)*6.241509324e11;
                    //The last number is the conversion from convensional units to cancel out to be
                    // unitless (as in g)

        }
    }
    //}

    if (free_e ==1)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            dv = (v[counter+1]-v[counter])/factor;
            df = (f0(energy[counter+1],e_f,T)-f0(energy[counter],e_f,T))/factor;

            if(energy[counter+1] == 0)
            	df = 0;
            	
            k_step = (k_grid[counter+1]-k_grid[counter])/factor;

            for (int i = 0;i<=factor-1;i++)
            {

		integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v[counter]+i*dv)*g_LO[counter]/E;
		// =1/E*int[g_LO(En)*DOS(En)*v(En)*dEn]
		integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy[counter],e_f,T)+i*df);
                // =int[f0(En)*DOS(En)*dEn]
            }
        }
    }
    else
    {
        for (int counter = 0;counter<=points-2;counter++)
        {
            de = (energy[counter+1]-energy[counter]);

            if(energy[counter+1] == 0)
            	de = 0;

            integral_numerator = integral_numerator + de*Ds[counter]/volume1*v[counter]*g_LO[counter]/E;
            // =1/E*int[g_LO(En)*DOS(En)*v(En)*dEn]

            integral_denominator = integral_denominator + de*Ds[counter]/volume1*f0(energy[counter],e_f,T);
                // =int[f0(En)*DOS(En)*dEn]
            
            /*            
            if(isnan(integral_numerator))
            {
            	cout<<"isnan for mu_po  counter   =    "<<counter<<"g_LO    =    "<<g_LO[counter]<<endl;
            	getchar();
            }
            
            if(isinf(integral_numerator))
            {
            	cout<<"isinf for mu_po  counter   =    "<<counter<<"g_LO    =    "<<g_LO[counter]<<endl;
            	getchar();
            }
            */
        }
    }

    double mobility_po;

    if(geometry==1)  // for 3D
    	mobility_po = 1/3.0*integral_numerator/integral_denominator;
    else if(geometry==2)   // for 2D
    	mobility_po = 1/2.0*integral_numerator/integral_denominator;
    
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm$
    return mobility_po;

}

