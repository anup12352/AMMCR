
#include"main.h"

double mu_elastic(double e_f,double T,double coefficients[5][7],double kindex[], double nu_elastic[],int points, int aa[], double energy[], double v[], double Ds[])
// It gives effect of an elastic scattering mechanism on mobility in units of (cm^2/V.s)
{
    // According to Equation
    // (A7) in	PRB(2011) doi: 10.1103/PhysRevB.84.075315  ,
    // but with calculated group velocity from band structure and DOS both calculated from DFT:

    double g_elastic[points]={0};
    double integral_numerator = 0;
    double integral_denominator = 0;
    int dos_intervals = points;
    int factor = 1;

    double k_step,dv,ddf,df,de,ds;
    
    
    if (free_e ==1)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            k_step = (k_grid[counter+1]-k_grid[counter])/factor;
            dv = (v[counter+1]-v[counter])/factor;


            ddf = (df0dk(k_grid[counter+1],T,e_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa))/factor;
	    df = (f0(energy[counter+1],e_f,T)-f0(energy[counter],e_f,T))/factor;
	    
	    if(energy[counter+1]==0)
	    	df = 0;
	    		
            for (int i=0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)+i*ddf)*1e-7;
                // The last number is the conversion from convensional units to cancel out to be unitless (as in g)

		integral_numerator = integral_numerator+ k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*
		(v[counter]+i*dv)*g_elastic[counter]/E;
	        // =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
		// units of group velocity (cm/s) and E (V/cm)
		
		integral_denominator = integral_denominator + 
		k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy[counter],e_f,T)+i*df);
		////// =int[f(En)*DOS(En)*dEn]

            }
        }
    }
    else
    {
        for (int counter = 0;counter<dos_intervals-1;counter++)
        {
            de = (energy[counter+1]-energy[counter])/factor;
            if(energy[counter+1]==0)
            	de = 0;
            	
            ds = (Ds[counter+1]-Ds[counter])/factor;
            dv = (v[counter+1]-v[counter])/factor;
            ddf = (df0dk(k_grid[counter+1],T,e_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa))/factor;

		df = (f0(energy[counter+1],e_f,T)-f0(energy[counter],e_f,T))/factor;
		if(energy[counter+1]==0)
			df = 0;
				
            for (int i = 0;i<=factor-1;i++)
            {
                g_elastic[counter] = (-1)*E/(h_bar*nu_elastic[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)+i*ddf)*1e-7;
                //// The last number is the conversion from convensional units to cancel out to be unitless (as in g)

		integral_numerator = integral_numerator+de*(Ds[counter]+i*ds)/volume1*(v[counter]+i*dv)*g_elastic[counter]/E;
		// =1/E*int[g_elastic(En)*DOS(En)*v(En)*dEn]
		// units of group velocity (cm/s) and E (V/cm)
		
		integral_denominator = integral_denominator+de*(Ds[counter]+i*ds)/volume1*(f0(energy[counter],e_f,T)+i*df);
		// =int[f0(En)*DOS(En)*dEn]
		// =int[f(En)*DOS(En)*dEn]
            }
        }
    }
    //cout<<"integral_numerator = "<<integral_numerator<<endl;
    //cout<<"integral_denominator = "<<integral_denominator<<endl;

    double mobility_elastic;
    
    if(geometry==1)  // for 3D
    	mobility_elastic = 1/3.0*integral_numerator/integral_denominator;
    else if(geometry==2)   // for 2D
    	mobility_elastic = 1/2.0*integral_numerator/integral_denominator;
    	
    // According to equation (46) in Rode's book; units of (cm^2/V.s)

    return mobility_elastic;
}

