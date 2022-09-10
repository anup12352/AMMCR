
#include"main.h"

double sigma_tensor(double e_f,double T,double coefficients[5][7],double kindex[],double g[], double nu_el[],int points, int aa[], double energy[], double v[], double Ds[], int jj)
// It gives the overall mobility in units of (cm^2/V.s)
{
    // According to Equation (46) in Semiconductors and Semimetals volume1 10
    // (Rode's chapter), but with calculated group velocity from band structure and DOS both calculated from DFT:

    double integral_numerator = 0;
    //double integral_denominator = 0;
    double df,k_step,de,dv,ds;
    int factor = 10;


	if (T < T_trans)
	{
		double beta1[points]={0}; 	
		// unitless

		for (int counter = 0;counter<points-1;counter++)
		{
			beta1[counter] = e*(v[counter]*0.01)*Bfield/((h_bar*e) * (k_grid[counter]*pow(10,9)) * nu_el[counter]);
			// unitless

			if(jj==1)
			{
				g[counter] = (-1)*E/(h_bar*nu_el[counter] * (1 + beta1[counter]*beta1[counter] ))  *
				 (df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)*1e-7);
			}
			else
			{
				g[counter] = beta1[counter]*E/(h_bar*nu_el[counter] * (1 + beta1[counter] * beta1[counter])) * 					(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)*1e-7);
			}
			// The last number is the conversion from convensional units to cancel out to be unitless (as in g and h)
			//end
		}
	}


	for (int counter = 0;counter<=points-2;counter++)
	{
		de = (energy[counter+1] - energy[counter]);
		if(energy[counter+1]==0)
			de = 0;
			
		integral_numerator = integral_numerator + de*(Ds[counter]/volume1)*v[counter]*g[counter];
		// =int[g(En)*DOS(En)*v(En)*dEn]/E
		// units of group velocity (cm/s) and E (V/cm)
	}
    
    //cout<<"integral_numerator   =   "<<integral_numerator<<endl;
    //getchar();
    
    
    double sig;
    
       
   if(geometry==1)  // 3D
   {
   	sig = e*integral_numerator/(3*E);	  
	// unit S/cm
   }   
   else if(geometry==2)  // for 2D
   {
   	sig = e*integral_numerator/(2*E);	  
	// unit S/cm
   }
   return sig;
}
