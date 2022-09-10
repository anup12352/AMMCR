
#include"main.h"

double sigma_freq(double e_f,double T,double coefficients[5][7],double kindex[],double fi[][limit9][limit10],int points, int aa[], double energy[], double v[], double Ds[], int ii, int mm)  // ii for freq variation and; mm for temp or doping variation
// It gives the overall mobility in units of (cm^2/V.s)
{
	// According to Equation (46) in Semiconductors and Semimetals volume1 10
	// (Rode's chapter), but with calculated group velocity from band structure and DOS both calculated from DFT:

	double integral_numerator = 0;
	//double integral_denominator = 0;
	double df,k_step,de,dv,ds;
	int factor = 10;


        for (int counter = 0;counter<=points-2;counter++)
        {
		de = (energy[counter+1] - energy[counter]);
		if(energy[counter+1]==0)
			de = 0;
			 
		integral_numerator = integral_numerator + 
		de*(Ds[counter]/volume1)*(v[counter]*v[counter])*fi[counter][ii][mm]*(-1*f0x1_f0[counter])/(k_B*T*e);
		// Ds[] is in desnity of states per unit eV and de is also in eV; 
		// =int[g(En)*DOS(En)*v(En)*dEn]/E
		// units of group velocity (cm/s) 

        }
    
	//cout<<"integral_numerator   =   "<<integral_numerator<<endl;
	//getchar();

	double sig;

	if(geometry==1)  
	{
		sig = e*e*integral_numerator/(3*E);
		// unit of E (V/cm)	  
		// unit S/cm
	}   
	else if(geometry==2)  
	{
		sig = e*e*integral_numerator/(2*E);	  
		// unit of E (V/cm)	  
		// unit S/cm
	}
	return sig;
}
