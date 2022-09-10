
#include"main.h"

double mu_overall(double e_f,double T,double coefficients[5][7],double kindex[], double g[],double nu_el[],int points,int aa[], double energy[], double v[], double Ds[])
// It gives the overall mobility in units of (cm^2/V.s)
{
    // According to Equation (46) in Semiconductors and Semimetals volume1 10
    // (Rode's chapter), but with calculated group velocity from band structure and DOS both calculated from DFT:
    double integral_numerator = 0;
    double integral_denominator = 0;
    double ddf,k_step,de,dv,ds,df;
    int factor = 1;


    if (T < T_trans)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            g[counter] = (-1)*E/(h_bar*nu_el[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)*1e-7);
            // hbar unit is eV-s so e in numerator is removed
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

            df = (f0(energy[counter+1],e_f,T) - f0(energy[counter],e_f,T))/factor;
	    if(energy[counter+1]==0)
	    	df = 0;
	    	 	
            for (int i = 0;i<=factor-1;i++)
            {
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v[counter]+i*dv)*g[counter]/E;
                        // =1/E*int[g(En)*DOS(En)*v(En)*dEn]
                integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy[counter],e_f,T)+i*df);
                    // =int[f(En)*DOS(En)*$
            }
        }
    }
    else
    {
        for (int counter = 0;counter<=points-2;counter++)
        {
		de = (energy[counter+1] - energy[counter]);
		if(energy[counter+1]==0)
			de = 0;
			
		integral_numerator = integral_numerator+de*(Ds[counter]/volume1)*v[counter]*g[counter]/E;
		// =1/E*int[g(En)*DOS(En)*v(En)*dEn]
		// units of group velocity (cm/s) and E (V/cm)

		integral_denominator = integral_denominator + de*(Ds[counter]/volume1)*f0(energy[counter],e_f,T);
		// =int[f(En)*DOS(En)*dEn]
		// =int[f(En)*DOS(En)*dEn]
		//for table form volume1 = 1
        }
    }
	 	
	/*
	cout<<endl<<"integral_numerator   =   "<<integral_numerator<<" 1/cm^3"<<endl;
	if(geometry==1)
		cout<<"integral_denominator   =   "<<integral_denominator<<" 1/cm^3"<<endl;
	else	
		cout<<"integral_denominator*thickness = "<<integral_denominator*thickness*100<<" 1/cm^2"<<endl;

	getchar();
	*/
	
	/*
	cout<<"mu_overall"<<endl;
	cout<<endl<<"integral_numerator   =   "<<integral_numerator<<" 1/cm^3"<<endl;
	cout<<"integral_denominator   =   "<<integral_denominator<<" 1/cm^3"<<endl;
	cout<<"integral_denominator*thickness*100   =   "<<integral_denominator*thickness*100<<" 1/cm^2"<<endl;
	getchar();
	*/
	
    double mobility_overall;

	
    if(geometry==1)  // for 3D   
    	mobility_overall = (1/3.0)*integral_numerator/integral_denominator;
    else if(geometry==2)   // for 2D		
	    mobility_overall = (1/2.0)*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)

    return mobility_overall;
}
