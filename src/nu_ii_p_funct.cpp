#include"main.h"


// ionized imprity scattering calculation //---- equation 3.29 of Ramu Thesis
void nu_ii_p_funct(int T_loop)
{
	
	double v;
	double k_dum;

	double A;
	double B;
	double iiB;
	double iiA;
	double epsilon_lf= epsilon_s[T_loop];
	
	double T= T_array[T_loop];
	
	double beta_const;
	double de_p, k_step_p;
    	double integral_1 = 0;
	
	for(int counter = 0;counter <= points-2; counter++)
	    {
		k_dum = k_grid[counter];
		de_p = energy_p[counter+1]-energy_p[counter];
		//cout<<"de_p= "<<de_p<<endl;
		k_step_p = k_grid[counter+1]-k_grid[counter];
		//cout<<"k_step = "<<k_step_p<<endl;
		
		if (free_e ==0)  // DOSCAR is used
		{
		   integral_1 = integral_1+ de_p*((1e6*Ds_p[counter])/volume1)*f0(energy_p[counter],efef_p,T)*(1-f0(energy_p[counter],efef_p,T));
		   // volume unit is cm^3
		   
		   //cout<< "integral doscar = "<<integral_1<<endl;
		   // unit is (1/m)^3; 1e6 is multiplied to convert volume from (1/cm)^3 to (1/m)^3
		}
		else
		{
		    integral_1 = integral_1+k_step_p*(k_dum/pi)*(k_dum/pi)*f0(energy_p[counter],efef_p,T)*(1-f0(energy_p[counter],efef_p,T))*1e27;
		    //cout<<"integral k_step = "<<integral_1<<endl;
		    // unit is (1/m)^3
		    // Part of equation (70) or Rode's book (book8)
		}
		
	     }
	    
    
	double beta = (e*e/(epsilon_s[T_loop]*epsilon_0*(k_B*e)*T)*integral_1);   // unit (1/m)^2		    								
	beta_const = pow(beta,0.5);   // unit (1/m) 
		    
	beta_const = beta_const*1e-9;
	
	//converted from 1/m to 1/nm
		   
	cout<<"screening length is ="<<beta_const<<" 1/nm"<<endl;
		    
	for (int counter = 0;counter<points;counter++)
	   {
		k_dum = k_grid[counter];
		v= v_p[counter]/100;   // divided with 100 to convert from cm/s to m/s
		
		A= 0.5*(2 + (beta_const*beta_const/(k_dum*k_dum)));		
		B= abs((A+1)/(A-1));
				
		//cout<<"A = "<<A<<endl;
		//cout<<"B = "<<B<<endl;
				

		// N_ii is convterd from 1/cm^3 to 1/m^3
		// k-dum is converted from 1/nm to 1/m
		// h_bar is multiplied with e to convert from eV-s to J-s
		iiA= (pow(e,4)* abs(N_ii*1e6))/(32*pi*(k_dum*k_dum*1e18)*epsilon_0*epsilon_0*epsilon_lf*epsilon_lf*v*pow(h_bar*e,2));
		
		iiB = abs((3*A -1)*(3*A -1)*log(B) - 18*A +12 -8/(A+1));
				
		nu_ionizedimpurity_p[counter][0][0]= iiA*iiB  ; 
		nu_ionizedimpurity[counter]= nu_ionizedimpurity_p[counter][0][0]; 
		/*
		cout<<"counter = "<<counter<<endl;	
		cout<<"v = "<<v<<endl;
		cout<<"nu_ionizedimpurity_p[counter][0][0]   =   "<<nu_ionizedimpurity_p[counter][0][0]<<endl;
		getchar();	
		*/		
	     }
				
	// saving results
	/*
	FILE *fid1;
	fid1= fopen("ii.dat","w");
	fprintf(fid1,"#energy              nu_ii_p\n");
	for(int i=0;i<points;i++)
	{
		fprintf(fid1,"%e     %e\n", k_grid[i], nu_ionizedimpurity_p[i][0][0]);
	}
	fclose(fid1);
	//*/	
}
