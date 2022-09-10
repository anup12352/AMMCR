#include"main.h"
// NPOP scattering calculation
void nu_npop_p_funct(double T)			
{
	double k_pp;
	double k_mm;
	double k_dum; 
	double arr[points];
	double v;
	double npop1,npop2,npop3,npop4;
	double N_op[limit5];
	double E_npop[limit5];
	double omega_npop[limit5];
	
	//cout<<"c_bar  =   "<<c_bar<<"  dyne/cm^2"<<endl;
	
	double kplus_grid_npop[limit2][npop_number], kminus_grid_npop[limit2][npop_number];
	int plus_index_npop[limit2][npop_number], minus_index_npop[limit2][npop_number];
	int minus_index, plus_index;
	
	for(int loop=0;loop<npop_number;loop++)
	{
		//cout<<"loop   =  "<<loop<<endl;
		E_npop[loop] = De_npop[loop]*e;    // e is multiplied to convert from ev to J  
		//cout<<"E_npop[loop]   = " <<E_npop[loop] <<endl;
		omega_npop[loop] = we_npop[loop]; 
		//cout<<"omega_npop[loop]   = " <<omega_npop[loop] <<endl;
		
		N_op[loop] = N_poph(omega_npop[loop],T);
		//cout<<"N_op[loop]   = " <<N_op[loop]<<endl;
	}
	
	
	for (int counter = 0;counter < points;counter++)
	{	
		for(int loop=0;loop<npop_number;loop++)
		{
			
			kplus_grid_npop[counter][loop] = kplus(counter,omega_npop[loop],points,energy_p);
			kminus_grid_npop[counter][loop] = kminus(counter,omega_npop[loop],points,energy_p);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kminus_grid_npop[counter][loop]);
			minus_index =FindMinInd(arr,points);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kplus_grid_npop[counter][loop]);
			plus_index =FindMinInd(arr,points);

			plus_index_npop[counter][loop] = plus_index;
			minus_index_npop[counter][loop] = minus_index;
		}
	}
	
		
	for (int i = 0;i<points;i++)
		nu_npop_total[i] = 0; 	
	
	for(int counter=0; counter<points;counter++)
	{	
		k_dum=k_grid[counter]*1e9;   // converted from 1/nm to 1/m --unit is 1/m
		v=v_p[counter]*1e-2;    // conveterd from unit  cm/s to m/s
		
		for(int loop=0;loop<npop_number;loop++)
		{	
			
			//cout<<"counter  = "<<counter<<endl;
			k_mm= kminus_grid_npop[counter][loop]*1e9;   // converted from 1/nm to 1/m --unit is 1/m
			k_pp= kplus_grid_npop[counter][loop]*1e9;   // converted from 1/nm to 1/m --unit is 1/m
			
			npop1= E_npop[loop] *E_npop[loop] *omega_npop[loop]*k_dum*k_pp/(2*pi*(h_bar*e)*(c_bar/10)*v);
			//cout<<"npop1 is = "<<npop1<<endl;
			npop2=  E_npop[loop]*E_npop[loop] *omega_npop[loop] *k_dum*k_mm/(2*pi*(h_bar*e)*(c_bar/10)*v);
			//cout<<"npop2 is = "<<npop2<<endl;

			npop3= (N_op[loop]*(1-f0(energy_p[plus_index_npop[counter][loop]],efef_p,T)) + (N_op[loop]+1)*f0(energy_p[plus_index_npop[counter][loop]],efef_p,T));
			npop4=	(N_op[loop]*f0(energy_p[minus_index_npop[counter][loop]],efef_p,T) + (N_op[loop]+1)*(1-f0(energy_p[minus_index_npop[counter][loop]],efef_p,T)));
			//cout<<"npop3 is = "<<npop3<<endl;
			//cout<<"npop4 is = "<<npop4<<endl;
			
			
			    if (energy_p[counter] < h_bar*omega_npop[loop])
			    {
				    nu_npop[loop][counter] = npop1*npop3;
			    }
			    else
			    {
			        
				nu_npop[loop][counter] = npop1*npop3 +npop2*npop4;
			    }
			    
			    
		 	nu_npop_total[counter] = nu_npop_total[counter] + nu_npop[loop][counter];  
		 	nu_npop_p[counter][0][0]= nu_npop_total[counter]; 
		 }	
		 
	 	/*
		cout<<"counter = "<<counter<<endl;	
		cout<<"v = "<<v<<endl;
		cout<<"nu_npop_p[counter][0][0] = "<<nu_npop_p[counter][0][0]<<endl;
	 	getchar();
	 	*/
	 }	
		
	/*	
	FILE *fid1;
	fid1 = fopen("npop_scattering_rate.dat","w");
	fprintf(fid1,"# energy           total	individual cattering rates \n");

	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e          %e  \t ", energy_p[i], 	nu_npop_p[i][0][0] );
		for(int j=0;j<npop_number;j++)
			fprintf(fid1,"%e \t", nu_npop_p_all[i][j]);
		
		fprintf(fid1,"\n "); 
	}
	fclose(fid1);
	//*/
	
			
}
