#include"main.h"


// POP scattering calculation
void nu_pop_p_funct(double T, int T_loop)
{	
	double k_pp;
	double k_mm;
	
	double c_minus[pop_number][points]={0};
	double c_plus[pop_number][points]={0};
	double B_minus[points]={0};
	double B_plus[points]={0};
	
	
	double k_dum; 
	double v;
	
	double arr[points];
	double pop1[points]={0},pop2[points]={0},pop3[points]={0};


	int minus_index, plus_index;

	double N_op[pop_number];	
//----------------------out scattering terms are calulated here ----------------------------------------------------------------
	for (int counter = 0;counter < points;counter++)
		nu_So_p[counter][0][0] = 0;
	
	
	for (int m3 = 0;m3 < pop_number;m3++)
	{	            
		N_op[m3] = N_poph(we_pop[m3],T);  // number of optical phonon at temperature T
		for (int counter = 0;counter < points;counter++)
		{
			kplus_grid_pop[m3][counter] = kplus(counter,we_pop[m3],points,energy_p);
			kminus_grid_pop[m3][counter] = kminus(counter,we_pop[m3],points,energy_p);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kminus_grid_pop[m3][counter]);
			minus_index =FindMinInd(arr,points);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kplus_grid_pop[m3][counter]);
			plus_index =FindMinInd(arr,points);

			plus_index_pop[m3][counter] = plus_index;
			minus_index_pop[m3][counter] = minus_index;

			//cout<<"counter = "<<counter<<endl;
			//cout<<"plus_index_pop[m3][counter1] = "<<plus_index_pop[m3][counter1]<<endl;
			//cout<<"minus_index_pop[m3][counter1] = "<<minus_index_pop[m3][counter1]<<endl;	
			//cout<<"kplus_grid_pop[m3][counter1] = "<<kplus_grid_pop[m3][counter1]<<endl;
			//cout<<"kminus_grid_pop[m3][counter1] = "<<kminus_grid_pop[m3][counter1]<<endl;
			//getchar();	
		}
	}		
	
	for (int m3 = 0;m3 < pop_number;m3++)
	{	            
		for(int counter=0;counter<points;counter++)
		{	
			k_dum=k_grid[counter];   // unit is 1/nm
			v=v_p[counter]*1e-2;    // conveterd from unit  cm/s to m/s
			k_mm= kminus_grid_pop[m3][counter];  // unit is 1/nm
			k_pp= kplus_grid_pop[m3][counter];   // unit is 1/nm
			c_minus[m3][counter]= (k_dum*k_dum+k_mm*k_mm)/(2*k_mm*k_dum);	// c_minus[m3][counter] from equation 3.22 of ramu thesis
	 		c_plus[m3][counter]= (k_dum*k_dum+k_pp*k_pp)/(2*k_pp*k_dum);	// c_minus[m3][counter] from equation 3.22 of ramu thesis	
		 	
		 	/*
			if((energy_p[counter] < h_bar*we_pop[m3])||(k_mm==k_dum))
				B_minus[counter] = 0;
			else
			//*/	 	 		
				B_minus[counter]= abs(((1+3.0*c_minus[m3][counter]*c_minus[m3][counter])/2.0)*log(abs((1+c_minus[m3][counter])/(1.0-c_minus[m3][counter]))) - 3.0*c_minus[m3][counter]);		
			//from equation 3.22 of ramu thesis
			/*
			if ((k_pp == k_dum))
				B_plus[counter] = 0;
			else
			//*/
				B_plus[counter] = abs(((1+3.0*c_plus[m3][counter]*c_plus[m3][counter])/2.0)*log(abs((1+c_plus[m3][counter])/(1-c_plus[m3][counter]))) - 3*c_plus[m3][counter]);		
			//from equation 3.22 of ramu thesis
			 
			pop1[counter]= (e*e*we_pop[m3]/(16*pi*h_bar*v*epsilon_0))*(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop]);
			pop1[counter]= pop1[counter]/e; // for unit conversion of h_bar
			pop2[counter]= B_plus[counter]*(N_op[m3]*(1-f0(energy_p[plus_index_pop[m3][counter]],efef_p,T)) + 
			(N_op[m3]+1)*f0(energy_p[plus_index_pop[m3][counter]],efef_p,T));
			pop3[counter]= B_minus[counter]*(N_op[m3]*f0(energy_p[minus_index_pop[m3][counter]],efef_p,T) + 
			(N_op[m3]+1)*(1-f0(energy_p[minus_index_pop[m3][counter]],efef_p,T)));
			
			
			// out scattering terms
			So_ab_pop[m3][counter] = pop1[counter]*(pop2[counter]);
			So_em_pop[m3][counter] = pop1[counter]*(pop3[counter]);
			
			nu_So_p[counter][0][0] = nu_So_p[counter][0][0] + So_ab_pop[m3][counter] + So_em_pop[m3][counter];
			
			So_pop[m3][counter] = So_ab_pop[m3][counter] + So_em_pop[m3][counter];
			
			/*
			cout<<"counter = "<<counter<<endl;	
			cout<<"v = "<<v<<endl;
			cout<<"nu_So_p[counter][0][0]   =   "<<nu_So_p[counter][0][0]<<endl;
			getchar();	
			*/	 
		}
	}

	for(int counter=0;counter<points;counter++)
	{	
		// Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
		if(nu_So_p[counter][0][0] == 0)  
			nu_So_p[counter][0][0] = 1; 
		
		
		// Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero	
		for (int m3 = 0;m3 < pop_number;m3++)
		{	            	
			if(So_pop[m3][counter] == 0) 	
				So_pop[m3][counter] = 1;	
		}
				
		nu_pop_total[counter] = nu_So_p[counter][0][0];
	}
	
	//------------------------- out scattering term for POP scattering rate calculated  --------------------------------
			
	FILE *fid1;
	
		
	/*	

	fid1 = fopen("B_plus.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e   \n", B_plus[i]);
	}
	fclose(fid1);

	fid1 = fopen("B_minus.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e     \n", B_minus[i]);
	}
	fclose(fid1);

	fid1 = fopen("c_plus.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e   \n", c_plus[i]);
	}
	fclose(fid1);


	fid1 = fopen("c_minus.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e     \n", c_minus[i]);
	}
	fclose(fid1);

	fid1 = fopen("pop1.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e         \n", pop1[i]);
	}
	fclose(fid1);


	fid1 = fopen("pop2.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e       \n", pop2[i]);
	}
	fclose(fid1);
	
	fid1 = fopen("pop3.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e       \n", pop3[i]);
	}
	fclose(fid1);

	//*/
	
	/*
	fid1 = fopen("pop.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e               %e  \n", k_grid[i], nu_So_p[i][0][0]);
	}
	fclose(fid1);
	//*/


//----------------------in scattering terms are calulated here ----------------------------------------------------------------
	for (int m3 = 0;m3 < pop_number;m3++)
	{	            		    
		// polar optical phonon scattering 
		N_poph_atT[m3] = N_poph(we_pop[m3],T);
		//cout<<"N_poph_atT[m3] = "<<N_poph_atT[m3]<<endl;
					
	//------------------------- In scattering term for POP scattering rate calculation  --------------------------------
		
		double const1, C_plus[points]={0}, C_minus[points]={0};
		double A, B, C, const2[points]={0}, const3;

		// Eq no. 20 of paper coupled band Ramu paper
		const1 = e*e*we_pop[m3]*(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop])*(1/epsilon_0)/(16*pi*(h_bar*e));

		for (int counter = 0;counter < points;counter++)
		{    
		    const2[counter]= const1/(v_p[counter]/100);		
		    
		    const3 = const2[counter];
		    
		    A = abs((1.0+c_plus[m3][counter])/(1.0-c_plus[m3][counter]));
		    B = (c_plus[m3][counter]+3.0*c_plus[m3][counter]*c_plus[m3][counter]*c_plus[m3][counter])/2.0;
		    C = 2.0 + 3.0 * (c_plus[m3][counter]*c_plus[m3][counter]);
		    
		    C_plus[counter] = abs(B*log(A) - C);    //
		     
		    A = abs((1.0+c_minus[m3][counter])/(1.0-c_minus[m3][counter]));
		    B = (c_minus[m3][counter]+3.0*c_minus[m3][counter]*c_minus[m3][counter]*c_minus[m3][counter])/2.0;
		    C = 2.0 + 3.0 * (c_minus[m3][counter]*c_minus[m3][counter]);

		    C_minus[counter] = abs(B*log(A) - C);   // 

		    lambda_i_plus_grid[m3][counter] = abs(const3*C_plus[counter]*((N_poph_atT[m3]+1)*(1 - f0(energy_p[counter],efef_p,T))
		    + (N_poph_atT[m3])*f0(energy_p[counter],efef_p,T)));

		    lambda_i_minus_grid[m3][counter] = abs(const3*C_minus[counter]*((N_poph_atT[m3])*(1 - f0(energy_p[counter],efef_p,T))
		    + (N_poph_atT[m3]+1)*f0(energy_p[counter],efef_p,T)));
		    
		    if(energy_p[counter]==0)
		    {
		    	lambda_i_minus_grid[m3][counter] = 0;
		    	lambda_i_plus_grid[m3][counter] = 0;
		    }
		    
		    // in scattering terms
		    Sa_pop[m3][counter] = lambda_i_minus_grid[m3][counter];
		    Se_pop[m3][counter] = lambda_i_plus_grid[m3][counter];	
		    
			Si_pop[m3][counter] = Sa_pop[m3][counter] + Se_pop[m3][counter];

			//---------------------------- code to debug -------------------------------------------------------------
			//cout<<"counter = "<<counter<<endl;

			//cout<<"lambda_i_plus_grid[counter] =  "<<lambda_i_plus_grid[counter]<<endl;
			//cout<<"lambda_i_minus_grid[counter] =  "<<lambda_i_minus_grid[counter]<<endl;
			//getchar();
		}
	} // 
	
	//------------------------- In scattering term for POP scattering rate calculated  --------------------------------

	/*
	fid1 = fopen("pop_scattering_rate.dat","w");	
	
	fprintf(fid1,"Energy (eV)   ab   em   total   \n ");
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		fprintf(fid1,"  %e	%e	%e  \n", So_ab_pop[m3][i], So_em_pop[m3][i], nu_So_p[i][0][0] );
	}		
	fclose(fid1);
	
	
	fid1 = fopen("pop_in_scattering_rate.dat","w");

	fprintf(fid1,"Energy (eV)   ab1   em1   \n");	
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		fprintf(fid1,"  %e   %e    \n", Sa_pop[m3][i], Se_pop[m3][i] );
	}
		
	fclose(fid1);
	
	
	//*/
	/*		
	fid1 = fopen("const2.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", const2[i]);	
	fclose(fid1);

	fid1 = fopen("lambda_i_plus_p.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", lambda_i_plus_grid[i]);	
	fclose(fid1);

	fid1 = fopen("lambda_i_minus_p.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", lambda_i_minus_grid[i]);	
	fclose(fid1);

	fid1 = fopen("C_plus_in_p.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", C_plus[i]);	
	fclose(fid1);

	fid1 = fopen("C_minus_in_p.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", C_minus[i]);	
	fclose(fid1);

	fid1 = fopen("k_minus.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", kminus_grid_pop[m3][i]);	
	fclose(fid1);

	fid1 = fopen("k_plus.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", kplus_grid_pop[m3][i]);	
	fclose(fid1);

	fid1 = fopen("plus_index.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%d \n", plus_index_pop[m3][i]);	
	fclose(fid1);

	fid1 = fopen("minus_index.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%d \n", minus_index_pop[m3][i]);	
	fclose(fid1);

	fid1 = fopen("denom.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", denom[i]);	
	fclose(fid1);
	*/
}
