#include"main.h"

void components_BTE(double T, int T_loop, double efefn, double efefp, int ii)   // ii for doping variation loop
{

	double integral_numerator = 0;
	double integral_denominator = 0;
	double k_dum;
	
//---------------------------- components for BTE ------------------------------------------------------------------
	if(geometry==1 && type == "n")
	{	
		    beta_constant = beta(T, T_loop);
		    // unit 1/nm

		    cout<< "Inverse screening length, beta = "<<beta_constant<<" (1/nm)"<<endl;

		    //--------------------------- df0dz_integral -----------------------------------------
		    int factr = 10;
		    for(int counter = 0;counter<=points-2;counter++)
		    {
		        double dk = (k_grid[counter+1]-k_grid[counter])/factr;
		        for (int ss = 0;ss<=factr-1;ss++)
		        {
		            integral_numerator = integral_numerator + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
		                                *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T))*energy_n[counter]/(k_B*T);
		            // Part of equation (54) of Rode's book
		            integral_denominator = integral_denominator + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
		            *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));
		            // Part of equation (54) of Rode's book
		        }
		    }

		    df0dz_integral = integral_numerator/integral_denominator;
		    //cout<<"df0dz_integral for n type = "<<df0dz_integral<<endl;
		    //--------------------------- df0dz_integral calculated -----------------------------------------
		    		    
//--------------------------------------common terms calculated -------------------------------------------------------

		for (int counter = 0;counter<points;counter++)
		{
			k_dum = k_grid[counter];
			// unit 1/nm
			//cout<<"counter+1 = "<<counter+1<<endl;
			//cout<<"k_dum = "<<k_dum<<endl;

		        df0dk_grid[counter] = df0dk(k_dum, T, efefn, coefficients_cond, kindex_cond, a11);
			 // unit (nm)

		        f_dist[counter] = f0(energy_n[counter],efefn,T);
		        //cout<<"In between "<<endl;
		        //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
		        //cout<<"efefn = "<<efefn<<endl;
		        //cout<<"T = "<<T<<endl;

		        thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, efefn, T, df0dz_integral,coefficients_cond, kindex_cond, a11);

		        f0x1_f0[counter] = f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));

		        electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
			// unit is 1/s , hbar unit is eV-s so e i s not multiplied in numerator

		        //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
		        //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
		        //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
		        //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
		        //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;
		}

//--------------------------------------common terms calculated completed  --------------------------------------------

//----------------------POP scattering rate --------------------------------------------------------------------
			    
		// pop scattering
		if(scattering_mechanisms[1]==1)
			nu_pop(T, efefn, ii, T_loop);
//-------------------------------------------------------- nu_pop calculated---------------------------------------

//----------------------POP scattering rate completed--------------------------------------------------------------------


//----------------------Ionized Impurity scattering rate --------------------------------------------------------------------
			        			
		// ionized impourity scattering
	        if (scattering_mechanisms[0]==1)
	        {
	        	nu_ii(epsilon_s[T_loop]);
	        }
//----------------------II scattering rate completed--------------------------------------------------------------------

//----------------------Acoustic scattering rate --------------------------------------------------------------------
		        
		// acoustic scattering 
	        if (scattering_mechanisms[3]==1)
	        {
            		nu_de(T);
		}	
//----------------------Acoustic scattering rate completed--------------------------------------------------------------------

//----------------------PZ scattering rate --------------------------------------------------------------------
		        
		// Piezoelectric scattering
	        if (scattering_mechanisms[4]==1)
	        {
         		nu_pe(T,P_piezo[T_loop],epsilon_s[T_loop]);
		}
//----------------------PZ scattering rate --------------------------------------------------------------------
			
//----------------------Disloaction scattering rate -------------------------------------------------------			
		// Dislocation scattering
	        if (scattering_mechanisms[6]==1)
	        {
			nu_dis(T,beta_constant,epsilon_s[T_loop]);
		}
//----------------------Disloaction scattering rate completed ------------------------------------------			

//----------------------Alloy scattering rate ----------------------------------------------------			
			
		// Alloy scattering
	        if (scattering_mechanisms[7]==1)
	        {
		        nu_alloy1();
		}
//----------------------Alloy scattering rate completed ------------------------------------------			

			
//---------------------- Neutral impurity scattering rate  ---- ------------------------------------------			

		// neutral impurity
	        if (scattering_mechanisms[9]==1)  
		{						
			nu_im(epsilon_s[T_loop],N_im[ii]);
		}
			
//---------------------- Neutral impurity scattering rate  completed ------------------------------------------			


//---------------------- Intervalley scattering rate  ---- ------------------------------------------			

		// intravalley scattering	
		if (scattering_mechanisms[8]==1)   
		{
			nu_iv_n(T);
		}
//------------------interavalley scattering rate completed ----------------------------
		
//----------------------  Interface roughness scattering ----------------------------------------------			

	 	// Interface roughness scattering
		if (scattering_mechanisms[11]==1)  
		{
			nu_ir(epsilon_s[T_loop]);						
		}
		
		
//--------------------------------------Interface roughness completed -------------------------------------------------------



		// total scattering rates calculated here
		for(int counter=0;counter<points;counter++)
		{
		        nu_el[counter] = nu_ionizedimpurity[counter]*scattering_mechanisms[0] + 
		        	nu_npop_total[counter] * scattering_mechanisms[2] + 
		        	nu_deformation[counter]*scattering_mechanisms[3] + 
		               nu_piezoelectric[counter]*scattering_mechanisms[4] + 
		         	nu_dislocation[counter] * scattering_mechanisms[6] + 
		         	nu_alloy[counter] * scattering_mechanisms[7] + 
		         	nu_iv_total[counter] * scattering_mechanisms[8] +
		         	nu_neutralimpurity[counter]*scattering_mechanisms[9] + 
		         	nu_irs[counter]*scattering_mechanisms[11];

			denom[counter] = (nu_pop_total[counter]*scattering_mechanisms[1] + nu_el[counter]);	
			//cout<<"nu_el[counter] = "<<nu_el[counter]<<endl;
		}			

//--------------------------------------Skew scattering -------------------------------------------------------

//----------------------  Skew scattering ----------------------------------------------			

	 	// Skew scattering
		if (scattering_mechanisms[12]==1)  
		{
			nu_skew();						
		}
		
//--------------------------------------skew scattering completed -------------------------------------------------------

	 	
//-----------------------------------------------------------------------------------------------------------------

				
	//------------------------------------ components for BTE END --------------------------------------------------------------
	// ----------------------------saving data ---------------------------------------------------
		    
		    /*		    
		    //FILE *fid1;
		    fid1 = fopen("nu_el.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_el[i]);
		fclose(fid1);

		    fid1 = fopen("df0dk_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, df0dk_grid[i]);
		fclose(fid1);

		    fid1 = fopen("f_dist.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, f_dist[i]);
		fclose(fid1);

		    fid1 = fopen("thermal_driving_force.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, thermal_driving_force[i]);
		fclose(fid1);

		    fid1 = fopen("f0x1_f0.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, f0x1_f0[i]);
		fclose(fid1);

		    fid1 = fopen("elecgtric_driving_force.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, electric_driving_force[i]);
		fclose(fid1);
		    */
		    //getchar();

	// ----------------------------saving data ---------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------
	}  // if condition for type 'n'
	else if(geometry==1 && type == "p")
	{
				
		//--------------------------- df0dz_integral -----------------------------------------
		int factr = 10;
		for(int counter = 0;counter<=points-2;counter++)
		{
			double dk = (k_grid[counter+1]-k_grid[counter])/factr;
			for (int ss = 0;ss<=factr-1;ss++)
			{
			    integral_numerator = integral_numerator + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
						*f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T))*energy_p[counter]/(k_B*T);
			    // Part of equation (54) of Rode's book
			    integral_denominator = integral_denominator + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
			    *f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T));
			    // Part of equation (54) of Rode's book
			}
		}

		df0dz_integral = integral_numerator/integral_denominator;
		//cout<<"df0dz_integral for p type = "<<df0dz_integral<<endl;
		//--------------------------- df0dz_integral calculated -----------------------------------------		    		    
//--------------------------------------common terms calculated -------------------------------------------------------

		for (int counter = 0;counter<points;counter++)
		{
			k_dum = k_grid[counter];
			// unit 1/nm
			//cout<<"counter+1 = "<<counter+1<<endl;
			//cout<<"k_dum = "<<k_dum<<endl;

		        df0dk_grid[counter] = df0dk(k_dum, T, efefp, coefficients_val, kindex_val, b11);
			 // unit (nm)

		        f_dist[counter] = f0(energy_p[counter],efefp,T);
		        //cout<<"In between "<<endl;
		        //cout<<"energy_p[counter]  = "<<energy_p[counter]<<endl;
		        //cout<<"efefp = "<<efefp<<endl;
		        //cout<<"T = "<<T<<endl;

		        thermal_driving_force[counter] = -1*v_p[counter]*df0dz(k_dum, efefp, T, df0dz_integral,coefficients_val, kindex_val, b11);

		        f0x1_f0[counter] = f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T));

		        electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
			// unit is 1/s , hbar unit is eV-s so e is not multiplied in numerator

		        //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
		        //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
		        //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
		        //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
		        //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;
		}
					
	        if (scattering_mechanisms[0]==1)
	        {
	        	nu_ii_p_funct(T_loop);
		}
		
		// POP scattering
	        if (scattering_mechanisms[1]==1)
		{
			nu_pop_p_funct(T, T_loop);
		}
			
		// npop scattering
		if (scattering_mechanisms[2]==1)
		{
			nu_npop_p_funct(T);
		}

		// acoustic scattering
		if (scattering_mechanisms[3]==1)
		{
			nu_de_p_funct(T_loop);
		}


		//*/

		// total scattering rates calculated here
		for(int counter=0;counter<points;counter++)
		{
		        nu_el[counter] = nu_ionizedimpurity_p[counter][0][0]*scattering_mechanisms[0] + 
		        	nu_npop_p[counter][0][0] * scattering_mechanisms[2] + 
		        	nu_deformation_p[counter][0][0]*scattering_mechanisms[3]; 
		       
			denom[counter] = (nu_So_p[counter][0][0]*scattering_mechanisms[1] + nu_el[counter]);	
			//cout<<"nu_el[counter] = "<<nu_el[counter]<<endl;
		}
					
	}
	else if(geometry==2)
	{
//--------------------------------------common terms calculated -------------------------------------------------------
		
		//cout<<"Components  result for 2D "<<endl;
		if(screening==1)
			polarizability(T, ii);
		
		for (int counter = 0;counter<points;counter++)
		{
			k_dum = k_grid[counter];
			// unit 1/nm
			//cout<<"counter+1 = "<<counter+1<<endl;
			//cout<<"k_dum = "<<k_dum<<endl;

		        df0dk_grid[counter] = df0dk(k_dum, T, efefn, coefficients_cond, kindex_cond, a11);
			 // unit (nm)

		        f_dist[counter] = f0(energy_n[counter],efefn,T);
		        //cout<<"In between "<<endl;
		        //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
		        //cout<<"efefn = "<<efefn<<endl;
		        //cout<<"T = "<<T<<endl;


		        f0x1_f0[counter] = f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));

		        electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
			// unit is 1/s , hbar unit is eV-s so e i s not multiplied in numerator
			
			
		        //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
		        //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
		        //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
		        //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
		        //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;
		}

//--------------------------------------common terms calculated completed  --------------------------------------------
		
		// remote impurity scattering
		if (scattering_mechanisms[0]==1)
		{
			nu_rim_2D(T);
		}

		// pop scattering
		if (scattering_mechanisms[1]==1)
		{
			nu_pop_2D(T, T_loop);
		}

		// npop scattering
		if (scattering_mechanisms[2]==1)
		{
			nu_npop_2D(T);
		}

		// acoustic scattering
		if (scattering_mechanisms[3]==1)
		{
			nu_de_2D(T);
		}
		
		
		// so pop scattering
		if (scattering_mechanisms[10]==1)
		{
			nu_so_pop_2D(T, T_loop);
		}

//----------------------  Interface roughness scattering ----------------------------------------------			

	 	// Interface roughness scattering
		if (scattering_mechanisms[11]==1)  
		{
			nu_ir(epsilon_s[T_loop]);						
		}		
		
//--------------------------------------Interface roughness completed -------------------------------------------------------

//----------------------  Skew scattering ----------------------------------------------			

	 	// Skew scattering
		if (scattering_mechanisms[12]==1)  
		{
			nu_skew();						
		}
		
//--------------------------------------skew scattering completed -------------------------------------------------------

		// total scattering rates calculated here
		for(int i=0;i<points;i++)
		{
		        nu_el[i] = nu_ionizedimpurity[i]*scattering_mechanisms[0] + 
		        	nu_npop_total[i] * scattering_mechanisms[2] + 
		        	nu_deformation[i]*scattering_mechanisms[3] + 
		               nu_piezoelectric[i]*scattering_mechanisms[4] +
				nu_irs[i]*scattering_mechanisms[11];

			denom[i] = (nu_pop_total[i]*scattering_mechanisms[1] + nu_el[i] 
			+ nu_so_pop_total[i]*scattering_mechanisms[10]);				
			//cout<<"nu_el[i] = "<<nu_el[i]<<endl;
			//cout<<"denom[i] = "<<denom[i]<<endl;
			
		}			

		//--------------------------- df0dz_integral -----------------------------------------
		int factr = 10;
		integral_numerator = 0;
		integral_denominator = 0;
		/*
		double de1;
		for(int counter = 0;counter<=points-2;counter++)
		{
			de1 = energy_n[counter+1]- energy_n[counter];
			if(de1<0)
				de1 = 0;			
			integral_numerator = integral_numerator + de1*(f0(energy_n[counter],efefn,T)*
			(1-f0(energy_n[counter],efefn,T)))*energy_n[counter]/(k_B*T);
			// Part of equation (54) of Rode's book

			integral_denominator = integral_denominator + de1*f0(energy_n[counter],efefn,T)*
			(1-f0(energy_n[counter],efefn,T));
			// Part of equation (54) of Rode's book
		}
		*/		
		double dk;
		for(int counter = 0;counter<=points-2;counter++)
		{
			dk = (k_grid[counter+1]-k_grid[counter])/factr;
			for (int ss = 0;ss<=factr-1;ss++)
			{
				integral_numerator = integral_numerator + dk*(((k_grid[counter]+ss*dk)/pi))
						*f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T))*energy_n[counter]/(k_B*T);
				// Part of equation (54) of Rode's book
				integral_denominator = integral_denominator + dk*((k_grid[counter]+ss*dk)/pi)
				*f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));
				// Part of equation (54) of Rode's book
			}
		}

		df0dz_integral = integral_numerator/integral_denominator;
		//cout<<"df0dz_integral for 2D geometry = "<<df0dz_integral<<endl;
		//--------------------------- df0dz_integral calculated -----------------------------------------
		
		//------------------- thermal driving force ----------------------------------------------------
		for (int counter = 0;counter<points;counter++)
		{
		        //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
		        //cout<<"efefn = "<<efefn<<endl;
		        //cout<<"T = "<<T<<endl;

		        thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, efefn, T, df0dz_integral,coefficients_cond,
		        kindex_cond, a11);           
		}		
		//--------------------------------------------------------------------------------------------
	}  // else if geometry==2;  2D material
	
	
	
	// Save componeents of scattering rate for diiferent scattering mechanism
		
	// save results
	/*
	FILE *fid1;
	fid1 = fopen("electric_driving_force.txt","w");
	
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", electric_driving_force[i]);	
	fclose(fid1);


	fid1 = fopen("df0dk.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", df0dk_grid[i]);	
	fclose(fid1);

	fid1 = fopen("thermal_driving_force.txt","w");		
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", thermal_driving_force[i]);	
	fclose(fid1);
	
	//*/
	
}


