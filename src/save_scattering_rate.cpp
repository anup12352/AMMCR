#include"main.h"

void save_scattering_rate()
{

	FILE *fid1,*fid2;
	
	double max_scattering_rate = 0;	

//------------------------------------Saving scattering rate--------------------------------------------------------
	if (ispin == 1)
	{			    
	    fid1 = fopen("scattering_rate.dat","w");
	    fid2 = fopen("mean_free_path.dat","w");
	}
	if (ispin == 2 && kk == 0)			    
	{
	    fid1 = fopen("scattering_rate_up_spin.dat","w");
	    fid2 = fopen("mean_free_path_up_spin.dat","w");
	}
	if (ispin == 2 && kk == 1)			    
	{
	   fid1 = fopen("scattering_rate_down_spin.dat","w");
	   fid2 = fopen("mean_free_path_down_spin.dat","w");
	}

	//3D bulk    	
	if(geometry==1)
	{		    
	    if(type=="n")
	    {

        	fprintf(fid1,"# energy           nu_im            nu_pop         nu_npop        nu_de            nu_pz           nu_dis        nu_alloy        nu_iv           nu_nim       nu_irs	nu_skew   nu_total (total relaxation time 1/second)\n");

        	fprintf(fid2,"# energy           Mean free path (nm)\n");

		    double nu_total[points]={0}, inv_nu_total[points]={0}, mfp[points]={0};	
		    
		    for (int i = 0; i < points; i++)
		    {
		        nu_total[i] = scattering_mechanisms[0] * nu_ionizedimpurity[i] + scattering_mechanisms[1] * nu_pop_total[i] + 
		        scattering_mechanisms[2] * nu_npop_total[i] + scattering_mechanisms[3] * nu_deformation[i] +
		        scattering_mechanisms[4] * nu_piezoelectric[i] + 
		        scattering_mechanisms[6] * nu_dislocation[i] + scattering_mechanisms[7] * nu_alloy[i] + 
		        scattering_mechanisms[8] * nu_iv_total[i] + scattering_mechanisms[9] * nu_neutralimpurity[i] 
			+ nu_irs[i]*scattering_mechanisms[11] + nu_skew_rate[i]*scattering_mechanisms[12];
		        
		        if(nu_total[i]!=0)
		        	inv_nu_total[i] = 1.0/nu_total[i];
		        
		        if(nu_total[i] > max_scattering_rate)
		        	max_scattering_rate = nu_total[i];
		        	
		        mfp[i] = v_n[i] * inv_nu_total[i];
		        
		        mfp[i] = mfp[i]/(1e-9);
		        
		        fprintf(fid1,"  %e    %e    %e    %e    %e    %e    %e    %e    %e 	 %e 	%e	%e     %e   \n", energy_n[i], nu_ionizedimpurity[i], 
			nu_pop_total[i], nu_npop_total[i], nu_deformation[i],nu_piezoelectric[i],nu_dislocation[i], 
			nu_alloy[i], nu_iv_total[i], nu_neutralimpurity[i], nu_irs[i], nu_skew_rate[i], nu_total[i] );

		        fprintf(fid2,"  %e    	%e   \n", energy_n[i], mfp[i] );
		    }
		    fclose(fid1);
		    fclose(fid2);
		    
		    cout<<"Maximum scattering rate for RTA = "<<max_scattering_rate<<endl;
	//-------------------------- save data ---------------------------------------------------


		    //cout<<"points = "<<points<<endl;
		    //getchar();
		    /*
		    fid1 = fopen("S_i_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, S_i_grid[i]);
		fclose(fid1);

		    fid1 = fopen("S_i_pop_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, S_i_pop_grid[i]);
		fclose(fid1);


		    fid1 = fopen("S_i_th_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, S_i_th_grid[i]);
		fclose(fid1);

		    fid1 = fopen("S_i_pop_th_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, S_i_pop_th_grid[i]);
		fclose(fid1);

		    fid1 = fopen("S_o_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, S_o_grid[i]);
		fclose(fid1);

		    fid1 = fopen("nu_pop_total.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_pop_total[i]);
		fclose(fid1);
		    */
		}
		else   // for 3D bulk p-type
		{

			fprintf(fid1,"# energy           nu_im            nu_pop         nu_npop        nu_de           nu_total (total relaxation time 1/second)\n");

			fprintf(fid2,"# energy           Mean free path (nm)\n");

			double nu_total[points]={0}, inv_nu_total[points]={0}, mfp[points]={0};	

			for (int i = 0; i < points; i++)
			{
				nu_total[i] = scattering_mechanisms[0] * nu_ionizedimpurity_p[i][0][0] + 
				scattering_mechanisms[1] * nu_So_p[i][0][0] + 
				scattering_mechanisms[2] * nu_npop_p[i][0][0] + scattering_mechanisms[3] * nu_deformation_p[i][0][0];

				if(nu_total[i]!=0)
					inv_nu_total[i] = 1.0/nu_total[i];

				if(nu_total[i] > max_scattering_rate)
					max_scattering_rate = nu_total[i];
					
				mfp[i] = v_n[i] * inv_nu_total[i];

				mfp[i] = mfp[i]/(1e-9);

				fprintf(fid1,"  %e    %e    %e    %e    %e    %e  \n", energy_n[i], nu_ionizedimpurity_p[i][0][0], 
				nu_So_p[i][0][0], nu_npop_p[i][0][0], nu_deformation_p[i][0][0], nu_total[i] );

				fprintf(fid2,"  %e    	%e   \n", energy_n[i], mfp[i] );
			}
			fclose(fid1);
			fclose(fid2);

			cout<<"Maximum scattering rate for RTA = "<<max_scattering_rate<<endl;
		
		} // end of else for type p	
	} // end of geometry==1 for 3D
	else if(geometry==2) // for 2D
	{

		fprintf(fid1,"# energy           nu_rim            nu_pop         nu_npop        nu_de      nu_pz    nu_so_pop   nu_irs 	nu_skew     nu_total (total relaxation time 1/second)\n");

		fprintf(fid2,"# energy           Mean free path (nm)\n");

		double inv_nu_total[points]={0}, mfp[points]={0};	

		for (int i = 0; i < points; i++)
		{
			if(denom[i]!=0)
				inv_nu_total[i] = 1.0/denom[i];

			if(denom[i] > max_scattering_rate)
				max_scattering_rate = denom[i];
				
			mfp[i] = v_n[i] * inv_nu_total[i];

			mfp[i] = mfp[i]/(1e-9);

		        fprintf(fid1,"  %e    %e    %e    %e    %e    %e    %e    %e   %e   %e  \n", energy_n[i], nu_ionizedimpurity[i], 
			nu_pop_total[i], nu_npop_total[i], nu_deformation[i], nu_piezoelectric[i], nu_so_pop_total[i],
			nu_irs[i], nu_skew_rate[i], denom[i] );

			fprintf(fid2,"  %e    	%e   \n", energy_n[i], mfp[i] );
		}
		fclose(fid1);
		fclose(fid2);

		cout<<"Maximum scattering rate for RTA = "<<max_scattering_rate<<endl;
		
	}
	
	// save individual scattering rates
//---------------------------------------------------------------------------------------------------------------
	save_scattering_rate_contribution();
//---------------------------------------------------------------------------------------------------------------
	
}	
