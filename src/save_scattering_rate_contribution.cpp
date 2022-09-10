#include"main.h"

void save_scattering_rate_contribution()
{			
	FILE *fid1;
	
	//cout<<"Reached here save_scattering_rate_contribution "<<endl;
	
	// pop scattering
	if(scattering_mechanisms[1]==1 && pop_number > 1 )
	{
	
		fid1 = fopen("pop_out_scattering_rate.dat","w");
			fprintf(fid1,"# energy           total	individual out scattering rates \n");

		for (int i = 0; i < points; i++)		
		{
			fprintf(fid1,"  %e    	%e   \t", energy_n[i], nu_pop_total[i]);
			
			for(int m3=0; m3<pop_number; m3++)
				fprintf(fid1,"  %e  \t", So_pop[m3][i] );
			
			fprintf(fid1,"  \n");

		}				
		fclose(fid1);
			
		fid1 = fopen("pop_in_scattering_rate.dat","w");
			fprintf(fid1,"# energy           individual in scattering rates \n");

		for (int i = 0; i < points; i++)		
		{
			fprintf(fid1,"  %e  \t", energy_n[i]);
			
			for(int m3=0; m3<pop_number; m3++)
				fprintf(fid1,"  %e  \t", Si_pop[m3][i] );
			
			fprintf(fid1,"  \n");

		}				
			
		fclose(fid1);	
			
	}	


	// npop scattering
	if (scattering_mechanisms[2]==1 && npop_number > 1)
	{
		fid1 = fopen("npop_scattering_rate.dat","w");
		fprintf(fid1,"# energy           total	individual scattering rates \n");

			for (int i = 0; i < points; i++)		
			{
				fprintf(fid1,"  %e    	%e   \t", energy_n[i], nu_npop_total[i]);
				
				for(int m3=0; m3<npop_number; m3++)
					fprintf(fid1,"  %e  \t", nu_npop[m3][i] );
				
				fprintf(fid1,"  \n");

			}
		fclose(fid1);	
	}



	// For acoustic scattering
	if(scattering_mechanisms[3]==1 && de_number > 1)   //  
	{
		fid1 = fopen("acoustic_scattering_rate.dat","w");
		
		if(de_number==2)
			fprintf(fid1,"# energy           	nu_LA 			nu_TA			total \n");
		else
			fprintf(fid1,"# energy           	nu_LA 			save_scattering_rate_contributionnu_TA		nu_ZA 			total  \n");
		
		if(de_number==2)
		{
			for (int i = 0; i < points; i++)		
				fprintf(fid1,"  %e    	%e    %e \n", energy_n[i], nu_def[0][i], nu_def[1][i] );
		}

		if(de_number==3)
		{
			for (int i = 0; i < points; i++)		
				fprintf(fid1,"  %e    	%e   	%e	%e \n", energy_n[i], nu_def[0][i], nu_def[1][i], nu_def[2][i] );
		}
		
		fclose(fid1);	
		
		/* 			
		fid1 = fopen("nu_deformation.txt","w");
		for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, nu_deformation[i]);
		fclose(fid1);
		*/
	}
	
	
	// only for 2D materials is required
	// Piezoelectric scattering
        if (scattering_mechanisms[4]==1 && geometry == 2 && de_number > 1)
        {
			fid1 = fopen("pz_scattering_rate.dat","w");
				
			if(de_number==2)
				fprintf(fid1,"# energy           nu_LA 	nu_TA	     total  \n");
			else
				fprintf(fid1,"# energy           nu_LA 	nu_TA		nu_ZA 		total  \n");
			

			if(de_number==2)
			{
				for (int i = 0; i < points; i++)		
					fprintf(fid1,"  %e    	%e    %e       %e \n", energy_n[i], nu_pz[0][i], 
					nu_pz[1][i],  nu_piezoelectric[i]);
			}

			if(de_number==3)
			{
				for (int i = 0; i < points; i++)		
					fprintf(fid1,"  %e    	%e   	%e	%e        %e \n", energy_n[i], nu_pz[0][i], nu_pz[1][i], 
					nu_pz[2][i], nu_piezoelectric[i]);
			}
			
			fclose(fid1);	

	}		

	// intravalley scattering	
	if (scattering_mechanisms[8]==1 && iv_number > 1)   
	{
		fid1 = fopen("intravalley_scattering_rate.dat","w");
		fprintf(fid1,"# energy           total	individual scattering rates \n");

			for (int i = 0; i < points; i++)		
			{
				fprintf(fid1,"  %e    	%e   \t", energy_n[i], nu_iv_total[i]);
				
				for(int m3=0; m3<iv_number; m3++)
					fprintf(fid1,"  %e  \t", nu_iv[m3][i] );
				
				fprintf(fid1,"  \n");

			}
		fclose(fid1);	
	
	}	
	
	// so pop scattering
	if (scattering_mechanisms[10]==1 && so_pop_number > 1)
	{
			
		fid1 = fopen("so_pop_out_scattering_rate.dat","w");
			fprintf(fid1,"# energy           total	individual out scattering rates \n");

		for (int i = 0; i < points; i++)		
		{
			fprintf(fid1,"  %e    	%e   \t", energy_n[i], nu_so_pop_total[i]);
			
			for(int m3=0; m3<pop_number; m3++)
				fprintf(fid1,"  %e  \t", So_so_pop[m3][i] );
			
			fprintf(fid1,"  \n");

		}				
		fclose(fid1);
			
		fid1 = fopen("so_pop_in_scattering_rate.dat","w");
			fprintf(fid1,"# energy           individual in scattering rates \n");

		for (int i = 0; i < points; i++)		
		{
			fprintf(fid1,"  %e  \t", energy_n[i]);
			
			for(int m3=0; m3<pop_number; m3++)
				fprintf(fid1,"  %e  \t", Si_so_pop[m3][i] );
			
			fprintf(fid1,"  \n");

		}				
		fclose(fid1);	
	}
					
	
	//cout<<"Reached here outside save_scattering_rate_contribution "<<endl;
	
}



