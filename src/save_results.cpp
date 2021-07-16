#include"main.h"

void save_results()
{

//------------------------------------------------------ saving results-------------------------------------------------------------
    if (variation==0)   // temperature variation
    {
    
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","a");
			for (int i = 0; i <count_t ; i++)
			{
				fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			}
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 0)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","a");

		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower_up_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_up_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_up_spin_hall.dat","a");

			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_down_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_down_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

	
    }
    else          // Doping variation
    {
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);
		
		if(Bfield!=0)
		{
			fid1 = fopen("mobility_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);
			
			fid1 = fopen("hall_factor.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 0)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);		
		}
	}


    	if (ispin == 2 && kk == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_down_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_down_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}
	
    }
//------------------------------------------------------ saving results completed ------------------------------------------------
}
