#include"main.h"

void save_results()
{

	FILE *fid1, *fid2;
		
	int count;
	if(variation==0)
		count = count_t;
	else
		count = count_d;


	// --------- saving mobility results ----------------------- 
	if (ispin == 1)
	{
		fid1 = fopen("mobility.dat","w");		
		if(Bfield!=0) // && type == "n")
			fid2 = fopen("mobility_hall.dat","w");				
	}
	else
	{	
		if(ispin==2 && kk==0)
		{
			fid1 = fopen("mobility_up_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("mobility_up_spin_hall.dat","w");				
		}
		
		if(ispin==2 && kk==1)		
		{
			fid1 = fopen("mobility_down_spin.dat","w");

			if(Bfield!=0) // && type == "n")
				fid2 = fopen("mobility_down_spin_hall.dat","w");				
		}
	}			

	if (variation==0)   // temperature variation
	{
		fprintf(fid1,"#Temperature(K)");
		if(Bfield!=0) // && type == "n")
			fprintf(fid2,"#Temperature(K)");		
	}
	else
	{
		if(geometry==1)
		{
			fprintf(fid1,"#Doping(cm^-3)");

			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-3)");						
		}
		else if(geometry==2)
		{
			fprintf(fid1,"#Doping(cm^-2)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-2)");
		}
	}
		
	
	if(geometry==1 && type=="n")
	{
		fprintf(fid1,"Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral   Mobility_ir  Mobility_skew   \n");

		for (int i = 0; i <count ; i++)
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e    %e    %e      %e \n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1], calc_mobility_npop[i][1], calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1], calc_mobility_ir[i][1], calc_mobility_skew[i][1]);


		if(Bfield!=0)
		{
						
		fprintf(fid2,"Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral   Mobility_ir  mobility_skew   \n");
			for (int i = 0; i <count ; i++)
			    fprintf(fid2,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e    %e    %e    %e  \n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1], 
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1], calc_mobility_hall_ir[i][1],  calc_mobility_hall_skew[i][1]);
			fclose(fid2);
		}
	}
	else if (geometry==1 && type=="p")
	{
		fprintf(fid1,"  Mobility(cm^2/V-s)    Mobility_rta     Mobility_remote_ii   Mobility_po  Mobility_npop   Mobility_de\n");

		for (int i = 0; i <count ; i++)
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    \n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], 
		calc_mobility_po[i][1], calc_mobility_npop[i][1], calc_mobility_de[i][1]);

		if(Bfield!=0)
		{				
			fprintf(fid2," Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii   Mobility_po   Mobility_npop	Mobility_de    \n");

			for (int i = 0; i <count ; i++)
			fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    \n",
			calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], 
			calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1], calc_mobility_hall_de[i][1]);

			fclose(fid2);
		}
	}   // end of geometry ==1 part if 
	else  // geometry==2
	{
		fprintf(fid1," Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop 	mobility_ir     mobility_skew    \n");

		for (int i = 0; i <count ; i++)
		{
		
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e \n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		calc_mobility_npop[i][1], calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_so_pop[i][1], 
		calc_mobility_ir[i][1],  calc_mobility_skew[i][1]);
		
		
		}

		if(Bfield!=0)
		{				
			fprintf(fid2," Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii   Mobility_po   Mobility_npop	Mobility_de    Mobility_pe 	Mobility_so_pop 	Mobility_ir    Mobility_skew  \n");
			for (int i = 0; i <count ; i++)
				fprintf(fid2,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e  \n",
				calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], 
				calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1],
				calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_so_po[i][1],
				calc_mobility_hall_ir[i][1], calc_mobility_hall_skew[i][1]);
			fclose(fid2);
		}

	}  // end of geometry ==2 part 
	
	fclose(fid1);
	// --------- saving mobility results completed -----------------------
	


	// --------- saving conductivity results ----------------------- 
	if (ispin == 1)
	{
		fid1 = fopen("conductivity.dat","w");		
		if(Bfield!=0) // && type == "n")
			fid2 = fopen("conductivity_hall.dat","w");				
	}
	else
	{	
		if (ispin == 2 && kk==0)
		{
			fid1 = fopen("conductivity_up_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("conductivity_up_spin_hall.dat","w");
		}
		else
		{
			fid1 = fopen("conductivity_down_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("conductivity_down_spin_hall.dat","w");
		}
	}			

	if (variation==0)   // temperature variation
	{
		fprintf(fid1,"#Temperature(K)");
		if(Bfield!=0) // && type == "n")
			fprintf(fid2,"#Temperature(K)");		
	}
	else
	{
		if(geometry==1)
		{
			fprintf(fid1,"#Doping(cm^-3)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-3)");							
		}
		else if(geometry==2)
		{
			fprintf(fid1,"#Doping(cm^-2)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-2)");		
		}
	}
	
	if(geometry==1)
	{	
		fprintf(fid1,"      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA) \n");
		for (int i = 0; i <count; i++)
		    fprintf(fid1," %e         %e              %e \n", calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		if(Bfield!=0) // && type == "n")
		{
			fprintf(fid2,"      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA) \n");
			for (int i = 0; i <count; i++)
			    fprintf(fid2," %e         %e              %e \n", 
			    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid2);
		}	
	}
	else // for 2D
	{

		fprintf(fid1,"      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA) \n");
		for (int i = 0; i <count; i++)
		    fprintf(fid1," %e         %e              %e \n", calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		if(Bfield!=0) // && type == "n")
		{
		fprintf(fid2,"        Conductivity_xx(S/cm)(Rode)  Conductivity_xx (S/cm)(RTA)     Conductivity_xy(S/cm)(Rode)  Conductivity_xy (S/cm)(RTA)\n");
		for (int i = 0; i <count; i++)
		    fprintf(fid2," %e         %e              %e           %e        %e  \n", calc_sigma[i][0], calc_sigma_xx[i][1],
		     calc_sigma_xx_rta[i][1], calc_sigma_xy[i][1], calc_sigma_xy_rta[i][1]);
		fclose(fid2);
		}					
	
	}

	// --------- saving thermopower results ----------------------- 
	if(geometry==1 && type == "n")
	{
		if (ispin == 1)
			fid1 = fopen("thermopower.dat","w");		
		else
		{	if (ispin == 2 && kk==0)
				fid1 = fopen("thermopower_up_spin.dat","w");
			else
				fid1 = fopen("thermopower_down_spin.dat","w");
		}			

		if (variation==0)   // temperature variation
			fprintf(fid1,"#Temperature(K)");
		else
			fprintf(fid1,"#Doping(cm^-3)");

		fprintf(fid1,"#	Thermopower(uV/K) \n");

		for (int i = 0; i <count ; i++)
		    fprintf(fid1," %e        %e\n", calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	
	}
	
	//-------------saving hall factor ----------------------------
	if(Bfield!=0)
	{
	
		if (ispin == 1)
			fid1 = fopen("hall_factor.dat","w");
		if (ispin == 2 && kk == 0)
			fid1 = fopen("hall_factor_up_spin.dat","w");
		if (ispin == 2 && kk == 1)
			fid1 = fopen("hall_factor_down_spin.dat","w");

		if (variation==0)   // temperature variation
			fprintf(fid1,"#Temperature(K)");
		else
		{
			if(geometry==1)
				fprintf(fid1,"#Doping(cm^-3)");
			else if(geometry==2)
				fprintf(fid1,"#Doping(cm^-2)");
		}
				
		if(geometry==1)
		{
			fprintf(fid1,"  	hall_factor	hall_factor_rta \n");

			for (int i = 0; i <count ; i++)
			{
				fprintf(fid1," %e         %e              %e \n",
				hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			}
			fclose(fid1);
		}
		else
		{
			fprintf(fid1,"  	hall_factor  		hall_factor_rta  	hall_coeff(Ohm-cm/G) 	hall_coeff_rta(Ohm-cm/G) long_restivity(cm-ohm) long_restivity_rta     Magneto_Resistance_Coefficient    Magneto_Resistance_Coefficient RTA \n");

			for (int i = 0; i <count ; i++)
			{
				fprintf(fid1," %e            %e           %e         %e           %e           %e         %e             %e           %e \n",
				hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1], 
				calc_hall_coeff[i][1], calc_hall_coeff_rta[i][1], 
				calc_long_restivity[i][1] , calc_long_restivity_rta[i][1], 
				calc_mg_resist[i][1], calc_mg_resist_rta[i][1]);
			}
			fclose(fid1);
		}
	}
	// ----------------------------- hall factor saved--------------------------------------------
	//----------------------------- saving freq variation conductivity -----------------------------

	if(freq_variation==1)
	{
	
		if (ispin == 1)
			fid1 = fopen("conductivity_freq.dat","w");
		if (ispin == 2 && kk == 0)
			fid1 = fopen("conductivity_freq_up_spin.dat","w");
		if (ispin == 2 && kk == 1)
			fid1 = fopen("conductivity_freq_down_spin.dat","w");


		fprintf(fid1,"Frequency (Hz)    Real_conductivity    Imag_conductivity (S/cm)    Real_Mobility	Imaginary_Mobility (cm^2/V-s) \n");
		
		for (int i = 0; i <count ; i++)
		{
			if(variation==0)
				fprintf(fid1,"# For temperature =  %e  K  \n",calc_mobility[i][0]);
			else
			{
				fprintf(fid1,"# For doping =  %e   ",calc_mobility[i][0]);
				if(geometry==1)
					fprintf(fid1," cm^-3   \n");
				else
					fprintf(fid1," cm^-2 \n");
			}		
			for (int j = 0; j <len_freq ; j++)
			{
				fprintf(fid1,"%e         %e              %e  	%e         %e    \n",
				freq[j], sigma_r[i][j],sigma_i[i][j], mobility_r[i][j], mobility_i[i][j]);
			}
			fprintf(fid1," \n ");
		}
		fclose(fid1);
	} // end of freq variation results saved
	//---------------------------- freq variation conductivity saved ----------------------------------------	
	//-----------------------------------------------------------------------------------------------
	//*/
	
	
	
	//---------------------------------- save acoustic mobility parts-----------------------------------------
	
	if(de_number>1 && scattering_mechanisms[3]==1)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("acoustic_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("acoustic_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("acoustic_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  Acoustic Mobility(cm^2/V-s)    individual mobility \n");
		
		for (int i = 0; i <count ; i++)
		{	fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_de[i][1]);
			
			for(int j=0;j<=de_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_de[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}
		fclose(fid1);
	} // end of if condition de_number > 1
	//---------------------------------------------------------------------------------------------------------------
	
	//---------------------------------- save PZ mobility parts only for 2D-----------------------------------------
	
	if(de_number>1  && scattering_mechanisms[4]==1 && geometry==2)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("piezoelectric_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("piezoelectric_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("piezoelectric_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  Piezoelectric Mobility (cm^2/V-s)      individual mobility   \n");
		
		for (int i = 0; i <count ; i++)
		{	fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_pe[i][1]);
			
			for(int j=0;j<=de_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_pe[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}

		fclose(fid1);
	} // end of if condition de_number > 1
	//---------------------------------------------------------------------------------------------------------------
	
	
	//---------------------------------- save npop mobility parts-----------------------------------------
	
	if(npop_number>1 && scattering_mechanisms[2]==1)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("npop_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("npop_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("npop_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  Npop Mobility(cm^2/V-s)    individual mobility \n");
		
		for (int i = 0; i <count ; i++)
		{	fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_npop[i][1]);
			
			for(int j=0;j<=npop_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_npop[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}
		fclose(fid1);
	} // end of if condition npop_number > 1
	//---------------------------------------------------------------------------------------------------------------
	
	//---------------------------------- save iv mobility parts-----------------------------------------
	
	if(iv_number>1 && scattering_mechanisms[2]==1)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("intravalley_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("intravalley_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("intravalley_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  Intravalley Mobility(cm^2/V-s)    individual mobility \n");
		
		for (int i = 0; i <count ; i++)
		{	fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_iv[i][1]);
			
			for(int j=0;j<=iv_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_iv[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}
		fclose(fid1);
	} // end of if condition iv_number > 1
	//---------------------------------------------------------------------------------------------------------------


	//---------------------------------- save pop mobility parts-----------------------------------------
	
	if(pop_number>1 && scattering_mechanisms[1]==1)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("pop_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("pop_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("pop_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  POP Mobility(cm^2/V-s)    individual mobility \n");
		
		for (int i = 0; i <count ; i++)
		{	fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_po[i][1]);
			
			for(int j=0;j<=pop_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_po[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}
		fclose(fid1);
	} // end of if condition pop_number > 1
	//---------------------------------------------------------------------------------------------------------------


	//---------------------------------- save so pop mobility parts-----------------------------------------
	
	if(so_pop_number>1 && scattering_mechanisms[10]==1)
	{	
		//cout<<"count = "<<count<<endl;
		if (ispin == 1)
		{
			fid1 = fopen("so_pop_mobility.dat","w");		
		}
		else
		{	
			if(ispin==2 && kk==0)
			{
				fid1 = fopen("so_pop_mobility_up_spin.dat","w");
			}
			
			if(ispin==2 && kk==1)		
			{
				fid1 = fopen("so_pop_mobility_down_spin.dat","w");
			}
		}			

		if (variation==0)   // temperature variation
		{
			fprintf(fid1,"#Temperature(K)");
		}
		else
		{
			if(geometry==1)
			{
				fprintf(fid1,"#Doping(cm^-3)");
			}
			else if(geometry==2)
			{
				fprintf(fid1,"#Doping(cm^-2)");
			}
		}

		fprintf(fid1,"  SO POP Mobility(cm^2/V-s)    individual mobility \n");
		
		for (int i = 0; i <count ; i++)
		{	
			fprintf(fid1,"%e     %e \t", calc_mobility[i][0], calc_mobility_so_pop[i][1]);
			
			for(int j=0;j<=so_pop_number-1;j++)
			{
				fprintf(fid1,"%e  \t ", calc_mobility_so_pop[i][j+2]);
			}
			fprintf(fid1,"\n");		
		}
		fclose(fid1);
	} // end of if condition so_pop_number > 1
	//---------------------------------------------------------------------------------------------------------------


	
}	

