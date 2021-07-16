#include"main.h"

void solve_g(double T)
{

//---------------------------------------- solve BTE for g ----------------------------------------------------------

            for (int i=0;i<points;i++)
            {
		S_o_grid[i] = 1;
		g[i]=0;

		g_rta[i]=0;

		g_old[i]=0;
		g_LO[i]=0;
		g_iv[i]=0;

		g_th[i]=0;
		g_th_old[i]=0;
		g_LO_th[i]=0;


		S_o_grid[i]=0;

		S_o_grid_total[i]=0;

		S_i_grid[i]=0;
		S_iLO_grid[i]=0;
		S_i_th_grid[i]=0;
		S_iLO_th_grid[i]=0;

		beta1[i]=0;
		gH[i]=0;
		hH[i]=0;

		gH_rta[i]=0;
		hH_rta[i]=0;

		gH_LO[i]=0;
		hH_LO[i]=0;

		S_i_grid_g[i]=0;
		S_i_grid_h[i]=0;

		S_iLO_grid_g[i]=0;
		S_iLO_grid_h[i]=0;

		S_o_gridH[i]=0;

		S_o_gridH[i] = 1;

		S_o_grid_totalH[i]=0;


	     }	
            for(int counter=0;counter<points;counter++)
            {
                result_g[counter][0] = k_grid[counter];
                result_g_LO[counter][0] = k_grid[counter];
                result_g_th[counter][0] = k_grid[counter];
                result_f[counter][0] = k_grid[counter];
            }


            for (int iteration = 0;iteration<iterations;iteration++)
            {

                double f_dist_temp[points],f_dist_temp_th[points];
                for (int counter1=0;counter1<points;counter1++)
                {
                    f_dist_temp[counter1] = f0(energy_n[counter1],E_F,T) + g[counter1];
                    f_dist_temp_th[counter1] = f0(energy_n[counter1],E_F,T) + g_th[counter1];
                    //cout<<"energy_n[counter1] = "<<energy_n[counter1]<<endl;
                    //cout<<"g[counter1] = "<<g[counter1]<<endl;
                    //cout<<"f_dist_temp[counter1] = "<<f_dist_temp[counter1]<<endl;
                    //cout<<"f_dist_temp_th[counter1] = "<<f_dist_temp_th[counter1]<<endl;
                    //getchar();
                }


                double sum=0;
                for (int counter1 = 0;counter1 < points;counter1++)
                    sum = sum +  S_o_grid[counter1];

                double average_dummy = sum/points;
                //cout<<"average dummy = "<<average_dummy<<endl;

                for (int counter1 = 0;counter1 < points;counter1++)
                {
                    double k_dum =	k_grid[counter1];
                    //f_dist(counter) = 1/(1+exp((energy(counter)-E_F)/(k_B*T)))+g(counter);

                    //f_dist_th(counter)= 1/(1+exp((energy(counter)-E_F)/(k_B*T)))+g_th(counter);

                    double arr[points];
                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kminus_grid[counter1]);
                    int minus_index =FindMinInd(arr,points);

                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kplus_grid[counter1]);
                    int plus_index =FindMinInd(arr,points);

                    // If POP scattering is included
                    if (scattering_mechanisms[1] == 1)
                    {
                        S_i_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g[plus_index];

                        S_iLO_grid[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_LO[plus_index];


                        S_i_th_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_th[minus_index]+
                        (N_poph_atT + 1 - f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_th[plus_index];


                        S_iLO_th_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_LO_th[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_LO_th[plus_index];


                        if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
                            S_o_grid[counter1] = average_dummy;
                            // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
                        else
                            S_o_grid[counter1] = (N_poph_atT+1-f_dist_temp[minus_index])*lambda_o_minus_grid[counter1] +
                            (N_poph_atT+f_dist_temp[plus_index])*lambda_o_plus_grid[counter1];
                            // Equation no. 116 of rode book
                            
                            /*
                            cout<<"counter1 = "<<counter1<<endl;
                            cout<<"S_i_grid[counter1] = "<<S_i_grid[counter1]<<endl;
                            cout<<"S_iLo_grid[counter1] = "<<S_iLO_grid[counter1]<<endl;
                            cout<<"S_i_th_grid[counter1] =    "<<S_i_th_grid[counter1]<<endl;
                            cout<<"S_iLO_th_grid[counter1] =   "<<S_iLO_th_grid[counter1]<<endl;
                            cout<<"S_o_grid[counter1] =    "<<S_o_grid[counter1]<<endl;
			     cout<<"N_poph_atT =    "<<N_poph_atT<<endl;
			     cout<<"f_dist_temp[minus_index]   "<<f_dist_temp[minus_index]<<endl;
			     cout<<"lambda_o_minus_grid[counter1] =    "<<lambda_o_minus_grid[counter1]<<endl;
			     cout<<"lambda_o_plus_grid[counter1]   "<<lambda_o_plus_grid[counter1]<<endl;
			     getchar();
			     //*/	                            
                    }

                }
                
                for (int counter1=0;counter1<points;counter1++)
                    S_o_grid_total[counter1] = S_o_grid[counter1];


                for (int counter1=0;counter1<points;counter1++)
                {
                    g[counter1] = (S_i_grid[counter1]+electric_driving_force[counter1])/(S_o_grid_total[counter1] + nu_el[counter1]);
                    g_th[counter1] = (S_i_th_grid[counter1] + thermal_driving_force[counter1])/(S_o_grid_total[counter1] + nu_el[counter1]);
                }

                // If POP scattering is included
                if (scattering_mechanisms[1] == 1)
                {
                    for (int counter1=0;counter1<points;counter1++)
                        g_LO[counter1] = (S_iLO_grid[counter1] + electric_driving_force[counter1])/S_o_grid[counter1];
                }

                for (int counter1=0;counter1<points;counter1++)
                {
                    result_g[counter1][iteration+1] = g[counter1];

                    result_g_th[counter1][iteration+1] = g_th[counter1];

                    result_f[counter1][iteration+1] = f_dist[counter1];

                //fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',iteration,T,sum(g-g_old)/points);

                    g_old[counter1] = g[counter1];
                    g_th_old[counter1] = g_th[counter1];

                    if (iteration==0)
                       g_rta[counter1] = g[counter1] ;
		     
		     /*
		     if(iteration == iterations-1)
		     {	
			     cout<<"counter1 =   "<<counter1<<endl;	
			     cout<<"g[counter1] =   "<<g[counter1]<<endl;
			     cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
			     cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
			     getchar();
		     }
		     //*/
                }
            }
		/*
                for (int counter1=0;counter1<points;counter1++)
                {
			     cout<<"counter1 =   "<<counter1<<endl;	
			     cout<<"g[counter1] =   "<<g[counter1]<<endl;
			     cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
			     cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
			     getchar();
		}
		*/
		
//---------------------------------------- solve BTE for g finished----------------------------------------------------------


//---------------------------------------- Solve BTE for g and h for hall mobility---------------------------------------------
            
        for (int i=0;i<points;i++)
            S_o_gridH[i] = 1;


	if(Bfield!=0)
	{
            for (int iteration = 0;iteration<iterations;iteration++)
            {

                double f_dist_temp[points];
                double beta1_LO[points];
                
                for (int counter1=0;counter1<points;counter1++)
                {
                    f_dist_temp[counter1] = f0(energy_n[counter1],E_F,T) + gH[counter1] + hH[counter1];
                    //cout<<"energy_n[counter1] = "<<energy_n[counter1]<<endl;
                    //cout<<"gH[counter1] = "<<gH[counter1]<<endl;
                    //cout<<"hH[counter1] = "<<hH[counter1]<<endl;
                    //cout<<"f_dist_temp[counter1] = "<<f_dist_temp[counter1]<<endl;
                    //getchar();
                }


                double sum1=0;
                for (int counter1 = 0;counter1 < points;counter1++)
                    sum1 = sum1 +  S_o_gridH[counter1];

                double average_dummy1 = sum1/points;
                //cout<<"average dummy1 = "<<average_dummy1<<endl;
		
                for (int counter1 = 0;counter1 < points;counter1++)
                {
                    double k_dum =	k_grid[counter1];

                    double arr[points];
                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kminus_grid[counter1]);
                    int minus_index =FindMinInd(arr,points);

                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kplus_grid[counter1]);
                    int plus_index =FindMinInd(arr,points);

                    // If POP scattering is included
                    if (scattering_mechanisms[1] == 1)
                    {
                        S_i_grid_g[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*gH[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*gH[plus_index];

                        S_iLO_grid_g[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*gH_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*gH_LO[plus_index];

                        S_i_grid_h[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*hH[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*hH[plus_index];

                        S_iLO_grid_h[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*hH_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*hH_LO[plus_index];

                        if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
                            S_o_gridH[counter1] = average_dummy1;
                            // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
                        else
                            S_o_gridH[counter1] = (N_poph_atT+1-f_dist_temp[minus_index])*lambda_o_minus_grid[counter1] +
                            (N_poph_atT+f_dist_temp[plus_index])*lambda_o_plus_grid[counter1];
                            // Equation no. 116 of rode book

                            /*		
                            cout<<"counter1 = "<<counter1<<endl;
                            cout<<"S_i_grid_g[counter1] = "<<S_i_grid_g[counter1]<<endl;
                            cout<<"S_iLo_grid_g[counter1] = "<<S_iLO_grid_g[counter1]<<endl;
                            cout<<"S_i_grid_h[counter1] = "<<S_i_grid_h[counter1]<<endl;
                            cout<<"S_iLo_grid_h[counter1] = "<<S_iLO_grid_h[counter1]<<endl;
                            cout<<"S_o_gridH[counter1] =    "<<S_o_gridH[counter1]<<endl;
			     cout<<"N_poph_atT =    "<<N_poph_atT<<endl;
			     cout<<"f_dist_temp[minus_index]   "<<f_dist_temp[minus_index]<<endl;
			     cout<<"lambda_o_minus_grid[counter1] =    "<<lambda_o_minus_grid[counter1]<<endl;
			     cout<<"lambda_o_plus_grid[counter1]   "<<lambda_o_plus_grid[counter1]<<endl;
			     getchar();
			     //*/	                            


                    }
                }
                
                for (int counter1=0;counter1<points;counter1++)
                    S_o_grid_totalH[counter1] = S_o_gridH[counter1];


                for (int counter1=0;counter1<points;counter1++)
                {

                	beta1[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(S_o_grid_totalH[counter1] + nu_el[counter1]));
                	
                	beta1_LO[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(S_o_grid_totalH[counter1]));

                    	gH[counter1] = (S_i_grid_g[counter1] + electric_driving_force[counter1] + beta1[counter1] * 			S_i_grid_h[counter1])/((S_o_grid_totalH[counter1] + nu_el[counter1])*(1 + beta1[counter1]*beta1[counter1]));
                    	
                    	hH[counter1] = (S_i_grid_h[counter1] - beta1[counter1] * electric_driving_force[counter1] -                    	 			beta1[counter1] * S_i_grid_g[counter1])/((S_o_grid_totalH[counter1] + nu_el[counter1]) * (1 + 				beta1[counter1]*beta1[counter1]));

                    	gH_LO[counter1] = (S_iLO_grid_g[counter1] + electric_driving_force[counter1] + beta1_LO[counter1] * 			S_iLO_grid_h[counter1])/((S_o_grid_totalH[counter1] )*(1 + beta1_LO[counter1]*beta1_LO[counter1]));
                    	
                    	hH_LO[counter1] = (S_iLO_grid_h[counter1] - beta1_LO[counter1] * electric_driving_force[counter1] 				- beta1_LO[counter1] * S_iLO_grid_g[counter1])/((S_o_grid_totalH[counter1]) * (1 + 				       beta1_LO[counter1]*beta1_LO[counter1]));
                }
                
                if (iteration==0)
	        {
		        for (int counter1=0;counter1<points;counter1++)
		        {
				gH_rta[counter1] = gH[counter1] ;
				hH_rta[counter1] = hH[counter1] ;
			}	        	
	        }

            }
	    	
 	}  // if condition !Bfield terminated
//---------------------------------------- Solve BTE for g and h for hall mobility finished---------------------------------------------


//-------------------------- save data ---------------------------------------------------	    
	    /*		
            FILE *fid1;
            fid1 = fopen("g_original.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g[i]);
	fclose(fid1);
		 
		 /* 
            fid1 = fopen("g_th.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_th[i]);
	fclose(fid1);

            fid1 = fopen("g_LO.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_LO[i]);
	fclose(fid1);

            fid1 = fopen("g_iv.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_iv[i]);
	fclose(fid1);

            fid1 = fopen("g_rta.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_rta[i]);
	fclose(fid1);
            */
            //cout<<"saved result check here";
            //getchar();
            //getchar();
//-----------------------------------------------------------------------------------------

}
