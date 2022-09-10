#include"main.h"

void conductivity_time(double T, int j)
{
	/*
        for (int counter1=0;counter1<points;counter1++)
	{
		denom1[counter1] = nu_pop_total[counter1] + nu_el[counter1] + omega_s;
		cout<<"counter1 = "<<counter1<<endl;
		cout<<"denom1[counter1] = "<<denom1[counter1]<<endl;
		getchar();
	}
	*/
	
	double g_result_time[points][10]={0}, electric_driving_force_new[points];	
	double denom1[points]={1};
	        
        for (int counter1 = 0;counter1 < points;counter1++)
	{
		electric_driving_force_new[counter1] = -(1*Efield_time[j]/h_bar)*df0dk_grid[counter1]*1e-7;
		// unit is 1/s , hbar unit is eV-s so e is not multiplied in numerator
		
		g_time[counter1] = 0;
		
		denom1[counter1] = nu_pop_total[counter1] + nu_el[counter1] + omega_s;			
	}
		
//----------------------------------S_O_grid calculated ---------------------------------------------------------------------

//------------------------------------ Iteration started ------------------------------------------------------------------
            for (int iteration = 0;iteration<iterations;iteration++)
            {

                for (int counter1 = 0;counter1 < points;counter1++)
			Si_grid[counter1] = 0;
			
                for (int counter1 = 0;counter1 < points;counter1++)
                {
                    double k_dum = k_grid[counter1];

		    // If POP scattering is included
	            if (scattering_mechanisms[1] == 1)
	            {
			for (int m3 = 0;m3 < pop_number;m3++)
			{	            
				double arr[points];
				for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - kminus_grid_pop[m3][counter1]);
				int minus_index =FindMinInd(arr,points);

				for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - kplus_grid_pop[m3][counter1]);
				int plus_index =FindMinInd(arr,points);

				Si_grid[counter1] = Si_grid[counter1] + 
				Sa_pop[m3][counter1]*g_time[minus_index] + Se_pop[m3][counter1]*g_time[plus_index];
					                    
				/*
				cout<<"counter1 = "<<counter1<<endl;
				cout<<"Si_grid[counter1] = "<<Si_grid[counter1]<<endl;
				cout<<"N_poph_atT[m3] =    "<<N_poph_atT[m3]<<endl;
				cout<<"Sa_pop[m3][counter1] =    "<<Sa_pop[m3][counter1]<<endl;
				cout<<"Se_pop[m3][counter1]   "<<Se_pop[m3][counter1]<<endl;
				getchar();
				//*/	                            
		            }
			}
                }
                

                for (int counter1=0;counter1<points;counter1++)
                {
                	if(j!=0)
			{
				g_time[counter1] = (Si_grid[counter1] + electric_driving_force_new[counter1] +
				 omega_s*g_time_old[counter1])/(denom1[counter1]);
			}
			else
			{
		    		g_time[counter1] = (Si_grid[counter1] + electric_driving_force_new[counter1])
		    		/denom1[counter1];
			
			}
			
               
                    /*
                    cout<<"counter1 = "<<counter1<<"    iteration = "<<iteration<<endl;
                    cout<<"g_time[counter1] = "<<g_time[counter1]<<endl;
		     cout<<"Si_grid[counter1] =    "<<Si_grid[counter1]<<endl;
		     cout<<"denom1[counter1]   "<<denom1[counter1]<<endl;
		     cout<<"electric_driving_force_new[counter1]   "<<electric_driving_force_new[counter1]<<endl;
		     cout<<" g_time_old[counter1] = "<<g_time_old[counter1]<<endl;
		     cout<<" omega_s*g_time_old[counter1] = "<<omega_s*g_time_old[counter1]<<endl;
		     getchar();
		     //*/
			     	                
			g_result_time[counter1][iteration]= g_time[counter1];
			            
                }
                   
            }  // uteration loop ends here
            
	for (int counter1=0;counter1<points;counter1++)
		g_time_old[counter1] = g_time[counter1];

	/*
		FILE *fid1;
		fid1 = fopen("g_time.txt","w");
		fprintf(fid1,"#index	1	2	3	4	5	6	7    8     9    10 \n" );
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d    ", i+1);
			for (int iteration = 0;iteration<iterations;iteration++)
			{
				fprintf(fid1,"%e   ", g_result_time[i][iteration]);
			}
			fprintf(fid1,"\n");
		}
		fclose(fid1);
		//cout<<"file saved"<<endl;
		//getchar();
	*/
	
	/*
		if(j==0)
		{
			FILE *fid1;
			fid1 = fopen("g_time_0.txt","w");
			fprintf(fid1,"#index	1	2	3	4	5	6	7    8     9    10 \n" );
			for (int i = 0; i < points; i++)
			{
				fprintf(fid1,"%d    ", i+1);
				for (int iteration = 0;iteration<iterations;iteration++)
				{
					fprintf(fid1,"%e   ", g_result_time[i][iteration]);
				}
				fprintf(fid1,"\n");
			}
			fclose(fid1);
		}
		else if(j==1)
		{
			FILE *fid1;
			fid1 = fopen("g_time_1.txt","w");

			fprintf(fid1,"#index	1	2	3	4	5	6	7    8     9    10 \n" );
			for (int i = 0; i < points; i++)
			{
				fprintf(fid1,"%d    ", i+1);
				for (int iteration = 0;iteration<iterations;iteration++)
				{
					fprintf(fid1,"%e   ", g_result_time[i][iteration]);
				}
				fprintf(fid1,"\n");
			}
			fclose(fid1);
		}            
		else if(j==10)
		{
			FILE *fid1;
			fid1 = fopen("g_time_10.txt","w");

			fprintf(fid1,"#index	1	2	3	4	5	6	7    8     9    10 \n" );
			for (int i = 0; i < points; i++)
			{
				fprintf(fid1,"%d    ", i+1);
				for (int iteration = 0;iteration<iterations;iteration++)
				{
					fprintf(fid1,"%e   ", g_result_time[i][iteration]);
				}
				fprintf(fid1,"\n");
			}
			fclose(fid1);
		}            

	*/
		    
                        
            mobility_time[j] = mu_overall_time(E_F,T,coefficients_cond,kindex_cond,g_time,nu_el,points,a11,j);
            //cout<<"j    =   "<<j<<"   mobility_time[j] = "<<mobility_time[j]<<endl;
            
            //unit cm^2/(V-s)
            
            sigma_time[j] = mobility_time[j] *  n0 * e;
            // unit S/cm 
            
            J_time[j] = sigma_time[j]*(Efield_time[j]);   
            // unit is A/cm^2      
	    	    
	    	    	    	            
}


