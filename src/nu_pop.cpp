#include"main.h"

void nu_pop(double T, double efef, int ii, int T_loop)
{


	int minus_index, plus_index;
	double arr[points];
	
	for(int m3=0; m3<pop_number; m3++)
	{			
		N_poph_atT[m3] = N_poph(we_pop[m3],T);
		//cout<<"N_poph_atT[m3] = "<<N_poph_atT[m3]<<endl;
		for (int counter1 = 0;counter1 < points;counter1++)
		{

			kplus_grid_pop[m3][counter1] = kplus(counter1,we_pop[m3],points,energy_n);
			kminus_grid_pop[m3][counter1] = kminus(counter1,we_pop[m3],points,energy_n);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kminus_grid_pop[m3][counter1]);
			minus_index =FindMinInd(arr,points);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kplus_grid_pop[m3][counter1]);
			plus_index =FindMinInd(arr,points);

			plus_index_pop[m3][counter1] = plus_index;
			minus_index_pop[m3][counter1] = minus_index;
		    
		    /*
		    cout<<"counter1 = "<<counter1<<endl;
		    cout<<"points = "<<points<<endl;
		    cout<<"plus_index_pop[m3][counter1] = "<<plus_index_pop[m3][counter1]<<endl;
		    cout<<"minus_index_pop[m3][counter1] = "<<minus_index_pop[m3][counter1]<<endl;	
		    cout<<"kplus_grid_pop[m3][counter1] = "<<kplus_grid_pop[m3][counter1]<<endl;
		    cout<<"kminus_grid_pop[m3][counter1] = "<<kminus_grid_pop[m3][counter1]<<endl;
		    getchar();	
		    */
		    
		}
	}
								
	for(int m3=0; m3<pop_number; m3++)
	{			
		for (int counter = 0;counter < points;counter++)
		{
		    Aminus_grid[m3][counter] = Aminus(counter,we_pop[m3],points,energy_n);
		    Aplus_grid[m3][counter] = Aplus(counter,we_pop[m3],points,energy_n);


		    betaplus_grid[m3][counter] = betaplus(counter,we_pop[m3],epsilon_s[T_loop],epsilon_inf[T_loop],points);
		    betaminus_grid[m3][counter] = betaminus(counter,we_pop[m3],epsilon_s[T_loop],epsilon_inf[T_loop],points);

		    lambda_i_plus_grid[m3][counter] = abs(lambda_i_plus(counter,we_pop[m3],Aplus_grid[m3][counter], epsilon_s[T_loop],epsilon_inf[T_loop], points));

		    lambda_i_minus_grid[m3][counter] = abs(lambda_i_minus(counter,we_pop[m3],Aminus_grid[m3][counter],                epsilon_s[T_loop],epsilon_inf[T_loop],points));

		    lambda_o_plus_grid[m3][counter] = abs(lambda_o_plus(counter,we_pop[m3],Aplus_grid[m3][counter],
		    epsilon_s[T_loop],epsilon_inf[T_loop],points));

		    lambda_o_minus_grid[m3][counter] = abs(lambda_o_minus(counter, we_pop[m3],Aminus_grid[m3][counter],
		    epsilon_s[T_loop],epsilon_inf[T_loop],points));

			//---------------------------- code to debug -------------------------------------------------------------
			//cout<<"counter = "<<counter<<endl;
			//cout<<"Aplus_grid[m3][counter] =   "<<Aplus_grid[m3][counter]<<endl;
			//cout<<"Aminus_grid[m3][counter] =   "<<Aminus_grid[m3][counter]<<endl;

			//cout<<"betAplus_grid[m3][counter] =  "<<betAplus_grid[m3][counter]<<endl;
			//cout<<"betAminus_grid[m3][counter] =  "<<betAminus_grid[m3][counter]<<endl;

			//cout<<"lambda_i_plus_grid[m3][counter] =  "<<lambda_i_plus_grid[m3][counter]<<endl;
			//cout<<"lambda_i_minus_grid[m3][counter] =  "<<lambda_i_minus_grid[m3][counter]<<endl;
			//cout<<"lambda_o_plus_grid[m3][counter] =  "<<lambda_o_plus_grid[m3][counter]<<endl;
			//cout<<"lambda_o_minus_grid[m3][counter] =  "<<lambda_o_minus_grid[m3][counter]<<endl;

			//getchar();
		}
	}	
	/*

	FILE *fid1;
	fid1 = fopen("Aplus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, Aplus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("Aminus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, Aminus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("betaplus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, betAplus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("betaminus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, betAminus_grid[m3][i]);
	fclose(fid1);


	fid1 = fopen("lambda_i_plus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_i_plus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("lambda_i_minus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_i_minus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("lambda_o_plus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_o_plus_grid[m3][i]);
	fclose(fid1);

	fid1 = fopen("lambda_o_minus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_o_minus_grid[m3][i]);
	fclose(fid1);
	*/
	

	for (int i=0;i<points;i++)
	{
		for(int m3=0; m3<pop_number; m3++)
			So_pop[m3][i] = 1;

		nu_pop_total[i]=0;
	}

//------------------------------------------ S_o_grid calculated ----------------------------------------------
	for(int m3=0; m3<pop_number; m3++)
	{	
		double sum=0;
		for (int counter1 = 0;counter1 < points;counter1++)
		    sum = sum +  So_pop[m3][counter1];

		double average_dummy = sum/points;
		//cout<<"average dummy = "<<average_dummy<<endl;
		//cout<<"N_poph_atT[m3]   =   "<<N_poph_atT[m3]<<endl;
		
		for (int counter1 = 0;counter1 < points;counter1++)
		{
			//cout<<"counter = "<<counter1<<endl;
			//cout<<"plus_index   =   "<<plus_index_pop[m3][counter1]<<endl;
			//cout<<"minus_index   =   "<<minus_index_pop[m3][counter1]<<endl;

			// out scattering terms
			if ((lambda_o_minus_grid[m3][counter1]==0) && (lambda_o_plus_grid[m3][counter1]==0))
			{
				So_ab_pop[m3][counter1] = 0;
				So_em_pop[m3][counter1] = 0;
				So_pop[m3][counter1] = average_dummy;
			    // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
			}
			else
			{
			    So_ab_pop[m3][counter1] = (N_poph_atT[m3]+f_dist[plus_index_pop[m3][counter1]])*lambda_o_plus_grid[m3][counter1];
			    So_em_pop[m3][counter1] = (N_poph_atT[m3]+1-f_dist[minus_index_pop[m3][counter1]])*lambda_o_minus_grid[m3][counter1];
			    So_pop[m3][counter1] =  So_ab_pop[m3][counter1] + So_em_pop[m3][counter1];
			    
			}
			
			// in scattering terms
			Sa_pop[m3][counter1] = (N_poph_atT[m3] + f_dist[counter1])*lambda_i_minus_grid[m3][counter1];
			Se_pop[m3][counter1] = (N_poph_atT[m3] + 1 - f_dist[counter1])*lambda_i_plus_grid[m3][counter1];
			
			Si_pop[m3][counter1] = Sa_pop[m3][counter1] + Se_pop[m3][counter1];

			/*
		cout<<"(N_poph_atT[m3]+1-f_dist[minus_index_pop[m3][counter1]])  = "<<(N_poph_atT[m3]+1-f_dist[minus_index_pop[m3][counter1]])<<endl;
		cout<<"(N_poph_atT[m3]+f_dist[plus_index_pop[m3][counter1]])   =   "<<(N_poph_atT[m3]+f_dist[plus_index_pop[m3][counter1]])<<endl;
			cout<<"lambda_o_plus_grid[m3][counter1]   = "<<lambda_o_plus_grid[m3][counter1]<<endl;
			cout<<"lambda_o_minus_grid[m3][counter1] = "<<lambda_o_minus_grid[m3][counter1]<<endl;
			cout<<"f_dist[plus_index_pop[m3][counter1]]  = "<<f_dist[plus_index_pop[m3][counter1]]<<endl;
			cout<<"f_dist[minus_index_pop[m3][counter1]]  = "<<f_dist[minus_index_pop[m3][counter1]]<<endl;
			cout<<"So_pop[m3][counter1]  =    "<<So_pop[m3][counter1]<<endl;
			getchar();
			*/
		}

		for (int counter1=0;counter1<points;counter1++)
		    nu_pop_total[counter1] = nu_pop_total[counter1] + So_pop[m3][counter1];
        }
            
	/*
	FILE *fid1;
	fid1 = fopen("pop_scattering_rate.dat","w");	
	
	fprintf(fid1,"Energy (eV)   ab   em   total   \n ");
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		fprintf(fid1,"  %e	%e	%e  \n", So_ab_pop[m3][i], So_em_pop[m3][i], So_pop[m3][i] );
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
            
}


