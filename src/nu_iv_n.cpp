#include"main.h"

// NPOP scattering calculation
void nu_iv_n(double T)			
{
	double k_dum;
	
	for (int aa=0;aa<iv_number;aa++)
	{
		N_e[aa] = N_poph(we[aa],T);
		//cout<<" N_e[aa] = "<<N_e[aa]<<endl;
	}

	//cout<<"For intervalley scattering   =   "<<endl;	
	for (int counter = 0;counter<points;counter++)
	{
	    k_dum = k_grid[counter];					
	    //cout<<"counter   =   "<<counter<<endl;
	    for (int aa = 0;aa<iv_number;aa++)
	    {
		lambda_e_plus_grid[counter][aa] = abs(lambda_e_plus(counter,we[aa],rho,De[aa],nfv[aa],points));
		lambda_e_minus_grid[counter][aa] = abs(lambda_e_minus(counter,we[aa],rho,De[aa],nfv[aa],points));
		//cout<<"counter = "<<counter<<endl;
		//cout<<"aa   =   "<<aa<<endl;
		//cout<<"lambda_e_plus_grid   =   "<<lambda_e_plus_grid[counter][aa]<<endl;
		//cout<<"lambda_e_minus_grid   =   "<<lambda_e_minus_grid[counter][aa]<<endl;
		
		// Equation number 129 of rode book
	    }
	    //getchar();
	}
	
	/*
	fid1 = fopen("lambda_e_plus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_e_plus_grid[i]);
	fclose(fid1);

	fid1 = fopen("lambda_e_minus.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, lambda_e_minus_grid[i]);
	fclose(fid1);
	*/

	int len = sizeof(nu_iv_total)/sizeof(nu_iv_total[0]);	
	for (int counter = 0;counter<len;counter++)
		nu_iv_total[counter] = 0;

	//cout<<endl<<"Inside intravalley scattering"<<endl;
	//cout<<" iv_number = "<<iv_number<<endl;
	for (int counter = 0;counter<points;counter++)
	{
	//cout<<"counter = "<<counter<<endl;
		for (int aa = 0;aa<iv_number;aa++)
		{
			//cout<<"aa = "<<aa<<endl;
			double k_minus = kminus(counter,we[aa],points,energy_n);
			//cout<<"k_minus = "<<k_minus<<endl;

			double arr[points];
			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_minus);
			int minus_index =FindMinInd(arr,points);


			double k_plus = kplus(counter,we[aa],points,energy_n);
			//cout<<"k_plus = "<<k_plus<<endl;

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_plus);
			int plus_index =FindMinInd(arr,points);

			//cout<<"plus_index = "<<plus_index<<endl;

			double f_negative = f0(energy_n[minus_index],efef_n,T);
			double f_positive =  f0(energy_n[plus_index],efef_n,T);

			//cout<<"f_negative = "<<f_negative<<endl;
			//cout<<"f_positive = "<<f_positive<<endl;

			if (energy_n[counter] < h_bar*we[aa])
			{
			    nu_iv[aa][counter] = (N_e[aa] + f_positive) *lambda_e_plus_grid[plus_index][aa];
			}
			else
			{    
			nu_iv[aa][counter] = (N_e[aa] + 1 - f_negative) * lambda_e_minus_grid[minus_index][aa] + 
						(N_e[aa] + f_positive)*lambda_e_plus_grid[plus_index][aa];
			}


			nu_iv_total[counter] = nu_iv_total[counter] + nu_iv[aa][counter];            

			//cout<<"nu_iv[aa][counter] = "<<nu_iv[aa][counter]<<endl;
			//cout<<"nu_iv_total[counter] = "<<nu_iv_total[counter]<<endl;
			//getchar(); 	

		}  
	}

			
}
