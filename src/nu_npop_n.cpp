#include"main.h"

// NPOP scattering calculation
void nu_npop_n(double T)			
{
	
	double k_dum;
	
	for (int aa=0;aa<npop_number;aa++)
	{
		N_npop[aa] = N_poph(we_npop[aa],T);
		//cout<<" N_npop[aa] = "<<N_npop[aa]<<endl;
	}   
	//getchar();	

	//cout<<"For npop scattering   =   "<<endl;	
	for (int counter = 0;counter<points;counter++)
	{
	    k_dum = k_grid[counter];					
    	
	    for (int aa = 0;aa<npop_number;aa++)
	    {
	        lambda_e_plus_grid_npop[counter][aa] = abs(lambda_e_plus(counter,we_npop[aa],rho,De_npop[aa],1,points));
	        lambda_e_minus_grid_npop[counter][aa] = abs(lambda_e_minus(counter,we_npop[aa],rho,De_npop[aa],1,points));
	        // Equation number 129 of rode book  for npop nfv = 1 always
	        //cout<<"aa   =   "<<aa<<endl;
	        //cout<<"lambda_e_plus_grid_npop[counter][aa]   =   "<<lambda_e_plus_grid_npop[counter][aa]<<endl;
	        //cout<<"lambda_e_minus_grid_npop[counter][aa]   =   "<<lambda_e_minus_grid_npop[counter][aa]<<endl; 
	    }
	    //getchar();
	}

	int len = sizeof(nu_npop_total)/sizeof(nu_npop_total[0]);	
	for (int counter = 0;counter<len;counter++)
		nu_npop_total[counter] = 0;

	//cout<<endl<<"Inside npop scattering"<<endl;
	//cout<<" npop_number = "<<npop_number<<endl;
	for (int counter = 0;counter<points;counter++)
	{
	//cout<<"counter = "<<counter<<endl;
		for (int aa = 0;aa<npop_number;aa++)
		{
		    //cout<<"aa = "<<aa<<endl;
		    double k_minus = kminus(counter,we_npop[aa],points,energy_n);
		    //cout<<"k_minus = "<<k_minus<<endl;

		    double arr[points];
		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_minus);
		    int minus_index =FindMinInd(arr,points);
		    //cout<<"minus_index = "<<minus_index<<endl;	

		    double k_plus = kplus(counter,we_npop[aa],points,energy_n);
		    //cout<<"k_plus = "<<k_plus<<endl;

		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_plus);
		    int plus_index =FindMinInd(arr,points);

		    //cout<<"plus_index = "<<plus_index<<endl;

		    double f_negative = f0(energy_n[minus_index],efef_n,T);
		    double f_positive =  f0(energy_n[plus_index],efef_n,T);

		    //cout<<"f_negative = "<<f_negative<<endl;
		    //cout<<"f_positive = "<<f_positive<<endl;

		    if (energy_n[counter] < h_bar*we_npop[aa])
		    {
			    nu_npop[aa][counter] = (N_npop[aa] + f_positive) *lambda_e_plus_grid_npop[plus_index][aa];
		    }
		    else
		    {    
			nu_npop[aa][counter] = (N_npop[aa] + 1 - f_negative) * lambda_e_minus_grid_npop[minus_index][aa] + 
						(N_npop[aa] + f_positive)*lambda_e_plus_grid_npop[plus_index][aa];
		    }
		    
		    
		    nu_npop_total[counter] = nu_npop_total[counter] + nu_npop[aa][counter];            
		    
		    //cout<<"nu_npop[aa][counter] = "<<nu_npop[aa][counter]<<endl;
		    //cout<<"nu_npop_total[counter] = "<<nu_npop_total[counter]<<endl;
		    //getchar(); 	
		}
	}
			
}
