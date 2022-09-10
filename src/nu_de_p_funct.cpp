#include"main.h"


// Acoustic scattering calculation  
void nu_de_p_funct(int T_loop)
{
	
	double T=T_array[T_loop];
	double v;
	double k_dum;
	
	//double nu[limit2][de_number]={0};

	for (int counter = 0;counter<points;counter++)
		nu_deformation[counter] = 0;

	//cout<< "temp  "<<T<< endl;
	//cout<<"C_long = "<<C_long<<" dyme/cm^2"<<endl;	
	for (int counter = 0;counter<points;counter++)
	    {
		k_dum = k_grid[counter]*1e9;
		v = v_p[counter]/100;   // unit m-s
			        
	        nu_def[0][counter] = ((k_B*e)*T*pow(E_deformation[0],2)*(k_dum*k_dum))/(2*pi*h_bar*h_bar*(C_long/10.0)*v); 
	        // equation taken from ramu thesis equation 3.32a 

		if(de_number==2)		
		        nu_def[1][counter] = ((k_B*e)*T*pow(E_deformation[1],2)*(k_dum*k_dum))/(2*pi*h_bar*h_bar*(C_trans/10.0)*v); 
		
		if(de_number==3)		
		{
		        nu_def[1][counter] = ((k_B*e)*T*pow(E_deformation[1],2)*(k_dum*k_dum))/(2*pi*h_bar*h_bar*(C_trans/10.0)*v); 
		        nu_def[2][counter] = ((k_B*e)*T*pow(E_deformation[2],2)*(k_dum*k_dum))/(2*pi*h_bar*h_bar*(C_za/10.0)*v); 
		}
		
		for(int i=0;i<de_number;i++)
			nu_deformation_p[counter][0][0] = nu_deformation_p[counter][0][0] + nu_def[i][counter];

		nu_deformation[counter]= nu_deformation_p[counter][0][0];
		
		/*
		cout<<"counter = "<<counter<<endl;
		cout<<"v = "<<v<<endl;
		cout<<"nu_deformation_p[counter][0][0] =  "<<nu_deformation_p[counter][0][0]<<endl;
		getchar();
		*/
	    }
	
	/*
	FILE *fid1;
	fid1 = fopen("acoustic_scattering_rate.dat","w");
	
	if(de_number==1)
		fprintf(fid1,"# energy           nu_LA \n");
	else if(de_number==2)
		fprintf(fid1,"# energy           nu_LA 	nu_TA	\n");
	else
		fprintf(fid1,"# energy           nu_LA 	nu_TA		nu_ZA \n");
	
	
	if(de_number==1)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e   \n", energy_p[i], nu[i][0] );
	}

	if(de_number==2)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e    %e \n", energy_p[i], nu[i][0], nu[i][1] );
	}

	if(de_number==3)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e   	%e	%e \n", energy_p[i], nu[i][0], nu[i][1], nu[i][2] );
	}
	fclose(fid1);
	
	//*/	
}
