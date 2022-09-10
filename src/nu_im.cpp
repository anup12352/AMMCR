#include"main.h"

void nu_im(double epsilon_s,double N_im)
// Neutral impurity scattering rate according to equation 4.25 from Garick Ng thesis this is write expression
{

    	double k, v;


	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];
		nu_neutralimpurity[counter] = (N_im*1e6)*(80*pi*epsilon_s*epsilon_0*(h_bar*e)*pow((v/100),2))/(e*e*pow((k*1e9),2));
		//cout<<"nu_neutralimpurity[counter] =  "<<nu_neutralimpurity[counter]<<endl;
	}

// 1e6 is to convert cm-3 to m-3
// e with hbar is to convert eV-s to J-s
// v(counter) is divided with 100 to convert cm/s to m/s
// 1e9 with k_dum to convert 1/nm to 1/m

	/*
	fid1 = fopen("nu_neutralimpurity.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, nu_neutralimpurity[i]);
	fclose(fid1);
	*/
	
}
