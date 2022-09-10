
#include"main.h"

void nu_alloy1()
{

	double k, v;

	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];

    		nu_alloy[counter] = (3*pi*pow((k*1e9),2)*pow(Uall,2)*(V0*1e-27)*xx*(1-xx))/(16*h_bar*h_bar*(v*1e-2) );
		// Equation No. 17 from Ramu paper 1
		// The coefficient of conversion is to get 1/s
		//cout<<"[counter] =  "<<counter<<endl;
		//cout<<"nu_alloy[counter] =  "<<nu_alloy[counter]<<endl;
		//getchar();
    }

	/*
	fid1 = fopen("nu_alloy.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, nu_alloy[i]);
	fclose(fid1);
	*/
	
}
