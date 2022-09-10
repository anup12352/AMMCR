
#include"main.h"

void nu_de(double T)
//deformation potential scattering rate according to equation (112) or Rode's book
{

	double k, v;
	

	for(int i=0;i<de_number;i++)
	{	
		for (int counter = 0;counter<points;counter++)
		{
			nu_def[i][counter] = 0;
		}
	}


	for (int counter = 0;counter<points;counter++)
		nu_deformation[counter] = 0;

	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];
    		nu_def[0][counter] = (k_B*T*pow(E_deformation[0],2)*k*k)/(3*pi*h_bar*h_bar*v*C_long)*(3-8*pow(c_n[counter],2)
+6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		// From equation (112) of Rode's book (book8):
		
		// 1e10 is coming from unit conversion (take a look at OneNote notes in Deformation potential //
		// section) and *1.60217657/1e8 is to get from cm/s (v) to hk/(md)


		if(de_number==2)		
		{
		nu_def[1][counter] = (k_B*T*pow(E_deformation[1],2)*k*k)/(3*pi*h_bar*h_bar*v*C_trans)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		}
		
		if(de_number==3)		
		{

		nu_def[1][counter] = (k_B*T*pow(E_deformation[1],2)*k*k)/(3*pi*h_bar*h_bar*v*C_trans)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;

		nu_def[2][counter] = (k_B*T*pow(E_deformation[2],2)*k*k)/(3*pi*h_bar*h_bar*v*C_za)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		}
		
		for(int i=0;i<de_number;i++)
			nu_deformation[counter] = nu_deformation[counter] + nu_def[i][counter];
		
		//cout<<"nu_deformation[counter] =  "<<nu_deformation[counter]<<endl;
	}
			
	//*/

//------------------------------ reading data -----------------------------------------------
	/*
	fid1 = fopen("nu_deformation.txt","r");
	for (int i = 0; i < points; i++)
	{
		fgets(line, 1000, fid1);
		sscanf(line, "%lf", &nu_deformation[i]);
	}
	fclose(fid1);
	*/
}



