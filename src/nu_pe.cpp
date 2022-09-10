#include"main.h"

void nu_pe(double T,double P_piezo,double epsilon_s)
// piezoelectric scattering rate according to equation (108) or Rode's book
{

    // From equation (108) of Rode's book (book8):

	double k, v;

	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];

    		nu_piezoelectric[counter] = (e*e*k_B*T*P_piezo*P_piezo)/(6*pi*h_bar*h_bar*epsilon_s*epsilon_0*v)*(3-6*pow(c_n[counter],2)+
            4*(pow(c_n[counter],4)))*100/(1.60217657e-19);

    // *100 is to convert 1/cm to 1/m and 1/e is to convert 1/ev to 1/(N.m)

    // 1e-4 is coming from unit conversion (take a look at OneNote notes in Piezoelectric Acoustic section in Book8 notes)
    // Please note that e^2 does not appear in this equation since P^2 has V^2
    // units and V^2 with e^2 is eV^2 canceled eV^2 of hbar^2 in the denominator
	}
	
	/*
	fid1 = fopen("nu_piezoelectric.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, nu_piezoelectric[i]);
	fclose(fid1);
	*/
	
	//------------------------------ reading data -----------------------------------------------
	/*
	fid1 = fopen("nu_piezoelectric.txt","r");
	for (int i = 0; i < points; i++)
	{
	fgets(line, 1000, fid1);
	sscanf(line, "%lf", &nu_piezoelectric[i]);
	}
	fclose(fid1);
	*/
	
}
