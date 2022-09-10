#include"main.h"
// expression from book of ferry
// Page no.138 Eq 3.167 of book Semiconductor transport
void nu_ir(double epsilon_s)
{
	double const1, k, v, const2;
	
	// converted to mks units 
	//fl_L = fl_L*1e-9;
	//fl_delta = fl_delta*1e-9;
	//ND = ND *1e4;
	//ns = ns*1e4;
	
	const1 = pow(e,4)*pow( (((ND + ns/2)*1e4) * (fl_delta*1e-9) * (fl_L*1e-9) ),2)/(2*pow((h_bar*e*epsilon_s*epsilon_0),2));
	
	int lt=2000;
	double intg, dtheta;
	double X[lt+1], theta[lt+1], denom[lt+1];    
	// lt is for theta variation
	
	//cout<<"const1  =  "<<const1<<endl;
	//cout<<"  =   "<<<<endl;
	//cout<<"  =   "<<<<endl;
	//cout<<"  =   "<<<<endl;
	//cout<<"  =   "<<<<endl;
	//cout<<"  =   "<<<<endl;
	
	dtheta = 2*pi/lt;
	theta[0] = 0;
	X[0] = cos(theta[0]/2)*cos(theta[0]/2);
	
	for(int i=1;i<=lt;i++)
	{
		theta[i] = theta[i-1] + dtheta;
		X[i] = cos(theta[i]/2)*cos(theta[i]/2);
	}
	
	for (int counter = 0;counter < points;counter++)
	{

		k = k_grid[counter]*1e9;
		// unit 1/m
		//k = k_dum*1e9;
		// converted to 1/m
		
		v = v_n[counter]*1e-2;
		// unit m/s - converted from cm/s to m/s
		
		const2= const1*k/v;

		//cout<<"counter =  "<<counter+1<<endl;
		//cout<<"k = "<<k<<endl;
		//cout<<"v = "<<v<<endl;
		
		intg = 0;
		
		for(int i=0;i<=lt;i++)
		{
			denom[i] = pow((1 + k*k*(fl_L*fl_L*1e-9*1e-9)*X[i]),1.5); 
			intg = intg + dtheta*1.0/denom[i];
		}
		
		nu_irs[counter] = const2*intg;
		// unit 1/second	

		//cout<<"const2 = "<<const2<<endl;
		//cout<<"intg = "<<intg<<endl;
		//cout<<"nu_irs[counter] = "<<nu_irs[counter]<<endl;
		//getchar();
	}
	
	/*
	FILE *fid1;
	fid1 = fopen("nu_irs.txt","w");
	for(int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e \n", i+1, nu_irs[i]);
	
	fclose(fid1);
	
	//*/
	
	/*
	fid1 = fopen("nu_irs.txt","r");
	for (int i = 0; i < points; i++)
	{
	fgets(line, 1000, fid1);
	sscanf(line, "%lf", &nu_irs[i]);
	}
	fclose(fid1);
	
	*/
}
