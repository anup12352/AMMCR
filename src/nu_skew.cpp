#include"main.h"

// link for bessel function
// https://www.bragitoff.com/2017/08/bessel-function-series-c-program/#:~:text=In%20this%20post%20we%20will,and%20thus%20find%20the%20values.

void nu_skew()
{
	double k,v;
	
	cout<<"Inside skew"<<endl;
	getchar();
	
	lambda_so =  5.3e-20;  // m^2 for GaAs
	double kf, eff;
	double theta[2000], dtheta;
	
	dtheta = -1*pi/(2*1000);
	
	theta[0] = -1*pi/2.0;
	
	for(int i=1;i<2000;i++)
	{
		theta[i] = theta[i-1]+ dtheta;
	}
	
	
	eff = abs(efef_n);
	
	double arr[points];
	for (int i=0;i<points;i++)
		arr[i] = abs(energy_n[i] - eff);
		
	int index = FindMinInd(arr,points);
	
	kf = k_grid[index];
			
	kf = kf*1e9;
	// converted from 1/nm to 1/m
	
	//cout<<"kf = "<<kf<<"  index = "<<index<<endl;
	//getchar();
	
	double a1,a2;
	
	double dd = 0;
	
	for(int i=0;i<2000;i++)
	{
		dd = dd + (1+cos(theta[i]*cos(theta[i])));
	}
	
	//cout<<"dd old = "<<dd<<endl;
	
	dd = dd*(kf*kf*kf*kf*lambda_so*lambda_so)/2000;
	
	//cout<<"dd new = "<<dd<<endl;
	//getchar();
	
	for (int i = 0;i<points;i++)
	{
	        k = k_grid[i]*1e9;	 // converted from 1/nm to 1/m				
    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s
   		
    		
		a1 = denom[i]*(1+2.0*kf*kf*kf*kf*lambda_so*lambda_so/3.0);
		
		a2 = denom[i]*dd/3.0;
		//a2 = 0;
		
		//cout<<"a1 = "<<a1<<"   a2 =  "<<a2<<"   i   = "<<i<<"   denom[i] =  "<<denom[i]<<endl;
		//getchar();
		nu_skew_rate[i] = a1 + a2 - denom[i];
		
		denom[i] = a1 + a2;
		
		nu_el[i] = nu_el[i] + nu_skew_rate[i];
		
	}	
	
}




