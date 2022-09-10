#include"main.h"

// link for bessel function
// https://www.bragitoff.com/2017/08/bessel-function-series-c-program/#:~:text=In%20this%20post%20we%20will,and%20thus%20find%20the%20values.

void nu_skew()
{
	double k,v;
	
	for (int i = 0;i<points;i++)
	{
	        k = k_grid[i]*1e9;	 // converted from 1/nm to 1/m				
    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s
    		
    		
    		
    		
		nu_skew_rate[i] = 1e10;
	}	
	
}




