
#include"main.h"

// 
void nu_rim_2D(double T)
//remote impurity scattering rate 
{
	
	double k, v;
			
	double intg, dtheta, rimp_const;
	double X[limit7+1], Y[limit7+1], Z[limit7+1], theta[limit7+1], aa[limit7+1], q1[limit7+1];    
	// limit7 is for theta variation
	
	rimp_const = (rimp*1e4)*e*e*e*e/(8*pi*h_bar*h_bar*e*e*eps_avg_low*epsilon_0*eps_avg_low*epsilon_0);
	
	//cout<<"rimp_const =  "<<rimp_const<<endl;
	
	dtheta = 2*pi/limit7;
	theta[0] = 0;
	X[0] = sin(theta[0]/2);
	Y[0] = 1 - cos(theta[0]);
	Z[0] = cos(theta[0]);
	
	for(int i=1;i<limit7+1;i++)
	{
		theta[i] = theta[i-1] + dtheta;
		X[i] = sin(theta[i]/2);
		Y[i] = 1 - cos(theta[i]);
		Z[i] = cos(theta[i]);
	}


	if(screening==0)
	{
		double pz1[limit7+1], ep, G;
		
		for (int i = 0;i<points;i++)
		{
			//cout<<"i = "<<i<<endl;
			//cout<<"k[i] = "<<k_grid[i]<<endl;
			

			k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
	    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	

			// change in wave vector for different values of theta
			for(int j=1;j<limit7;j++)			
			{
				q1[j] = 2*k*X[j];
			}
			// to remove zero elements
			q1[0] = q1[1];
			q1[limit7] = q1[limit7-1];
			
			// overlap integral caculated here
			if(overlap==1)
			{			
				for(int j=0;j<limit7+1;j++)			
				{
					aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
					aa[j] = aa[j] * aa[j];
				}
			}
			else
			{
				for(int j=0;j<limit7+1;j++)			
					aa[j] = 1;
			}					
								
			int index, last_index=0;
			double min=1e20;
			
			// polarizibility is calculated for different values of theta 			
			for(int j=0;j<=limit7;j++)
			{				
						
				min=1e15;
								
				for(int l=last_index;l<=limit6;)   // limit6+1 size of polarizability array
				{
					if(min > abs(q1[j]-q[l]))
					{
						index = l;
						min = abs(q1[j]-q[l]);		
					}
					else
					{
						last_index = index;
						break;
					}

					if (q1[j] > q1[j-1])
						l++;						
					else
						l--;					
				}
				pz1[j] = pz[index];  // polarizibility for different theta				
					
			} // loop for theta, polarizibility is calculated for all values of theta of a particular k	
			
			intg = 0;
			// scattering rate is calculated 			
			for(int j=0;j<=limit7;j++)
				intg = intg + aa[j]*exp(-2*q1[j]*dist)*dtheta;	
					
	    		nu_ionizedimpurity[i] = rimp_const*intg/(2*k*v);

			//cout<<"nu_ionizedimpurity[i] =  "<<nu_ionizedimpurity[i]<<endl;
											
		} // loop for different values of k		
	}
	else
	{	
		double pz1[limit7+1], ep, G;
		
		for (int i = 0;i<points;i++)
		{
			//cout<<"i = "<<i<<endl;
			//cout<<"k[i] = "<<k_grid[i]<<endl;
			

			k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
	    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	

			// change in wave vector for different values of theta
			for(int j=1;j<limit7;j++)			
			{
				q1[j] = 2*k*X[j];
			}
			// to remove zero elements
			q1[0] = q1[1];
			q1[limit7] = q1[limit7-1];
			
			// overlap integral caculated here
			if(overlap==1)
			{						
				for(int j=0;j<limit7+1;j++)			
				{
					aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
					aa[j] = aa[j] * aa[j];
				}
			}
			else
			{
				for(int j=0;j<limit7+1;j++)			
					aa[j] = 1;
			}					
						
			int index, last_index=0;
			double min=1e20;
			
			// polarizibility is calculated for different values of theta 			
			for(int j=0;j<=limit7;j++)
			{				
						
				min=1e15;
								
				for(int l=last_index;l<=limit6;)   // limit6+1 size of polarizability array
				{
					if(min > abs(q1[j]-q[l]))
					{
						index = l;
						min = abs(q1[j]-q[l]);		
					}
					else
					{
						last_index = index;
						break;
					}

					if (q1[j] > q1[j-1])
						l++;						
					else
						l--;					
				}
				pz1[j] = pz[index];  // polarizibility for different theta				
					
			} // loop for theta, polarizibility is calculated for all values of theta of a particular k	
			
			intg = 0;
			// scattering rate is calculated 			
			for(int j=0;j<=limit7;j++)
			{		
				G = 1/(epsilon_0*(eps_sub_low + eps_up_low)*q1[j]);
				ep = 1 + e*e*G*pz1[j];  // screening factor
				intg = intg + aa[j]*exp(-2*q1[j]*dist)/(ep*ep)*dtheta;	
			}
			
			
	    		nu_ionizedimpurity[i] = rimp_const*intg/(2*k*v);

			//cout<<"nu_ionizedimpurity[i] =  "<<nu_ionizedimpurity[i]<<endl;
											
		} // loop for different values of k		
	} // else part with screening completed
	
		
	/*		
	
	FILE *fid1;
	if(screening==0)
		fid1 = fopen("remote_impurity_scattering_rate.dat","w");
	else
		fid1 = fopen("remote_impurity_scattering_rate_screening.dat","w");
	
		
	fprintf(fid1,"# energy           nu_LA \n");
	for (int i = 0; i < points; i++)		
		fprintf(fid1,"  %e    	%e   \n", energy_n[i], nu_ionizedimpurity[i] );

	fclose(fid1);
	//*/
}


