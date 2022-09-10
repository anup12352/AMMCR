
#include"main.h"


void nu_de_2D(double T)
//deformation potential scattering rate 
{
	
	double k, v;
	
	for(int i=0;i<de_number;i++)
	{	
		for (int counter = 0;counter<points;counter++)
		{
			nu_def[i][counter] = 0;
		}
	}

	/*
	cout<<"E_deformation[0] = "<<E_deformation[0]<<endl;
	cout<<"C_long = "<<C_long<<endl;
	cout<<"C_trans  = "<<C_trans<<endl;
	cout<<"de_number =  "<<de_number<<endl;
	cout<<"T =  "<<T<<endl;
	*/
		
	//cout<<"de_number =  "<<de_number<<endl;
	//cout<<"T =  "<<T<<endl;

	for (int counter = 0;counter<points;counter++)
	{
		nu_deformation[counter] = 0;
		nu_piezoelectric[counter] = 0;
	}
	
	if(screening==0)
	{
		double aa[limit2];
		for (int counter = 0;counter<points;counter++)
		{
			k = k_grid[counter]*1e9;   // converted from 1/nm to 1/m					
	    		v = v_n[counter]*1e-2;	      // converted from cm/s to m/s	
			
			//aa[counter] = pi*(2*pow(a_n[counter],4) + pow(c_n[counter],4) - 2*pow(a_n[counter],2)*pow(c_n[counter],2));
			
			aa[counter] = 2*pi;
			
	    		nu_def[0][counter] = ((k_B*e)*T*pow(E_deformation[0],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_long*1e-3))*
	    		aa[counter];
			// E_deformation[0] is in eV and h_bar is in eV-s, so e^2 is both numerator and denominator cancelled 
			// C_long is in dyne/cm so converted into mks N/m

			if(de_number==2)		
			{
	    		nu_def[1][counter] = ((k_B*e)*T*pow(E_deformation[1],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_trans*1e-3))*
	    		aa[counter];
			}
			
			if(de_number==3)		
			{

	    		nu_def[1][counter] = ((k_B*e)*T*pow(E_deformation[1],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_trans*1e-3))*
	    		aa[counter];

	    		nu_def[2][counter] = ((k_B*e)*T*pow(E_deformation[2],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_za*1e-3))*
	    		aa[counter];

			}
			
		}
	}
	else
	{	
		double pz1[limit7+1], intg, dtheta, ep, G;
		double X[limit7+1], Y[limit7+1], Z[limit7+1], theta[limit7+1], aa[limit7+1], q1[limit7+1];    
		// limit7 is for theta variation
		
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
			
			for(int j=0;j<limit7+1;j++)			
			{
				aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
				aa[j] = aa[j] * aa[j];
			}
								
			int index, last_index=0;
			double min=1e20;
			
			// polarizibility is calculated for different values of theta 			
			for(int j=0;j<=limit7;j++)
			{				
						
				min=1e15;
				/*	
				if(j!=0)
				{
					if(q[j]==q[j-1] )
					{
						index = last_index;	
						pz1[j] = pz[index];
						continue;
					}
				}
				*/
								
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
				intg = intg + aa[j]*Y[j]/(ep*ep)*dtheta;	
			}
			
			
	    		nu_def[0][i] = ((k_B*e)*T*pow(E_deformation[0],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_long*1e-3))*intg;
			// E_deformation[0] is in eV and h_bar is in eV-s, so e^2 is both numerator and denominator cancelled 
			// C_long is in dyne/cm so converted into mks

			if(de_number==2)		
	    			nu_def[1][i] = ((k_B*e)*T*pow(E_deformation[1],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_trans*1e-3))*intg;
			
			
			if(de_number==3)		
			{

	    		nu_def[1][i] = ((k_B*e)*T*pow(E_deformation[1],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_trans*1e-3))*intg;

	    		nu_def[2][i] = ((k_B*e)*T*pow(E_deformation[2],2)*k)/(2*pi*(h_bar*h_bar)*v*(C_za*1e-3))*intg;

			}
						
											
		} // loop for different values of k		
	} // else part with screening completed
	
	for(int counter=0;counter<points;counter++)
	{
		for(int i=0;i<de_number;i++)
		{	nu_deformation[counter] = nu_deformation[counter] + nu_def[i][counter];
			//cout<<"nu_deformation[counter] =  "<<nu_deformation[counter]<<endl;
		}
	}
		
//----------------------------------------------- acoustic scattering rate completed ----------------------------------


//----------------------------------------------- pz scattering rate  ----------------------------------
	
	// 
	if(scattering_mechanisms[4]==1)   // pz scattering
	{
		
		//cout<<"e11 C/cm   = "<<P_piezo[0]<<" C/cm"<<endl;
		
		double const_pz[de_number];
		
		for(int i=0;i<de_number;i++)
		{
			const_pz[i] = ((P_piezo[0]*1e2)/(epsilon_0*E_deformation[i]));
			const_pz[i] = const_pz[i] * const_pz[i]/2.0;
			
			for (int counter = 0;counter<points;counter++)
				nu_pz[i][counter] = nu_def[i][counter]*const_pz[i];
		}
		
		for (int counter = 0;counter<points;counter++)
		{
			for(int i=0;i<de_number;i++)
				nu_piezoelectric[counter] = nu_piezoelectric[counter] + nu_pz[i][counter];
		}
	}
	
	//*/
}


