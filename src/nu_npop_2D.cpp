
#include"main.h"

// 
void nu_npop_2D(double T) 
{

	for (int i = 0;i<points;i++)
		nu_npop_total[i] = 0; 	

	double k, v;
	
	double X[limit7+1], Y[limit7+1], Z[limit7+1], theta[limit7+1], aa[limit7+1], q_ab[limit7+1], q_em[limit7+1];    
	double arr[points], dtheta;
	int minus_index, plus_index;
	double E_plus, E_minus, k_plus, k_minus, f, fp, fm, E1, Nv, npopconst, J_plus, J_minus;
	double I_plus, I_minus, v_plus, v_minus;
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

	if(screening==0)
	{
		
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
		

		for(int m3=0; m3<npop_number; m3++)
		{
			Nv = 1/(exp((h_bar*we_npop[m3])/(k_B*T))-1);
			npopconst = (De_npop[m3]*De_npop[m3]*1e10*1e10*e*e)*gd[m3]/(2*pi*(h_bar*e)*2*rho*we_npop[m3]);
			
			//cout<<"Nv = "<<Nv<<endl;
			//cout<<"npopconst = "<<npopconst<<endl;
			//getchar();
	
			for (int i = 0;i<points;i++)
			{
				//cout<<"i = "<<i<<endl;
				//cout<<"k[i] = "<<k_grid[i]<<endl;
				

				k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
		    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	
				
						
				//-------- overlap integral calculated for different theta------------------------------------
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
					
				//-------- overlap integral calculated completed for different theta------------------------------------

				//--------- q_ab and q_em is calculated for different values of theta ----------------------
				E1 = energy_n[i];
				
				E_plus =  E1 + h_bar*we_npop[m3];
				if(E_plus > energy_n[points-1])
					E_plus =  0;
					
				E_minus = E1- h_bar*we_npop[m3];
				if(E_minus < 0)
					E_minus = 0;
												
				f = 1/(1+exp(((E1 - efef_n)/(k_B*T))));
				fp = 1/(1+exp((E_plus - efef_n)/(k_B*T)));
				fm = 1/(1+exp((E_minus - efef_n)/(k_B*T)));

				for (int i1=0;i1<points;i1++)
					arr[i1] = abs(energy_n[i1] - E_plus);

				plus_index =FindMinInd(arr,points);

				for (int i1=0;i1<points;i1++)
					arr[i1] = abs(energy_n[i1] - E_minus);
				minus_index =FindMinInd(arr,points);					
				
				k_plus = k_grid[plus_index]*1e9;
				k_minus = k_grid[minus_index]*1e9;
		    		v_plus = v_n[plus_index]*1e-2;	      // converted from cm/s to m/s	
		    		v_minus = v_n[minus_index]*1e-2;	      // converted from cm/s to m/s	
				
				/*
				cout<<"i = "<<i<<endl;
				cout<<"E1 = "<<E1<<endl;
				cout<<"E_plus = "<<E_plus<<endl;
				cout<<"E_minus = "<<E_minus<<endl;
				cout<<"h_bar*we_npop[m3] ="<<h_bar*we_npop[m3]<<endl;
				cout<<"k_plus = "<<k_plus<<endl;
				cout<<"k_minus = "<<k_minus<<endl;
				cout<<"plus_index = "<<plus_index<<endl;
				cout<<"minus_index = "<<minus_index<<endl;
				getchar();
				//*/
								
				// ---------integral is calculated next for scattering rate for different kpoints----	
				// --- for absorption --
				I_plus = 0;
				J_plus = 0; 			
				for(int j=0;j<=limit7;j++)
				{							
					I_plus = I_plus + dtheta*aa[j];
					J_plus = J_plus + dtheta*aa[j]*Z[j];	
				}
				
				// --- for emission ----  

				I_minus = 0;
				J_minus = 0;
				// scattering rate is calculated 			
				for(int j=0;j<=limit7;j++)
				{							
					I_minus = I_minus + dtheta*aa[j];
					J_minus = J_minus + dtheta*aa[j]*Z[j];	
				}
	    			//
	    			
	    			/*
				cout<<"J_plus =  "<<J_plus<<endl;
				cout<<"J_minus =  "<<J_minus<<endl;
				cout<<"I_plus =  "<<I_plus<<endl;
				cout<<"I_minus =  "<<I_minus<<endl;
				*/
	    			// ----------------- out scattering rate is calculated here -------------------------
	    			if(E_plus < energy_n[points-1])
					So_ab_npop[m3][i] = npopconst*(Nv)*(1 - fp)*I_plus*k_plus/(v_plus*(1-f));
				else
					So_ab_npop[m3][i] = 0;
					
				if(E1 > h_bar*we_npop[m3])
					So_em_npop[m3][i] = npopconst*(Nv+1)*(1 - fm)*I_minus*k_minus/(v_minus*(1-f));
				else
					So_em_npop[m3][i] = 0;
				
				nu_npop[m3][i] = So_ab_npop[m3][i] + So_em_npop[m3][i];
				
				//--------- in scattering terms calculated here ----------------------------------
				if(E1 > h_bar*we_npop[m3])
				Sa_npop[m3][i] = npopconst*(Nv+1)*(1 - fm)*J_minus*sqrt(E_minus/E1)*k_minus/(v_minus*(1-f));
				else
					Sa_npop[m3][i] = 0;
					
				if(E_plus < energy_n[points-1])
				Se_npop[m3][i] = npopconst*(Nv)*(1 - fp)*J_plus*sqrt(E_plus/E1)*k_plus/(v_plus*(1-f));	
				else
					Se_npop[m3][i] = 0;
				
				/*
				cout<<"So_ab_npop[m3][i] = "<<So_ab_npop[m3][i]<<endl;
				cout<<"So_em_npop[m3][i] = "<<So_em_npop[m3][i]<<endl;
				cout<<"nu_npop[m3][i] = "<<nu_npop[m3][i]<<endl;
				cout<<"Sa_npop[m3][i] = "<<Sa_npop[m3][i]<<endl;
				cout<<"Se_npop[m3][i] = "<<Se_npop[m3][i]<<endl;

				getchar();
				//*/
				
				//--------------------------------------------------------------------------------
			} // loop for different k points finished here			
		} // loop for different npop constants finished here
		
		/*	
		for(int m3=0; m3<npop_number; m3++)
		{
			// to remove nan due to E1=0
			Sa_npop[m3][0] = Sa_npop[m3][1];

			Se_npop[m3][0] = Se_npop[m3][1];
		}
		*/
		
		for (int i = 0;i<points;i++)
		{
			for(int m3=0; m3<npop_number; m3++)
				nu_npop_total[i] = nu_npop_total[i] + nu_npop[m3][i];
			
			/*
			cout<<"nu_npop[m3][i] = "<<nu_npop[0][i]<<endl;

 			cout<<"nu_npop_total[i] = "<<nu_npop_total[i]<<endl;
			getchar();
			*/
		}									
	}
	else
	{
		double pz_ab[limit7+1], pz_em[limit7+1], G_em, G_ab, ep_ab, ep_em;		

		for(int m3=0; m3<npop_number; m3++)
		{
			Nv = 1/(exp((h_bar*we_npop[m3])/(k_B*T))-1);
			npopconst = (De_npop[m3]*De_npop[m3]*1e10*1e10*e*e)*gd[m3]/(2*pi*(h_bar*e)*2*rho*we_npop[m3]);
			
			/*
			cout<<"Nv = "<<Nv<<endl;
			cout<<"npopconst = "<<npopconst<<endl;
			getchar();
			//*/
			
			for (int i = 0;i<points;i++)
			{
				/*
				cout<<"i = "<<i<<endl;
				cout<<"k[i] = "<<k_grid[i]<<endl;
				getchar();
				//*/
				
				k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
		    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	
				
						
				//-------- overlap integral calculated for different theta------------------------------------
				if(overlap==1)
				{
					for(int j=0;j<limit7+1;j++)			
					{
						aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
						aa[j] = aa[j] * aa[j];
					}
				}
				//-------- overlap integral calculated completed for different theta------------------------------------
				else
				{
					for(int j=0;j<limit7+1;j++)			
						aa[j] = 1;
				}					
				//--------- q_ab and q_em is calculated for different values of theta ----------------------
				E1 = energy_n[i];
				
				E_plus =  E1 + h_bar*we_npop[m3];
				if(E_plus > energy_n[points-1])
					E_plus =  0;
					
				E_minus = E1- h_bar*we_npop[m3];
				if(E_minus < 0)
					E_minus = 0;
												
				f = 1/(1+exp(((E1 - efef_n)/(k_B*T))));
				fp = 1/(1+exp((E_plus - efef_n)/(k_B*T)));
				fm = 1/(1+exp((E_minus - efef_n)/(k_B*T)));

				for (int i1=0;i1<points;i1++)
					arr[i1] = abs(energy_n[i1] - E_plus);
				plus_index =FindMinInd(arr,points);

				for (int i1=0;i1<points;i1++)
					arr[i1] = abs(energy_n[i1] - E_minus);
				minus_index =FindMinInd(arr,points);					
				
				k_plus = k_grid[plus_index]*1e9;
				k_minus = k_grid[minus_index]*1e9;
		    		v_plus = v_n[plus_index]*1e-2;	      // converted from cm/s to m/s	
		    		v_minus = v_n[minus_index]*1e-2;	      // converted from cm/s to m/s	
				
				/*
				cout<<"i = "<<i<<endl;
				cout<<"E1 = "<<E1<<endl;
				cout<<"E_plus = "<<E_plus<<endl;
				cout<<"E_minus = "<<E_minus<<endl;
				cout<<"h_bar*we_npop[m3] ="<<h_bar*we_npop[m3]<<endl;
				cout<<"k_plus = "<<k_plus<<endl;
				cout<<"k_minus = "<<k_minus<<endl;
				cout<<"plus_index = "<<plus_index<<endl;
				cout<<"minus_index = "<<minus_index<<endl;
				getchar();
				//*/
				
				// change in wave vector for different values of theta
				for(int j=0;j<limit7+1;j++)			
				{				
					q_ab[j] = sqrt(k*k + k_plus*k_plus - 2*k*k_plus*Z[j]);            
					q_em[j] = sqrt(k*k + k_minus*k_minus - 2*k*k_minus*Z[j]);				
	    			}
				//--------- q_ab and q_em is calculated completed for different values of theta -----		
				
				//------------ polarizibility is calculated for different values of theta --------		
				
				//-- for absorption first ---			
				int index, last_index=0;
				double min=1e20;
				for(int j=0;j<=limit7;j++)
				{				
					// ---- search for q_ab[j]------------------------------		
					min=1e15;
									
					for(int l=last_index;l<=limit6;)   // limit6+1 size of polarizability array
					{
						if(min > abs(q_ab[j]-q[l]))
						{
							index = l;
							min = abs(q_ab[j]-q[l]);		
						}
						else
						{
							last_index = index;
							break;
						}

						if (q_ab[j] > q_ab[j-1])
							l++;						
						else
							l--;					
					}
					pz_ab[j] = pz[index];  // polarizibility for different theta				
					// ---- search for q_ab[j] finished and pz_ab calculated -----------------------------		
						
				} // loop for theta for polarizibility ab finished here
				
				// --for emission --			
				last_index=0;
				
				for(int j=0;j<=limit7;j++)
				{				
					// ---- search for q_em[j]------------------------------		
					min=1e15;
									
					for(int l=last_index;l<=limit6;)   // limit6+1 size of polarizability array
					{
						if(min > abs(q_em[j]-q[l]))
						{
							index = l;
							min = abs(q_em[j]-q[l]);		
						}
						else
						{
							last_index = index;
							break;
						}

						if (q_em[j] > q_em[j-1])
							l++;						
						else
							l--;					
					}
					pz_em[j] = pz[index];  // polarizibility for different theta				
					// ---- search for q_em[j] finished and pz_em calculated ---------		
						
				} // loop for theta for polarizibility em finished here
				//-- for emission  ----			
				
				//------------ completed polarizibility is calculated for different values of theta ------------------

				// ---------integral is calculated next for scattering rate for different kpoints----	
				// --- for absorption --
				I_plus = 0;
				J_plus = 0; 			
				for(int j=0;j<=limit7;j++)
				{		
					G_ab = 1/(epsilon_0*(eps_sub_low + eps_up_low)*q_ab[j]);
					ep_ab = 1 + e*e*G_ab*pz_ab[j];  // screening factor
					
					I_plus = I_plus + dtheta*aa[j]/(ep_ab*ep_ab);
					J_plus = J_plus + dtheta*aa[j]*Z[j]/(ep_ab*ep_ab);	
				}
				
				// --- for emission ----  
				I_minus = 0;
				J_minus = 0;
				// scattering rate is calculated 			
				for(int j=0;j<=limit7;j++)
				{		
					G_em = 1/(epsilon_0*(eps_sub_low + eps_up_low)*q_em[j]);
					ep_em = 1 + e*e*G_em*pz_em[j];  // screening factor
					
					I_minus = I_minus + dtheta*aa[j]/(ep_em*ep_em);
					J_minus = J_minus + dtheta*aa[j]*Z[j]/(ep_em*ep_em);	
				}
	    			//
	    			
	    			/*
				cout<<"J_plus =  "<<J_plus<<endl;
				cout<<"J_minus =  "<<J_minus<<endl;
				cout<<"I_plus =  "<<I_plus<<endl;
				cout<<"I_minus =  "<<I_minus<<endl;
				//*/
				
	    			// ----------------- out scattering rate is calculated here -------------------------
	    			if(E_plus < energy_n[points-1])
					So_ab_npop[m3][i] = npopconst*(Nv)*(1 - fp)*I_plus*k_plus/(v_plus*(1-f));
				else
					So_ab_npop[m3][i] = 0;
					
				if(E1 > h_bar*we_npop[m3])
					So_em_npop[m3][i] = npopconst*(Nv+1)*(1 - fm)*I_minus*k_minus/(v_minus*(1-f));
				else
					So_em_npop[m3][i] = 0;
				
				nu_npop[m3][i] = So_ab_npop[m3][i] + So_em_npop[m3][i];
				
				//--------- in scattering terms calculated here ----------------------------------
				if(E1 > h_bar*we_npop[m3])
				Sa_npop[m3][i] = npopconst*(Nv+1)*(1 - fm)*J_minus*sqrt(E_minus/E1)*k_minus/(v_minus*(1-f));
				else
					Sa_npop[m3][i] = 0;
					
				if(E_plus < energy_n[points-1])
				Se_npop[m3][i] = npopconst*(Nv)*(1 - fp)*J_plus*sqrt(E_plus/E1)*k_plus/(v_plus*(1-f));	
				else
					Se_npop[m3][i] = 0;
				/*
				cout<<"So_ab_npop[m3][i] = "<<So_ab_npop[m3][i]<<endl;
				cout<<"So_em_npop[m3][i] = "<<So_em_npop[m3][i]<<endl;
				cout<<"nu_npop[m3][i] = "<<nu_npop[m3][i]<<endl;
				cout<<"Sa_npop[m3][i] = "<<Sa_npop[m3][i]<<endl;
				cout<<"Se_npop[m3][i] = "<<Se_npop[m3][i]<<endl;
				getchar();
				//*/
				
				//--------------------------------------------------------------------------------
			} // loop for different k points finished here			
		} // loop for different npop constants finished here
		
		/*	
		for(int m3=0; m3<npop_number; m3++)
		{
			// to remove nan due to E1=0
			Sa_npop[m3][0] = Sa_npop[m3][1];
			Se_npop[m3][0] = Se_npop[m3][1];
		}
		*/
		
		for (int i = 0;i<points;i++)
		{
			for(int m3=0; m3<npop_number; m3++)
				nu_npop_total[i] = nu_npop_total[i] + nu_npop[m3][i];
			
			/*
			cout<<"nu_npop[m3][i] = "<<nu_npop[0][i]<<endl;
 			cout<<"nu_npop_total[i] = "<<nu_npop_total[i]<<endl;
			getchar();
			*/
		}									
	}  // else condiction for sceerning finshed here			
	
	/*
	FILE *fid1;
	
	if(screening==0)
		fid1 = fopen("npop_scattering_rate.dat","w");
	else
		fid1 = fopen("npop_scattering_rate_screening.dat","w");
		
	
	fprintf(fid1,"Energy (eV)   scattering rates            total (last) \n ");
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		for(int m3=0; m3<npop_number; m3++)
			fprintf(fid1,"  %e  \t", nu_npop[m3][i] );

		fprintf(fid1,"  %e  \n", nu_npop_total[i] );
	}
		
	fclose(fid1);
	
	if(screening==0)		
		fid1 = fopen("npop_in_scattering_rate.dat","w");
	else
		fid1 = fopen("npop_in_scattering_rate_screening.dat","w");
	
	fprintf(fid1,"Energy (eV)   ab1   em1   ab2   em2 ..... \n");	
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		for(int m3=0; m3<npop_number; m3++)
			fprintf(fid1,"  %e   %e    \t", Sa_npop[m3][i], Se_npop[m3][i] );

		fprintf(fid1," \n" );

	}
		
	fclose(fid1);
	//*/
}


