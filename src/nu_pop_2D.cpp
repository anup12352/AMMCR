
#include"main.h"

void nu_pop_2D(double T, int T_loop) 
{
	//cout<<"nu_pop called"<<endl;
	//getchar();
	
	for (int i = 0;i<points;i++)
		nu_pop_total[i] = 0; 	

	double k, v;
	double X[limit7+1], Y[limit7+1], Z[limit7+1], theta[limit7+1], aa[limit7+1], q_ab[limit7+1], q_em[limit7+1], dtheta;    
	double arr[points];
	int minus_index, plus_index;
	double E_plus, E_minus, k_plus, k_minus, f, fp, fm, fm1, fpm1, fmm1, E1, Nv, popconst, J_plus, J_minus;
	double I_plus, I_minus, v_plus, v_minus;
	// limit7 is for theta variation
	
	dtheta = 2*pi/limit7;
	theta[0] = 0;
	X[0] = sin(theta[0]/2);
	Y[0] = 1 - cos(theta[0]);
	Z[0] = cos(theta[0]);
	
	int a,b;
	
	for(int i=1;i<limit7;i++)
	{
		theta[i] = theta[i-1] + dtheta;
		X[i] = sin(theta[i]/2);
		Y[i] = 1 - cos(theta[i]);
		Z[i] = cos(theta[i]);
		/*
		if(Z[i]==1 || Z[i]>1)
		{
			cout<<"Z[]   =  "<<Z[i]<<"   i   = "<<i<<endl;
			getchar();
		}
		*/
	}
	
	for(int m3=0; m3<pop_number; m3++)
	{			
		for (int i = 0;i<points;i++)
		{
			kplus_grid_pop[m3][i] = kplus(i,we_pop[m3],points,energy_n);
			kminus_grid_pop[m3][i] = kminus(i,we_pop[m3],points,energy_n);

			for (int ii=0;ii<points;ii++)
				arr[ii] = abs(k_grid[ii] - kminus_grid_pop[m3][i]);
			minus_index =FindMinInd(arr,points);

			for (int ii=0;ii<points;ii++)
				arr[ii] = abs(k_grid[ii] - kplus_grid_pop[m3][i]);
			plus_index =FindMinInd(arr,points);

			plus_index_pop[m3][i] = plus_index;
			minus_index_pop[m3][i] = minus_index;
		}
	}

	if(screening==0)
	{
		for(int m3=0; m3<pop_number; m3++)
		{
		
			Nv = 1/(exp((h_bar*we_pop[m3])/(k_B*T))-1);			
			popconst = 
			abs((e*e*we_pop[m3])*(1/(epsilon_inf[T_loop]*epsilon_0) - 1/(epsilon_s[T_loop]*epsilon_0))/(8.0*pi*h_bar*e));
			
			//cout<<" we_pop[m3] = "<<we_pop[m3]<<endl;
			//cout<<"  (1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop]) = "<<(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop])<<endl;
			//cout<<"Nv = "<<Nv<<endl;
			//cout<<"popconst = "<<popconst<<endl;
			//getchar();
	
			for (int i = 0;i<points;i++)
			{

				//cout<<"i = "<<i<<endl;
				//cout<<"k[i] = "<<k_grid[i]<<endl;				

				k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
		    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	
										

				//--------- q_ab and q_em is calculated for different values of theta ----------------------


				E1 = energy_n[i];
				E_plus =  energy_n[plus_index_pop[m3][i]];					
				E_minus = energy_n[minus_index_pop[m3][i]];
												
				f = 1/(1+exp(((E1 - efef_n)/(k_B*T))));
				fp = 1/(1+exp((E_plus - efef_n)/(k_B*T)));
				fm = 1/(1+exp((E_minus - efef_n)/(k_B*T)));

				fm1 = exp(((E1 - efef_n)/(k_B*T)))/(1+exp(((E1 - efef_n)/(k_B*T))));    // 1-f
				fpm1 = exp((E_plus - efef_n))/(1+exp((E_plus - efef_n)/(k_B*T)));       //1-fp
				fmm1 = exp((E_minus - efef_n))/(1+exp((E_minus - efef_n)/(k_B*T)));     // 1-fm
				
				k_plus = k_grid[plus_index_pop[m3][i]]*1e9;
				k_minus = k_grid[minus_index_pop[m3][i]]*1e9;
		    		v_plus = v_n[plus_index_pop[m3][i]]*1e-2;	      // converted from cm/s to m/s	
		    		v_minus = v_n[minus_index_pop[m3][i]]*1e-2;	      // converted from cm/s to m/s	
				
				if(k_plus==k || k_minus==k)
				{
					a=limit7;
					b=1;
				}
				else
				{
					a=limit7;
					b = 0;
				}				
				
				/*
				cout<<"i = "<<i<<endl;
				cout<<"E1 = "<<E1<<endl;
				cout<<"E_plus = "<<E_plus<<endl;
				cout<<"E_minus = "<<E_minus<<endl;
				cout<<"h_bar*we_pop[m3] ="<<h_bar*we_pop[m3]<<endl;
				cout<<"k_plus = "<<k_plus<<endl;
				cout<<"k_minus = "<<k_minus<<endl;
				cout<<"plus_index = "<<plus_index<<endl;
				cout<<"minus_index = "<<minus_index<<endl;
				getchar();
				//*/
				
				/*
				if(f==1 || fp==1)
				{
					cout<<"i = "<<i<<endl;
					cout<<"E1 = "<<E1<<endl;
					cout<<"E_plus = "<<E_plus<<endl;
					cout<<"E_minus = "<<E_minus<<endl;
					cout<<"h_bar*we_pop[m3] ="<<h_bar*we_pop[m3]<<endl;
					cout<<"k_plus = "<<k_plus<<endl;
					cout<<"k_minus = "<<k_minus<<endl;
					cout<<"plus_index = "<<plus_index<<endl;
					cout<<"minus_index = "<<minus_index<<endl;
					cout<<"efef_n    =  "<<efef_n<<endl;
					//cout<<"(E1 - efef_n)/(k_B*T))  =  "<<((E1 - efef_n)/(k_B*T))<<endl;
					//cout<<"((Eplus - efef_n)/(k_B*T))  =  "<<((E_plus - efef_n)/(k_B*T))<<endl;
					//cout<<"((Eminus - efef_n)/(k_B*T))  =  "<<((E_minus - efef_n)/(k_B*T))<<endl;
					
					getchar();
				}
				//*/
				//-------- overlap integral calculated for different theta------------------------------------
				if(overlap==1)
				{
					for(int j=b;j<a;j++)			
					{
						aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
						aa[j] = aa[j] * aa[j];
					}
				}
				else // this is by default run always
				{
					for(int j=b;j<a;j++)			
						aa[j] = 1;
				}					
				//-------- overlap integral calculated completed for different theta------------------------------------


				// change in wave vector for different values of theta
				for(int j=b;j<a;j++)			
				{	
					//if(k_plus!=k)			
					q_ab[j] = sqrt(k*k + k_plus*k_plus - 2*k*k_plus*Z[j]);            
					q_em[j] = sqrt(k*k + k_minus*k_minus - 2*k*k_minus*Z[j]);
					
					/*
					if(q_ab[j]==0 || q_em[j]==0)	
					{
						cout<<"j     =    "<<j<<"  q_ab[j] = "<<q_ab[j]<<endl;
						cout<<"j     =    "<<j<<"  q_em[j] = "<<q_em[j]<<endl;
						cout<<"k   =  "<<k<<"   k_minus   = "<<k_minus<<"  k_plus   =   "<<k_plus<<"   Z[j]  =  "<<Z[j]<<endl;
						getchar();
					}	
					//*/
													
	    			}
				//--------- q_ab and q_em is calculated completed for different values of theta -----		
				
				// ---------integral is calculated next for scattering rate for different kpoints----	
				// --- for absorption --
				I_plus = 0;
				J_plus = 0; 			
				for(int j=b;j<a;j++)
				{							
					I_plus = I_plus + dtheta*aa[j]/(q_ab[j]);
					J_plus = J_plus + dtheta*aa[j]*Z[j]/(q_ab[j]);	
					
					/*
					if(isnan(J_plus))
					{
						cout<<"isnan J_plus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_ab[j]   = "<<q_ab[j]<<endl;
						getchar();
					}
					if(isinf(J_plus))
					{
						cout<<"isinf J_plus =    i =  "<<i<<"  m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_ab[j]   = "<<q_ab[j]<<endl;
						getchar();
					}

					if(isnan(J_plus))
					{
						cout<<"isnan J_plus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_ab[j]   = "<<q_ab[j]<<endl;
						getchar();
					}
					if(isnan(J_plus))
					{
						cout<<"isnan J_plus =    i =  "<<i<<"  m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_ab[j]   = "<<q_ab[j]<<endl;
						getchar();
					}

					//*/

				}
				
				// --- for emission ----  
				I_minus = 0;
				J_minus = 0;
				// scattering rate is calculated 			
				for(int j=b;j<a;j++)
				{							
					I_minus = I_minus + dtheta*aa[j]/(q_em[j]);
					J_minus = J_minus + dtheta*aa[j]*Z[j]/(q_em[j]);	
					
					/*
					if(isnan(J_minus))
					{
						cout<<"isnan J_minus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"aa[j] =    i =  "<<i<<"  m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_em[j]   = "<<q_em[j]<<endl;
						getchar();
					}
					if(isinf(J_minus))
					{
						cout<<"isinf J_minus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_em[j]   = "<<q_em[j]<<endl;
						getchar();
					}
					if(isnan(I_minus))
					{
						cout<<"isnan I_minus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"aa[j] =    i =  "<<i<<"  m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_em[j]   = "<<q_em[j]<<endl;
						getchar();
					}
					if(isinf(I_minus))
					{
						cout<<"isinf I_minus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
						cout<<"j   = "<<j<<"   Z[j]   =   "<<Z[j]<<"     q_em[j]   = "<<q_em[j]<<endl;
						getchar();
					}
					//*/
					
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
					So_ab_pop[m3][i] = popconst*(Nv)*fpm1*I_plus*k_plus/(v_plus*(fm1));
				else
					So_ab_pop[m3][i] = 0;
					
				if(E1 > h_bar*we_pop[m3])
					So_em_pop[m3][i] = popconst*(Nv+1)*(1 - fm)*I_minus*k_minus/(v_minus*(fm1));
				else
					So_em_pop[m3][i] = 0;
					
				/*
				if(isinf(So_ab_pop[m3][i]))
				{
					cout<<"i   =   "<<i<<endl;
					cout<<"isinf So_ab_pop numer   = "<<popconst*(Nv)*fpm1*I_plus*k_plus<<"   denom = "<<(v_plus*(fm1))<<endl;
					cout<<"I_plus = "<<I_plus<<endl;
					cout<<"k_plus = "<<k_plus<<endl;
					cout<<"v_plus = "<<v_plus<<endl;
					cout<<"fm1 = "<<(fm1)<<endl;    // 1-f
					cout<<"fpm1 = "<<fpm1<<endl;    // 1-fp
					
					getchar();
				}
					
				if(isnan(So_ab_pop[m3][i]))
				{
					cout<<"i   =   "<<i<<endl;
					cout<<"isnan So_ab_pop numer   = "<<popconst*(Nv)*fpm1*I_plus*k_plus<<"   denom = "<<(v_plus*(fm1))<<endl;
					cout<<"I_plus = "<<I_plus<<endl;
					cout<<"k_plus = "<<k_plus<<endl;
					cout<<"v_plus = "<<v_plus<<endl;
					cout<<"fm1 = "<<(fm1)<<endl;     // 1-f
					cout<<"fpm1 = "<<fpm1<<endl;     // 1-fp
					
					getchar();
				}

				//*/				
				//--------- in scattering terms calculated here ----------------------------------
				
				// mrethod=1 always fixed				
				if(method==1)   // rode direct method
				{
					//cout<<"Direct rode method is selected "<<endl;				
					if(E1 > h_bar*we_pop[m3])
						Sa_pop[m3][i] = popconst*(Nv+1)*(f)*J_minus*k_minus/(v_minus*(fm));		
					else
						Sa_pop[m3][i] = 0;
						
					if(E_plus < energy_n[points-1])
						Se_pop[m3][i] = popconst*(Nv)*(f)*J_plus*k_plus/(v_plus*(fp));		
					else	
						Se_pop[m3][i] = 0;
				}
				else   // alternate rode direct method
				{
				//cout<<"ALternate rode method is selected "<<endl;
				if(E1 > h_bar*we_pop[m3])
				Sa_pop[m3][i] = popconst*(Nv+1)*(1 - fm)*J_minus*sqrt(E_minus/E1)*k_minus/(v_minus*(fm1));		
				else
					Sa_pop[m3][i] = 0;
					
				if(E_plus < energy_n[points-1])
				Se_pop[m3][i] = popconst*(Nv)*fpm1*J_plus*sqrt(E_plus/E1)*k_plus/(v_plus*(fm1));		
				else	
					Se_pop[m3][i] = 0;
				
				}
				
				
				/*
				if(isnan(J_minus))
				{
					cout<<"isnan J_minus =    i =  "<<i<<"  m3  =  "<<m3<<endl;
					getchar();
				}
				if(isinf(J_minus))
				{
					cout<<"isinf J_minus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
					getchar();
				}

				if(isnan(J_plus))
				{
					cout<<"isnan J_plus =    i =  "<<i<<"   m3  =  "<<m3<<endl;
					getchar();
				}
				if(isinf(J_plus))
				{
					cout<<"isinf J_plus =    i =  "<<i<<"  m3  =  "<<m3<<endl;
					getchar();
				}
				//*/
									
				/*
				cout<<"So_ab_pop[m3][i] = "<<So_ab_pop[m3][i]<<endl;
				cout<<"So_em_pop[m3][i] = "<<So_em_pop[m3][i]<<endl;

				cout<<"So_pop[m3][i] = "<<So_pop[m3][i]<<endl;
				cout<<"Sa_pop[m3][i] = "<<Sa_pop[m3][i]<<endl;
				cout<<"Se_pop[m3][i] = "<<Se_pop[m3][i]<<endl;
				getchar();
				//*/
				
				//--------------------------------------------------------------------------------
			} // loop for different k points finished here			
		} // loop for different pop constants finished here		
	}
	else    // with screening
	{		
		double pz_ab[limit7+1], pz_em[limit7+1], G_em, G_ab, ep_ab, ep_em;
		for(int m3=0; m3<pop_number; m3++)
		{
		
			Nv = 1/(exp((h_bar*we_pop[m3])/(k_B*T))-1);			
        	popconst = abs((e*e*we_pop[m3])*(1/(epsilon_inf[T_loop]*epsilon_0) - 1/(epsilon_s[T_loop]*epsilon_0))/(8.0*pi*h_bar*e));
			
			//cout<<" we_pop[m3] = "<<we_pop[m3]<<endl;
			//cout<<"  (1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop]) = "<<(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop])<<endl;
			//cout<<"Nv = "<<Nv<<endl;
			//cout<<"popconst = "<<popconst<<endl;
			//getchar();
	
			for (int i = 0;i<points;i++)
			{
				//cout<<"i = "<<i<<endl;
				//cout<<"k[i] = "<<k_grid[i]<<endl;
				

				k = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
		    		v = v_n[i]*1e-2;	      // converted from cm/s to m/s	
								
				//--------- q_ab and q_em is calculated for different values of theta ----------------------
				E1 = energy_n[i];
				E_plus =  energy_n[plus_index_pop[m3][i]];					
				E_minus = energy_n[minus_index_pop[m3][i]];
																
				f = 1/(1+exp(((E1 - efef_n)/(k_B*T))));
				fp = 1/(1+exp((E_plus - efef_n)/(k_B*T)));
				fm = 1/(1+exp((E_minus - efef_n)/(k_B*T)));
				
				fm1 = exp(((E1 - efef_n)/(k_B*T)))/(1+exp(((E1 - efef_n)/(k_B*T))));    // 1-f
				fpm1 = exp((E_plus - efef_n))/(1+exp((E_plus - efef_n)/(k_B*T)));       //1-fp
				fmm1 = exp((E_minus - efef_n))/(1+exp((E_minus - efef_n)/(k_B*T)));     // 1-fm

				k_plus = k_grid[plus_index_pop[m3][i]]*1e9;
				k_minus = k_grid[minus_index_pop[m3][i]]*1e9;
		    		v_plus = v_n[plus_index_pop[m3][i]]*1e-2;	      // converted from cm/s to m/s	
		    		v_minus = v_n[minus_index_pop[m3][i]]*1e-2;	      // converted from cm/s to m/s	
				

				/*
				cout<<"i = "<<i<<endl;
				cout<<"E1 = "<<E1<<endl;
				cout<<"E_plus = "<<E_plus<<endl;
				cout<<"E_minus = "<<E_minus<<endl;
				cout<<"h_bar*we_pop[m3] ="<<h_bar*we_pop[m3]<<endl;
				cout<<"k_plus = "<<k_plus<<endl;
				cout<<"k_minus = "<<k_minus<<endl;
				cout<<"plus_index = "<<plus_index<<endl;
				cout<<"minus_index = "<<minus_index<<endl;
				getchar();
				//*/
				
				if(k_plus==k || k_minus==k)
				{
					a=limit7;
					b=1;
				}
				else
				{
					a=limit7;
					b = 0;
				}				

				//-------- overlap integral calculated for different theta------------------------------------
				if(overlap==1)
				{
					for(int j=b;j<a;j++)			
					{
						aa[j] = pow(a_n[i],2) + pow(c_n[i],2)*Z[j] ;
						aa[j] = aa[j] * aa[j];
					}
					//-------- overlap integral calculated completed for different theta --------------------
				}
				else
				{
					for(int j=b;j<a;j++)			
						aa[j] = 1;
				}					
				


				// change in wave vector for different values of theta
				for(int j=b;j<a;j++)			
				{				
					q_ab[j] = sqrt(k*k + k_plus*k_plus - 2*k*k_plus*Z[j]);            
					q_em[j] = sqrt(k*k + k_minus*k_minus - 2*k*k_minus*Z[j]);				
	    			}
				//--------- q_ab and q_em is calculated completed for different values of theta -----		
				
				//------------ polarizibility is calculated for different values of theta --------		
				
				//-- for absorption first ---			
				int index, last_index=0;
				double min=1e20;
				for(int j=b;j<a;j++)
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
					// ---- search for q_ab[j] finished and pz_ab calculated -----------------	
						
				} // loop for theta for polarizibility ab finished here
				
				// --for emission --			
				last_index=0;
				
				for(int j=b;j<a;j++)
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
				for(int j=b;j<a;j++)
				{		
					G_ab = 1/(epsilon_0*(eps_sub_low + eps_up_low)*q_ab[j]);
					ep_ab = 1 + e*e*G_ab*pz_ab[j];  // screening factor
					
					I_plus = I_plus + dtheta*aa[j]/(ep_ab*ep_ab*q_ab[j]);
					J_plus = J_plus + dtheta*aa[j]*Z[j]/(ep_ab*ep_ab*q_ab[j]);	
				}
				
				// --- for emission ----  
				I_minus = 0;
				J_minus = 0;
				// scattering rate is calculated 			
				for(int j=b;j<a;j++)
				{		
					G_em = 1/(epsilon_0*(eps_sub_low + eps_up_low)*q_em[j]);
					ep_em = 1 + e*e*G_em*pz_em[j];  // screening factor
					
					I_minus = I_minus + dtheta*aa[j]/(ep_em*ep_em*q_em[j]);
					J_minus = J_minus + dtheta*aa[j]*Z[j]/(ep_em*ep_em*q_em[j]);	
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
					So_ab_pop[m3][i] = popconst*(Nv)*fpm1*I_plus*k_plus/(v_plus*(fm1));
				else
					So_ab_pop[m3][i] = 0;
					
				if(E1 > h_bar*we_pop[m3])
					So_em_pop[m3][i] = popconst*(Nv+1)*(1 - fm)*I_minus*k_minus/(v_minus*(fm1));
				else
					So_em_pop[m3][i] = 0;
								
				//--------- in scattering terms calculated here ----------------------------------
								
				if(E1 > h_bar*we_pop[m3])
				Sa_pop[m3][i] = popconst*(Nv+1)*(1 - fm)*J_minus*sqrt(E_minus/E1)*k_minus/(v_minus*(fm1));		
				else
					Sa_pop[m3][i] = 0;
					
				if(E_plus < energy_n[points-1])
				Se_pop[m3][i] = popconst*(Nv)*fpm1*J_plus*sqrt(E_plus/E1)*k_plus/(v_plus*(fm1));		
				else	
					Se_pop[m3][i] = 0;
					
				/*
				cout<<"So_ab_pop[m3][i] = "<<So_ab_pop[m3][i]<<endl;
				cout<<"So_em_pop[m3][i] = "<<So_em_pop[m3][i]<<endl;
				cout<<"So_pop[m3][i] = "<<So_pop[m3][i]<<endl;
				cout<<"Sa_pop[m3][i] = "<<Sa_pop[m3][i]<<endl;
				cout<<"Se_pop[m3][i] = "<<Se_pop[m3][i]<<endl;
				getchar();
				//*/
				
				//--------------------------------------------------------------------------------
			} // loop for different k points finished here			
		} // loop for different pop constants finished here
				
	}  // else condiction for sceerning finshed here			

	/*
	for(int m3=0; m3<pop_number; m3++)
	{
		// to remove nan due to E1=0
		Sa_pop[m3][0] = Sa_pop[m3][1];
		Se_pop[m3][points-1] = Se_pop[m3][points-2];
		
		So_ab_pop[m3][points-1] = So_ab_pop[m3][points-2];
		So_em_pop[m3][0] = So_em_pop[m3][1];
		
	}
	//*/

	
	for (int i = 0;i<points;i++)
	{
		for(int m3=0; m3<pop_number; m3++)
		{
			So_pop[m3][i] = So_ab_pop[m3][i] + So_em_pop[m3][i];
			Si_pop[m3][i] = Sa_pop[m3][i] + Se_pop[m3][i];
		}
		
		for(int m3=0; m3<pop_number; m3++)
			nu_pop_total[i] = nu_pop_total[i] + So_pop[m3][i];
	}
	
	/*
	FILE *fid1;
	fid1 = fopen("in_scattering.dat","w");
	fprintf(fid1,"# points       Se_pop[0][i], Se_pop[1][i], Se_pop[2][i], Sa_pop[0][i], Sa_pop[1][i], Sa_pop[2][i] \n");
	for (int i = 0;i<points;i++)
	{
		fprintf(fid1," %d   %e    %e    %e    %e    %e    %e   \n", i, Se_pop[0][i], Se_pop[1][i], Se_pop[2][i], 
		Sa_pop[0][i], Sa_pop[1][i], Sa_pop[2][i] );
	}
	fclose(fid1);
	//*/
	
	/*
	for (int i = 0;i<points;i++)
	{
		if(isinf(nu_pop_total[i]))
		{
					
			cout<<" i  =  "<<i<<endl;
			for(int m3=0; m3<pop_number; m3++)
				cout<<"m3   =   "<<m3<<"    So_pop[m3][i] = "<<So_pop[m3][i]<<"  So_ab_pop[m3][i]   =  "<<So_ab_pop[m3][i]<<"  So_em_pop[m3][i]   =  "<<So_em_pop[m3][i]<<endl;
				
			cout<<"nu_pop_total[i] = "<<nu_pop_total[i]<<endl;
			
						
			nu_pop_total[i] = nu_pop_total[i-1];
						
			//getchar();
		}
	}
	//*/
			
	/*
	for(int m3=0; m3<pop_number; m3++)
	{
		if(isnan(Se_pop[m3][i]))
		{
			cout<<"isnan for Se    i   =  "<<i<<"   m3  =  "<<m3<<endl;
			getchar();
		}
		if(isnan(Sa_pop[m3][i]))
		{
			cout<<"isnan for Sa    i   =  "<<i<<"   m3  =  "<<m3<<endl;
			getchar();
		}
	
		if(isinf(Se_pop[m3][i]))
		{
			cout<<"isinf for Se    i   =  "<<i<<"   m3  =  "<<m3<<endl;
			getchar();
		}
		if(isinf(Sa_pop[m3][i]))
		{
			cout<<"isinf for Sa    i   =  "<<i<<"   m3  =  "<<m3<<endl;
			getchar();
		}
	
	}
	//*/
		
										
	
	/*
	FILE *fid1;
	if(screening==0)
		fid1 = fopen("pop_scattering_rate.dat","w");
	else
		fid1 = fopen("pop_scattering_rate_screening.dat","w");
	
	
	fprintf(fid1,"Energy (eV)   ab   em   total   ....  scattering rates ....      total (last) \n ");
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		for(int m3=0; m3<pop_number; m3++)
			fprintf(fid1,"  %e	%e	%e  \t", So_ab_pop[m3][i], So_em_pop[m3][i], So_pop[m3][i] );

		fprintf(fid1,"  %e  \n", nu_pop_total[i] );
	}
		
	fclose(fid1);
	
	if(screening==0)		
		fid1 = fopen("pop_in_scattering_rate.dat","w");
	else
		fid1 = fopen("pop_in_scattering_rate_screening.dat","w");

	fprintf(fid1,"Energy (eV)   ab1   em1   ab2   em2 ..... \n");	
	for (int i = 0; i < points; i++)		
	{
		fprintf(fid1,"  %e \t", energy_n[i]);
		for(int m3=0; m3<pop_number; m3++)
			fprintf(fid1,"  %e   %e    \t", Sa_pop[m3][i], Se_pop[m3][i] );

		fprintf(fid1," \n" );

	}
		
	fclose(fid1);
	//*/
}


