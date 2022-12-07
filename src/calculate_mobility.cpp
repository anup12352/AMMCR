#include"main.h"

void calculate_mobility(double T, int T_loop, int d_loop)
{

	//cout<<"Reached here calculate_mobility "<<endl;
//----------------------------------calculate drift mobility -------------------------------------------------------------
            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);

            cout<<endl;
		cc = cc+1;

	double dummy[points]={0};
	
	FILE *fid;
	char line[1000];
	
	if(geometry==1)
	{
		if(type=="n")
		{	
		//--------------------------- conductivity with time ------------------------------------------------------------------------
			//cout<<"iterations = "<<iterations<<endl;
			
			if(time_variation==1)
			{
				
				double t, dt;

				//tau = 75e-15;
				//initial = +0.8e-12;
				dt = 1/omega_s; 
				/*
				fid = fopen("Efield_time.dat", "r");
				if (fid==NULL)
				{
				    if (fid==NULL)
				    {
					cout<<"Efield_time.dat file is not present. Exit from program";
					exit(EXIT_FAILURE);
				    }
				}
				fgets(line, 1000, (FILE*)fid);   // pass first line
				double dummy;
				*/
				
				for(int i=0;i<time_limit;i++)
				{
					//t = initial + i*dt + dt;
					//t = initial + i*dt ;
					cout<<"i = "<<i<<endl;
					Efield_time[i] = 1000;   // unit V/cm
					//Efield_time[i] = -100*exp(-t*t/(tau*tau))*(2*t*t/(tau*tau) - 1);   // unit V/cm
							
					//fgets(line, 1000, (FILE*)fid);   
					//sscanf(line, "%lf %lf ", &dummy, &Efield_time[i]);  
					
					conductivity_time(T,i);
				}
				
				
				//conductivity_freq();
				
				FILE *fid1;
				fid1 = fopen("time_variation.dat","w");
				fprintf(fid1,"#index	time		mobility(cm^2/(V-s))  sigma(S/cm)  J(A/cm^2)      Efield (V/cm)\n");	
				{	
				for (int i = 0; i < time_limit; i++)
					fprintf(fid1,"%d	%e	%e	%e	%e	%e \n", i+1, ((i+1)/omega_s), mobility_time[i], sigma_time[i], J_time[i], Efield_time[i]);
				fclose(fid1);

				}	
			}   // end of time variation
			
			
		//--------------------------- conductivity with time completed ----------------------------------------------------------
			/*
			if(freq_variation==1)
				conductivity_with_freq(T);
			*/

			    // ionized impurity scattering
			    if (scattering_mechanisms[0]==1)
			    {
				mobility_ii = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_ii = "<<mobility_ii<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_ii = 1e10;




			    // POP scattering
			    if (scattering_mechanisms[1]==1)
			    {
				mobility_po = mu_po(E_F,T,coefficients_cond,kindex_cond,g_pop,g,nu_el,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_po = "<<mobility_po<<" cm^2/(V-s)"<<endl;
				
				//calc_mobility_po[cc][1] = mobility_po;
				
				if(pop_number>1) 
				{
					
					for(int aa=0;aa<=pop_number-1;aa++)	
					{	
						for (int counter = 0;counter<points;counter++)
							dummy[counter] = g_pop_parts[aa][counter];

						calc_mobility_po[cc][aa+2] = mu_po(E_F,T,coefficients_cond,kindex_cond,dummy,g,nu_el,points,a11,energy_n,v_n,Ds_n);;
					}
				}	
				
			    }
			    else
				mobility_po=1e10;



			    // npop scattering
			    if (scattering_mechanisms[2]==1)
			    {
				mobility_npop = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_npop_total,points,a11, energy_n, v_n, Ds_n);
				cout<<"mobility_npop = "<<mobility_npop<<" cm^2/(V-s)"<<endl;
				
				//calc_mobility_npop[cc][1] = mobility_npop;
				
				if(npop_number>1) 
				{
					
					for(int aa=0;aa<=npop_number-1;aa++)	
					{	
						for (int counter = 0;counter<points;counter++)
							dummy[counter] = nu_npop[aa][counter] ;

						calc_mobility_npop[cc][aa+2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
					}
				}	

			    }
			    else
				mobility_npop = 1e10;




			    // Acoustic deformation scattering
			    if (scattering_mechanisms[3]==1)
			    {
				mobility_de = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_deformation,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_de = "<<mobility_de<<" cm^2/(V-s)"<<endl;
				
				//calc_mobility_de[cc][1] = mobility_de;
				
				if(de_number>1) 
				{
							
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[0][counter] ;

					calc_mobility_de[cc][2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
					cout<<"mobility_de for LA = "<<calc_mobility_de[cc][2]<<" cm^2/(V-s)"<<endl;
					
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[1][counter] ;

					calc_mobility_de[cc][3]  = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
					cout<<"mobility_de for TA = "<<calc_mobility_de[cc][3]<<" cm^2/(V-s)"<<endl;

					if(de_number==3)
					{
						for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[2][counter] ;

						calc_mobility_de[cc][4]  = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
						cout<<"mobility_de for ZA = "<<calc_mobility_de[cc][4]<<" cm^2/(V-s)"<<endl;
					}
				}	
			    }
			    else
				mobility_de=1e10;
			    
			    	

			    // piezoelectric scattering
			    if (scattering_mechanisms[4]==1)
			    {
				mobility_pe = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_pe = "<<mobility_pe<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_pe=1e10;

			   //calc_mobility_pe[cc][1] = mobility_e;
			   

			   // Transverse optical POP scattering 
			    if (scattering_mechanisms[5]==1)
			    {
				//C_long = C_long/10; // To convert dyn/cm2 back to Pa
				double zz = C_long/(10*rho);
				double vs = pow(zz,0.5); // speed of sound (longitudinal here [m/s])
				// C_long is divided with 10 in pow to convert dyn/cm2 back to Pa
				cout<<"\nCalculated speed of sound (longitudinal) is  "<<vs<<" m/s\n" <<endl;
				double a;
				if (E_F>0)
				    a =E_F;
				else
				    a=0;
				mobility_to = ( pow(2,0.5)*pi*e*pow(h_bar,3)*rho*(exp(h_bar*omega_TO/(k_B*T))-1)*vs*vs ) / ( pow(m,2.5)*pow(m_e,2.5)*omega_TO*E_deformation[0]*E_deformation[0]*
				                    pow((a+h_bar*omega_TO),0.5) ) * 1e4*pow(e,0.5);
				                    //Last coeff. is to convert to cm2/V.s
				//nu_to(:) = e/(m*m_e*mobility_to)*1e4;
				// Last part is to convert units to [1/s]

				cout<<"mobility_to = "<<mobility_to<<" cm^2/(V-s)"<<endl;
			    }
			    



			    // dislocation scattering
			    if (scattering_mechanisms[6]==1)
			    {
				mobility_dis = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_dislocation,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_dis = "<<mobility_dis<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_dis=1e10;

			    


			    // alloy scattering
			    if (scattering_mechanisms[7]==1)
			    {
				mobility_alloy = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_alloy,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_alloy = "<<mobility_alloy<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_alloy=1e10;




			    // inter-valley scattering
			    if (scattering_mechanisms[8]==1)
			    {
				mobility_iv = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_iv_total,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_iv  = "<<mobility_iv<<" cm^2/(V-s)"<<endl;

				calc_mobility_iv[cc][1] = mobility_iv;
				
				if(iv_number>1) 
				{
					
					for(int aa=0;aa<=iv_number-1;aa++)	
					{	
						for (int counter = 0;counter<points;counter++)
							dummy[counter] = nu_iv[aa][counter] ;

						calc_mobility_iv[cc][aa+2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
					}
				}	


			    }
			    else
				mobility_iv=1e10;



			    // neutral impurity scattering
			    if (scattering_mechanisms[9]==1)
			    {
			     mobility_neutral = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_neutralimpurity,points,a11,energy_n,v_n,Ds_n);

				cout<<"mobility_neutral = "<<mobility_neutral<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_neutral=1e10;

			    double nu_to[points]={0};

			    			    
			    			    
			    // Interface roughness scattering
			    if (scattering_mechanisms[11]==1)
			    {
				mobility_ir = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_irs,points,a11, energy_n, v_n, Ds_n);
				cout<<"mobility_ir = "<<mobility_ir<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_ir = 1e10;
				
			
			    // skew scattering
			    if (scattering_mechanisms[12]==1)
			    {
				mobility_skew = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_skew_rate,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_skew = "<<mobility_skew<<" cm^2/(V-s)"<<endl;
			    }
			    else
			    {
			
				mobility_skew=1e10;
			    }	

			    mobility_all[0] = mobility_ii;
			    mobility_all[1] = mobility_po;
			    mobility_all[2] = mobility_npop;
			    mobility_all[3] = mobility_de;
			    mobility_all[4] = mobility_pe;
			    mobility_all[5] = mobility_to;
			    mobility_all[6] = mobility_dis;
			    mobility_all[7] = mobility_alloy;
			    mobility_all[8] = mobility_iv;
			    mobility_all[9] = mobility_neutral;
			    mobility_all[11] = mobility_ir;
			    mobility_all[12] = mobility_skew;


			    //scattering_mechanisms
			    double sum =0 ;
			    for (int i=0;i<15;i++)
			    {
				if (scattering_mechanisms[i]!=0 && mobility_all[i]!=0)
				sum = sum + 1/ (mobility_all[i]*scattering_mechanisms[i]);
			    }

			    //cout<<"sum = "<<sum<<endl;

			    mobility_avg = 1/sum;

			    mobility = mu_overall(E_F,T,coefficients_cond,kindex_cond,g,nu_el,points,a11,energy_n,v_n,Ds_n);
			    // unit cm^2/(V-s)
			    
			    mobility_rta = mu_overall(E_F,T,coefficients_cond,kindex_cond,g_rta,nu_el,points,a11,energy_n,v_n,Ds_n);
			     // unit cm^2/(V-s)
			     	

			    if (omega_TO > 0.0)
				mobility = 1 / (1/mobility + 1/mobility_to);
				// unit cm^2/(V-s)

			    if (mobility < 0)
				    mobility = mobility_avg;
				// unit cm^2/(V-s)
					
				
		//----------------------------------calculate mobility completed -----------------------------
		
			    sigma = mobility *  n0 * e;
			    // unit S/cm 	
				
			    sigma_rta = mobility_rta * n0 * e;
			    // unit S/cm
			    	
			    if (n0 == 0)
			    {
			    	sigma = mobility *  abs(n_e) * e;
			    	sigma_rta = mobility_rta *  abs(n_e) * e;
			    	// unit S/cm
			    }
			    thermopower = -k_B*(df0dz_integral- E_F /(k_B*T))*1e6 + (J(T,g_th,points,v_n)/sigma)/dTdz*1e6;
			    // Equation No. 52 of rode book   micro V/K

		///---------------------------------------------------------------------------------------------------------------------
				

		// -------------------------------- calculate hall mobility -------------------------------------------------------------
			if(Bfield!=0)
			{

			    cout.precision(6);                        //set precision
			    cout.setf(ios::scientific);

			    cout<<endl;
			    // ionized impurity scattering
			    if (scattering_mechanisms[0]==1)
			    {
				mobility_hall_ii =
				 mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_ii = "<<mobility_hall_ii<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_ii = 1e10;


			    // POP scattering
			    if (scattering_mechanisms[1]==1)
			    {
				mobility_hall_po =
				mu_poH(E_F,T,coefficients_cond,kindex_cond,gH_pop,hH_pop,nu_el,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_po = "<<mobility_hall_po<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_po=1e10;

			    // npop scattering
			    if (scattering_mechanisms[2]==1)
			    {
				mobility_hall_npop = 
				mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_npop_total,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_npop = "<<mobility_hall_npop<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_npop = 1e10;


			    // Acoustic deformation scattering
			    if (scattering_mechanisms[3]==1)
			    {
				mobility_hall_de = 
				mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_deformation,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_de = "<<mobility_hall_de<<" cm^2/(V-s)"<<endl;
								
			    }
			    else
				mobility_hall_de=1e10;

			    // piezoelectric scattering
			    if (scattering_mechanisms[4]==1)
			    {
				mobility_hall_pe =
				mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_pe = "<<mobility_hall_pe<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_pe=1e10;

			   // Transverse optical POP scattering 
			    if (scattering_mechanisms[5]==1)
			    {
				//C_long = C_long/10; // To convert dyn/cm2 back to Pa
				double zz = C_long/(10*rho);
				double vs = pow(zz,0.5); // speed of sound (longitudinal here [m/s])
				// C_long is divided with 10 in pow to convert dyn/cm2 back to Pa
				cout<<"\nCalculated speed of sound (longitudinal) is  "<<vs<<" m/s\n" <<endl;
				double a;
				if (E_F>0)
				    a =E_F;
				else
				    a=0;
				mobility_hall_to = ( pow(2,0.5)*pi*e*pow(h_bar,3)*rho*(exp(h_bar*omega_TO/(k_B*T))-1)*vs*vs ) / ( pow(m,2.5)*pow(m_e,2.5)*omega_TO*E_deformation[0]*E_deformation[0]*
				                    pow((a+h_bar*omega_TO),0.5) ) * 1e4*pow(e,0.5);
				                    //Last coeff. is to convert to cm2/V.s
				//nu_to(:) = e/(m*m_e*mobility_to)*1e4;
				// Last part is to convert units to [1/s]

				cout<<"mobility_hall_to = "<<mobility_hall_to<<" cm^2/(V-s)"<<endl;
			    }

			    // dislocation scattering
			    if (scattering_mechanisms[6]==1)
			    {
				mobility_hall_dis =
				 mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_dislocation,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_dis = "<<mobility_hall_dis<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_dis=1e10;



			    // alloy scattering
			    if (scattering_mechanisms[7]==1)
			    {
				mobility_hall_alloy = 
				mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_alloy,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_alloy = "<<mobility_hall_alloy<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_alloy=1e10;


			    // inter-valley scattering
			    if (scattering_mechanisms[8]==1)
			    {
				mobility_hall_iv = 
				mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_iv_total,points,a11,energy_n,v_n,Ds_n);
				cout<<"mobility_hall_iv  = "<<mobility_hall_iv<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_iv=1e10;

			    // neutral impurity scattering
			    if (scattering_mechanisms[9]==1)
			    {
				mobility_hall_neutral = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,
				    nu_neutralimpurity,points,a11,energy_n,v_n,Ds_n);

				cout<<"mobility_hall_neutral = "<<mobility_hall_neutral<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_neutral=1e10;

			    // interface roughness scattering
			    if (scattering_mechanisms[11]==1)
			    {
				mobility_hall_ir = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,
				    nu_irs,points,a11,energy_n,v_n,Ds_n);

				cout<<"mobility_hall_ir = "<<mobility_hall_ir<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_ir=1e10;


			    // skew scattering
			    if (scattering_mechanisms[12]==1)
			    {
				mobility_hall_skew = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,
				    nu_skew_rate,points,a11,energy_n,v_n,Ds_n);

				cout<<"mobility_hall_skew = "<<mobility_hall_skew<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_skew=1e10;

			    double nu_hall_to[points]={0};
			    
			    
			    mobility_hall_all[0] = mobility_hall_ii;
			    mobility_hall_all[1] = mobility_hall_po;
			    mobility_hall_all[2] = mobility_hall_npop;
			    mobility_hall_all[3] = mobility_hall_de;
			    mobility_hall_all[4] = mobility_hall_pe;
			    mobility_hall_all[5] = mobility_hall_to;
			    mobility_hall_all[6] = mobility_hall_dis;
			    mobility_hall_all[7] = mobility_hall_alloy;
			    mobility_hall_all[8] = mobility_hall_iv;
			    mobility_hall_all[9] = mobility_hall_neutral;
			    mobility_hall_all[11] = mobility_hall_ir;
			    mobility_hall_all[12] = mobility_hall_skew;
			    
			    //scattering_mechanisms
			    
			    double sum_hall =0 ;
			    for (int i=0;i<15;i++)
			    {
				if (scattering_mechanisms[i]!=0 && mobility_hall_all[i]!=0)
				sum_hall = sum_hall + 1/ (mobility_hall_all[i]*scattering_mechanisms[i]);
			    }

			    //cout<<"sum_hall = "<<sum_hall<<endl;

			    mobility_hall_avg = 1/sum_hall;

			    mobility_hall = mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH,hH,nu_el,points,a11,energy_n,v_n,Ds_n);

			    mobility_hall_rta = 
			    mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH_rta, hH_rta, nu_el,points,a11,energy_n,v_n,Ds_n);

			    if (omega_TO > 0.0)
				mobility_hall = 1 / (1/mobility_hall + 1/mobility_hall_to);


			    if (mobility_hall < 0)
				    mobility_hall = mobility_hall_avg;

			     hall_factor1 = mobility_hall/mobility;
			     hall_factor_rta1 = mobility_hall_rta/mobility_rta;
			     
			    sigma_hall = mobility_hall *  n0 * e;

			    sigma_hall_rta = mobility_hall_rta * n0 * e;

			    if (n0 == 0)
			    {
			    	sigma_hall = mobility_hall *  abs(n_e) * e;
			    	sigma_hall_rta = mobility_hall_rta *  abs(n_e) * e;
			    }

			}	     

		//----------------------------------calculated hall mobility -------------------------------------------------------------
			
			calc_mobility[cc][1] = mobility;
			calc_mobility_rta[cc][1] = mobility_rta;
			calc_thermopower[cc][1] = thermopower;  // unit micro V/K
			calc_sigma[cc][1] = sigma;
			calc_sigma_rta[cc][1] = sigma_rta;


			calc_mobility_ii[cc][1] = mobility_ii;
			calc_mobility_po[cc][1] = mobility_po;
			calc_mobility_npop[cc][1] = mobility_npop;
			calc_mobility_de[cc][1] = mobility_de;
			calc_mobility_pe[cc][1] = mobility_pe;
			calc_mobility_to[cc][1] = mobility_to;
			calc_mobility_dis[cc][1] = mobility_dis;
			calc_mobility_alloy[cc][1] = mobility_alloy;
              	        calc_mobility_iv[cc][1] = mobility_iv;
			calc_mobility_neutral[cc][1] = mobility_neutral;
			calc_mobility_ir[cc][1] = mobility_ir;
			calc_mobility_skew[cc][1] = mobility_skew;

			if(Bfield!=0)
			{

			    calc_mobility_hall[cc][1] = mobility_hall;
			    calc_mobility_hall_rta[cc][1] = mobility_hall_rta;
			    calc_sigma_hall[cc][1] = sigma_hall;
			    calc_sigma_hall_rta[cc][1] = sigma_hall_rta;


			    calc_mobility_hall_ii[cc][1] = mobility_hall_ii;
			    calc_mobility_hall_po[cc][1] = mobility_hall_po;
			    calc_mobility_hall_npop[cc][1] = mobility_hall_npop;
			    calc_mobility_hall_de[cc][1] = mobility_hall_de;
			    calc_mobility_hall_pe[cc][1] = mobility_hall_pe;
			    calc_mobility_hall_to[cc][1] = mobility_hall_to;
			    calc_mobility_hall_dis[cc][1] = mobility_hall_dis;
			    calc_mobility_hall_alloy[cc][1] = mobility_hall_alloy;
			    calc_mobility_hall_iv[cc][1] = mobility_hall_iv;
			    calc_mobility_hall_neutral[cc][1] = mobility_hall_neutral;
			    calc_mobility_hall_ir[cc][1] = mobility_hall_ir;
			    calc_mobility_hall_skew[cc][1] = mobility_hall_skew;
			    
			    hall_factor[cc][1] = hall_factor1; 
			    hall_factor_rta[cc][1] = hall_factor_rta1; 
			}
			
			/*
			cout<<"cc = "<<cc<<endl;
			cout<<"calc_mobility[cc][0] = "<<calc_mobility[cc][0]<<endl;
			cout<<"calc_mobility[cc][1] =  "<<calc_mobility[cc][1]<<endl;
			cout<<"calc_mobility_rta[cc][1] =  "<<calc_mobility_rta[cc][1]<<endl;
			cout<<"calc_mobility_ii[cc][1] =  "<<calc_mobility_ii[cc][1]<<endl;
			cout<<"calc_mobility_po[cc][1] =  "<<calc_mobility_po[cc][1]<<endl;
			cout<<"calc_mobility_npop[cc][1] =  "<<calc_mobility_npop[cc][1]<<endl;
			cout<<"calc_mobility_de[cc][1] =  "<<calc_mobility_de[cc][1]<<endl;
			cout<<"calc_mobility_pe[cc][1] =  "<<calc_mobility_pe[cc][1]<<endl;
			getchar();
			//*/

			
		} // if condiction for type n completed 
		else   // for p type
		{
			// ii impurity scattering
			if (scattering_mechanisms[0]==1)
			{

				for (int counter = 0;counter<points;counter++)
					dummy[counter] = nu_ionizedimpurity_p[counter][0][0];		    		

				mobility_ii = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_ii = "<<mobility_ii<<" cm^2/(V-s)"<<endl;
				
				calc_mobility_ii[cc][1] = mobility_ii;
			}
			else
				mobility_ii = 1e10;



			// POP scattering
			if (scattering_mechanisms[1]==1)
			{
				mobility_po = mu_po(E_F,T,coefficients_val,kindex_val,g_pop,g,nu_el,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_po = "<<mobility_po<<" cm^2/(V-s)"<<endl;

				calc_mobility_po[cc][1] = mobility_po;
				
				if(pop_number>1) 
				{
					
					for(int aa=0;aa<=pop_number-1;aa++)	
					{	
						for (int counter = 0;counter<points;counter++)
							dummy[counter] = g_pop_parts[aa][counter];

						calc_mobility_po[cc][aa+2] = mu_po(E_F,T,coefficients_val,kindex_val,dummy,g,nu_el,points,b11,energy_p,v_p,Ds_p);
					}
				}	



			}
			else
				mobility_po=1e10;

			// npop scattering
			if (scattering_mechanisms[2]==1)
			{
				for (int counter = 0;counter<points;counter++)
					dummy[counter] = nu_npop_p[counter][0][0];		    		

				mobility_npop = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_npop = "<<mobility_npop<<" cm^2/(V-s)"<<endl;


				calc_mobility_npop[cc][1] = mobility_npop;
				
				if(npop_number>1) 
				{
					
					for(int aa=0;aa<=npop_number-1;aa++)	
					{	
						for (int counter = 0;counter<points;counter++)
							dummy[counter] = nu_npop[aa][counter] ;

						calc_mobility_npop[cc][aa+2] = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
					}
				}	
			}
			else
				mobility_npop = 1e10;


			// Acoustic deformation scattering
			if (scattering_mechanisms[3]==1)
			{	 
				for (int counter = 0;counter<points;counter++)
					dummy[counter] = nu_deformation_p[counter][0][0];
							    		
				mobility_de = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_de = "<<mobility_de<<" cm^2/(V-s)"<<endl;
				
				calc_mobility_de[cc][1] = mobility_de;

				if(de_number>1) 
				{
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[counter][0] ;

					calc_mobility_de[cc][2] = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
					cout<<"mobility_de for LA = "<<calc_mobility_de[cc][2]<<" cm^2/(V-s)"<<endl;
					
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[counter][1] ;

					calc_mobility_de[cc][3]  = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
					cout<<"mobility_de for TA = "<<calc_mobility_de[cc][3]<<" cm^2/(V-s)"<<endl;

					if(de_number==3)
					{
						for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[counter][2] ;

						calc_mobility_de[cc][4]  = mu_elastic(E_F,T,coefficients_val,kindex_val,dummy,points,b11,energy_p,v_p,Ds_p);
						cout<<"mobility_de for ZA = "<<calc_mobility_de[cc][4]<<" cm^2/(V-s)"<<endl;
					}
				}	
			}
			else
				mobility_de=1e10;
		
		
			mobility_all[0] = mobility_ii;
			mobility_all[1] = mobility_po;
			mobility_all[2] = mobility_npop;
			mobility_all[3] = mobility_de;
			    
			//scattering_mechanisms
			double sum =0 ;
			for (int i=0;i<=14;i++)
			{
				if (scattering_mechanisms[i]!=0 && mobility_all[i]!=0)
					sum = sum + 1/ (mobility_all[i]*scattering_mechanisms[i]);
			}

			//cout<<"sum = "<<sum<<endl;

			mobility_avg = 1/sum;

			mobility = mu_overall(E_F,T,coefficients_val,kindex_val,g,nu_el,points,b11,energy_p,v_p,Ds_p);
			// unit cm^2/(V-s)

			mobility_rta = mu_overall(E_F,T,coefficients_val,kindex_val,g_rta,nu_el,points,b11,energy_p,v_p,Ds_p);
			// unit cm^2/(V-s)

			if (omega_TO > 0.0)
				mobility = 1 / (1/mobility + 1/mobility_to);
			// unit cm^2/(V-s)

			if (mobility < 0)
				mobility = mobility_avg;
			// unit cm^2/(V-s)

				
		//----------------------------------calculated mobility ----------------------------------------------------------

			sigma = mobility *  n0 * e;
			// unit S/cm 	

			sigma_rta = mobility_rta * n0 * e;
			// unit S/cm

			if (n0 == 0)
			{
			sigma = mobility *  abs(n_h) * e;
			sigma_rta = mobility_rta *  abs(n_h) * e;
			// unit S/cm
			}

			thermopower = -k_B*(df0dz_integral- E_F /(k_B*T))*1e6 + (J(T,g_th,points,v_p)/sigma)/dTdz*1e6;
			// Equation No. 52 of rode book   micro V/K

		// -------------------------------- calculate hall mobility -------------------------------------------------------------
			if(Bfield!=0)   // for ptype
			{
				//cout<<"reached isnide"<<endl;
				//getchar();
			    cout.precision(6);                        //set precision
			    cout.setf(ios::scientific);

			    cout<<endl;
			    
			    // ionized impurity scattering
			    if (scattering_mechanisms[0]==1)
			    {
				mobility_hall_ii =
				 mu_elasticH(E_F,T,coefficients_val,kindex_val,
				 nu_ionizedimpurity,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_hall_ii = "<<mobility_hall_ii<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_ii = 1e10;


			    // POP scattering
			    if (scattering_mechanisms[1]==1)
			    {
				mobility_hall_po =
				mu_poH(E_F,T,coefficients_val,kindex_val,gH_pop,hH_pop,nu_el,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_hall_po = "<<mobility_hall_po<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_po=1e10;

			    // npop scattering
			    if (scattering_mechanisms[2]==1)
			    {
				mobility_hall_npop = 
				mu_elasticH(E_F,T,coefficients_val,kindex_val,nu_npop_total,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_hall_npop = "<<mobility_hall_npop<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_npop = 1e10;


			    // Acoustic deformation scattering
			    if (scattering_mechanisms[3]==1)
			    {
				mobility_hall_de = 
				mu_elasticH(E_F,T,coefficients_val,kindex_val,nu_deformation,points,b11,energy_p,v_p,Ds_p);
				cout<<"mobility_hall_de = "<<mobility_hall_de<<" cm^2/(V-s)"<<endl;
			    }
			    else
				mobility_hall_de=1e10;
						    
			    mobility_hall_all[0] = mobility_hall_ii;
			    mobility_hall_all[1] = mobility_hall_po;
			    mobility_hall_all[2] = mobility_hall_npop;
			    mobility_hall_all[3] = mobility_hall_de;
			    
			    //scattering_mechanisms
			    
			    double sum_hall =0 ;
			    for (int i=0;i<15;i++)
			    {
				if (scattering_mechanisms[i]!=0 && mobility_hall_all[i]!=0)
				sum_hall = sum_hall + 1/ (mobility_hall_all[i]*scattering_mechanisms[i]);
			    }

			    //cout<<"sum_hall = "<<sum_hall<<endl;

			    mobility_hall_avg = 1/sum_hall;

			    mobility_hall = mu_overallH(E_F,T,coefficients_val,kindex_val,gH,hH,nu_el,points,b11,energy_p,v_p,Ds_p);

			    mobility_hall_rta = 
			    mu_overallH(E_F,T,coefficients_val,kindex_val,gH_rta, hH_rta, nu_el,points,b11,energy_p,v_p,Ds_p);

			    if (omega_TO > 0.0)
				mobility_hall = 1 / (1/mobility_hall + 1/mobility_hall_to);


			    if (mobility_hall < 0)
				    mobility_hall = mobility_hall_avg;

			     hall_factor1 = mobility_hall/mobility;
			     hall_factor_rta1 = mobility_hall_rta/mobility_rta;
			     
			    sigma_hall = mobility_hall *  n0 * e;

			    sigma_hall_rta = mobility_hall_rta * n0 * e;

			    if (n0 == 0)
			    {
			    	sigma_hall = mobility_hall *  abs(n_e) * e;
			    	sigma_hall_rta = mobility_hall_rta *  abs(n_e) * e;
			    }

			}  // bfield condition finished	     

		//----------------------------------calculated hall mobility ----------------------------------------

			calc_mobility[cc][1] = mobility;
			calc_mobility_rta[cc][1] = mobility_rta;
			calc_sigma[cc][1] = sigma;
			calc_sigma_rta[cc][1] = sigma_rta;


			if(Bfield!=0)
			{

			    calc_mobility_hall[cc][1] = mobility_hall;
			    calc_mobility_hall_rta[cc][1] = mobility_hall_rta;
			    calc_sigma_hall[cc][1] = sigma_hall;
			    calc_sigma_hall_rta[cc][1] = sigma_hall_rta;


			    calc_mobility_hall_ii[cc][1] = mobility_hall_ii;
			    calc_mobility_hall_po[cc][1] = mobility_hall_po;
			    calc_mobility_hall_npop[cc][1] = mobility_hall_npop;
			    calc_mobility_hall_de[cc][1] = mobility_hall_de;
		    
			    hall_factor[cc][1] = hall_factor1; 
			    hall_factor_rta[cc][1] = hall_factor_rta1; 
			}

		///---------------------------------------------------------------------------------------------------------------------
		} // else condition for type p completed	

		cout.precision(6);                        //set precision
		cout.setf(ios::scientific);
		cout<<endl<<"Drift mobility results"<<endl;
		cout<<"Temperature = "<<T<<" K"<<endl;
		cout<<"doping = "<<n_array[d_loop]<<" cm^-3"<<endl;
		cout<<"mobility = "<<mobility<<" cm^2/(V-s)"<<endl;
		cout<<"mobility_rta = "<<mobility_rta<<" cm^2/(V-s)"<<endl;
		cout<<"mobility_avg = "<<mobility_avg<<" cm^2/(V-s)"<<endl;
		cout<<"sigma = "<<sigma<<" S/cm "<<endl;
		cout<<"sigma_rta = "<<sigma_rta<<" S/cm "<<endl;
		cout<<"thermopower = "<<thermopower<<" micro V/K"<<endl<<endl;

		if(Bfield!=0)
		{

		    cout<<"Hall mobility results"<<endl;
		    cout<<"mobility_hall = "<<mobility_hall<<" cm^2/(V-s)"<<endl;
		    cout<<"mobility_hall_rta = "<<mobility_hall_rta<<" cm^2/(V-s)"<<endl;
		    cout<<"mobility_hall_avg = "<<mobility_hall_avg<<" cm^2/(V-s)"<<endl;
		    cout<<"Hall factor = "<<hall_factor1<<endl;
		    cout<<"Hall factor RTA = "<<hall_factor_rta1<<endl;
		    cout<<"sigma_hall = "<<sigma_hall<<" S/cm "<<endl;
		    cout<<"sigma_hall_rta = "<<sigma_hall_rta<<" S/cm "<<endl<<endl<<endl;

		}
	} // end of if part for goemetry == 1
	else if(geometry==2) // for 2D
	{
		    //cout<<"calculating result for 2D "<<endl;
		    	
		    // ionized impurity scattering
		    if (scattering_mechanisms[0]==1)
		    {
			mobility_ii = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_ii = "<<mobility_ii<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_ii = 1e10;


		    // POP scattering
		    if (scattering_mechanisms[1]==1)
		    {
			mobility_po = mu_po(E_F,T,coefficients_cond,kindex_cond,g_pop,g,nu_el,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_po = "<<mobility_po<<" cm^2/(V-s)"<<endl;

			calc_mobility_po[cc][1] = mobility_po;
			
			if(pop_number>1) 
			{
				
				for(int aa=0;aa<=pop_number-1;aa++)	
				{	
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = g_pop_parts[aa][counter];

					calc_mobility_po[cc][aa+2] = mu_po(E_F,T,coefficients_cond,kindex_cond,dummy,g,nu_el,points,a11,energy_n,v_n,Ds_n);;
				}
			}	

		    }
		    else
			mobility_po=1e10;

		    // npop scattering
		    if (scattering_mechanisms[2]==1)
		    {
			mobility_npop = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_npop_total,points,a11, energy_n, v_n, Ds_n);
			cout<<"mobility_npop = "<<mobility_npop<<" cm^2/(V-s)"<<endl;
			

			calc_mobility_npop[cc][1] = mobility_npop;
				
			if(npop_number>1) 
			{
				
				for(int aa=0;aa<=npop_number-1;aa++)	
				{	
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_npop[aa][counter] ;

					calc_mobility_npop[cc][aa+2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
				}
			}				
		    }
		    else
			mobility_npop = 1e10;


		    // Acoustic deformation scattering
		    if (scattering_mechanisms[3]==1)
		    {
			mobility_de = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_deformation,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_de = "<<mobility_de<<" cm^2/(V-s)"<<endl;
			
			calc_mobility_de[cc][1] = mobility_de;
			
			if(de_number>1)
			{
						
				for (int counter = 0;counter<points;counter++)
					dummy[counter] = nu_def[0][counter] ;

				calc_mobility_de[cc][2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
				//cout<<"mobility_de for LA = "<<calc_mobility_de[cc][2]<<" cm^2/(V-s)"<<endl;
				
				for (int counter = 0;counter<points;counter++)
					dummy[counter] = nu_def[1][counter] ;

				calc_mobility_de[cc][3] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
				//cout<<"mobility_de for TA = "<<calc_mobility_de[cc][3]<<" cm^2/(V-s)"<<endl;

				if(de_number==3)
				{
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_def[2][counter];

					calc_mobility_de[cc][4] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
					//cout<<"mobility_de for ZA = "<<calc_mobility_de[cc][4]<<" cm^2/(V-s)"<<endl;
				}
		         }
		    }
		    else
			mobility_de=1e10;

		    // piezoelectric scattering
		    if (scattering_mechanisms[4]==1)
		    {
			mobility_pe = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_pe = "<<mobility_pe<<" cm^2/(V-s)"<<endl;

			calc_mobility_pe[cc][1] = mobility_pe;
			
			if(de_number>1) 
			{
				for(int aa=0;aa<=de_number-1;aa++)	
				{	
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = nu_pz[aa][counter] ;

					calc_mobility_pe[cc][aa+2] = mu_elastic(E_F,T,coefficients_cond,kindex_cond,dummy,points,a11,energy_n,v_n,Ds_n);
				}
			}	
		    }
		    else
			mobility_pe=1e10;
		
		
		
		    // SO POP scattering
		    if (scattering_mechanisms[10]==1)
		    {
			mobility_so_pop = mu_po(E_F,T,coefficients_cond,kindex_cond,g_so_pop,g,nu_el,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_so_pop = "<<mobility_so_pop<<" cm^2/(V-s)"<<endl;
			
			calc_mobility_so_pop[cc][1] = mobility_so_pop;
			
			if(so_pop_number>1) 
			{
				
				for(int aa=0;aa<=so_pop_number-1;aa++)	
				{	
					for (int counter = 0;counter<points;counter++)
						dummy[counter] = g_so_pop_parts[aa][counter];

					calc_mobility_so_pop[cc][aa+2] = mu_po(E_F,T,coefficients_cond,kindex_cond,dummy,g,nu_el,points,a11,energy_n,v_n,Ds_n);;
				}
			}	

		    }
		    else
			mobility_so_pop = 1e10;


		    // Interface roughness scattering
		    if (scattering_mechanisms[11]==1)
		    {
			mobility_ir = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_irs,points,a11, energy_n, v_n, Ds_n);
			cout<<"mobility_ir = "<<mobility_ir<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_ir = 1e10;


		    // skew scattering
		    if (scattering_mechanisms[12]==1)
		    {
			mobility_skew = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_skew_rate,points,a11, energy_n, v_n, Ds_n);
			cout<<"mobility_skew = "<<mobility_skew<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_skew = 1e10;


		    mobility_all[0] = mobility_ii;
		    mobility_all[1] = mobility_po;
		    mobility_all[2] = mobility_npop;
		    mobility_all[3] = mobility_de;
		    mobility_all[4] = mobility_pe;
		    mobility_all[10] = mobility_so_pop;
		    mobility_all[11] = mobility_ir;
		    mobility_all[12] = mobility_skew;

		    //scattering_mechanisms
		    double sum =0 ;
		    for (int i=0;i<15;i++)
		    {
			if (scattering_mechanisms[i]!=0 && mobility_all[i]!=0)
			sum = sum + 1/ (mobility_all[i]*scattering_mechanisms[i]);
		    }

		    //cout<<"sum = "<<sum<<endl;

		    mobility_avg = 1/sum;

		    mobility = mu_overall(E_F,T,coefficients_cond,kindex_cond,g,nu_el,points,a11,energy_n,v_n,Ds_n);
		    // unit cm^2/(V-s)
		    
		    cout<<"mobility  = "<<mobility<<endl;
		    
		    mobility_rta = mu_overall(E_F,T,coefficients_cond,kindex_cond,g_rta,nu_el,points,a11,energy_n,v_n,Ds_n);
		     // unit cm^2/(V-s)
		     	
		    if (mobility < 0)
			    mobility = mobility_avg;
			// unit cm^2/(V-s)

		//----------------------------------calculate mobility -------------------------------------------------------------

		    sigma = mobility *  n0 * e;
		    // unit S 	
		    //cout<<"sigma  = "<<sigma<<endl;
		
			
		    sigma_rta = mobility_rta * n0 * e;
		    // unit S, since n0 is per cm^2
		    	
		    if (n0 == 0)
		    {
		    	sigma = mobility *  abs(n_h) * e;
		    	sigma_rta = mobility_rta *  abs(n_h) * e;
		    	// unit S
		    }
		
		sigma = sigma/(thickness*100);
		sigma_rta = sigma_rta/(thickness*100);
		// unit S/cm
		//cout<<"sigma after thickness = "<<sigma<<endl;
		
		// other transport coefficients
		// according to paper Thermoelectric transport coefficients in monolayer Mos2 and Wes2 Role of substarte, interface 
		// phonons plasmon and dynamic screening
		
		/*
		sx = electric_current_density(E_F, T, coefficients_cond, kindex_cond, g, nu_el, points, a11, energy_n, v_n, 
		Ds_n)/(E*1e2);
		// unit S/m
		 

		//cout<<"sigma/thickness = "<<sigma/thickness<<"  S/m"<<endl; 
		//cout<<"sx = "<<sx<<"  S/m"<<endl;			
		//getchar();

		px = heat_current_density(E_F, T, coefficients_cond, kindex_cond, g, nu_el, points, a11, energy_n, v_n, Ds_n)
		/(E*1e2);
			
		// px unit is Joule/V-s-m
		// E unit is V/cm 	
		
		Bx = electric_current_density(E_F, T, coefficients_cond, kindex_cond, g_th, nu_el, points, a11, energy_n, v_n, Ds_n)
		*T*T/(dTdz*100);
		// unit A-K/m ampere Kelvin
		
		Kx = heat_current_density(E_F, T, coefficients_cond, kindex_cond, g_th, nu_el, points, a11, energy_n, v_n, Ds_n)
		*T*T/(dTdz*100);  
		// unit J-K/(m-s)  or N-K/s
		// heat current denisty unit is J/sm or W/m   dtdz unit is K/cm
		
		thermopower = -k_B*(df0dz_integral- E_F /(k_B*T))*1e6 + (J(T,g_th,points,v_p)*1e4)*1e6/(sx*(dTdz*100));
		// mks unit micro V/K
		
		//thermopower = Bx/(sx*T*T)*1e6;   // mks unit micro V/K
		peltier = (px)/sx;  // mks units volt or Joule/Coulomb     
		thermal_conductivity = (Kx - (px*Bx/sx))/(T*T);   // mks unit J/(K-s-m) or W/(K-m)
		
		cout<<"T    = "<<T<<endl;		
		cout<<"ef    = "<<E_F<<endl;
		cout<<"sigma   =   "<<sigma<<endl;
		cout<<"sx   =   "<<sx<<endl;
		cout<<"px   =   "<<px<<endl;
		cout<<"Bx   =   "<<Bx<<endl;
		cout<<"Kx   =   "<<Kx<<endl;
		cout<<"thermopower   =   "<<thermopower<<"   micro V/K "<<endl;
		getchar();
		
		FILE *fid1;
		fid1 = fopen("data.dat","a");	
		fprintf(fid1,"# Temperature(K)       E_F 	sx	px	Bx	Kx \n");
		fprintf(fid1," %e	%e        %e		%e  	%e	%e \n",T, E_F, sx, px, Bx, Kx);
		
		fclose(fid1);	
		*/
		

		calc_mobility[cc][1] = mobility;
		calc_mobility_rta[cc][1] = mobility_rta;
		
		//calc_thermopower[cc][1] = thermopower;  //converted from V/K to micro V/K 
		//calc_peltier[cc][1] = peltier;
		//calc_thermal_conductivity[cc][1] = thermal_conductivity;
		
		calc_sigma[cc][1] = sigma;
		calc_sigma_rta[cc][1] = sigma_rta;


		calc_mobility_ii[cc][1] = mobility_ii;
		calc_mobility_po[cc][1] = mobility_po;
		calc_mobility_so_pop[cc][1] = mobility_so_pop;
		calc_mobility_ir[cc][1] = mobility_ir;
		
		/*
		cout<<"cc = "<<cc<<endl;
		cout<<"calc_mobility[cc][0] = "<<calc_mobility[cc][0]<<endl;
		cout<<"calc_mobility[cc][1] =  "<<calc_mobility[cc][1]<<endl;
		cout<<"calc_mobility_rta[cc][1] =  "<<calc_mobility_rta[cc][1]<<endl;
		cout<<"calc_mobility_ii[cc][1] =  "<<calc_mobility_ii[cc][1]<<endl;
		cout<<"calc_mobility_po[cc][1] =  "<<calc_mobility_po[cc][1]<<endl;
		cout<<"calc_mobility_npop[cc][1] =  "<<calc_mobility_npop[cc][1]<<endl;
		cout<<"calc_mobility_de[cc][1] =  "<<calc_mobility_de[cc][1]<<endl;
		cout<<"calc_mobility_pe[cc][1] =  "<<calc_mobility_pe[cc][1]<<endl;
		cout<<"calc_mobility_so_pop[cc][1] =  "<<calc_mobility_so_pop[cc][1]<<endl;
		getchar();
		//*/

		if(Bfield!=0)
		{

		    cout.precision(6);                        //set precision
		    cout.setf(ios::scientific);

		    cout<<endl;

		    // ionized impurity scattering
		    if (scattering_mechanisms[0]==1)
		    {
			mobility_hall_ii =
			mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_ii = "<<mobility_hall_ii<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_ii = 1e10;


		    // POP scattering
		    if (scattering_mechanisms[1]==1)
		    {
			mobility_hall_po = 
			mu_poH(E_F,T,coefficients_cond,kindex_cond,gH_pop,hH_pop,nu_el,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_po = "<<mobility_hall_po<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_po=1e10;

		    // npop scattering
		    if (scattering_mechanisms[2]==1)
		    {
			mobility_hall_npop = 
			mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_npop_total,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_npop = "<<mobility_hall_npop<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_npop = 1e10;


		    // Acoustic deformation scattering
		    if (scattering_mechanisms[3]==1)
		    {
			mobility_hall_de = 
			mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_deformation,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_de = "<<mobility_hall_de<<" cm^2/(V-s)"<<endl;

		    }
		    else
			mobility_hall_de=1e10;

		    // piezoelectric scattering
		    if (scattering_mechanisms[4]==1)
		    {
			mobility_hall_pe = 
			mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_pe = "<<mobility_hall_pe<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_pe=1e10;

		    // SO POP scattering
		    if (scattering_mechanisms[10]==1)
		    {
			mobility_hall_so_po = 
			mu_poH(E_F,T,coefficients_cond,kindex_cond,gH_so_pop,hH_so_pop,nu_el,points,a11,energy_n,v_n,Ds_n);
			cout<<"mobility_hall_so_po = "<<mobility_hall_so_po<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_so_po=1e10;


		    // interface roughness scattering
		    if (scattering_mechanisms[11]==1)
		    {
			mobility_hall_ir = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_irs,points,a11,energy_n,v_n,Ds_n);

			cout<<"mobility_hall_ir = "<<mobility_hall_ir<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_ir=1e10;


		    // skew scattering
		    if (scattering_mechanisms[12]==1)
		    {
			mobility_hall_skew = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_skew_rate,points,a11,energy_n,v_n,Ds_n);

			cout<<"mobility_hall_skew = "<<mobility_hall_skew<<" cm^2/(V-s)"<<endl;
		    }
		    else
			mobility_hall_skew=1e10;


		    mobility_hall_all[0] = mobility_hall_ii;
		    mobility_hall_all[1] = mobility_hall_po;
		    mobility_hall_all[2] = mobility_hall_npop;
		    mobility_hall_all[3] = mobility_hall_de;
		    mobility_hall_all[4] = mobility_hall_pe;
		    mobility_hall_all[10] = mobility_hall_so_po;
		    mobility_hall_all[11] = mobility_hall_ir;
		    mobility_hall_all[12] = mobility_hall_skew;
	
		    double sum_hall =0 ;
		    for (int i=0;i<15;i++)
		    {
			if (scattering_mechanisms[i]!=0 && mobility_hall_all[i]!=0)
			sum_hall = sum_hall + 1/ (mobility_hall_all[i]*scattering_mechanisms[i]);
		    }

		    //cout<<"sum_hall = "<<sum_hall<<endl;

		    mobility_hall_avg = 1/sum_hall;    
		    
		    /*	
		    mobility_hall = mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH,hH,nu_el,points,a11,energy_n,v_n,Ds_n);
		    		    
		    if(mobility_hall < 0 )
		    {
		    	cout<<"mobility_hall is negative = "<<mobility_hall<<endl;
			getchar();		    	
		    }	
		    
		    	    	
		    mobility_hall_rta = 
		    mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH_rta, hH_rta, nu_el,points,a11,energy_n,v_n,Ds_n);

		    sigma_hall = mobility_hall *  n0 * e;
		    // unit S	
		    sigma_hall_rta = mobility_hall_rta * n0 * e;
		    // unit S
		    	
		    if (n0 == 0)
		    {
		    	sigma_hall = mobility_hall *  abs(n_e) * e;
		    	sigma_hall_rta = mobility_hall_rta *  abs(n_e) * e;
		    	// unit S
		    }
		    
		    sigma_hall = sigma_hall/(thickness*100);
		    sigma_hall_rta = sigma_hall_rta/(thickness*100);
		    // unit S/cm

		    */
		    	
		    sigma_xx = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,gH,nu_el,points,a11,energy_n,v_n,Ds_n,1);
		    // unit S/cm
		    
		    sigma_xy = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,hH,nu_el,points,a11,energy_n,v_n,Ds_n,2);
		    // unit S/cm
		    
		    sigma_xy = abs(sigma_xy);
		    //cout<<"sigma_xx = "<<sigma_xx<<" S/cm "<<endl;
		    //cout<<"sigma_xy = "<<sigma_xy<<" S/cm "<<endl;
		    //getchar();
		    
		    sigma_xx_rta = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,gH_rta,nu_el,points,a11,energy_n,v_n,Ds_n,1);
		    // unit S/cm
		    
		    sigma_xy_rta = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,hH_rta,nu_el,points,a11,energy_n,v_n,Ds_n,2);
		    // unit S/cm
		    
		    sigma_xy_rta = abs(sigma_xy_rta);
		    	
		    double dd,dd_rta;
		    dd = (sigma_xx*sigma_xx + sigma_xy*sigma_xy);
		    // unit S^2/cm^2
		    
		    dd_rta = (sigma_xx_rta*sigma_xx_rta + sigma_xy_rta*sigma_xy_rta);
		    // unit S^2/cm^2
		    
		    hall_coeff = -sigma_xy/(Bfield*0.0001*dd);
                    //Ohm-cm/G  
                    
		    hall_coeff_rta = -sigma_xy_rta/(Bfield*0.0001*dd_rta);
                    //Ohm-cm/G  
                    
                    // Bfield is multplied with 0.0001 to convert into cgs unit
		    // 1 gauss = 1e-4 Tesla; 1T = 1 kg/(A s^2) have unit of inverse of mobility
		    
		    //The two most widely used units for the Hall coefficients are SI units, m3/A-sec = m3/C, and 
		    //the hybrid unit Ohm-cm/G (which combines the practical quantities volt and amp with the 
		    // cgs quantities centimeter and Gauss).
		    
		    mobility_hall = sigma*abs(hall_coeff);			    
		    // unit cm^2/V-s
		    
		    mobility_hall_rta = sigma_rta*abs(hall_coeff_rta);			    
		    // unit cm^2/V-s

		    
		    long_restivity = abs(sigma_xx)/(sigma_xx*sigma_xx + sigma_xy*sigma_xy);
		    // unit cm-ohm
		    
		    long_restivity_rta = abs(sigma_xx_rta)/(sigma_xx_rta*sigma_xx_rta + sigma_xy_rta*sigma_xy_rta);
		    // unit cm-ohm

		    mg_resist = (sigma_xx*sigma/(dd) - 1);
		    // unit less

		    mg_resist_rta = (sigma_xx_rta*sigma_rta/(dd_rta) - 1);
		    // unit less

		    if (mobility_hall < 0)
			    mobility_hall = mobility_hall_avg;

		     hall_factor1 = mobility_hall/mobility;
		     hall_factor_rta1 = mobility_hall_rta/mobility_rta;

		    calc_mobility_hall_ii[cc][1] = mobility_hall_ii;
		    calc_mobility_hall_po[cc][1] = mobility_hall_po;
		    calc_mobility_hall_npop[cc][1] = mobility_hall_npop;
		    calc_mobility_hall_de[cc][1] = mobility_hall_de;

		    calc_mobility_hall_pe[cc][1] = mobility_hall_pe;
		    calc_mobility_hall_so_po[cc][1] = mobility_hall_so_po;
		    calc_mobility_hall_ir[cc][1] = mobility_hall_ir;

		    calc_mobility_hall[cc][1] = mobility_hall;
		    calc_mobility_hall_rta[cc][1] = mobility_hall_rta;

		    hall_factor[cc][1] = hall_factor1; 
		    calc_sigma_xx[cc][1] = sigma_xx;
		    calc_sigma_xy[cc][1] = sigma_xy;
		    calc_hall_coeff[cc][1] = hall_coeff;
		    calc_long_restivity[cc][1] = long_restivity;
		    calc_mg_resist[cc][1] = mg_resist;
		    
		    hall_factor_rta[cc][1] = hall_factor_rta1; 
		    calc_sigma_xx_rta[cc][1] = sigma_xx_rta;
		    calc_sigma_xy_rta[cc][1] = sigma_xy_rta;
		    calc_hall_coeff_rta[cc][1] = hall_coeff_rta;
		    calc_long_restivity_rta[cc][1] = long_restivity_rta;
		    calc_mg_resist_rta[cc][1] = mg_resist_rta;

		}  // if condition for bfield finished	     
		
		
		
		// optical conductivity part
	//-------------------------------------------------------------------------------------------------------------------------------------------------------
		// optical conductivity part
		if(freq_variation==1)
		{
		    cout.precision(6);                        //set precision
		    cout.setf(ios::scientific);

		    cout<<endl;
			
			int mm;
			if(variation==0) // temp variation			
				mm = T_loop;
			else		  // doping variation
				mm = d_loop;
		    
		    for(int i1=0;i1<len_freq;i1++) 	
		    {
		    	sigma_real = sigma_freq(E_F,T,coefficients_cond,kindex_cond,fqr,points,a11,energy_n,v_n,Ds_n,i1, mm);
		    	// unit S/cm
		    
		    	sigma_img = sigma_freq(E_F,T,coefficients_cond,kindex_cond,fqi,points,a11,energy_n,v_n,Ds_n,i1, mm);
		    	// unit S/cm
		    	
		    	sigma_real = abs(sigma_real);
		    	sigma_img = abs(sigma_img);
		    	
			//sigma_img = abs(sigma_img);
			//cout<<"sigma_real = "<<sigma_real<<" S/cm "<<endl;
			//cout<<"sigma_img = "<<sigma_img<<" S/cm "<<endl;
			//getchar();

		    	sigma_r[cc][i1] = sigma_real;
		    	sigma_i[cc][i1] = sigma_img;
		    	// unit S/cm				    		    
		    	
		    	if(n0!=0)
		    	{
		    		mobility_r[cc][i1] = sigma_real*(thickness*100)/(e*n0);
		    		mobility_i[cc][i1] = sigma_img*(thickness*100)/(e*n0);
			}
			else
			{
		    		mobility_r[cc][i1] = sigma_real*(thickness*100)/(e*abs(n_e));
		    		mobility_i[cc][i1] = sigma_img*(thickness*100)/(e*abs(n_e));			
			}
			// unit cm^2/V-s
			
			cout<<"Temperature = "<<T<<" K"<<endl;
			cout<<"freq = "<<freq[i1]<<" Hz "<<endl;
			cout<<"Sigma real  =  "<<sigma_r[cc][i1]<<"  S/cm "<<endl;
			cout<<"Sigma imaginary = "<<sigma_i[cc][i1]<<"  S/cm "<<endl;
			cout<<"Mobility real =  "<<mobility_r[cc][i1]<<"  cm^2/V-s "<<endl;
			cout<<"Mobility imaginary  =  "<<mobility_i[cc][i1]<<"  cm^2/V-s "<<endl;
			cout<<endl;
			//getchar();
			//sigma_r_rta = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,fqr_rta,points,a11,energy_n,v_n,Ds_n);
			// unit S/cm

			//sigma_i_rta = sigma_tensor(E_F,T,coefficients_cond,kindex_cond,fqi_rta,points,a11,energy_n,v_n,Ds_n);
			// unit S/cm

			//sigma_i_rta = abs(sigma_i_rta);
		    }
		
		} // end of freq variation condition or optical conductivity
	//-------------------------------------------------------------------------------------------------------------------------------------------------------	
		
		
		cout.setf(ios::scientific);
		cout<<"Temperature = "<<T<<" K"<<endl;
		cout<<"doping = "<<n_array[d_loop]<<" cm^-3"<<endl;
		cout<<"mobility = "<<mobility<<" cm^2/(V-s)"<<endl;
		cout<<"mobility_rta = "<<mobility_rta<<" cm^2/(V-s)"<<endl;
		cout<<"mobility_avg = "<<mobility_avg<<" cm^2/(V-s)"<<endl;
		cout<<"sigma = "<<sigma<<" S/cm "<<endl;
		cout<<"sigma_rta = "<<sigma_rta<<" S/cm "<<endl;
		
		if(Bfield!=0)
		{
		    cout<<endl;	
		    cout<<"Hall mobility results"<<endl;
		    cout<<"mobility_hall = "<<mobility_hall<<" cm^2/(V-s)"<<endl;
		    cout<<"mobility_hall_rta = "<<mobility_hall_rta<<" cm^2/(V-s)"<<endl;
		    cout<<"mobility_hall_avg = "<<mobility_hall_avg<<" cm^2/(V-s)"<<endl;
		    //cout<<"sigma_hall = "<<sigma_hall<<" S/cm "<<endl;
		    //cout<<"sigma_hall_rta = "<<sigma_hall_rta<<" S/cm "<<endl<<endl<<endl;
		    cout<<"Hall factor = "<<hall_factor1<<endl;
		    cout<<"Hall factor RTA = "<<hall_factor_rta1<<endl;
		    cout<<"sigma_xx = "<<sigma_xx<<" S/cm "<<endl;
		    cout<<"sigma_xx_rta = "<<sigma_xx_rta<<" S/cm "<<endl;
		    cout<<"sigma_xy = "<<sigma_xy<<" S/cm "<<endl;
		    cout<<"sigma_xy_rta = "<<sigma_xy_rta<<" S/cm "<<endl;
		    cout<<"hall_coeff = "<<hall_coeff<<"  Ohm-cm/G "<<endl;
		    cout<<"hall_coeff_rta = "<<hall_coeff_rta<<"  Ohm-cm/G "<<endl;
		    cout<<"long_restivity = "<<long_restivity<<"  cm-ohm "<<endl;
		    cout<<"long_restivity_rta = "<<long_restivity_rta<<"  cm-ohm "<<endl;
		    cout<<"Magneto resistance coefficient = "<<mg_resist<<endl;
		    cout<<"Magneto resistance coefficient RTA = "<<mg_resist_rta<<endl;


		    //getchar();
                  		    	
		}

		//cout<<"sx = "<<sx<<" S/m "<<endl;
		//cout<<"thermopower = "<<thermopower<<" micro V/K"<<endl;
		//cout<<"peltier = "<<peltier<<" Volt "<<endl;
		//cout<<"thermal_conductivity = "<<thermal_conductivity<<" W/K-m"<<endl<<endl;
							
	}  // end of else part for goemetry == 2, 
	
	//cout<<"Reached here calculate_mobility outside"<<endl;
	
}
