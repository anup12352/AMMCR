#include"main.h"

void components_BTE(double T, int T_loop, double efef, int ii)
{
	
//---------------------------- components for BTE ------------------------------------------------------------------
            beta_constant = beta(T, T_loop);
            // unit 1/nm

            cout<< "Inverse screening length, beta = "<<beta_constant<<" (1/nm)"<<endl;
            double integral_numerator_n = 0;
            double integral_denominator_n = 0;

            int factr = 10;
            for(int counter = 0;counter<=points-2;counter++)
            {
                double dk = (k_grid[counter+1]-k_grid[counter])/factr;
                for (int ss = 0;ss<=factr-1;ss++)
                {
                    integral_numerator_n = integral_numerator_n + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
                                        *f0(energy_n[counter],efef,T)*(1-f0(energy_n[counter],efef,T))*energy_n[counter]/(k_B*T);
                    // Part of equation (54) of Rode's book
                    integral_denominator_n = integral_denominator_n + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
                    *f0(energy_n[counter],efef,T)*(1-f0(energy_n[counter],efef,T));
                    // Part of equation (54) of Rode's book
                }
            }

            df0dz_integral_n = integral_numerator_n/integral_denominator_n;
            //cout<<"df0dz_integral_n = "<<df0dz_integral_n<<endl;

            N_poph_atT = N_poph(omega_LO,T);

            //cout<<"N_poph_atT = "<<N_poph_atT<<endl;
            

            if (scattering_mechanisms[8]==1)   // intravalley scattering
            {
                for (int aa=0;aa<iv_number;aa++)
                    N_e[aa] = N_poph(we[aa],T);
            }


	    double k_dum;	
            for (int counter = 0;counter<points;counter++)
            {
                k_dum = k_grid[counter];
                // unit 1/nm
                //cout<<"counter+1 = "<<counter+1<<endl;
                //cout<<"k_dum = "<<k_dum<<endl;
                
		 // polar optical phonon scattering 
                if (scattering_mechanisms[1]==1)
                {
                    kplus_grid[counter] = kplus(counter,omega_LO,points);
                    kminus_grid[counter] = kminus(counter,omega_LO,points);
                    //cout<<"kplus_grid[counter] = "<<kplus_grid[counter]<<endl;
                    //cout<<"kminus_grid[counter] = "<<kminus_grid[counter]<<endl;
                }
		
		// ionized impourity scattering
                if (scattering_mechanisms[0]==1)
                {
                    //cout<<"k_dum = "<<k_dum<<endl;
                    //cout<<"beta_constant =  "<<beta_constant<<endl;
                    //cout<<"c_n[counter] =  "<<c_n[counter]<<endl;
                    //cout<<"counter =  "<<counter+1<<endl;

                    B_ii = (4*(k_dum*k_dum)/(beta_constant*beta_constant))/(1+4*k_dum*k_dum/(beta_constant*beta_constant))
                    +8*(beta_constant*beta_constant+2*k_dum*k_dum)/(beta_constant*beta_constant+4*k_dum*k_dum)*(pow(c_n[counter],2))
                    +(3*pow(beta_constant,4)+
                      6*pow(beta_constant,2)*k_dum*k_dum-8*k_dum*k_dum*k_dum*k_dum)/((beta_constant*beta_constant+4*k_dum*k_dum)*k_dum*k_dum)*pow(c_n[counter],4);
                    // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)


                    D_ii = 1+(2*pow(beta_constant,2)*(pow(c_n[counter],2)/pow(k_dum,2))+
                              (3*pow(beta_constant,4)*(pow(c_n[counter],4))/(4*pow(k_dum,4))));
                              // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)

                    nu_ionizedimpurity[counter] = nu_ii(k_dum,counter,beta_constant,v_n[counter],epsilon_s[T_loop]);
		     // unit 1/second	
                    //cout<<"B_ii = "<<B_ii<<endl;
                    //cout<<"D_ii = "<<D_ii<<endl;
                    //cout<<"N_ii = "<<N_ii<<endl;
                    //cout<<"nu_ionizedimpurity[counter] = "<<nu_ionizedimpurity[counter]<<endl;
                }

		// POP scattering
                if (scattering_mechanisms[1]==1)
                {
                    Aminus_grid[counter] = Aminus(counter,omega_LO,points);
                    Aplus_grid[counter] = Aplus(counter,omega_LO,points);


                    betaplus_grid[counter] = betaplus(counter,omega_LO,epsilon_s[T_loop],epsilon_inf[T_loop],points);
                    betaminus_grid[counter] = betaminus(counter,omega_LO,epsilon_s[T_loop],epsilon_inf[T_loop],points);

                    lambda_i_plus_grid[counter] = abs(lambda_i_plus(counter,omega_LO,Aplus_grid[counter], epsilon_s[T_loop],epsilon_inf[T_loop], points));

                    lambda_i_minus_grid[counter] = abs(lambda_i_minus(counter,omega_LO,Aminus_grid[counter],                epsilon_s[T_loop],epsilon_inf[T_loop],points));

                    lambda_o_plus_grid[counter] = abs(lambda_o_plus(counter,omega_LO,Aplus_grid[counter],
                    epsilon_s[T_loop],epsilon_inf[T_loop],points));

                    lambda_o_minus_grid[counter] = abs(lambda_o_minus(counter, omega_LO,Aminus_grid[counter],
                    epsilon_s[T_loop],epsilon_inf[T_loop],points));
                }
                
		// acoustic scattering 
                if (scattering_mechanisms[3]==1)
                    nu_deformation[counter] = nu_de(k_dum,counter,T,v_n[counter]);
                /*
                if (T==200)
                {
                    cout<<"counter = "<<counter+1<<endl;
                    cout<<"nu_deformation[counter] = "<<nu_deformation[counter]<<endl;
                    getchar();
                }
                */
                
		// Piezoelectric scattering
                if (scattering_mechanisms[4]==1)

                    nu_piezoelectric[counter] = nu_pe(k_dum,counter,T,P_piezo[T_loop],epsilon_s[T_loop],v_n[counter]);

		// Dislocation scattering
                if (scattering_mechanisms[6]==1)
                    nu_dislocation[counter] = nu_dis(k_dum,counter,T,beta_constant,epsilon_s[T_loop], v_n[counter]);

		// Alloy scattering
                if (scattering_mechanisms[7]==1)
                    nu_alloy[counter] = nu_alloy1(k_dum, v_n[counter]);


                if (scattering_mechanisms[8]==1)  // intravalley scattering
                {
                    for (int aa = 0;aa<iv_number;aa++)
                    {
                        lambda_e_plus_grid[counter][aa] = abs(lambda_e_plus(counter,we[aa],rho,De[aa],nfv[aa],points));
                        lambda_e_minus_grid[counter][aa] = abs(lambda_e_minus(counter,we[aa],rho,De[aa],nfv[aa],points));
                        // Equation number 129 of rode book
                    }
                }

                if (scattering_mechanisms[9]==1)  // neutral impurity
                    {			
			nu_neutralimpurity[counter] = nu_im(k_dum,counter,epsilon_s[T_loop],N_im[ii],v_n[counter]);
		    }

                df0dk_grid[counter] = df0dk(k_dum, T, E_F, coefficients_cond, kindex_cond, a11);

                f_dist[counter] = f0(energy_n[counter],E_F,T);
                //cout<<"In between "<<endl;
                //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
                //cout<<"E_F = "<<E_F<<endl;
                //cout<<"T = "<<T<<endl;

                thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, E_F, T, df0dz_integral_n,coefficients_cond, kindex_cond, a11);

                f0x1_f0[counter] = f0(energy_n[counter],E_F,T)*(1-f0(energy_n[counter],E_F,T));

                nu_el[counter] = nu_deformation[counter] + nu_piezoelectric[counter] + nu_ionizedimpurity[counter]
                + nu_dislocation[counter] + nu_alloy[counter] + nu_neutralimpurity[counter];

                electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
		// unit is 1/s 
		
		//---------------------------- code to debug -------------------------------------------------------------
                //cout<<"counter+1 = "<<counter+1<<endl;
                //cout<<"Aplus_grid[counter] =   "<<Aplus_grid[counter]<<endl;
                //cout<<"Aminus_grid[counter] =   "<<Aminus_grid[counter]<<endl;

                //cout<<"betaplus_grid[counter] =  "<<betaplus_grid[counter]<<endl;
                //cout<<"betaminus_grid[counter] =  "<<betaminus_grid[counter]<<endl;

                //cout<<"lambda_i_plus_grid[counter] =  "<<lambda_i_plus_grid[counter]<<endl;
                //cout<<"lambda_i_minus_grid[counter] =  "<<lambda_i_minus_grid[counter]<<endl;
                //cout<<"lambda_o_plus_grid[counter] =  "<<lambda_o_plus_grid[counter]<<endl;
                //cout<<"lambda_o_minus_grid[counter] =  "<<lambda_o_minus_grid[counter]<<endl;

                //cout<<"lambda_e_plus_grid[counter][aa] =  "<<lambda_e_plus_grid[counter][0]<<endl;
                //cout<<"lambda_e_minus_grid[counter][aa] =  "<<lambda_e_minus_grid[counter][0]<<endl;

                //cout<<"nu_deformation[counter] =  "<<nu_deformation[counter]<<endl;
                //cout<<"nu_piezoelectric[counter] =  "<<nu_piezoelectric[counter]<<endl;
                //cout<<"nu_ionizedimpurity[counter] =  "<<nu_ionizedimpurity[counter]<<endl;
                //cout<<"nu_dislocation[counter] =  "<<nu_dislocation[counter]<<endl;
                //cout<<"nu_alloy[counter] =  "<<nu_alloy[counter]<<endl;
                //cout<<"nu_neutralimpurity[counter] =  "<<nu_neutralimpurity[counter]<<endl;
                //cout<<"nu_el[counter] = "<<nu_el[counter]<<endl;
                //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
                //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
                //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
                //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
                //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;

                //if (counter==100||counter==200||counter==300||counter==400||counter==500)
                //    getchar();
		//---------------------------- code to debug -------------------------------------------------------------

            }
 
        if (scattering_mechanisms[8]==1)  // intravalley scattering
        {
		int len = sizeof(nu_iv_total)/sizeof(nu_iv_total[0]);	
		for (int counter = 0;counter<len;counter++)
			nu_iv_total[counter] = 0;

		//cout<<endl<<"Inside"<<endl;
		//cout<<" iv_number = "<<iv_number<<endl;
            for (int counter = 0;counter<points;counter++)
            {
            	//cout<<"counter = "<<counter<<endl;
		for (int aa = 0;aa<iv_number;aa++)
		{
		    //cout<<"aa = "<<aa<<endl;
		    double k_minus = kminus(counter,we[aa],points);
		    //cout<<"k_minus = "<<k_minus<<endl;

		    double arr[points];
		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_minus);
		    int minus_index =FindMinInd(arr,points);


		    double k_plus = kplus(counter,we[aa],points);
		    //cout<<"k_plus = "<<k_plus<<endl;

		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_plus);
		    int plus_index =FindMinInd(arr,points);

		    //cout<<"plus_index = "<<plus_index<<endl;

		    double f_negative = f0(energy_n[minus_index],efef,T);
		    double f_positive =  f0(energy_n[plus_index],efef,T);

		    //cout<<"f_negative = "<<f_negative<<endl;
		    //cout<<"f_positive = "<<f_positive<<endl;

		    if (energy_n[counter] < h_bar*we[aa])
		    {
			    nu_iv[counter][aa] = (N_e[aa] + f_positive) *lambda_e_plus_grid[plus_index][aa];
		    }
		    else
		    {    
			nu_iv[counter][aa] = (N_e[aa] + 1 - f_negative) * lambda_e_minus_grid[minus_index][aa] + 
						(N_e[aa] + f_positive)*lambda_e_plus_grid[plus_index][aa];
		    }
		    
		    nu_iv_total[counter] = nu_iv_total[counter] + nu_iv[counter][aa];            
		    
		    //cout<<"nu_iv[counter][aa] = "<<nu_iv[counter][aa]<<endl;
		    //cout<<"nu_iv_total[counter] = "<<nu_iv_total[counter]<<endl;
		    //getchar(); 	
		     	
		}
		nu_el[counter] = nu_el[counter] + nu_iv_total[counter];  
            }
        }

//------------------------------------ components for BTE END --------------------------------------------------------------
// ----------------------------saving data ---------------------------------------------------
            /*
            FILE *fid1;
            fid1 = fopen("Aplus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, Aplus_grid[i]);
	fclose(fid1);

            fid1 = fopen("Aminus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, Aminus_grid[i]);
	fclose(fid1);

            fid1 = fopen("betaplus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, betaplus_grid[i]);
	fclose(fid1);

            fid1 = fopen("betaminus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, betaminus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_i_plus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_i_plus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_i_minus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_i_minus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_o_plus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_o_plus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_o_minus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_o_minus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_e_plus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_e_plus_grid[i]);
	fclose(fid1);

            fid1 = fopen("lambda_e_minus.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, lambda_e_minus_grid[i]);
	fclose(fid1);
            */
            
	    /* 			
            fid1 = fopen("nu_deformation.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_deformation[i]);
	fclose(fid1);

            fid1 = fopen("nu_piezoelectric.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_piezoelectric[i]);
	fclose(fid1);

            fid1 = fopen("nu_ionizedimpurity.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_ionizedimpurity[i]);
	fclose(fid1);
            
            fid1 = fopen("nu_dislocation.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_dislocation[i]);
	fclose(fid1);    

            fid1 = fopen("nu_alloy.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_alloy[i]);
	fclose(fid1);

            fid1 = fopen("nu_neutralimpurity.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_neutralimpurity[i]);
	fclose(fid1);
	    
            //FILE *fid1;
            fid1 = fopen("nu_el.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, nu_el[i]);
	fclose(fid1);

            fid1 = fopen("df0dk_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, df0dk_grid[i]);
	fclose(fid1);

            fid1 = fopen("f_dist.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, f_dist[i]);
	fclose(fid1);

            fid1 = fopen("thermal_driving_force.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, thermal_driving_force[i]);
	fclose(fid1);

            fid1 = fopen("f0x1_f0.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, f0x1_f0[i]);
	fclose(fid1);

            fid1 = fopen("elecgtric_driving_force.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, electric_driving_force[i]);
	fclose(fid1);
            */
            //getchar();

// ----------------------------saving data ---------------------------------------------------
//------------------------------ reading data -----------------------------------------------
            /*
            fid1 = fopen("nu_deformation.txt","r");
            for (int i = 0; i < points; i++)
            {
                fgets(line, 1000, fid1);
                sscanf(line, "%lf", &nu_deformation[i]);
            }
	fclose(fid1);

           fid1 = fopen("nu_ionizedimpurity.txt","r");
            for (int i = 0; i < points; i++)
            {
                fgets(line, 1000, fid1);
                sscanf(line, "%lf", &nu_ionizedimpurity[i]);
            }
	fclose(fid1);

           fid1 = fopen("nu_piezoelectric.txt","r");
            for (int i = 0; i < points; i++)
            {
                fgets(line, 1000, fid1);
                sscanf(line, "%lf", &nu_piezoelectric[i]);
            }
	fclose(fid1);
            */

//------------------------------------------------------------------------------------------------------------------
}


