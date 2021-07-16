#include"main.h"

void calculate_mobility(double T, int ii)
{

//----------------------------------calculate drift mobility -------------------------------------------------------------
            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);

            cout<<endl;
            // ionized impurity scattering
            if (scattering_mechanisms[0]==1)
            {
                mobility_ii = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,g,points,a11);
                cout<<"mobility_ii = "<<mobility_ii<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_ii = 1e10;


            // POP scattering
            if (scattering_mechanisms[1]==1)
            {
                mobility_po = mu_po(E_F,T,coefficients_cond,kindex_cond,g_LO,g,nu_el,points,a11);
                cout<<"mobility_po = "<<mobility_po<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_po=1e10;

            // Acoustic deformation scattering
            if (scattering_mechanisms[3]==1)
            {
                mobility_de = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_deformation,g,points,a11);
                cout<<"mobility_de = "<<mobility_de<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_de=1e10;

            // piezoelectric scattering
            if (scattering_mechanisms[4]==1)
            {
                mobility_pe = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,g,points,a11);
                cout<<"mobility_pe = "<<mobility_pe<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_pe=1e10;

            // dislocation scattering
            if (scattering_mechanisms[6]==1)
            {
                mobility_dis = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_dislocation,g,points,a11);
                cout<<"mobility_dis = "<<mobility_dis<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_dis=1e10;



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
                mobility_to = ( pow(2,0.5)*pi*e*pow(h_bar,3)*rho*(exp(h_bar*omega_TO/(k_B*T))-1)*vs*vs ) / ( pow(m,2.5)*pow(m_e,2.5)*omega_TO*E_deformation_n*E_deformation_n*
                                    pow((a+h_bar*omega_TO),0.5) ) * 1e4*pow(e,0.5);
                                    //Last coeff. is to convert to cm2/V.s
                //nu_to(:) = e/(m*m_e*mobility_to)*1e4;
                // Last part is to convert units to [1/s]

                cout<<"mobility_to = "<<mobility_to<<" cm^2/(V-s)"<<endl;
            }

            // alloy scattering
            if (scattering_mechanisms[7]==1)
            {
                mobility_alloy = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_alloy,g,points,a11);
                cout<<"mobility_alloy = "<<mobility_alloy<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_alloy=1e10;

            // inter-valley scattering
            if (scattering_mechanisms[8]==1)
            {
                mobility_iv = mu_elastic(E_F,T,coefficients_cond,kindex_cond,nu_iv_total,g,points,a11);
                cout<<"mobility_iv  = "<<mobility_iv<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_iv=1e10;

            // neutral impurity scattering
            if (scattering_mechanisms[9]==1)
            {
                mobility_neutral = mu_elastic(E_F,T,coefficients_cond,kindex_cond, nu_neutralimpurity,g,points,a11);

                cout<<"mobility_neutral = "<<mobility_neutral<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_neutral=1e10;

            double nu_to[points]={0};
	    
	    
            mobility_all[0] = mobility_ii;
            mobility_all[1] = mobility_po;
            mobility_all[2] = mobility_npop;
            mobility_all[3] = mobility_de;
            mobility_all[4] = mobility_pe;
            mobility_all[5] = mobility_dis;
            mobility_all[6] = mobility_to;
            mobility_all[7] = mobility_alloy;
            mobility_all[8] = mobility_iv;
            mobility_all[9] = mobility_neutral;
            //scattering_mechanisms
            double sum =0 ;
            for (int i=0;i<10;i++)
            {
                if (scattering_mechanisms[i]!=0 && mobility_all[i]!=0)
                sum = sum + 1/ (mobility_all[i]*scattering_mechanisms[i]);
            }

            //cout<<"sum = "<<sum<<endl;

            mobility_avg = 1/sum;

            mobility = mu_overall(E_F,T,coefficients_cond,kindex_cond,g,nu_el,points,a11);
	    // unit cm^2/(V-s)
	    
            mobility_rta = mu_overall(E_F,T,coefficients_cond,kindex_cond,g_rta,nu_el,points,a11);
	     // unit cm^2/(V-s)
	     	
            if (omega_TO > 0.0)
                mobility = 1 / (1/mobility + 1/mobility_to);
		// unit cm^2/(V-s)

            if (mobility < 0)
                    mobility = mobility_avg;
		// unit cm^2/(V-s)
		
		
//----------------------------------calculate mobility -------------------------------------------------------------

            sigma = mobility *  n0 * e;

            sigma_rta = mobility_rta * n0 * e;

	    if (n0 == 0)
	    {
	    	sigma = mobility *  abs(n_e) * e;
	    	sigma_rta = mobility_rta *  abs(n_e) * e;
	    }
            thermopower = -k_B*(df0dz_integral_n- E_F /(k_B*T))*1e6 + (J(T,m,g_th,points)/sigma)/dTdz*1e6;
            // Equation No. 52 of rode book


//----------------------------------------------------------------------------------------------------------------------



// -------------------------------- calculate hall mobility -------------------------------------------------------------
	if(Bfield!=0)
	{

            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);

            cout<<endl;
            // ionized impurity scattering
            if (scattering_mechanisms[0]==1)
            {
                mobility_hall_ii = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_ionizedimpurity,points,a11);
                cout<<"mobility_hall_ii = "<<mobility_hall_ii<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_ii = 1e10;


            // POP scattering
            if (scattering_mechanisms[1]==1)
            {
                mobility_hall_po = mu_poH(E_F,T,coefficients_cond,kindex_cond,gH_LO,hH_LO,nu_el,points,a11);
                cout<<"mobility_hall_po = "<<mobility_hall_po<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_po=1e10;

            // Acoustic deformation scattering
            if (scattering_mechanisms[3]==1)
            {
                mobility_hall_de = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_deformation,points,a11);
                cout<<"mobility_hall_de = "<<mobility_hall_de<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_de=1e10;

            // piezoelectric scattering
            if (scattering_mechanisms[4]==1)
            {
                mobility_hall_pe = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_piezoelectric,points,a11);
                cout<<"mobility_hall_pe = "<<mobility_hall_pe<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_pe=1e10;

            // dislocation scattering
            if (scattering_mechanisms[6]==1)
            {
                mobility_hall_dis = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_dislocation,points,a11);
                cout<<"mobility_hall_dis = "<<mobility_hall_dis<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_dis=1e10;



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
                mobility_hall_to = ( pow(2,0.5)*pi*e*pow(h_bar,3)*rho*(exp(h_bar*omega_TO/(k_B*T))-1)*vs*vs ) / ( pow(m,2.5)*pow(m_e,2.5)*omega_TO*E_deformation_n*E_deformation_n*
                                    pow((a+h_bar*omega_TO),0.5) ) * 1e4*pow(e,0.5);
                                    //Last coeff. is to convert to cm2/V.s
                //nu_to(:) = e/(m*m_e*mobility_to)*1e4;
                // Last part is to convert units to [1/s]

                cout<<"mobility_hall_to = "<<mobility_hall_to<<" cm^2/(V-s)"<<endl;
            }

            // alloy scattering
            if (scattering_mechanisms[7]==1)
            {
                mobility_hall_alloy = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_alloy,points,a11);
                cout<<"mobility_hall_alloy = "<<mobility_hall_alloy<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_alloy=1e10;

            // inter-valley scattering
            if (scattering_mechanisms[8]==1)
            {
                mobility_hall_iv = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,nu_iv_total,points,a11);
                cout<<"mobility_hall_iv  = "<<mobility_hall_iv<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_iv=1e10;

            // neutral impurity scattering
            if (scattering_mechanisms[9]==1)
            {
                mobility_hall_neutral = mu_elasticH(E_F,T,coefficients_cond,kindex_cond,
                    nu_neutralimpurity,points,a11);

                cout<<"mobility_hall_neutral = "<<mobility_hall_neutral<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_neutral=1e10;

            double nu_hall_to[points]={0};
	    
	    
            mobility_hall_all[0] = mobility_hall_ii;
            mobility_hall_all[1] = mobility_hall_po;
            mobility_hall_all[2] = mobility_hall_npop;
            mobility_hall_all[3] = mobility_hall_de;
            mobility_hall_all[4] = mobility_hall_pe;
            mobility_hall_all[5] = mobility_hall_dis;
            mobility_hall_all[6] = mobility_hall_to;
            mobility_hall_all[7] = mobility_hall_alloy;
            mobility_hall_all[8] = mobility_hall_iv;
            mobility_hall_all[9] = mobility_hall_neutral;
            //scattering_mechanisms
            
            double sum_hall =0 ;
            for (int i=0;i<10;i++)
            {
                if (scattering_mechanisms[i]!=0 && mobility_hall_all[i]!=0)
                sum_hall = sum_hall + 1/ (mobility_hall_all[i]*scattering_mechanisms[i]);
            }

            //cout<<"sum_hall = "<<sum_hall<<endl;

            mobility_hall_avg = 1/sum_hall;

            mobility_hall = mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH,hH,nu_el,points,a11);

            mobility_hall_rta = mu_overallH(E_F,T,coefficients_cond,kindex_cond,gH_rta, hH_rta, nu_el,points,a11);

            if (omega_TO > 0.0)
                mobility_hall = 1 / (1/mobility_hall + 1/mobility_hall_to);


            if (mobility_hall < 0)
                    mobility_hall = mobility_hall_avg;

	     hall_factor1 = mobility_hall/mobility;
	     hall_factor_rta1 = mobility_hall_rta/mobility_rta;
	     
	}	     

//----------------------------------calculated hall mobility -------------------------------------------------------------
	if(Bfield!=0)
	{


            sigma_hall = mobility_hall *  n0 * e;

            sigma_hall_rta = mobility_hall_rta * n0 * e;

	    if (n0 == 0)
	    {
	    	sigma_hall = mobility_hall *  abs(n_e) * e;
	    	sigma_hall_rta = mobility_hall_rta *  abs(n_e) * e;
	    }

	}
//--------------------------------------------------------------------------------------------------------------------------------



            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);
            cout<<endl<<"Drift mobility results"<<endl;
            cout<<"Temperature = "<<T<<" K"<<endl;
            cout<<"doping = "<<n_array[ii]<<" cm^-3"<<endl;
            cout<<"mobility = "<<mobility<<" cm^2/(V-s)"<<endl;
            cout<<"mobility_rta = "<<mobility_rta<<" cm^2/(V-s)"<<endl;
            cout<<"mobility_avg = "<<mobility_avg<<" cm^2/(V-s)"<<endl;
            cout<<"sigma = "<<sigma<<" S/cm "<<endl;
            cout<<"sigma_rta = "<<sigma_rta<<" S/cm "<<endl;
            cout<<"thermopower = "<<thermopower<<" V/K"<<endl<<endl;

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
	
            cc = cc+1;

            calc_mobility[cc][1] = mobility;
            calc_mobility_rta[cc][1] = mobility_rta;
            calc_thermopower[cc][1] = thermopower;
            calc_sigma[cc][1] = sigma;
            calc_sigma_rta[cc][1] = sigma_rta;


            calc_mobility_ii[cc][1] = mobility_ii;
            calc_mobility_po[cc][1] = mobility_po;
            calc_mobility_de[cc][1] = mobility_de;
            calc_mobility_pe[cc][1] = mobility_pe;
            calc_mobility_dis[cc][1] = mobility_dis;
            calc_mobility_to[cc][1] = mobility_to;
            calc_mobility_alloy[cc][1] = mobility_alloy;
            calc_mobility_iv[cc][1] = mobility_iv;
            calc_mobility_neutral[cc][1] = mobility_neutral;

		if(Bfield!=0)
		{

		    calc_mobility_hall[cc][1] = mobility_hall;
		    calc_mobility_rta[cc][1] = mobility_hall_rta;
		    calc_sigma_hall[cc][1] = sigma_hall;
		    calc_sigma_hall_rta[cc][1] = sigma_hall_rta;


		    calc_mobility_hall_ii[cc][1] = mobility_hall_ii;
		    calc_mobility_hall_po[cc][1] = mobility_hall_po;
		    calc_mobility_hall_de[cc][1] = mobility_hall_de;
		    calc_mobility_hall_pe[cc][1] = mobility_hall_pe;
		    calc_mobility_hall_dis[cc][1] = mobility_hall_dis;
		    calc_mobility_hall_to[cc][1] = mobility_hall_to;
		    calc_mobility_hall_alloy[cc][1] = mobility_hall_alloy;
		    calc_mobility_hall_iv[cc][1] = mobility_hall_iv;
		    calc_mobility_hall_neutral[cc][1] = mobility_hall_neutral;
		    hall_factor[cc][1] = hall_factor1; 
		    hall_factor_rta[cc][1] = hall_factor_rta1; 
		}
	
}
