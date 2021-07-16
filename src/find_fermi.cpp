
#include"main.h"


void find_fermi(double n, double T, int ii)
{

    if (n > 5e20)
        n = 5e20;

    //cout<<"n = "<<n<<endl;
    //cout<<" T=  "<<T<<endl;
    //cout<<"points = "<<points<<endl;
    //cout<<" volume1  "<<volume1<<endl;

    double e_f,E1,E2,E_mid,n1,n2,n_mid,temp,E_mid_old,integral_n,integral_p,x = 1.0;
    double de_n,de_p,dk;

    for(int i=1;i<=3;i++)
    {
        if (i==1)
        {
            E1 = x;
            e_f = E1;
        }
        else if (i==2)
        {
            E2 = -1*(Bgap[ii] + x);
            e_f = E2;
        }
        else
        {
            E_mid = (E1+E2)/2;
            e_f = E_mid;
        }

        integral_n = 0;
        integral_p = 0;

        if (free_e==0)
        {
            for(int counter=0;counter<=points-2;counter++)
            {
                de_n = energy_n[counter+1]-energy_n[counter];
                de_p = energy_p[counter+1]-energy_p[counter];

                integral_n = integral_n+de_n*Ds_n[counter]/volume1*1000*f0(energy_n[counter],e_f,T)*N_cb;
                integral_p = integral_p+de_p*Ds_p[counter]/volume1*1000*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;

                //cout<<"counter = "<<counter+1<<endl;
                //cout<<"de_n = "<<de_n<<endl;
                //cout<<"de_p = "<<de_p<<endl;
                //cout<<"N_cb = "<<N_cb<<endl;
                //cout<<"integral_n = "<<integral_n<<endl;
                //cout<<"integral_p = "<<integral_p<<endl;
                //getchar();
            }
        }
        else
        {
            for(int counter = 0;counter<=points-2;counter++)
            {
                dk = (k_grid[counter+1]-k_grid[counter]);
                integral_n = integral_n+dk*pow((k_grid[counter]/pi),2)*f0(energy_n[counter],e_f,T)*N_cb;
                integral_p = integral_p+dk*pow((k_grid[counter]/pi),2)*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
            }
        }

        temp = (integral_n - integral_p)*1e21;  // to have n in units of [1/cm^3]

        if (i==1)
            n1 = temp; //cout<<"n1 = "<<n1<<endl;  }
        else if (i ==2)
            n2 = temp; //cout<<"n2 = "<<n2<<endl;  }
        else
            n_mid = temp; //cout<<"n_mid = "<<n_mid<<endl;  }

    }
    //getchar();

    E_mid_old = -20;

    while (abs(abs(n_mid)/n-1) > 0.001)
    {
        if (n1 > n && n2 < n)
        {
            if (n_mid < n)
            {
                E2 = E_mid;
                E_mid = (E1+E2)/2;
                e_f = E_mid;
            }
            if (n_mid > n)
            {
                E1=E_mid;
                E_mid=(E1+E2)/2;
                e_f = E_mid;
            }
        }

        integral_n = 0;
        integral_p = 0;

        if (free_e ==0)
        {
            for (int counter=0;counter<=points-2;counter++)
            {
                de_n = energy_n[counter+1]-energy_n[counter];
                de_p = energy_p[counter+1]-energy_p[counter];
                integral_n = integral_n+de_n*Ds_n[counter]/volume1*1000*f0(energy_n[counter],e_f,T)*N_cb;
                integral_p = integral_p+de_p*Ds_p[counter]/volume1*1000*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
            }
        }
        else
        {
            for (int counter = 1;counter<=points-2;counter++)
            {
                dk = (k_grid[counter+1]-k_grid[counter]);
                integral_n = integral_n+dk*pow((k_grid[counter]/pi),2)*f0(energy_n[counter],e_f,T)*N_cb;
                integral_p = integral_p+dk*pow((k_grid[counter]/pi),2)*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
                // integral_p = integral_p+dk*(k_grid(counter)/pi)^2*(1-f0(energy_p(counter),-(e_f+Bgap[ii]),T))*N_vb;
            }
        }

        temp = (integral_n - integral_p)*1e21;  // to have n in units of [1/cm^3]
        n_mid = temp;

        if (E1 == E2 || E_mid_old == E_mid)
        {
            //cout<<"Calculated concentration is not so accurate, it may lead to wrong answer"<<endl;
            break;
        }
        E_mid_old = E_mid;
    }

    E_F = e_f;
    n_e = integral_n * 1e21;    // to have n in units of [1/cm^3]
    n_h = integral_p * 1e21;    // to have n in units of [1/cm^3]


        cout<<"E_F = "<<E_F<<" eV"<<endl;
        cout<<"n_e = "<<n_e<<" cm^-3"<<endl;
        cout<<"n_h = "<<n_h<<" cm^-3"<<endl;
	
	//cout<<"n = "<<n<<endl;
		
	//cout<<"(n_e-n_h)/n = "<<(n_e-n_h)/n<<endl;
	
        if ((n_e-n_h)/n < 0.5 && (n_e-n_h)/n > 1.5)
        {
            cout<<"Calculated concentration is not within 50% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if ((n_e-n_h)/n < 0.6 && (n_e-n_h)/n > 1.4)
        {
            cout<<"Calculated concentration is not within 60% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if ((n_e-n_h)/n < 0.8 && (n_e-n_h)/n > 1.2)
        {
            cout<<"Calculated concentration is not within 80% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if ((n_e-n_h)/n < 0.9 && (n_e-n_h)/n > 1.1)
        {
            cout<<"Calculated concentration is not within 90% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

    
            
            if (isnan(n_e)!=0)
                n_e = 0;
            if (isnan(n_h)!=0)
                n_h = 0;

            if (De_ionization ==0)
            {
                if (n_e < Nd1)
                {
                    n_e = Nd1;
                    //cout<<"Modified ionized Donors =  "<<n_e<<endl;
                }

                if (n_h < Na1)
                {
                    n_h = Na1;
                    //cout<<"Modified ionized Acceptor = "<<n_h<<endl;
                }
            }

	    double a1;
	    if (scattering_mechanisms[6]==1)   // disloaction scattering
		a1 = abs(N_dis/c_lattice*1e7);
	    else
		a1=0;
	    				
            N_ii = (n_e + n_h) + a1;
	    cout<<"Net ionized donors concentration =  "<<N_ii<<"  cm^(-3)"<<endl;

            // net neutral impurity
            cout<<"Net neutral impurity =  "<<N_im[ii]<<" cm^(-3)"<<endl;

}


