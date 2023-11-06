
#include"main.h"


void find_fermi(double n, double T, int ii)
{
		
    if (n > 5e20)
        n = 5e20;

    //cout<<"n = "<<n<<endl;
    //cout<<" T=  "<<T<<endl;
    //cout<<"points = "<<points<<endl;
    //cout<<" volume1  "<<volume1<<endl;
    // unit 1/cm^3 
    //cout<<" thickness  "<<thickness<<"   m "<<endl;
    //getchar();
    
    
    double e_f,E1,E2,E_mid,n1,n2,n_mid,temp,E_mid_old,integral_n,integral_p,x = 1.0;
    double de_n,de_p,dk;
	
	int flag=0;

	while(flag==0)
	{		
		e_f = x;
		integral_n = 0;
		integral_p = 0;
		if (free_e==0)
		{
		    for(int counter=0;counter<=points-2;counter++)
		    {
		        de_n = energy_n[counter+1]-energy_n[counter];
		        de_p = energy_p[counter+1]-energy_p[counter];

			if(de_n<0)
				de_n = 0;
				
			if(de_p<0)
				de_p = 0;
			
			// Ds_n unit is (per eV per cm^3) for 3D and 2D for VASP inputs
			// Ds_n unit is (per eV per cm^2) for VASP =2  inputs tight binding and thickness = 0.01, volume1 = 1

			if(type == "n")
			{
				integral_n = integral_n+
				de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],e_f,T)*N_cb; 
				// 100 is multiplied to convert thicknesss from m to cm
				// volume1 unit cm^3 	
				
				integral_p = integral_p+
				de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
				// 100 is multiplied to convert thicknesss from m to cm
				// volume1 unit cm^3 	

				// integral_n unit is 1/cm^3 ; for 3D
				// integral_n unit is 1/cm^2 ; for 2D
				 
			}
			else
			{
				integral_n = integral_n + 
				de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
				// 100 is multiplied to convert thicknesss from m to cm
				// volume1 unit cm^3 	
				
				integral_p = integral_p + 
				de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],e_f,T)*N_vb;
				// 100 is multiplied to convert thicknesss from m to cm
				// volume1 unit cm^3 	

				// integral_n unit is 1/cm^3 ; for 3D
				// integral_n unit is 1/cm^2 ; for 2D
			}
			
			/*
		        cout<<"counter = "<<counter+1<<endl;
		        cout<<"de_n = "<<de_n<<endl;
		        cout<<"de_p = "<<de_p<<endl;
		        cout<<"N_cb = "<<N_cb<<endl;
		        cout<<"integral_n = "<<integral_n<<endl;
		        cout<<"integral_p = "<<integral_p<<endl;
		        getchar();
		        //*/
		    }

		}   // if free==0 condition 
		else
		{
		    for(int counter = 0;counter<=points-2;counter++)
		    {
		        dk = (k_grid[counter+1]-k_grid[counter]);
			if(type == "n")
			{
				integral_n = integral_n + 
				(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],e_f,T)*N_cb;
				integral_p = integral_p + 
				(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
				// multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
				// integral_n unit is 1/cm^3 ; for 3D
				// integral_n unit is 1/cm^2 ; for 2D
			}
			else
			{
				integral_n = integral_n+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
				integral_p = integral_p+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],e_f,T)*N_vb;
				// multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
				// integral_n unit is 1/cm^3 ; for 3D
				// integral_n unit is 1/cm^2 ; for 2D
			}	
		    }
		}  // if else free==0 condition finished here

		if(type == "n")
			temp = (integral_n - integral_p);  // to have n in units of [1/cm^3] for 3D and [1/cm^2] for 2D
		else
			temp = (integral_p - integral_n);  // to have n in units of [1/cm^3] for 3D and [1/cm^2] for 2D

		/*
		cout<<"temp = "<<temp<<endl;
		cout<<"n = "<<n<<endl;
		cout<<"x = "<<x<<endl;
        	cout<<"integral_n = "<<integral_n<<endl;  
        	cout<<"integral_p = "<<integral_p<<endl;
		getchar();
		//*/
				
		if(temp < n)
			x = x + 1;
		else
			flag=1;
		
		if(x==6)
		{
			cout<<"Program cannot simulate for such concnetration. Exit from program"<<endl;
	               exit(EXIT_FAILURE);
			break;
		}
		/*			
		cout<<"x = "<<x<<endl;
		cout<<"flag = "<<flag<<endl;
		getchar();
		*/
		
	}  // end of while loop flag==0 
	
	//cout<<"x = "<<x<<endl;
	//getchar();
	
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

		if(de_n<0)
			de_n = 0;
			
		if(de_p<0)
			de_p = 0;
		
		if(type == "n")
		{
		        integral_n = integral_n + 
		        de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],e_f,T)*N_cb;
			// volume1 unit 1/cm^3 	

		        
		        integral_p = integral_p + 
		        de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
			// volume1 unit 1/cm^3 	

			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
			 
		}
		else
		{
		        integral_n = integral_n + 
		        de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
		        integral_p = integral_p + 
		        de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],e_f,T)*N_vb;
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
		
		/*
                cout<<"counter = "<<counter+1<<endl;
                cout<<"de_n = "<<de_n<<endl;
                cout<<"de_p = "<<de_p<<endl;
                cout<<"N_cb = "<<N_cb<<endl;
                cout<<"integral_n = "<<integral_n<<endl;
                cout<<"integral_p = "<<integral_p<<endl;
                getchar();
                //*/
            }
        }
        else
        {
            for(int counter = 0;counter<=points-2;counter++)
            {
                dk = (k_grid[counter+1]-k_grid[counter]);
		if(type == "n")
		{
		        integral_n = integral_n+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],e_f,T)*N_cb;
		        integral_p = integral_p+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
		        // multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
		else
		{
		        integral_n = integral_n+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
		        integral_p = integral_p+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],e_f,T)*N_vb;
		        // multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}	
            }
        }

	if(type == "n")
	        temp = (integral_n - integral_p);  // to have n in units of [1/cm^3] for 3D and [1/cm^2] for 2D
	else
	        temp = (integral_p - integral_n);  // to have n in units of [1/cm^3] for 3D and [1/cm^2] for 2D
	
	
        if (i==1)
        {
        	n1 = temp; 
        		
        	/*
        	cout<<endl<<"n1 = "<<n1<<endl;  
        	cout<<"integral_n = "<<integral_n<<endl;  
        	cout<<"integral_p = "<<integral_p<<endl;
        	cout<<"e_f = "<<e_f<<endl;  
        	getchar();
        	//*/
        }
        else if (i ==2)
        {
            n2 = temp; 
            
            /*
            cout<<endl<<"n2 = "<<n2<<endl;  
		cout<<"integral_n = "<<integral_n<<endl;  
		cout<<"integral_p = "<<integral_p<<endl; 
        	cout<<"e_f = "<<e_f<<endl;
        	getchar();  
		//*/
		 
        }
        else
        {
        	n_mid = temp; 
        	/*
        	cout<<endl<<"n_mid = "<<n_mid<<endl;  
        	cout<<"integral_n = "<<integral_n<<endl;  
        	cout<<"integral_p = "<<integral_p<<endl;  
        	cout<<"e_f = "<<e_f<<endl; 
        	getchar(); 
        	//*/
        }

    }
    //getchar();

    E_mid_old = -20;

    while (abs(abs(n_mid)/n-1) > 0.001)
    {
    	/*
    	cout<<"Before E1 = "<<E1<<endl;
    	cout<<"E2 = "<<E2<<endl;
    	cout<<"E_mid = "<<E_mid<<endl;
    	cout<<"n1 = "<<n1<<endl;
    	cout<<"n2 = "<<n2<<endl;
    	cout<<"n_mid = "<<n_mid<<endl;
    	cout<<"e_f = "<<e_f<<endl;
	//*/
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
	/*
    	cout<<endl<<"After E1 = "<<E1<<endl;
    	cout<<"E2 = "<<E2<<endl;
    	cout<<"E_mid = "<<E_mid<<endl;
    	cout<<"e_f = "<<e_f<<endl;
	//*/
	integral_n = 0;
        integral_p = 0;

        if (free_e ==0)
        {
            for (int counter=0;counter<=points-2;counter++)
            {
                de_n = energy_n[counter+1]-energy_n[counter];
                de_p = energy_p[counter+1]-energy_p[counter];

		if(de_n<0)
			de_n = 0;
			
		if(de_p<0)
			de_p = 0;
                
		if(type == "n")
		{
			integral_n = integral_n+de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],e_f,T)*N_cb;
			integral_p = integral_p+de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
		else
		{
			integral_n = integral_n+de_n*(Ds_n[counter]*thickness*100)/volume1*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
			integral_p = integral_p+de_p*(Ds_p[counter]*thickness*100)/volume1*f0(energy_p[counter],e_f,T)*N_vb;
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
            }
        }
        else
        {
            for (int counter = 1;counter<=points-2;counter++)
            {
                dk = (k_grid[counter+1]-k_grid[counter]);

		if(type == "n")
		{
		        integral_n = integral_n+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],e_f,T)*N_cb;
		        integral_p = integral_p+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],-(e_f+Bgap[ii]),T)*N_vb;
		        // multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
		else
		{
		        integral_n = integral_n+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_n[counter],-(e_f+Bgap[ii]),T)*N_cb;
		        integral_p = integral_p+(dk*(pow((k_grid[counter]/pi),2))*1e21)*f0(energy_p[counter],e_f,T)*N_vb;
		        // multipiled with 1e21 to convert (1/nm)^3 to (1/cm)^3  
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D
		}
               
            }
        }
	if(type == "n")
	        temp = (integral_n - integral_p);  
        else
	        temp = (integral_p - integral_n);  
        	
        n_mid = temp;
	/*
	cout<<"n_mid = "<<n_mid<<endl;
	cout<<"integral_n = "<<integral_n<<endl;
	cout<<"integral_p = "<<integral_p<<endl;
	getchar();
	//*/
        if (E1 == E2 || E_mid_old == E_mid)
        {
        	if(n==0)
        	{
        		if(n_mid > 1)
        		{
        			cout<<"Calculated concentration is not so accurate, it may lead to wrong answer"<<endl;
        		}
        	}
        	else
        	{
        		cout<<"Calculated concentration is not so accurate, it may lead to wrong answer"<<endl;
        	}


            //getchar();
            /*
        	cout<<"n_mid = "<<n_mid<<endl;
        	cout<<"integral_n = "<<integral_n<<endl;
        	cout<<"integral_p = "<<integral_p<<endl;
        	getchar();
        	//*/
            break;
        }
        E_mid_old = E_mid;
    }

    E_F = e_f;
    n_e = integral_n;    // to have n in units of [1/cm^3]
    n_h = integral_p;    // to have n in units of [1/cm^3]
			// integral_n unit is 1/cm^3 ; for 3D
			// integral_n unit is 1/cm^2 ; for 2D

	if(geometry==1)  // for 3D
	{
		cout<<"E_F = "<<E_F<<" eV"<<endl;
		cout<<"n_e = "<<n_e<<" cm^-3"<<endl;
		cout<<"n_h = "<<n_h<<" cm^-3"<<endl;
	}
	else if(geometry==2)  // for 2D
	{
		cout<<"E_F = "<<E_F<<" eV"<<endl;
		cout<<"n_e = "<<n_e<<" cm^-2"<<endl;
		cout<<"n_h = "<<n_h<<" cm^-2"<<endl;
	}
	//cout<<"n = "<<n<<endl;
		
	//cout<<"(n_e-n_h)/n = "<<(n_e-n_h)/n<<endl;
	
        if (abs(n_e-n_h)/n < 0.5 && abs(n_e-n_h)/n > 1.5)
        {
            cout<<"Calculated concentration is not within 50% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if (abs(n_e-n_h)/n < 0.6 && abs(n_e-n_h)/n > 1.4)
        {
            cout<<"Calculated concentration is not within 60% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if (abs(n_e-n_h)/n < 0.8 && abs(n_e-n_h)/n > 1.2)
        {
            cout<<"Calculated concentration is not within 80% of given concentrtion. Programs end here"<<endl;
            exit(EXIT_FAILURE);
        }

        if (abs(n_e-n_h)/n < 0.9 && abs(n_e-n_h)/n > 1.1)
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

	    if(geometry==1)  // for 3D
	    {	
		    double a1;
		    if (scattering_mechanisms[6]==1)   // disloaction scattering
			a1 = abs(N_dis/c_lattice*1e7);
		    else
			a1=0;

		    N_ii = (n_e + n_h) + a1;

		    if(scattering_mechanisms[0] == 1)
		    {
		    	cout<<"Net ionized donors concentration =  "<<N_ii<<"  cm^(-3)"<<endl;
		    }

		    if(scattering_mechanisms[9] == 1)	// Neutral impurity scattering  // for 3D
			{
		    	cout<<"Net neutral impurity =  "<<N_im[ii]<<" cm^(-3)"<<endl;
			}
			    // net neutral impurity
	    }	
}


