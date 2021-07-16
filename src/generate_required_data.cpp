#include"main.h"
double E_F,n_e,n_h;

void generate_required_data(double T)
{
	FILE *fid1;	

	    k_min = k_min0;
	    k_step_fine = k_step_fine0;

	    if (T <= 35)
	    {
		k_trans = k_trans0;
		k_min = k_min0/10;
		k_step_fine = k_min;
	    }
	    else if ((T>35)&&(T<100))
	    {
		k_trans = k_trans0;
		k_min = k_min0/5;
	    }
	    else
		k_trans = k_trans0;

	    k_step = k_step0;
	    points1 = floor((k_trans-k_min)/k_step_fine)+1;
	    points2 = ceil((k_max-k_trans)/k_step)+1;
	    points = points1+points2;

	    //cout<<"points = "<<points<<endl;
	    //getchar();
	    if(points>2000)
	    {
	    	cout<<"Error points size is larger than maximum size of arrays given "<<endl;
	    	exit(EXIT_FAILURE);
	    }	
	    		

	    double k_dum;

	    for (int counter=0;counter<points;counter++)
	    {
		if ((k_min+(counter)*k_step_fine)<k_trans)
		    k_dum = k_min+(counter)*k_step_fine;
		else
		    k_dum = k_trans+(counter-points1+1)*k_step;

		k_grid[counter] = k_dum;
		//cout<<"k_grid[counter] = "<<k_grid[counter]<<endl;

		energy_n[counter] = conduction_dispersion(k_dum, coefficients_cond, kindex_cond, a11);
		energy_p[counter] = conduction_dispersion(k_dum, coefficients_val, kindex_val, b11);

		//cout<<"energy_n[counter] = "<<energy_n[counter]<<endl;
		//cout<<"energy_p[counter] = "<<energy_p[counter]<<endl;


		v_n[counter] = abs(dedk(k_dum,coefficients_cond,kindex_cond,a11)/h_bar*1e-7);
		// group velocity in cm/s
		//cout<<"v_n[counter] = "<<v_n[counter]<<endl;
		v_p[counter] = abs(dedk(k_dum,coefficients_val,kindex_val,b11)/h_bar*1e-7);
		//cout<<"v_p[counter] = "<<v_p[counter]<<endl;
		//getchar();
		
		if (count_orbital!=0)
		{
		    a_n[counter] = admixture_value(k_dum,2);
		    c_n[counter] = admixture_value(k_dum,3);
		}
		else
		{
		    a_n[counter] = 1;
		    c_n[counter] = 0;
		}


		if (count_orbital_p!=0)
		{
		    a_p[counter] = admixture_value_p(k_dum,2);
		    c_p[counter] = admixture_value_p(k_dum,3);
		}
		else
		{
		    a_p[counter] = 0;
		    c_p[counter] = 1;
		}

		//cout<<"a_n[counter] = "<<a_n[counter]<<endl;
		//cout<<"c_n[counter] = "<<c_n[counter]<<endl;
		//getchar();

		//cout<<"a_p[counter] = "<<a_p[counter]<<endl;
		//cout<<"c_p[counter] = "<<c_p[counter]<<endl;
		//getchar();

		if (free_e==0)
		{
		    //if (counter==699)
		    //    Ds_n[counter] = DOS_value1(energy_n[counter],1);   // 1 for n means conduction band
		    //else
		    Ds_n[counter] = DOS_value(energy_n[counter],1);   // 1 for n means conduction band

		    Ds_p[counter] = DOS_value(energy_p[counter],2);   // 2 for p means valence band
		    //cout<<"Ds_n[counter] = "<<Ds_n[counter]<<endl;
		    //cout<<"Ds_p[counter] = "<<Ds_p[counter]<<endl;
		}
	    }

	//-------------------------------------------------saving conduction and valence band --------------------------------
	    
	    //cout<<"saving conduction band"<<endl;
	    if(ispin==1)	  	    	          
	    	fid1 = fopen("conduction_band.dat","w");
	    if(ispin==2 && kk == 0)
	    	fid1 = fopen("conduction_band_up_spin.dat","w");
	    if(ispin==2 && kk == 1)
	    	fid1 = fopen("conduction_band_down_spin.dat","w");
	    
	    fprintf(fid1,"# Sr.No.    k (1/nm)             Energy(eV) \n");
						
	    for (int i = 0; i < points; i++)
		fprintf(fid1,"  %d         %e         %e   \n", i+1, k_grid[i], energy_n[i]);
	fclose(fid1);


	    if(ispin==1)	  	    	          
		    fid1 = fopen("valence_band.dat","w");
	    if(ispin==2 && kk == 0)
		    fid1 = fopen("valence_band_up_spin.dat","w");
	    if(ispin==2 && kk == 1)
		    fid1 = fopen("valence_band_down_spin.dat","w");

	    fprintf(fid1,"# Sr.No.    k (1/nm)            Energy(eV) \n");

	    for (int i = 0; i < points; i++)
		fprintf(fid1,"  %d         %e        %e   \n", i+1, k_grid[i], energy_p[i] );
	fclose(fid1);

	//-------------------------------------------------saving data --------------------------------
	    /*	
	    fid1 = fopen("a_n.txt","w");
	    for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, a_n[i]);
	fclose(fid1);

	    fid1 = fopen("c_n.txt","w");
	    for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, c_n[i]);
	fclose(fid1);

	    fid1 = fopen("Ds_n.txt","w");
	    for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, Ds_n[i]);
	fclose(fid1);

	    fid1 = fopen("Ds_p.txt","w");
	    for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, Ds_p[i]);
	fclose(fid1);
	    */

//--------------------------------------------------------------------------------------------------


            double dum,dum1,dum2;
            dum1 = energy_n[0];
            dum2 = energy_p[0];
            for (int counter=0;counter<points;counter++)
            {
                energy_n[counter] = energy_n[counter] - dum1;
                energy_p[counter] = energy_p[counter] - dum2;

                if (energy_n[counter] < 0)
                    energy_n[counter] = 0;

                if (energy_p[counter] < 0)
                    energy_p[counter] = 0;

                if (a_n[counter] < 0)
                    a_n[counter] = 1e-6;
                if (c_n[counter] < 0)
                    c_n[counter] = 1e-6;

                dum = pow(a_n[counter],2) + pow(c_n[counter],2);
                a_n[counter] = a_n[counter]/(pow(dum,0.5));
                c_n[counter] = c_n[counter]/(pow(dum,0.5));

                if (Ds_n[counter] < 0)
                    Ds_n[counter] = 1e-10;
                if (Ds_p[counter] < 0)
                    Ds_p[counter] = 1e-10;
            }
//-----------------------------------------saving data ----------------------------------------------------------------
            /*
            fid1 = fopen("energy_n.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, energy_n[i]);
	fclose(fid1);

            fid1 = fopen("energy_p.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, energy_p[i]);
	fclose(fid1);

            fid1 = fopen("v_n.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, v_n[i]);
	fclose(fid1);


            fid1 = fopen("a_n.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, a_n[i]);
	fclose(fid1);

            fid1 = fopen("c_n.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, c_n[i]);
	fclose(fid1);

            fid1 = fopen("Ds_n.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, Ds_n[i]);
	fclose(fid1);

            fid1 = fopen("Ds_p.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, Ds_p[i]);
	fclose(fid1);

            */
            /*
            fid1 = fopen("a_n.txt","r");
            for (int i = 0; i < points; i++)
            {
                fgets(line, 1000, fid1);
                sscanf(line, "%lf", &a_n[i]);
            }
	fclose(fid1);

            fid1 = fopen("c_n.txt","r");
            for (int i = 0; i < points; i++)
            {
                fgets(line, 1000, fid1);
                sscanf(line, "%lf", &c_n[i]);
            }
	fclose(fid1);
            */
//-----------------------------------------------------------------------------------------------------------------------

//-------------------------------------For debugging -------------------------------------------------------------------------
            /*
            for (int counter=0;counter<points;counter++)
            {
                cout<<"i+1 = "<<counter+1<<"   k_grid[counter] =   "<<k_grid[counter]<<endl;
                cout<<"i+1 = "<<counter+1<<"   energy_n[counter] =   "<<energy_n[counter]<<endl;
                cout<<"i+1 = "<<counter+1<<"   energy_p[counter] =   "<<energy_p[counter]<<endl;
                //cout<<"i+1 = "<<counter+1<<"   Ds_n[counter] =   "<<Ds_n[counter]<<endl;
                //cout<<"i+1 = "<<counter+1<<"   Ds_p[counter] =   "<<Ds_p[counter]<<endl;
                if (counter==99 || counter == 199 || counter == 299 || counter == 399 || counter == 499 || counter == 599)
                    getchar();
            }
            getchar();
	    */
}



