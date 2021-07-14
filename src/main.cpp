# include "main.h"
int a[2],b[2];
double h_bar,CBM,VBM;
int count1,count2;
int count_orbital, count_orbital_p;

double mobility_all[10]={0};
double calc_mobility[30][2] = {0};
double calc_mobility_rta[30][2] = {0};
double calc_thermopower[30][2] = {0};
double calc_sigma[30][2] = {0};
double calc_sigma_rta[30][2] = {0};

double calc_mobility_pe[30][1] = {0};
double calc_mobility_de[30][1] = {0};
double calc_mobility_dis[30][1] = {0};
double calc_mobility_ii[30][1] = {0};
double calc_mobility_po[30][1] = {0};
double calc_mobility_to[30][1] = {0};
double calc_mobility_alloy[30][1] = {0};
double calc_mobility_iv[30][1] = {0};
double calc_mobility_neutral[30][1] = {0};

double mobility_hall_all[10]={0};
double calc_mobility_hall[30][2] = {0};
double calc_mobility_hall_rta[30][2] = {0};
double calc_sigma_hall[30][2] = {0};
double calc_sigma_hall_rta[30][2] = {0};
double hall_factor[30][2] = {0};
double hall_factor_rta[30][2] = {0};

double calc_mobility_hall_pe[30][1] = {0};
double calc_mobility_hall_de[30][1] = {0};
double calc_mobility_hall_dis[30][1] = {0};
double calc_mobility_hall_ii[30][1] = {0};
double calc_mobility_hall_po[30][1] = {0};
double calc_mobility_hall_to[30][1] = {0};
double calc_mobility_hall_alloy[30][1] = {0};
double calc_mobility_hall_iv[30][1] = {0};
double calc_mobility_hall_neutral[30][1] = {0};


double N_im_de, Nd_plus,N_im_modified;

int kk, ispin;

bool FileExists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main()
{
		
    cout.precision(15);
    FILE *fid1;
    FILE *fid2;
	char line[1000];
	h_bar =  h_planck / (2 * pi);    // Planck constant divided by 2pi[eV.s]
    //cout<<" hbar = "<<h_bar<<endl;
    copyright();	
    read_input();	
	
	N_cb = 1;
	N_vb = 1;

    cout<<endl;

    cout<<"Inside main function"<<endl;	
    //cout<<"iterations = "<<iterations<<endl;

    if (scattering_mechanisms[1] == 0)
    {
        if (iterations != 1 )
        {
            iterations = 1;
            cout<<"Since polar optical phonon scattering is not used, so iterations is set to 1"<<endl;
        }
        T_trans = 0;
    }

    if (De_ionization == 1)
    {	// neutral impurity scattering
        if (scattering_mechanisms[9] == 0)
        {
            scattering_mechanisms[9] = 1;
            cout<<"De-ionization is to be included so neutral impurity scattering is to included in simulation";
        }
    }
	
	
	
    if (free_e == 0 )
        cout<<"free_e = false"<<endl;
    else
        cout<<"free_e = true"<<endl;


    printf("\n C_long  =   %e dyne/cm^2 \n", C_long);
    printf("\n C_trans =   %e dyne/cm^2 \n", C_trans);
    printf("\n c_bar   =   %e dyne/cm^2 \n", c_bar);
    //printf("\n P_piezo =   %e \n",P_piezo);

    read_OUTCAR();

    //cout<<"number = "<<number<<endl;

    if (number==1)
    {
        printf("\nion_mass =  %lf \n", ion_mass1[0] );   // unit amu
        printf("\nion_numbers =  %d \n", ion_numbers1[0] );
    }
    else if (number==2)
    {
        printf("\nion_mass =  %lf %lf \n", ion_mass1[0], ion_mass1[1]);   // unit amu
        printf("\nion_numbers =  %d %d\n", ion_numbers1[0], ion_numbers1[1]);
    }
    else if (number==3)
    {
        printf("\nion_mass =  %lf %lf %lf \n", ion_mass1[0], ion_mass1[1],ion_mass1[2]);   // unit amu
        printf("\nion_numbers =  %d %d %d \n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2]);
    }
    else if (number==4)
    {
        printf("\nion_mass =  %lf %lf %lf %lf\n", ion_mass1[0], ion_mass1[1],ion_mass1[2], ion_mass1[3]);   // unit amu
        printf("\nion_numbers =  %d %d %d %d\n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2], ion_numbers1[3]);
    }
    else
    {
        printf("\nion_mass =  %lf %lf %lf %lf\n", ion_mass1[0], ion_mass1[1],ion_mass1[2], ion_mass1[3]);     // unit amu
        printf("\nion_numbers =  %d %d %d %d\n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2], ion_numbers1[3]);
    }
	printf("\nvolume =  %lf A^(-3)\n", volume1);
	// unit Angestron^3  volume = a^3/4   (a -> lattice constant)

    //cout<<"spin_orbit_coupling = "<<spin_orbit_coupling<<endl;
    if(spin_orbit_coupling == 0)
        cout<<endl<<"spin_orbit_coupling = false"<<endl;
    if(spin_orbit_coupling == 1)
        cout<<endl<<"spin_orbit_coupling = true"<<endl;

	cout<<"Lattice matrix : "<<endl;

	for (int i = 0; i < 3; i++)
        {
		for (int j = 3; j < 6; j++)
        	{
			printf(" %lf", 10 * lm[i][j]);
		}
		printf("\n");
	}

	
	if (rho==0)
    	{
		double sum=0;
		for (int i = 0; i < number; i++)
		{
		    sum += (ion_mass1[i] * ion_numbers1[i]);
			//cout<<"sum = "<<sum<<endl;
		}
		rho = (sum*1.67377e-27) / (volume1*1e-30);    // unit kg/m^3
	    	cout<<endl<<"Calculated density = "<<rho/1000<<" g/cm^3"<<endl;       
	}
	else	
	{
		rho = rho*1000;   // to convert g/cm^3 to kg/cm^3
				  // user has given density into g/cm^3 in calculation kg/m^3 is used	
	}
	if (scattering_mechanisms[6]==1)   // dislocation scattering 
	{
		cout<<endl<<"c_lattice constant = "<<c_lattice<<" nm"<<endl;       		
	}	
	//	
	//return 0;

	//getchar();

    double kcbm[3],kvbm[3];

    read_ispin();
    	
    int jj;	
    cout<<"ispin = "<<ispin<<endl;
    //cout<<"Next "<<endl;
    if (ispin == 1)    	
    {
	cout<<"Material is non magnetic"<<endl;
	jj = 0;
    }
    else
    {
    	cout<<"Material is magnetic"<<endl;
	jj = 1;
    }


//----------------------------------------------generating output file -----------------------------------------------------------------

    if (variation==0)   // temperature variation
    {
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","w");
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		fid1 = fopen("thermopower.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);
		
		if(Bfield!=0)
		{
			fid1 = fopen("mobility_hall.dat","w");
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}
	
    	if (ispin == 2 )
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","w");
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		fid1 = fopen("thermopower_up_spin.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);

		fid1 = fopen("mobility_down_spin.dat","w");
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		fid1 = fopen("thermopower_down_spin.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);
		
		fid1 = fopen("mobility_hall_up_spin.dat","w");
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		if(Bfield!=0)
		{
			fid1 = fopen("conductivity_hall_up_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);


			fid1 = fopen("mobility_hall_down_spin.dat","w");
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			fclose(fid1);

			fid1 = fopen("conductivity_hall_down_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}		
	}

    }
    else          // Doping variation
    {
	if (ispin == 1)
	{
		FILE *fid1;

		fid1 = fopen("mobility.dat","w");
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		fid1 = fopen("thermopower.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_hall.dat","w");
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","w");
			fprintf(fid1,"# Doping(cm^-3)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}


	if (ispin == 2)
	{
		FILE *fid1;

		fid1 = fopen("mobility_up_spin.dat","w");
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		fid1 = fopen("thermopower_up_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);

		fid1 = fopen("mobility_down_spin.dat","w");
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		fid1 = fopen("thermopower_down_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_hall_up_spin.dat","w");
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			fclose(fid1);

			fid1 = fopen("conductivity_hall_up_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);
			
			fid1 = fopen("hall_factor_up_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);

			
			fid1 = fopen("mobility_hall_down_spin.dat","w");
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			fclose(fid1);

			fid1 = fopen("conductivity_hall_down_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)     Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}

    }

//---------------------------------------------saved output file ---------------------------------------------------------

	
    for (kk=0; kk<= jj ; kk++)
    {

    if (ispin == 2 && kk==0)    	
	cout<<"-------------------- For up spin ----------------------------------------- "<<endl;

    if (ispin == 2 && kk==1)
    	cout<<"-------------------- For down spin ---------------------------------------"<<endl;

	
    find_cbm_vbm(1,spin_orbit_coupling);   // 1 means for CB and 2 means for VB
    	
    kcbm[0]=kcbm1[0];
    kcbm[1]=kcbm1[1];
    kcbm[2]=kcbm1[2];

    find_cbm_vbm(2,spin_orbit_coupling);   // 2 means for VB
        	
    kvbm[0]=kvbm1[0];
    kvbm[1]=kvbm1[1];
    kvbm[2]=kvbm1[2];

    cout<<endl<<"kcbm = "<<"   "<<kcbm[0]<<"   "<<kcbm[1]<<"   "<<kcbm[2]<<endl;
    cout<<endl<<"kvbm = "<<"   "<<kvbm[0]<<"   "<<kvbm[1]<<"   "<<kvbm[2]<<endl;
    
    		

//------------------------------------------- Reading Procar file not for conduction band dependent on ispin----------------------
    //cout<<"fitting_2 = "<<fitting_2<<endl;
    if (fitting_2 == 0)
    {
        count_orbital = decompose_wave();
        //cout<<"count_orbital = "<<count_orbital<<endl;
    }
    else     // this part is for debugging code only
    {
        int i;
        ifstream cin("orbital_decomposedd.txt");

        i=0;
        while(!cin.eof())
        {
            cin>>orbital_decomposedd[i][0]>>orbital_decomposedd[i][1]>>orbital_decomposedd[i][2]>>orbital_decomposedd[i][3];
            /*

            cout<<orbital_decomposedd[i][0]<<"    "<<orbital_decomposedd[i][1]<<"    "<<
            orbital_decomposedd[i][2]<<"    "<<orbital_decomposedd[i][3]<<endl;
            getchar();
            */
            i++;
        }
        count_orbital = i;
    }

    /*
	if(fitting_2 == 0)
	{
	    fid1 = fopen("orbital_decomposedd.txt","w");

	    for (int i = 0; i < count_orbital; i++)
		fprintf(fid1,"%d    %e    %e    %e    %e\n", i+1, orbital_decomposedd[i][0],
		orbital_decomposedd[i][1], orbital_decomposedd[i][2],orbital_decomposedd[i][3]);
	    fclose(fid1);

	}
    */
//------------------------------------------- Reading Procar file not for valence band dependent on ispin----------------------
	
    //cout<<"fitting_2 = "<<fitting_2<<endl;
    if (fitting_2 == 0)
    {
        count_orbital_p = decompose_wave_p();
        //cout<<"count_orbital_p = "<<count_orbital_p<<endl;
    }
    else     // this part is for debugging code only
    {
        int i;
        ifstream cin("orbital_decomposedd_p.txt");

        i=0;
        while(!cin.eof())
        {
            cin>>orbital_decomposedd_p[i][0]>>orbital_decomposedd_p[i][1]>>orbital_decomposedd_p[i][2]>>orbital_decomposedd_p[i][3];
            /*
            cout<<orbital_decomposedd_p[i][0]<<"    "<<orbital_decomposedd_p[i][1]<<"    "<<
            orbital_decomposedd_p[i][2]<<"    "<<orbital_decomposedd_p[i][3]<<endl;
            getchar();
            */
            i++;
        }
        count_orbital_p = i;
    }


//------------------------------------------------------------------------------------------------------------------------------------


    make_band(1); // for conduction band
    count1 = countx;   // length of conduction band
    //cout<<"count1 = "<<count1<<endl;


    make_band(2); // for valence band
    count2 = countx;   // valence band elements number
    //cout<<"count2 = "<<count2<<endl;
    /*
    for (int i=0;i<count1;i++)
    {
        for(int j=0;j<2;j++)
            cout<<"i =  "<<i<<"  "<<cond_band[i][j]<<"   ";
        cout<<endl;
        //if (i>2500)
        getchar();
    }
    //*/
    /*
    for (int i=0;i<count2;i++)
    {
        for(int j=0;j<2;j++)
            cout<<"i =  "<<i<<"  "<<val_band[i][j]<<"   ";
        cout<<endl;
        //if (i>2500)
        //getchar();
    }

    getchar();
    getchar();
    */
    
    cout<<"In main.cpp"<<endl;	
    cout<<endl<<"CBM = "<<CBM<<" eV"<<endl;
    cout<<endl<<"VBM = "<<VBM<<" eV"<<endl;
    if (Bgap[0]==0)
    {
        int l = (sizeof(Bgap)/sizeof(*Bgap));
        for (int i=0;i<l;i++)
            Bgap[i] = abs(CBM - VBM);
        cout<<endl<<"Bgap = "<<Bgap[0]<<" eV"<<endl;
    }


    if (fitting_1==0)
    {
	//cout<<"Reached here"<<endl; getchar();
        cout<<endl<<"Analytical fitting of conduction band"<<endl;
        int *c = analytical_fitting(cond_band,count1,1);    // coloumn - a[0]+1   row - a[1]+1    of coefficients
        a[0] = c[0];
        a[1] = c[1];

        //cout<<"a[0] = "<<a[0]<<endl;
        //cout<<"a[1] = "<<a[1]<<endl;
        //getchar();
        // a[0] --contain mamimum degree   coloumn - a[0]+1
        // a[1] = contain length_dd        row - a[1]+1

//--------------------------- save coefficient and kindex -----------------------------
	/*
        fid1 = fopen("coefficients_conduction_band.txt","w");
        for (int i = 0; i <a[1]+1 ; i++)
        {
            for (int j = 0; j<a[0]+1 ; j++)
                fprintf(fid1,"%e    ", coefficients_cond[i][j]);
            fprintf(fid1,"\n");
        }
	fclose(fid1);
	
        fid1 = fopen("kindex_conduction_band.txt","w");
        for (int i = 0; i <a[1] ; i++)
        {
            fprintf(fid1,"%e \n", kindex_cond[i]);

        }
	fclose(fid1);
	*/
//--------------------------- ---------------------------- -----------------------------
    }
    else    // this part is for debugging code only
    {
        fid1 = fopen("coefficients_conduction_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"coefficients_conduction_band.txt `f\ECle is not present exit from program"<<endl;
            return 0;
        }
        //a[2] and b[2] contains size of coefficient of conduction and valence band
        a[0] = 6;  // a[0]+1 coloumn of conduction band
        a[1] = 4;  // a[1]+1 row of conduction band   and a[1] number of elements in k_index

        for (int i = 0; i < 5; i++)
        {
            fgets(line, 1000, fid1);

                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &coefficients_cond[i][0], &coefficients_cond[i][1], &coefficients_cond[i][2],
                       &coefficients_cond[i][3],&coefficients_cond[i][4],&coefficients_cond[i][5], &coefficients_cond[i][6]);
        }
	fclose(fid1);

        fid1 = fopen("kindex_conduction_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"kindex_conduction_band.txt `f\ECle is not present exit from program"<<endl;
            return 0;
        }
        for (int i = 0; i < 4; i++)
        {
            fgets(line, 1000, fid1);

                sscanf(line, "%lf", &kindex_cond[i]);
        }
	fclose(fid1);

    }

    if (fitting_1==0)
    {
        cout<<endl<<"Analytical fitting of valence band"<<endl;
        int *c = analytical_fitting(val_band,count2,2);
        cout<<endl;
        b[0] = c[0];
        b[1] = c[1];

        //cout<<"b[0] = "<<b[0]<<endl;
        //cout<<"b[1] = "<<b[1]<<endl;
        //getchar();
        // b[0] --contain mamimum degree   coloumn - b[0]+1
        // b[1] = contain length_dd         row - b[1]+1

//--------------------------- save coefficient and kindex -----------------------------
	/*
        fid1 = fopen("coefficients_valence_band.txt","w");
        for (int i = 0; i <b[1]+1 ; i++)
        {
            for (int j = 0; j<b[0]+1 ; j++)
                fprintf(fid1,"%e    ", coefficients_val[i][j]);
            fprintf(fid1,"\n");
        }
	fclose(fid1);

        fid1 = fopen("kindex_valence_band.txt","w");
        for (int i = 0; i <b[1] ; i++)
        {
            fprintf(fid1,"%e \n", kindex_val[i]);
        }
	fclose(fid1);
	*/
//----------------------------------------------------------- -----------------------------
    }
    else                 // this part is for debugging code only
    {
        b[0] = 6;  // b[0]+1 coloumn of conduction band
        b[1] = 4;  // b[1]+1 row of conduction band  and b[1] number of elements in k_index
        fid1 = fopen("coefficients_valence_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"coefficients_valence_band.txt `f\ECle is not present exit from program"<<endl;
            return 0;
        }
        for (int i = 0; i < 5; i++)
        {
            fgets(line, 1000, fid1);

                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &coefficients_val[i][0], &coefficients_val[i][1], &coefficients_val[i][2],
                       &coefficients_val[i][3],&coefficients_val[i][4],&coefficients_val[i][5], &coefficients_val[i][6]);
        }
	fclose(fid1);

        fid1 = fopen("kindex_valence_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"kindex_valence__band.txt `f\ECle is not present exit from program"<<endl;
            return 0;
        }

        for (int i = 0; i < 4; i++)
        {
            fgets(line, 1000, fid1);
            sscanf(line, "%lf", &kindex_val[i]);
        }
	fclose(fid1);
    }
//-------------------------------------------Save coefficient result in txt file -----------------------

//------------------------------- print coefficient -and kindex --------------------------
    cout<<"coefficient of conduction band"<<endl;

    for (int i = 0; i < a[1]+1; i++)
	{
		for (int j = 0; j < a[0]+1; j++)
			cout<<coefficients_cond[i][j]<<"     " ;
        cout<<endl;
	}
	cout<<endl;
    cout<<"kindex of conduction band"<<endl;
	for (int i = 0; i < a[1]; i++)
    {
			cout<<kindex_cond[i]<<endl;
    }
    cout<<endl;


    cout<<"coefficient of valence band"<<endl;
    for (int i = 0; i < b[0]+1; i++)
	{
		for (int j = 0; j < b[1]+1; j++)
			cout<<coefficients_val[i][j]<<"     ";
        cout<<endl;
    }
    cout<<endl;

    cout<<"kindex of valence band"<<endl;
	for (int i = 0; i < b[1]; i++)
        cout<<kindex_val[i]<<endl;
    cout<<endl;

//---------------------------------------------------------DOSCAR --------------------------------------------------
    if (free_e==0)
    {
        n_DOS();
    }

    // Electron effective mass calculation 
    double qq = coefficients_cond[0][a[0]-2]*2*1e-18/(pow(h_bar,2)*e);
    m = pow(qq,-1)/m_e;
    //cout<<"m = "<<m<<endl;
    
    // hole effective mass calculation
    qq = coefficients_val[0][b[0]-2]*2*1e-18/(pow(h_bar,2)*e);
    m_h = pow(qq,-1)/m_e;
    
    cout<<"Average electron effective mass, m*_e/m =  "<<m<<endl;
    cout<<"Average hole effective mass, m*_h/m =  "<<m_h<<endl;


    int count_d=0,count_t=0;
    int ss=0;
    while(n_array[ss]!=0)
        ss++;

    count_d = ss;
	
	ss=0;

    while(T_array[ss]!=0)
        ss++;

    count_t = ss;

    int cc = -1;
    
    if(count_d == 0)
    	count_d = count_d+1;
    	
    //cout<<"Count_d = "<<count_d<<endl;
    //cout<<"Count_t = "<<count_t<<endl;
    
    
    if (variation==1)
    {
        for (int i=0;i<count_d;i++)
        {
            calc_mobility[i][0] = n_array[i];
            calc_mobility_rta[i][0] = n_array[i];
            calc_thermopower[i][0] = n_array[i];
            calc_sigma[i][0] = n_array[i];
            calc_sigma_rta[i][0] = n_array[i];

            calc_mobility_ii[i][0] = n_array[i];
            calc_mobility_po[i][0] = n_array[i];
            calc_mobility_de[i][0] = n_array[i];
            calc_mobility_pe[i][0] = n_array[i];
            calc_mobility_dis[i][0] = n_array[i];
            calc_mobility_to[i][0] = n_array[i];
            calc_mobility_alloy[i][0] = n_array[i];
            calc_mobility_iv[i][0] = n_array[i];
            calc_mobility_neutral[i][0] = n_array[i];


            calc_mobility_hall[i][0] = n_array[i];
            calc_mobility_hall_rta[i][0] = n_array[i];
            calc_sigma_hall[i][0] = n_array[i];
            calc_sigma_hall_rta[i][0] = n_array[i];

            calc_mobility_hall_ii[i][0] = n_array[i];
            calc_mobility_hall_po[i][0] = n_array[i];
            calc_mobility_hall_de[i][0] = n_array[i];
            calc_mobility_hall_pe[i][0] = n_array[i];
            calc_mobility_hall_dis[i][0] = n_array[i];
            calc_mobility_hall_to[i][0] = n_array[i];
            calc_mobility_hall_alloy[i][0] = n_array[i];
            calc_mobility_hall_iv[i][0] = n_array[i];
            calc_mobility_hall_neutral[i][0] = n_array[i];

            hall_factor[i][0] = n_array[i];
            hall_factor_rta[i][0] = n_array[i];
        }
    }
    else
    {
        for (int i=0;i<count_t;i++)
        {
            calc_mobility[i][0] = T_array[i];
            calc_mobility_rta[i][0] = T_array[i];
            calc_thermopower[i][0] = T_array[i];
            calc_sigma[i][0] = T_array[i];
            calc_sigma_rta[i][0] = T_array[i];

            calc_mobility_ii[i][0] = T_array[i];
            calc_mobility_po[i][0] = T_array[i];
            calc_mobility_de[i][0] = T_array[i];
            calc_mobility_pe[i][0] = T_array[i];
            calc_mobility_dis[i][0] = T_array[i];
            calc_mobility_to[i][0] = T_array[i];
            calc_mobility_alloy[i][0] = T_array[i];
            calc_mobility_iv[i][0] = T_array[i];
            calc_mobility_neutral[i][0] = T_array[i];


            calc_mobility_hall[i][0] = T_array[i];
            calc_mobility_hall_rta[i][0] = T_array[i];
            calc_sigma_hall[i][0] = T_array[i];
            calc_sigma_hall_rta[i][0] = T_array[i];

            calc_mobility_hall_ii[i][0] = T_array[i];
            calc_mobility_hall_po[i][0] = T_array[i];
            calc_mobility_hall_de[i][0] = T_array[i];
            calc_mobility_hall_pe[i][0] = T_array[i];
            calc_mobility_hall_dis[i][0] = T_array[i];
            calc_mobility_hall_to[i][0] = T_array[i];
            calc_mobility_hall_alloy[i][0] = T_array[i];
            calc_mobility_hall_iv[i][0] = T_array[i];
            calc_mobility_hall_neutral[i][0] = T_array[i];
            hall_factor[i][0] = T_array[i];
            hall_factor_rta[i][0] = T_array[i];

        }
    }

    double mobility_ii = 1e10;
    double mobility_po = 1e10;
    double mobility_to = 1e10;
    double mobility_npop = 1e10;
    double mobility_de = 1e10;
    double mobility_pe = 1e10;
    double mobility_dis = 1e10;
    double mobility_alloy =1e10;
    double mobility_iv =1e10;
    double mobility_neutral = 1e10;

    double mobility_avg = 1e10;
    double mobility = 1e10;



    double mobility_hall_ii = 1e10;
    double mobility_hall_po = 1e10;
    double mobility_hall_to = 1e10;
    double mobility_hall_npop = 1e10;
    double mobility_hall_de = 1e10;
    double mobility_hall_pe = 1e10;
    double mobility_hall_dis = 1e10;
    double mobility_hall_alloy =1e10;
    double mobility_hall_iv =1e10;
    double mobility_hall_neutral = 1e10;

    double mobility_hall_avg = 1e10;
    double mobility_hall = 1e10;
    double mobility_hall_rta = 1e10;
    double hall_factor1;
    double hall_factor_rta1;

    double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;


    double n0,Nd1,Na1;
    double T,epsilon_s1,epsilon_inf1,P_piezo1;
    double k_min, k_trans, k_step_fine, k_step;
    int points, points1, points2;


    for(int ii=0;ii<count_d;ii++)  // For doping variation
    {
        n0 = n_array[ii];
        Nd1 = Nd[ii];
        Na1 = Na[ii];
	cout<<"Nd1 = "<<Nd1<<endl;
	cout<<"Na1 = "<<Na1<<endl;

        for(int T_loop=0;T_loop<count_t;T_loop++)
        {
             T = T_array[T_loop];

             printf("\nFor Temperature = %f ",T);
             printf("\nFor Doping = %e \n",n_array[ii]);
             epsilon_s1 = epsilon_s[T_loop];
             epsilon_inf1 = epsilon_inf[T_loop];
             P_piezo1 = P_piezo[T_loop];

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

            double k_grid[points]={0};

            double v_n[points]={0};	// Group velocity of electrons over k_grid
            double v_p[points]={0};	// Group velocity of electrons over k_grid

            double energy_ref[points]={0};
            double energy_n[points]={0};
            double energy_p[points]={0};
            double a_n[points]={0};
            double c_n[points]={0};

            double a_p[points]={0};
            double c_p[points]={0};

            double Ds[points]={0};
            double Ds_n[points]={0};
            double Ds_p[points]={0};

            double B_ii = 0;
            double D_ii = 0;

            double electric_driving_force_n[points]={0};
            double thermal_driving_force_n[points]={0};

            double A_ii = 0;

            double k_dum;
            for (int counter=0;counter<points;counter++)
            {
                if ((k_min+(counter)*k_step_fine)<k_trans)
                    k_dum = k_min+(counter)*k_step_fine;
                else
                    k_dum = k_trans+(counter-points1+1)*k_step;

                k_grid[counter] = k_dum;
                //cout<<"k_grid[counter] = "<<k_grid[counter]<<endl;

                energy_n[counter] = conduction_dispersion(k_dum, coefficients_cond, kindex_cond, a);
                energy_p[counter] = conduction_dispersion(k_dum,coefficients_val, kindex_val, b);

                //cout<<"energy_n[counter] = "<<energy_n[counter]<<endl;
                //cout<<"energy_p[counter] = "<<energy_p[counter]<<endl;

                //if (counter==99 || counter == 199 || counter == 299 || counter == 399 || counter == 499 || counter == 599)
                //    getchar();

                v_n[counter] = abs(dedk(k_dum,coefficients_cond,kindex_cond,a)/h_bar*1e-7);
                // group velocity in cm/s
                //cout<<"v_n[counter] = "<<v_n[counter]<<endl;
		
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
//---------------------------------------------------------------------------------------------------------------

            if (energy_n[points-1]<6*k_B*T)
                cout<<"The selected value of k_max is too low"<<endl;
            //cout<<"energy_n[points-1] = "<<energy_n[points-1]<<endl;
                
            if (n_array[ii]>5e20)
                n0 = 5e20;
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
            /*
            FILE *fid1;
            fid1 = fopen("energy_n.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d   %e    %e \n", counter+1, k_grid[counter], energy_n[counter]);
	fclose(fid1);

            fid1 = fopen("energy_p.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d   %e    %e \n", counter+1, k_grid[counter], energy_p[counter]);
	fclose(fid1);

            fid1 = fopen("Ds_n.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d    %e    %e \n", counter+1, k_grid[counter], Ds_n[counter]);
	fclose(fid1);

            fid1 = fopen("Ds_p.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d    %e    %e \n", counter+1, k_grid[counter], Ds_p[counter]);
	fclose(fid1);

            fid1 = fopen("a_n.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d    %e    %e \n", counter+1, k_grid[counter], a_n[counter]);
	fclose(fid1);

            fid1 = fopen("c_n.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d    %e    %e \n", counter+1, k_grid[counter], c_n[counter]);
	fclose(fid1);

            fid1 = fopen("v_n.txt","w");
            for (int counter=0;counter<points;counter++)
                fprintf(fid1,"%d    %e    %e \n", counter+1, k_grid[counter], v_n[counter]);
	fclose(fid1);
            */
	
	//cout<<"Reached here"<<endl;
		
            if (De_ionization==0)
            {
                find_fermi(k_grid, n0, T, T_loop, energy_n, energy_p, Ds_n, Ds_p, points);

                cout<<"E_F = "<<E_F<<" eV"<<endl;
                cout<<"n_e = "<<n_e<<" cm^-3"<<endl;
                cout<<"n_h = "<<n_h<<" cm^-3"<<endl;
		
		//cout<<"n0 = "<<n0<<endl;
			
			
                N_im_de = 0;
                Nd_plus = Nd1;

		//cout<<"(n_e-n_h)/n0 = "<<(n_e-n_h)/n0<<endl;
		
                if ((n_e-n_h)/n0 < 0.5 && (n_e-n_h)/n0 > 1.5)
                {
                    cout<<"Calculated concentration is not within 50% of given concentrtion. Programs end here"<<endl;
                    return 0;
                }

                if ((n_e-n_h)/n0 < 0.6 && (n_e-n_h)/n0 > 1.4)
                {
                    cout<<"Calculated concentration is not within 60% of given concentrtion. Programs end here"<<endl;
                    return 0;
                }

                if ((n_e-n_h)/n0 < 0.8 && (n_e-n_h)/n0 > 1.2)
                {
                    cout<<"Calculated concentration is not within 80% of given concentrtion. Programs end here"<<endl;
                    return 0;
                }

                if ((n_e-n_h)/n0 < 0.9 && (n_e-n_h)/n0 > 1.1)
                {
                    cout<<"Calculated concentration is not within 90% of given concentrtion. Programs end here"<<endl;
                    return 0;
                }

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
	    
				
            double N_ii = (n_e + n_h) + a1;
	    cout<<"Net ionized donors concentration =  "<<N_ii<<"  cm^(-3)"<<endl;

            N_im_modified = N_im[ii] + N_im_de;    // net neutral impurity
            cout<<"Net neutral impurity =  "<<N_im_modified<<" cm^(-3)"<<endl;

            double efef_n = 0; // effective ef (fermi energy)
            double efef_p = 0; // effective ef (fermi energy)

            efef_n = E_F;
            efef_p = -(E_F+Bgap[T_loop]);

            //E_F =  0.002111407852173;

            double beta_constant = beta(T, Ds_n, energy_n, k_grid, v_n, points, T_loop);
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
                                        *f0(energy_n[counter],efef_n,T)*(1-f0(energy_n[counter],efef_n,T))*energy_n[counter]/(k_B*T);
                    // Part of equation (54) of Rode's book
                    integral_denominator_n = integral_denominator_n + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
                    *f0(energy_n[counter],efef_n,T)*(1-f0(energy_n[counter],efef_n,T));
                    // Part of equation (54) of Rode's book
                }
            }

            double df0dz_integral_n = integral_numerator_n/integral_denominator_n;
            //cout<<"df0dz_integral_n = "<<df0dz_integral_n<<endl;

            double N_poph_atT = N_poph(omega_LO,T);

            //cout<<"N_poph_atT = "<<N_poph_atT<<endl;
            double N_e[iv_number]={0};

            if (scattering_mechanisms[8]==1)   // intravalley scattering
            {
                for (int aa=0;aa<iv_number;aa++)
                    N_e[aa] = N_poph(we[aa],T);
            }

//---------------------------- Tools for BTE ------------------------------------------------------------------

            double df0dk_grid[points]={0};
            double f0x1_f0[points]={0};
            double thermal_driving_force[points]={0};

            double kplus_grid[points]={0};
            double kminus_grid[points]={0};

            double betaplus_grid[points]={0};
            double betaminus_grid[points]={0};

            double Aminus_grid[points]={0};
            double Aplus_grid[points]={0};

            double lambda_i_plus_grid[points]={0};
            double lambda_o_plus_grid[points]={0};
            double lambda_i_minus_grid[points]={0};
            double lambda_o_minus_grid[points]={0};

            double lambda_e_plus_grid[points][iv_number]={0};
            double lambda_e_minus_grid[points][iv_number]={0};

            double nu_deformation[points]={0};
            double nu_piezoelectric[points]={0};
            double nu_ionizedimpurity[points]={0};
            double nu_dislocation[points]={0};
            double nu_alloy[points]={0};
            double nu_neutralimpurity[points]={0};
	    double nu_iv[points][iv_number]={0};
	    double nu_iv_total[points]={0};	
            double nu_el[points]={0};

            double electric_driving_force[points]={0};
            double f_dist[points]={0};

            for (int counter = 0;counter<points;counter++)
            {
                k_dum = k_grid[counter];
                // unit 1/nm
                //cout<<"counter+1 = "<<counter+1<<endl;
                //cout<<"k_dum = "<<k_dum<<endl;
                
		 // polar optical phonon scattering 
                if (scattering_mechanisms[1]==1)
                {
                    kplus_grid[counter] = kplus(counter,k_grid,omega_LO,energy_n,points);
                    kminus_grid[counter] = kminus(counter,k_grid,omega_LO,energy_n,points);
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

                    //B_ii = (4*k_dum^2/beta_constant^2)/(1+4*k_dum^2/beta_constant^2)+8*(beta_constant^2+2*k_dum^2)/(beta_constant^2+4*k_dum^2)*(c(counter))^2
                     //       +(3*beta_constant^4+6*beta_constant^2*k_dum^2-8*k_dum^4)/((beta_constant^2+4*k_dum^2)*k_dum^2)*(c(counter))^4;  % According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)


                    D_ii = 1+(2*pow(beta_constant,2)*(pow(c_n[counter],2)/pow(k_dum,2))+
                              (3*pow(beta_constant,4)*(pow(c_n[counter],4))/(4*pow(k_dum,4))));
                              // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)
                //D_ii = 1+(2*beta_constant^2*(c(counter))^2/k_dum^2)+
                //(3*beta_constant^4*(c(counter))^4/(4*k_dum^4));  % According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)

                    nu_ionizedimpurity[counter] = nu_ii(k_dum,counter,B_ii,D_ii,beta_constant,N_ii,v_n[counter],epsilon_s1);
		     // unit 1/second	
                    //cout<<"B_ii = "<<B_ii<<endl;
                    //cout<<"D_ii = "<<D_ii<<endl;
                    //cout<<"N_ii = "<<N_ii<<endl;
                    //cout<<"nu_ionizedimpurity[counter] = "<<nu_ionizedimpurity[counter]<<endl;
                }

		// POP scattering
                if (scattering_mechanisms[1]==1)
                {
                    Aminus_grid[counter] = Aminus(counter,k_grid,omega_LO,a_n,c_n,energy_n,points);
                    Aplus_grid[counter] = Aplus(counter,k_grid,omega_LO,a_n,c_n,energy_n,points);


                    betaplus_grid[counter] = betaplus(counter,k_grid,omega_LO,epsilon_s1,epsilon_inf1,energy_n,v_n,points);
                    betaminus_grid[counter] = betaminus(counter,k_grid,omega_LO,epsilon_s1,epsilon_inf1,energy_n,v_n,points);

                    lambda_i_plus_grid[counter] = abs(lambda_i_plus(counter,k_grid,omega_LO,Aplus_grid[counter],c_n,
                                                                epsilon_s1,epsilon_inf1,energy_n,v_n,points));

                    lambda_i_minus_grid[counter] = abs(lambda_i_minus(counter,k_grid,omega_LO,Aminus_grid[counter],c_n,
                    epsilon_s1,epsilon_inf1,energy_n,v_n,points));

                    lambda_o_plus_grid[counter] = abs(lambda_o_plus(counter,k_grid,omega_LO,Aplus_grid[counter],a_n,c_n,
                    epsilon_s1,epsilon_inf1,energy_n,v_n,points));

                    lambda_o_minus_grid[counter] = abs(lambda_o_minus(counter,k_grid,omega_LO,Aminus_grid[counter],a_n,c_n,
                    epsilon_s1,epsilon_inf1,energy_n,v_n,points));
                }
                
		// acoustic scattering 
                if (scattering_mechanisms[3]==1)
                    nu_deformation[counter] = nu_de(k_dum,counter,T,c_n,v_n[counter]);
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

                    nu_piezoelectric[counter] = nu_pe(k_dum,counter,T,c_n,P_piezo[T_loop],epsilon_s1,v_n[counter]);

		// Dislocation scattering
                if (scattering_mechanisms[6]==1)
                    nu_dislocation[counter] = nu_dis(k_dum,counter,T,beta_constant,epsilon_s1,v_n[counter]);

		// Alloy scattering
                if (scattering_mechanisms[7]==1)
                    nu_alloy[counter] = nu_alloy1(k_dum, v_n[counter]);


                if (scattering_mechanisms[8]==1)  // intravalley scattering
                {
                    for (int aa = 0;aa<iv_number;aa++)
                    {
                        lambda_e_plus_grid[counter][aa] = abs(lambda_e_plus(counter,k_grid,we[aa],rho,De[aa],nfv[aa],energy_n,v_n,points));
                        lambda_e_minus_grid[counter][aa] = abs(lambda_e_minus(counter,k_grid,we[aa],rho,De[aa],nfv[aa],energy_n,v_n,points));
                        // Equation number 129 of rode book
                    }
                }

                if (scattering_mechanisms[9]==1)  // neutral impurity
                    {			
			nu_neutralimpurity[counter] = nu_im(k_dum,counter,epsilon_s1,N_im[ii],v_n[counter]);
		    }

                df0dk_grid[counter] = df0dk(k_grid, k_dum, T, E_F, coefficients_cond, kindex_cond, a);

                f_dist[counter] = f0(energy_n[counter],E_F,T);
                //cout<<"In between "<<endl;
                //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
                //cout<<"E_F = "<<E_F<<endl;
                //cout<<"T = "<<T<<endl;

                thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, E_F, T, df0dz_integral_n,coefficients_cond, kindex_cond, a);;

                f0x1_f0[counter] = f0(energy_n[counter],E_F,T)*(1-f0(energy_n[counter],E_F,T));

                nu_el[counter] = nu_deformation[counter] + nu_piezoelectric[counter] + nu_ionizedimpurity[counter]
                + nu_dislocation[counter] + nu_alloy[counter] + nu_neutralimpurity[counter];

                electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
		// unit is 1/s 
		
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


            }
 
        if (scattering_mechanisms[8]==1)  // intravalley scattering
        {
            for (int counter = 0;counter<points;counter++)
            {
		for (int aa = 0;aa<iv_number;aa++)
		{
		    double k_minus = kminus(counter,k_grid,we[aa],energy_n,points);

		    //cout<<endl<<"Inside"<<endl;
		    //cout<<"k_minus = "<<k_minus<<endl;

		    double arr[points];
		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_minus);
		    int minus_index =FindMinInd(arr,points);


		    double k_plus = kplus(counter,k_grid,we[aa],energy_n,points);
		    //cout<<"k_plus = "<<k_plus<<endl;

		    for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - k_plus);
		    int plus_index =FindMinInd(arr,points);

		    //cout<<"plus_index = "<<plus_index<<endl;

		    double f_negative = f0(energy_n[minus_index],efef_n,T);
		    double f_positive =  f0(energy_n[plus_index],efef_n,T);


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
		     	
		}
		nu_el[counter] = nu_el[counter] + nu_iv_total[counter];  
            }
        }

//------------------------------------ tools for BTE END --------------------------------------------------------------
// ----------------------------saving data ---------------------------------------------------
            /*
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
 
//---------------------------------------- solve BTE for g ----------------------------------------------------------

            double g[points]={0};
	    	
            double g_rta[points]={0};

            double g_old[points]={0};
            double g_LO[points]={0};
            double g_iv[points]={0};

            double g_th[points]={0};
            double g_th_old[points]={0};
            double g_LO_th[points]={0};


            double S_o_grid[points]={0};
            for (int i=0;i<points;i++)
                S_o_grid[i] = 1;

            double S_o_grid_total[points]={0};

            double S_i_grid[points]={0};
            double S_iLO_grid[points]={0};
            double S_i_th_grid[points]={0};
            double S_iLO_th_grid[points]={0};


            double result_g[points][iterations+1]={0};
            double result_g_LO[points][iterations+1]={0};
            double result_g_th[points][iterations+1]={0};
            double result_f[points][iterations+1]={0};


            for(int counter=0;counter<points;counter++)
            {
                result_g[counter][0] = k_grid[counter];
                result_g_LO[counter][0] = k_grid[counter];
                result_g_th[counter][0] = k_grid[counter];
                result_f[counter][0] = k_grid[counter];
            }


            for (int iteration = 0;iteration<iterations;iteration++)
            {

                double f_dist_temp[points],f_dist_temp_th[points];
                for (int counter1=0;counter1<points;counter1++)
                {
                    f_dist_temp[counter1] = f0(energy_n[counter1],E_F,T) + g[counter1];
                    f_dist_temp_th[counter1] = f0(energy_n[counter1],E_F,T) + g_th[counter1];
                    //cout<<"energy_n[counter1] = "<<energy_n[counter1]<<endl;
                    //cout<<"g[counter1] = "<<g[counter1]<<endl;
                    //cout<<"f_dist_temp[counter1] = "<<f_dist_temp[counter1]<<endl;
                    //cout<<"f_dist_temp_th[counter1] = "<<f_dist_temp_th[counter1]<<endl;
                    //getchar();
                }


                double sum=0;
                for (int counter1 = 0;counter1 < points;counter1++)
                    sum = sum +  S_o_grid[counter1];

                double average_dummy = sum/points;
                //cout<<"average dummy = "<<average_dummy<<endl;

                for (int counter1 = 0;counter1 < points;counter1++)
                {
                    k_dum =	k_grid[counter1];
                    //f_dist(counter) = 1/(1+exp((energy(counter)-E_F)/(k_B*T)))+g(counter);

                    //f_dist_th(counter)= 1/(1+exp((energy(counter)-E_F)/(k_B*T)))+g_th(counter);

                    double arr[points];
                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kminus_grid[counter1]);
                    int minus_index =FindMinInd(arr,points);

                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kplus_grid[counter1]);
                    int plus_index =FindMinInd(arr,points);

                    // If POP scattering is included
                    if (scattering_mechanisms[1] == 1)
                    {
                        S_i_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g[plus_index];

                        S_iLO_grid[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_LO[plus_index];


                        S_i_th_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_th[minus_index]+
                        (N_poph_atT + 1 - f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_th[plus_index];


                        S_iLO_th_grid[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*g_LO_th[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*g_LO_th[plus_index];


                        if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
                            S_o_grid[counter1] = average_dummy;
                            // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
                        else
                            S_o_grid[counter1] = (N_poph_atT+1-f_dist_temp[minus_index])*lambda_o_minus_grid[counter1] +
                            (N_poph_atT+f_dist_temp[plus_index])*lambda_o_plus_grid[counter1];
                            // Equation no. 116 of rode book
                    }

                }
                
                for (int counter1=0;counter1<points;counter1++)
                    S_o_grid_total[counter1] = S_o_grid[counter1];


                for (int counter1=0;counter1<points;counter1++)
                {
                    g[counter1] = (S_i_grid[counter1]+electric_driving_force[counter1])/(S_o_grid_total[counter1] + nu_el[counter1]);
                    g_th[counter1] = (S_i_th_grid[counter1] + thermal_driving_force[counter1])/(S_o_grid_total[counter1] + nu_el[counter1]);
                }

                // If POP scattering is included
                if (scattering_mechanisms[1] == 1)
                {
                    for (int counter1=0;counter1<points;counter1++)
                        g_LO[counter1] = (S_iLO_grid[counter1] + electric_driving_force[counter1])/S_o_grid[counter1];
                }

                for (int counter1=0;counter1<points;counter1++)
                {
                    result_g[counter1][iteration+1] = g[counter1];

                    result_g_th[counter1][iteration+1] = g_th[counter1];

                    result_f[counter1][iteration+1] = f_dist[counter1];

                //fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',iteration,T,sum(g-g_old)/points);

                    g_old[counter1] = g[counter1];
                    g_th_old[counter1] = g_th[counter1];

                    if (iteration==0)
                       g_rta[counter1] = g[counter1] ;
                }
            }
//---------------------------------------- solve BTE for g finished----------------------------------------------------------





//---------------------------------------- Solve BTE for g and h for hall mobility---------------------------------------------
	    double beta1[points]={0};
	    double gH[points]={0};
	    double hH[points]={0};
	    
	    double gH_rta[points]={0};
	    double hH_rta[points]={0};

	    double gH_LO[points]={0};
	    double hH_LO[points]={0};

	    double S_i_grid_g[points]={0};
	    double S_i_grid_h[points]={0};

            double S_iLO_grid_g[points]={0};
            double S_iLO_grid_h[points]={0};

            double S_o_gridH[points]={0};
            
            for (int i=0;i<points;i++)
                S_o_gridH[i] = 1;

            double S_o_grid_totalH[points]={0};


	if(Bfield!=0)
	{
            for (int iteration = 0;iteration<iterations;iteration++)
            {

                double f_dist_temp[points];
                double beta1_LO[points];
                
                for (int counter1=0;counter1<points;counter1++)
                {
                    f_dist_temp[counter1] = f0(energy_n[counter1],E_F,T) + gH[counter1] + hH[counter1];
                    //cout<<"energy_n[counter1] = "<<energy_n[counter1]<<endl;
                    //cout<<"gH[counter1] = "<<gH[counter1]<<endl;
                    //cout<<"hH[counter1] = "<<hH[counter1]<<endl;
                    //cout<<"f_dist_temp[counter1] = "<<f_dist_temp[counter1]<<endl;
                    //getchar();
                }


                double sum1=0;
                for (int counter1 = 0;counter1 < points;counter1++)
                    sum1 = sum1 +  S_o_gridH[counter1];

                double average_dummy1 = sum1/points;
                //cout<<"average dummy1 = "<<average_dummy1<<endl;

                for (int counter1 = 0;counter1 < points;counter1++)
                {
                    k_dum =	k_grid[counter1];

                    double arr[points];
                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kminus_grid[counter1]);
                    int minus_index =FindMinInd(arr,points);

                    for (int i=0;i<points;i++)
                        arr[i] = abs(k_grid[i] - kplus_grid[counter1]);
                    int plus_index =FindMinInd(arr,points);

                    // If POP scattering is included
                    if (scattering_mechanisms[1] == 1)
                    {
                        S_i_grid_g[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*gH[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*gH[plus_index];

                        S_iLO_grid_g[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*gH_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*gH_LO[plus_index];

                        S_i_grid_h[counter1] = (N_poph_atT + f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*hH[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*hH[plus_index];

                        S_iLO_grid_h[counter1] = (N_poph_atT+f_dist_temp[counter1])*lambda_i_minus_grid[counter1]*hH_LO[minus_index]+
                        (N_poph_atT+1-f_dist_temp[counter1])*lambda_i_plus_grid[counter1]*hH_LO[plus_index];

                        if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
                            S_o_gridH[counter1] = average_dummy1;
                            // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
                        else
                            S_o_gridH[counter1] = (N_poph_atT+1-f_dist_temp[minus_index])*lambda_o_minus_grid[counter1] +
                            (N_poph_atT+f_dist_temp[plus_index])*lambda_o_plus_grid[counter1];
                            // Equation no. 116 of rode book
                    }
                }
                
                for (int counter1=0;counter1<points;counter1++)
                    S_o_grid_totalH[counter1] = S_o_gridH[counter1];


                for (int counter1=0;counter1<points;counter1++)
                {

                	beta1[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(S_o_grid_totalH[counter1] + nu_el[counter1]));
                	
                	beta1_LO[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(S_o_grid_totalH[counter1]));

                    	gH[counter1] = (S_i_grid_g[counter1] + electric_driving_force[counter1] + beta1[counter1] * 			S_i_grid_h[counter1])/((S_o_grid_totalH[counter1] + nu_el[counter1])*(1 + beta1[counter1]*beta1[counter1]));
                    	
                    	hH[counter1] = (S_i_grid_h[counter1] - beta1[counter1] * electric_driving_force[counter1] -                    	 			beta1[counter1] * S_i_grid_g[counter1])/((S_o_grid_totalH[counter1] + nu_el[counter1]) * (1 + 				beta1[counter1]*beta1[counter1]));

                    	gH_LO[counter1] = (S_iLO_grid_g[counter1] + electric_driving_force[counter1] + beta1_LO[counter1] * 			S_iLO_grid_h[counter1])/((S_o_grid_totalH[counter1] )*(1 + beta1_LO[counter1]*beta1_LO[counter1]));
                    	
                    	hH_LO[counter1] = (S_iLO_grid_h[counter1] - beta1_LO[counter1] * electric_driving_force[counter1] 				- beta1_LO[counter1] * S_iLO_grid_g[counter1])/((S_o_grid_totalH[counter1]) * (1 + 				       beta1_LO[counter1]*beta1_LO[counter1]));
                }
                
                if (iteration==0)
	        {
		        for (int counter1=0;counter1<points;counter1++)
		        {
				gH_rta[counter1] = gH[counter1] ;
				hH_rta[counter1] = hH[counter1] ;
			}	        	
	        }

            }
	    	
 	
//-------------------------------------------------------------------------------------------------------------------------------
	}

//------------------------------------Saving scattering rate--------------------------------------------------------
            if (ispin == 1)
            {			    
	            fid1 = fopen("scattering_rate.dat","w");
	            fid2 = fopen("mean_free_path.dat","w");
            }
            if (ispin == 2 && kk == 0)			    
	    {
	            fid1 = fopen("scattering_rate_up_spin.dat","w");
	            fid2 = fopen("mean_free_path_up_spin.dat","w");
            }
            if (ispin == 2 && kk == 1)			    
	    {
	           fid1 = fopen("scattering_rate_down_spin.dat","w");
	           fid2 = fopen("mean_free_path_down_spin.dat","w");
	    }

	    	
        	fprintf(fid1,"# energy           nu_im            nu_po          nu_de            nu_pe           nu_dis        nu_alloy        nu_iv           nu_nim       nu_total (total relaxation time 1/second)\n");

        	fprintf(fid2,"# energy           Mean free path (nm)\n");
	    
	    double nu_total[points]={0}, inv_nu_total[points]={0}, mfp[points]={0};	
	    
            for (int i = 0; i < points; i++)
            {
                nu_total[i] = scattering_mechanisms[0] * nu_ionizedimpurity[i] + scattering_mechanisms[1] * S_o_grid_total[i] + 
                scattering_mechanisms[3] * nu_deformation[i] + scattering_mechanisms[4] * nu_piezoelectric[i] + 
                scattering_mechanisms[6] * nu_dislocation[i] + scattering_mechanisms[7] * nu_alloy[i] + 
                scattering_mechanisms[8] * nu_iv_total[i] + scattering_mechanisms[9] * nu_neutralimpurity[i];
                
                if(nu_total[i]!=0)
                	inv_nu_total[i] = 1.0/nu_total[i];
                
                mfp[i] = v_n[i] * inv_nu_total[i];
                
                mfp[i] = mfp[i]/(1e-9);
                
                fprintf(fid1,"  %e    %e    %e    %e    %e    %e    %e    %e    %e 		%e \n", energy_n[i], nu_ionizedimpurity[i], 
		S_o_grid_total[i],nu_deformation[i],nu_piezoelectric[i],nu_dislocation[i], 
		nu_alloy[i], nu_iv_total[i], nu_neutralimpurity[i], nu_total[i] );

                fprintf(fid2,"  %e    	%e   \n", energy_n[i], mfp[i] );
	    }
	    fclose(fid1);
	    fclose(fid2);
		
//----------------------------------save perturbation g(E)--------------------------------------------------------
            if (ispin == 1 )			    
	            fid1 = fopen("g.dat","w");

            if (ispin == 2 && kk == 0)			    
	            fid1 = fopen("g_up_spin.dat","w");

            if (ispin == 2 && kk == 1)			    
	            fid1 = fopen("g_down_spin.dat","w");

       	    fprintf(fid1,"# S.No.    Energy (eV)       g(E)\n");

            for (int i = 0; i < points; i++)
                {
			//cout<<"i+1 = "<<i+1<<"    g[i] = "<<g[i]<<endl;			
			fprintf(fid1,"  %d        %e      %e  \n", i+1, energy_n[i], g[i]);
			//getchar();
		}
		fclose(fid1);
	    


//------------------------------------Save gH(E) and hH(E)-----------------------------------------------------
	if(Bfield!=0)
	{

           if (ispin == 1 )
           {			    
	            fid1 = fopen("gH.dat","w");
	            fid2 = fopen("hH.dat","w");
	   }	
           if (ispin == 2 && kk == 0)			    
	   {
	            fid1 = fopen("gH_up_spin.dat","w");
	            fid2 = fopen("hH_up_spin.dat","w");
	   }	
            if (ispin == 2 && kk == 1)			    
	     {
	     	    fid1 = fopen("gH_down_spin.dat","w");
	            fid2 = fopen("hH_down_spin.dat","w");
	     }	
       	    fprintf(fid1,"# S.No.    Energy (eV)       gH(E)\n");
       	    fprintf(fid2,"# S.No.    Energy (eV)       hH(E)\n");


            for (int i = 0; i < points; i++)
                {
			//cout<<"i+1 = "<<i+1<<"    gH[i] = "<<gH[i]<<endl;			
			//cout<<"i+1 = "<<i+1<<"    hH[i] = "<<hH[i]<<endl;			
			fprintf(fid1,"  %d        %e      %e  \n", i+1, energy_n[i], gH[i]);
			fprintf(fid2,"  %d        %e      %e  \n", i+1, energy_n[i], hH[i]);
			//getchar();
		}
		fclose(fid1);
	    	fclose(fid2);
	}
//-------------------------- save data ---------------------------------------------------


            //cout<<"points = "<<points<<endl;
            //getchar();
            /*
            fid1 = fopen("S_i_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_i_grid[i]);
	fclose(fid1);

            fid1 = fopen("S_iLO_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_iLO_grid[i]);
	fclose(fid1);


            fid1 = fopen("S_i_th_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_i_th_grid[i]);
	fclose(fid1);

            fid1 = fopen("S_iLO_th_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_iLO_th_grid[i]);
	fclose(fid1);

            fid1 = fopen("S_o_grid.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_o_grid[i]);
	fclose(fid1);

            fid1 = fopen("S_o_grid_total.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, S_o_grid_total[i]);
	fclose(fid1);
            */
	    
			
            /*
            fid1 = fopen("g_th.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_th[i]);
	fclose(fid1);

            fid1 = fopen("g_LO.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_LO[i]);
	fclose(fid1);

            fid1 = fopen("g_iv.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_iv[i]);
	fclose(fid1);

            fid1 = fopen("g_rta.txt","w");
            for (int i = 0; i < points; i++)
                fprintf(fid1,"%d    %e\n", i+1, g_rta[i]);
	fclose(fid1);
            */
            //cout<<"saved result check here";
            //getchar();
            //getchar();
//-----------------------------------------------------------------------------------------

//----------------------------------calculate drift mobility -------------------------------------------------------------
            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);


            double mobility_npop = 1e10;
            double mobility_to = 1e10;

            cout<<endl;
            // ionized impurity scattering
            if (scattering_mechanisms[0]==1)
            {
                mobility_ii = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_ionizedimpurity,g,points,a);
                cout<<"mobility_ii = "<<mobility_ii<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_ii = 1e10;


            // POP scattering
            if (scattering_mechanisms[1]==1)
            {
                mobility_po = mu_po(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,g_LO,g,nu_el,points,a);
                cout<<"mobility_po = "<<mobility_po<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_po=1e10;

            // Acoustic deformation scattering
            if (scattering_mechanisms[3]==1)
            {
                mobility_de = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_deformation,g,points,a);
                cout<<"mobility_de = "<<mobility_de<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_de=1e10;

            // piezoelectric scattering
            if (scattering_mechanisms[4]==1)
            {
                mobility_pe = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_piezoelectric,g,points,a);
                cout<<"mobility_pe = "<<mobility_pe<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_pe=1e10;

            // dislocation scattering
            if (scattering_mechanisms[6]==1)
            {
                mobility_dis = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_dislocation,g,points,a);
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
                mobility_alloy = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_alloy,g,points,a);
                cout<<"mobility_alloy = "<<mobility_alloy<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_alloy=1e10;

            // inter-valley scattering
            if (scattering_mechanisms[8]==1)
            {
                mobility_iv = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_iv_total,g,points,a);
                cout<<"mobility_iv  = "<<mobility_iv<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_iv=1e10;

            // neutral impurity scattering
            if (scattering_mechanisms[9]==1)
            {
                mobility_neutral = mu_elastic(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,
                    Ds_n,v_n,nu_neutralimpurity,g,points,a);

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

            mobility = mu_overall(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,g,nu_el,points,a);
	    // unit cm^2/(V-s)
	    
            double mobility_rta = mu_overall(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,g_rta,nu_el,points,a);
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
            thermopower = -k_B*(df0dz_integral_n- E_F /(k_B*T))*1e6 + (J(k_grid,T,m,g_th,Ds_n,energy_n,v_n,points)/sigma)/dTdz*1e6;
            // Equation No. 52 of rode book


//----------------------------------------------------------------------------------------------------------------------



// -------------------------------- calculate hall mobility -------------------------------------------------------------
	if(Bfield!=0)
	{

            cout.precision(6);                        //set precision
            cout.setf(ios::scientific);

            double mobility_hall_npop = 1e10;
            double mobility_hall_to = 1e10;

            cout<<endl;
            // ionized impurity scattering
            if (scattering_mechanisms[0]==1)
            {
                mobility_hall_ii = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_ionizedimpurity,points,a);
                cout<<"mobility_hall_ii = "<<mobility_hall_ii<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_ii = 1e10;


            // POP scattering
            if (scattering_mechanisms[1]==1)
            {
                mobility_hall_po = mu_poH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,gH_LO,hH_LO,nu_el,points,a);
                cout<<"mobility_hall_po = "<<mobility_hall_po<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_po=1e10;

            // Acoustic deformation scattering
            if (scattering_mechanisms[3]==1)
            {
                mobility_hall_de = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_deformation,points,a);
                cout<<"mobility_hall_de = "<<mobility_hall_de<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_de=1e10;

            // piezoelectric scattering
            if (scattering_mechanisms[4]==1)
            {
                mobility_hall_pe = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_piezoelectric,points,a);
                cout<<"mobility_hall_pe = "<<mobility_hall_pe<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_pe=1e10;

            // dislocation scattering
            if (scattering_mechanisms[6]==1)
            {
                mobility_hall_dis = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_dislocation,points,a);
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
                mobility_hall_alloy = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_alloy,points,a);
                cout<<"mobility_hall_alloy = "<<mobility_hall_alloy<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_alloy=1e10;

            // inter-valley scattering
            if (scattering_mechanisms[8]==1)
            {
                mobility_hall_iv = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,nu_iv_total,points,a);
                cout<<"mobility_hall_iv  = "<<mobility_hall_iv<<" cm^2/(V-s)"<<endl;
            }
            else
                mobility_hall_iv=1e10;

            // neutral impurity scattering
            if (scattering_mechanisms[9]==1)
            {
                mobility_hall_neutral = mu_elasticH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,
                    Ds_n,v_n,nu_neutralimpurity,points,a);

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

            mobility_hall = mu_overallH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,gH,hH,nu_el,points,a);

            mobility_hall_rta = mu_overallH(k_grid,energy_n,E_F,T,coefficients_cond,kindex_cond,Ds_n,v_n,gH_rta, hH_rta, nu_el,points,a);

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
    }

    if (variation==0)   // temperature variation
    {
    
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","a");
			for (int i = 0; i <count_t ; i++)
			{
				fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			}
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 0)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","a");

		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower_up_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_up_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_up_spin_hall.dat","a");

			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);

		fid1 = fopen("thermopower_down_spin.dat","a");
		for (int i = 0; i <count_t ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_down_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_down_spin_hall.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","a");
			for (int i = 0; i <count_t ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

	
    }
    else          // Doping variation
    {
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);
		
		if(Bfield!=0)
		{
			fid1 = fopen("mobility_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);
			
			fid1 = fopen("hall_factor.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}

    	if (ispin == 2 && kk == 0)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower_up_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);		
		}
	}


    	if (ispin == 2 && kk == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
		            calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		            calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		            calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e         %e              %e \n",
		            calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
		fclose(fid1);
		
		fid1 = fopen("thermopower_down_spin.dat","a");
		for (int i = 0; i <count_d ; i++)
		    fprintf(fid1," %e        %e\n",
		            calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);

		if(Bfield!=0)
		{

			fid1 = fopen("mobility_down_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1],
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid1);

			fid1 = fopen("conductivity_down_spin_hall.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","a");
			for (int i = 0; i <count_d ; i++)
			    fprintf(fid1," %e         %e              %e \n",
				    hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
			fclose(fid1);
		}
	}
	
    }

    
    }	// end of loop ue to ispin

    return 0;
}
