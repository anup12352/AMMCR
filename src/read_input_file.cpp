#include"main.h"

double n0, Nd1,Na1, efef_n, efef_p, N_ii;
int cc=-1, count_d, count_t;
double mobility_ii, mobility_po, mobility_to, mobility_npop, mobility_de, mobility_pe, mobility_dis;
double mobility_alloy, mobility_iv, mobility_neutral, mobility_avg, mobility, mobility_rta;
double mobility_hall_ii, mobility_hall_po, mobility_hall_to, mobility_hall_npop, mobility_hall_de;
double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv;
double mobility_hall_neutral, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;

double g[2000], g_rta[2000], g_old[2000], g_LO[2000], g_iv[2000], g_th[2000], g_th_old[2000], g_LO_th[2000];
double S_o_grid[2000], S_o_grid_total[2000], S_i_grid[2000], S_iLO_grid[2000], S_i_th_grid[2000], S_iLO_th_grid[2000];
double result_g[2000][15+1], result_g_LO[2000][15+1], result_g_th[2000][15+1], result_f[2000][15+1];

double N_poph_atT, df0dz_integral_n, N_e[10], beta_constant; 

double k_min, k_trans, k_step_fine, k_step;
int points, points1, points2;
double df0dk_grid[2000], f0x1_f0[2000], electric_driving_force[2000], thermal_driving_force[2000], f_dist[2000];

double kplus_grid[2000], kminus_grid[2000], betaplus_grid[2000], betaminus_grid[2000];

double  Aminus_grid[2000], Aplus_grid[2000], lambda_i_plus_grid[2000], lambda_o_plus_grid[2000];
double lambda_i_minus_grid[2000], lambda_o_minus_grid[2000], lambda_e_plus_grid[2000][5], lambda_e_minus_grid[2000][5];

double nu_deformation[2000], nu_piezoelectric[2000], nu_ionizedimpurity[2000], nu_dislocation[2000], nu_alloy[2000];
double nu_neutralimpurity[2000], nu_iv[2000][5], nu_iv_total[2000], nu_el[2000];

double Ed;
int a11[2],b11[2];
double h_bar,CBM,VBM;
int count1,count2;
int count_orbital, count_orbital_p;

double beta1[2000], gH[2000], hH[2000], gH_rta[2000], hH_rta[2000];
double gH_LO[2000], hH_LO[2000], S_i_grid_g[2000], S_i_grid_h[2000];
double S_iLO_grid_g[2000], S_iLO_grid_h[2000], S_o_gridH[2000], S_o_grid_totalH[2000];

double mobility_all[10]={0} , calc_mobility[30][2] = {0}, calc_mobility_rta[30][2] = {0};
double calc_thermopower[30][2] = {0}, calc_sigma[30][2] = {0}, calc_sigma_rta[30][2] = {0};

double calc_mobility_pe[30][1] = {0}, calc_mobility_de[30][1] = {0}, calc_mobility_dis[30][1] = {0}, calc_mobility_ii[30][1] = {0};
double calc_mobility_po[30][1] = {0}, calc_mobility_to[30][1] = {0}, calc_mobility_alloy[30][1] = {0}, calc_mobility_iv[30][1] = {0};
double calc_mobility_neutral[30][1] = {0};

double mobility_hall_all[10]={0}, calc_mobility_hall[30][2] = {0}, calc_mobility_hall_rta[30][2] = {0};
double calc_sigma_hall[30][2] = {0}, calc_sigma_hall_rta[30][2] = {0}, hall_factor[30][2] = {0}, hall_factor_rta[30][2] = {0};

double calc_mobility_hall_pe[30][1] = {0}, calc_mobility_hall_de[30][1] = {0}, calc_mobility_hall_dis[30][1] = {0};
double calc_mobility_hall_ii[30][1] = {0}, calc_mobility_hall_iv[30][1] = {0}, calc_mobility_hall_neutral[30][1] = {0};
double calc_mobility_hall_po[30][1] = {0}, calc_mobility_hall_to[30][1] = {0}, calc_mobility_hall_alloy[30][1] = {0}; 
double kcbm[3],kvbm[3];

double k_grid[2000]={0}, v_n[2000]={0}, v_p[2000]={0};	
double energy_n[2000]={0}, energy_p[2000]={0}, a_n[2000]={0}, c_n[2000]={0};

double a_p[2000]={0}, c_p[2000]={0}, Ds_n[2000]={0}, Ds_p[2000]={0};

double B_ii = 0, D_ii = 0, A_ii = 0;

double N_im_de, Nd_plus,N_im_modified;

int kk, ispin;

double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];

double we[10],De[10];
int nfv[10];

int degree1;
double fraction[4];
int length_fraction;

int De_ionization,N_cb, N_vb, iterations, variation, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

double rho, k_max, N_dis, omega_LO, omega_TO, E_deformation_n, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

double Bfield;

int free_e,len_T,len_n;


void read_input_file()
{

	cout.precision(6);                        //set precision
	cout.setf(ios::scientific);
	cout<<"Data from input.dat file"<<endl;
//-------------------------------------------------------	
	C_11=0;
	C_12=0;
	C_44=0;
//-----------------------------------------------------
 	for (int i=0;i<30;i++)
	{
	T_array[i]=0; epsilon_s[i]=0; epsilon_inf[i]=0; Bgap[i]=0; P_piezo[i] = 0;
	C_piezo_h14[i]=0; n_array[i]=0;  Nd[i]=0; Na[i]=0; N_im[i]=0;
	}

	for (int i=0;i<10;i++)
	{
	we[i]=0; nfv[i]=0; De[i]=0;
	}

	De_ionization = 0;
	double Ed = 0.0092; // in eV   ionization energy  For ag 1 paper

	if (De_ionization==0)
	Ed = 0; // in eV   ionization energy

	vector<string> all_data;
	ifstream fe("input.dat");


	// Read all data line by line in vector all_data
	string line;
	while(getline(fe, line))
	all_data.push_back(line);

	/*
	for(int i=0 ; i < all_data.size(); ++i)
	cout << all_data[i] << '\n';
	*/


//-------------------------------------------------------------------------
	// Read temperature array
	stringstream ss(all_data[0]);

	int count = 0,count1;int i=0; 
	while (ss >> T_array[i]) 
	{ count++; i++; } 

	//cout<<"Number of temp array = "<<count<<endl;

	len_T = count;
	
	
	// To display temperature array
	cout<<"Temp array"<<endl;
	for(int i=0;i<len_T;i++)
	{
	  cout << T_array[i]<<"K    ";
	}		
	cout<<endl<<endl;


//-------------------------------------------------------------------------
	// Read donor doping array
	// Count No. of doping points
	stringstream ss1(all_data[1]);

	count = 0; i=0;
	while (ss1 >> Nd[i]) 
	{ count++;  i++;  }

	//cout<<"Number of donor doping array = "<<count<<endl;

	len_n = count;
	// To display donor doping array
	cout<<"donor doping array "<<endl;
	for(int i=0;i<len_n;i++)
	{
	  cout << Nd[i]<<" per cm^3    ";
	}		
	cout<<endl<<endl;
	

//-------------------------------------------------------------------------
	// Read acceptor doping array
	// Count No. of doping points
	stringstream ss2(all_data[2]);

	count = 0; i=0;
	while (ss2 >> Na[i]) 
	{ count++;  i++;  }
	//cout<<"count    = "<<count<<endl;
	//cout<<"Number of acceptor doping array = "<<count<<endl;
		
	// To display acceptor doping array
	cout<<"acceptor doping array "<<endl;
	for(int i=0;i<len_n;i++)
	{
	  cout << Na[i]<<" per cm^3    ";
	}		
	cout<<endl<<endl;

//-----------------------------------------------------------------------

if (count != len_n )
{
	cout<<"Error Donor array length and acceptor array length must be same. Exit from program "<<endl;
	exit (EXIT_FAILURE);
}
else
{
	cout<<"Net doping array"<<endl;
	for(int j=0;j<len_n;j++)
		{n_array[j] = Nd[j] - Na[j];
		cout<<n_array[j]<<"   "<<"cm^(-3)";}
		cout<<endl;		

}
//-------------------------------------------------------------------------
	// Read neutral impurity doping array
	// Count No. of neutral impurity doping points

	stringstream ss3(all_data[3]);

	count = 0; i=0;
	while (ss3 >> N_im[i]) 
	{  count++;   i++;   }
	
	//cout<<"Number of neutral impurity array = "<<count<<endl;


	if (len_n != count)
	{
		cout<<"Error both doping array and neutral impurity array should have same length, exit from program"<<endl;
		exit (EXIT_FAILURE);
	}
		
	
	// To display neutral impurity doping array
	cout<<"Neutral impurity doping array "<<endl;
	for(int i=0;i<len_n;i++)
	{
	  cout << N_im[i]<<" per cm^3    ";
	}		
	cout<<endl<<endl;
	

//-------------------------------------------------------------------------
	// Read scattering array
	
	stringstream  ss10(all_data[10]);
	int a1[9];
	for(int i=0;i<9;i++)
	  ss10 >> a1[i];
	
	scattering_mechanisms[0] = a1[0];       // Ionized imourity 
	scattering_mechanisms[1] = a1[1];     	// Polar Optical phonon scattering due to longitudinal phonon
	scattering_mechanisms[3] = a1[2];	// Acoustic deformation scattering
	scattering_mechanisms[4] = a1[3];	// Piezoelectric scattering
	//scattering_mechanisms[5] = a1[4];     // TO phonon
	scattering_mechanisms[6] = a1[4];	// Dislocation scattering 
	scattering_mechanisms[7] = a1[5];	// Alloy scattering
	scattering_mechanisms[8] = a1[6];	// Intra-valley scattering
	scattering_mechanisms[9] = a1[7];	// Neutral impurity scattering
	
	/*
	// To display a1 array
	cout<<"a1 array "<<endl;
	for(int i=0;i<9;i++)
	{
	  cout << a1[i]<<"    ";
	}		
	cout<<endl;
	*/
//-----------------------------------------------------------------------------------

	// read intra_valley constants
	// Read we means frequency of phonon
	stringstream ss20(all_data[20]);

	count1 = 0; i=0;
	while (ss20 >> we[i]) 
	{  count1++;   i++;   }
	
	// Read De coupling constant
	stringstream ss21(all_data[21]);

	int count2 = 0; i=0;
	while (ss21 >> De[i]) 
	{  count2++;   i++;   }

	// Read nfv number of final valley
	stringstream ss22(all_data[22]);

	int count3 = 0; i=0;
	while (ss22 >> nfv[i]) 
	{  count3++;   i++;   }



	if (count1 != count2)
	{
		cout<<"Error same number of intra valley constants are required"<<endl;
		exit (EXIT_FAILURE);
	}
		
	if (count1 != count3)
	{
		cout<<"Error same number of intra valley constants are required"<<endl;
		exit (EXIT_FAILURE);
	}
	
	iv_number = count1; 
	// To display intra-valley constants
	for(int i=0;i<iv_number;i++)
	{
		cout<<"frequency["<<i<<"] = "<<we[i]<<"THz  "<<endl;
		we[i] = we[i]*2*pi*1e12;

		cout<<"Coupling constant["<<i<<"] = "<<De[i]<<"e8 eV/cm  "<<endl;

		cout<<"number of final valleys["<<i<<"] = "<<nfv[i]<<endl<<endl;
	}


//-------------------------------------------------------------------------
	FILE *fid;
	fid = fopen("input.dat", "r");
	if (fid==NULL)
	{
		cout<<"input.dat fìle is not present exit from program"<<endl;
		exit;
	}
	char line1[1000];
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	
	
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &epsilon_s[0]);
	cout<<"epsilon_s = "<<epsilon_s[0]<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &epsilon_inf[0]);
	cout<<"epsilon_inf = "<<epsilon_inf[0]<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &Bgap[0]);
	cout<<"Bgap = "<<Bgap[0]<<" eV"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%d", &N_vb);
	cout<<"N_vb = "<<N_vb<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%d", &N_cb);
	cout<<"N_cb = "<<N_cb<<endl;
	/*
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &c_lattice);
	cout<<"c_lattice = "<<c_lattice<<" nm"<<endl;
	*/
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &rho);
	cout<<"rho = "<<rho<<" g/cm^3"<<endl;


	fgets(line1, 1000, fid);

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &N_dis);
	cout<<"N_dis= "<<N_dis<<" /cm^2"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &omega_LO);
	cout<<"freqency of Logitudinal Optical phonon = "<<omega_LO<<"THz"<<endl;	
	omega_LO = omega_LO*2*pi*1e12;

	/*
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &omega_TO);
	cout<<"omega_TO = "<<omega_TO<<endl;
	omega_TO = omega_TO*2*pi*1e12;	
	*/

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &E_deformation_n);
	cout<<"E_deformation_n = "<<E_deformation_n<<" eV"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &P_piezo[0]);
	cout<<"P_piezo =  "<<P_piezo[0]<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &C_long);
	cout<<" C_long =  "<<C_long<<" dyne/cm^2"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &C_trans);
	cout<<"C_trans =  "<<C_trans<<" dyne/cm^2"<<endl;
	
	/*
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &C_piezo_c14);
	cout<<"C_piezo_c14 =  "<<C_piezo_c14<<endl;
	*/

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &Uall);
	cout<<"Uall (Alloy Potential) = "<<Uall<<" eV"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &V0);
	cout<<"V0 (Voulme of the primitive cell)= "<<V0<<" nm^3"<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &xx);
	cout<<"xx (fraction)= "<<xx<<endl;
	

	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	
	fgets(line1, 1000, fid);
	sscanf(line1, "%d", &free_e);
	cout<<"free_e = "<<free_e<<endl;

	fgets(line1, 1000, fid);
	sscanf(line1, "%d", &iterations);
	cout<<"iterations = "<<iterations<<endl;
	
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	fgets(line1, 1000, fid);
	sscanf(line1, "%lf", &Bfield);
	cout<<"Magnetic Field = "<<Bfield<<endl;

//--------------------Read fraction of band--------------------------------------------------

	// Read fraction of band 
	stringstream ss26(all_data[26]);
	
	fraction[0] = 0;
	fraction[1] = 0;
	fraction[2] = 0;
	fraction[3] = 0; 

	int count26 = 0; i=0;
	while (ss26 >> fraction[i]) 
	{  count26++;   i++;   }

	// To display fraction array
	cout<<"Fraction array"<<endl;
	for(int i=0;i<count26;i++)
	{
	  cout << fraction[i]<<"    "<<endl;
	}

	length_fraction = count26;
	//cout<<"Length of fraction array given by user is  "<<length_fraction<<endl;
	if (count26==1 && fraction[0]==0)
	{
		//cout<<"Code will automatically calculate points for division of wavevector for minimum discontinuity"<<endl;
	}	
	else if ((count26) > 4 || (count26) < 3)
	{
		cout<<"Three or four points should be given for division of k segment. Exit from program "<<endl;
		exit (EXIT_FAILURE);		
	}		
	cout<<endl;

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------------
	if (len_T > 1)
 	{
		for (int i=1;i<len_T;i++)
		{
			epsilon_s[i]=epsilon_s[0]; 
			epsilon_inf[i]=epsilon_inf[0]; 
			Bgap[i]=Bgap[0]; 
		}
	}

	for(int i=0;i<len_T;i++)
	C_piezo_h14[i] = C_piezo_c14/(epsilon_s[i]*epsilon_0);   // unit - N/C

	if (C_trans == 0)
	C_trans =  (C_11 - C_12 + 3 * C_44)/5;  // in dyne/cm2   Equation 99 from rode book

	if (C_long == 0)
	C_long =  (3*C_11 + 2 * C_12 + 4 * C_44)/5;         // in dyne/cm2   Equation 100 from rode book

	c_bar = (1.0/3.0)*C_long + (2.0/3)*C_trans;   // in dyne/cm^2

	if (P_piezo[0] == 0)
	{
		for(int i=0;i<len_T;i++)
		{
		    P_piezo[i] = (pow(C_piezo_h14[i],2)*epsilon_0*epsilon_s[i]*(12/C_long+16/C_trans)/35*1e1 );
		    P_piezo[i] = pow(P_piezo[i],0.5);
		// P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
		}
	}
	else
	{
		if (len_T>1)		
		for(int i=1;i<len_T;i++)
		{
		    P_piezo[i] = P_piezo[0];
		// P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
		}
		
	}


	//cout<<"len_n = "<<len_n<<endl;
	//cout<<"len_T = "<<len_T<<endl;
	
	if (len_n==1 & len_T==1)
		variation = 0;    // Any variation can be taken temperature variation
	else if (len_n==1)
		variation = 0;   // temp variation
	else if (len_T==1)
		variation=1;   // Doping variation
	else
	{	
		cout<<"Error, out of doping and temperature array one should be one length exit from the program"<<endl;
		exit (EXIT_FAILURE);
	}	

	k_max = 6;
	T_trans = 40;
	
	cout<<"End of Reading input.dat file"<<endl;
	
	//c_lattice = 0.591;  // in nm

	fclose(fid);


	count_d=0;
	count_t=0;

	int s1=0;
	while(n_array[s1]!=0)
	s1++;

	count_d = s1;

	s1=0;

	while(T_array[s1]!=0)
	s1++;

	count_t = s1;

	if(count_d == 0)
	count_d = count_d+1;

	N_cb = 1;
	N_vb = 1;

	//cout<<"Count_d = "<<count_d<<endl;
	//cout<<"Count_t = "<<count_t<<endl;

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

	h_bar =  h_planck / (2 * pi);    // Planck constant divided by 2pi[eV.s]

	fitting_1 = 0;  // FOR BAND
	fitting_2 = 0;  // FOR PROCAR
	fitting_3 = 1;  // FOR DOSCAR

	return ;
}

