#include"main.h"

//double k_min0,k_trans0,k_step_fine0,k_step0;

double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];

double we[10],De[10];
int nfv[10];

int degree1;
double fraction[4];
int length_fraction;

double kcbm[3],kvbm[3];

int De_ionization,N_cb, N_vb, iterations, variation, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

double rho, k_max, N_dis, omega_LO, omega_TO, E_deformation_n, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

double Bfield;

string type;
int free_e,len_T,len_n;
char line[1000];

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
	type = "n";
	T_trans = 40;
	
	cout<<"End of Reading input.dat file"<<endl;
	
	//c_lattice = 0.591;  // in nm

	fclose(fid);

	return ;
}

