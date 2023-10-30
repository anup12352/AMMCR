#include"main.h"

int flag[50]={0};

double energies[limit3][2]={0}, kpoints[limit3][4]={0}, temp1[limit3][3]={0};
int vbm_index=0, cbm_index=0;

double ecbm,evbm,kcbm1[3],kvbm1[3];
int NKPTS,NBVAL,NBTOT;

int spin_orbit_coupling;
double ion_mass1[5];
double lm[10][10];
int ion_numbers1[5];
double volume1;
double c_lattice=0;
int number;
string TB_material;

//---------------------------------------------------------------------------------------------------------------
double vf, vf_cb, vf_vb;
int linear_fit=0, SORT=0, save_data=0;
//---------------------------------------------------------------------------------------------------------------
double lambda_so;  // effective spin orbit coupling

//---------------------------------------------------------------------------------------------------------------
double Efield_time[limit1];
double omega_s, initial, freq[limit9];
double J_time[limit1]={0}, sigma_time[limit1]={0}, mobility_time[limit1]={0};
int time_variation;
int show_acoustic_contribution=0;
int time_limit, len_freq;
double g_time[limit2]={0},g_time_old[limit2]={0};
double omegam[limit10];
double J_freqr[limit10]={0}, sigma_freqr[limit10]={0}, Efield_freqr[limit10]={0}, J_freqi[limit10]={0}, sigma_freqi[limit10]={0}, Efield_freqi[limit10]={0};
double mobility_freqr2[limit10]={0}, mobility_freqi2[limit10]={0}, sigma_freqr2[limit10]={0}, sigma_freqi2[limit10]={0};
double mobility_drude_freqr2[limit10]={0}, mobility_drude_freqi2[limit10]={0}, sigma_drude_freqr2[limit10]={0}, sigma_drude_freqi2[limit10]={0};
int freq_variation=0;
double tau;

//-------------------------------- 2D variables -------------------------------------------------------------
double q[limit6+1], pz[limit6+1], X[limit7+1], Y[limit7+1], Z[limit7+1], theta[limit7+1], thickness=0.01;
double So_ab_npop[limit5][limit2], So_em_npop[limit5][limit2], Se_npop[limit5][limit2], Sa_npop[limit5][limit2];
double So_ab_pop[limit5][limit2], So_em_pop[limit5][limit2], Se_pop[limit5][limit2], Sa_pop[limit5][limit2];
double So_ab_so_pop[limit5][limit2], So_em_so_pop[limit5][limit2], Se_so_pop[limit5][limit2], Sa_so_pop[limit5][limit2];
double So_pop[limit5][limit2], Si_pop[limit5][limit2], So_so_pop[limit5][limit2], Si_so_pop[limit5][limit2];
int gd[limit5], method=1, overlap = 0;
double eps_sub_low, eps_sub_high, eps_up_low, eps_up_high, eps_avg_low, screening=0;
int pop_number, so_pop_number;
double we_pop[limit5], we_to[limit5];
double nu_pop_total[limit2], dist=0, nu_so_pop_total[limit2];
double rimp;

double fl_delta=0, fl_L=0, ns=0, ND=0, ni=0;

double beta1[limit2], gH[limit2], hH[limit2], gH_rta[limit2], hH_rta[limit2];
double gH_pop[limit2], hH_pop[limit2], gH_so_pop[limit2], hH_so_pop[limit2], Si_grid_g[limit2], Si_grid_h[limit2];
double Si_pop_grid_g[limit2], Si_pop_grid_h[limit2], Si_so_pop_grid_g[limit2], Si_so_pop_grid_h[limit2];

double Si_pop_parts[limit5][limit2], Si_so_pop_parts[limit5][limit2],  g_pop_parts[limit5][limit2]={0}, g_so_pop_parts[limit5][limit2]={0};

double g[limit2], g_rta[limit2], g_old[limit2], g_iv[limit2], g_so_pop[limit2], g_pop[limit2];
double g_th[limit2], g_th_old[limit2], g_th_pop[limit2], g_th_so_pop[limit2];

double result_g[limit2][limit8+1]={0}, result_g_pop[limit2][limit8+1]={0}, result_g_so_pop[limit2][limit8+1]={0}; 
double result_g_th[limit2][limit8+1]={0}, result_g_th_pop[limit2][limit8+1]={0}, result_g_th_so_pop[limit2][limit8+1]={0};
double result_gH_pop[limit2][limit8+1]={0}, result_hH_pop[limit2][limit8+1]={0}; 
double result_gH_so_pop[limit2][limit8+1]={0}, result_hH_so_pop[limit2][limit8+1]={0}; 

double Si_grid[limit2], Si_pop_grid[limit2], Si_so_pop_grid[limit2];
double Si_th_grid[limit2], Si_th_pop_grid[limit2], Si_th_so_pop_grid[limit2];

double result_fqr[limit2][limit8+1]={0}, result_fqi[limit2][limit8+1]={0}; 
double result_fqr_pop[limit2][limit8+1]={0}, result_fqi_pop[limit2][limit8+1]={0}; 
double result_fqr_so_pop[limit2][limit8+1]={0}, result_fqi_so_pop[limit2][limit8+1]={0}; 

double fqr[limit2][limit9][limit10]={0}, fqi[limit2][limit9][limit10]={0};
double fqr_rta[limit2][limit9][limit10]={0}, fqi_rta[limit2][limit9][limit10]={0};
double fqr_pop[limit2][limit9][limit10]={0}, fqi_pop[limit2][limit9][limit10]={0};
double fqr_so_pop[limit2][limit9][limit10]={0}, fqi_so_pop[limit2][limit9][limit10]={0};
double Si_grid_fqr[limit2], Si_grid_fqi[limit2], Si_pop_grid_fqr[limit2], Si_pop_grid_fqi[limit2];
double Si_so_pop_grid_fqr[limit2], Si_so_pop_grid_fqi[limit2];
double sigma_r[limit10][limit9], sigma_i[limit10][limit9], mobility_r[limit10][limit9], mobility_i[limit10][limit9];   
double sigma_real, sigma_img;

//---------------------------------------------------------------------------------------------------------------

double nu_deformation_p[limit2][2][2]={0}, nu_ionizedimpurity_p[limit2][2][2]={0}, nu_el_p[limit2][2][2]={0};
double nu_npop_p[limit2][2][2]={0};
double nu_So_p[limit2][2][2]={0};
//---------------------------------------------------------------------------------------------------------------

double nu_deformation[limit2]={0}, nu_piezoelectric[limit2]={0}, nu_def[3][limit2]={0}, nu_pz[3][limit2]={0};
double nu_ionizedimpurity[limit2]={0}, nu_dislocation[limit2]={0}, nu_alloy[limit2]={0};
double nu_neutralimpurity[limit2]={0}, nu_npop[limit5][limit2]={0}, nu_npop_total[limit2]={0}, nu_iv[limit4][limit2]={0}, nu_iv_total[limit2]={0}, nu_el[limit2]={0}, nu_irs[limit2]={0}, nu_skew_rate[limit2]={0};

//---------------------------------------------------------------------------------------------------------------

string  type="n";

double n0, Nd1,Na1, efef_n, efef_p, N_ii;
int cc=-1, count_d, count_t, VASP=1, geometry=1;
double mobility_ii, mobility_po, mobility_to, mobility_pe, mobility_dis, mobility_so_pop;
double mobility_de, mobility_de1, mobility_de2, mobility_de3;
double mobility_alloy, mobility_iv, mobility_neutral, mobility_npop, mobility_avg, mobility, mobility_rta, mobility_ir, mobility_skew;
double mobility_hall_ii, mobility_hall_po, mobility_hall_so_po, mobility_hall_to;
double mobility_hall_de, mobility_hall_de1, mobility_hall_de2, mobility_hall_de3;
double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv, mobility_hall_ir, mobility_hall_skew;
double mobility_hall_neutral, mobility_hall_npop, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta, peltier, thermal_conductivity;

double calc_sigma_xx[limit10][2], calc_sigma_xy[limit10][2], calc_hall_coeff[limit10][2], calc_long_restivity[limit10][2], calc_mg_resist[limit10][2];
double sigma_xx, sigma_xy, hall_coeff, long_restivity, mg_resist;

double calc_sigma_xx_rta[limit10][2], calc_sigma_xy_rta[limit10][2], calc_hall_coeff_rta[limit10][2];
double calc_long_restivity_rta[limit10][2], calc_mg_resist_rta[limit10][2];
double sigma_xx_rta, sigma_xy_rta, hall_coeff_rta, long_restivity_rta, mg_resist_rta;


double denom[limit2];
int plus_index_pop[limit5][limit2], minus_index_pop[limit5][limit2];
int plus_index_so_pop[limit5][limit2], minus_index_so_pop[limit5][limit2];

double N_poph_atT[limit5], df0dz_integral, N_e[limit4], beta_constant; 

double k_min, k_trans, k_step_fine, k_step;
int points, points1, points2;
double df0dk_grid[limit2], f0x1_f0[limit2], electric_driving_force[limit2], thermal_driving_force[limit2], f_dist[limit2];

double kplus_grid_pop[limit5][limit2], kminus_grid_pop[limit5][limit2];
double kplus_grid_so_pop[limit5][limit2], kminus_grid_so_pop[limit5][limit2];

double betaplus_grid[limit5][limit2], betaminus_grid[limit5][limit2];
double Aminus_grid[limit5][limit2], Aplus_grid[limit5][limit2];
double lambda_i_plus_grid[limit5][limit2], lambda_o_plus_grid[limit5][limit2];
double lambda_i_minus_grid[limit5][limit2], lambda_o_minus_grid[limit5][limit2];
double lambda_e_plus_grid[limit2][limit4], lambda_e_minus_grid[limit2][limit4];
double lambda_e_plus_grid_npop[limit2][limit5], lambda_e_minus_grid_npop[limit2][limit5], N_npop[limit5];
int npop_number;

double Ed, Emax=0.5;
int a11[2],b11[2];
double h_bar,CBM,VBM;
int count1,count2;
int count_orbital, count_orbital_p;
double px, Bx, Kx, sx;
double result_gH[limit2][limit8+1], result_hH[limit2][limit8+1]; 

double S_o_gridH[limit2]={0}, S_o_grid_totalH[limit2]={0};

double mobility_all[13]={0} , calc_mobility[limit10][2] = {0}, calc_mobility_rta[limit10][2] = {0};
double calc_thermopower[limit10][2] = {0}, calc_sigma[limit10][2] = {0}, calc_sigma_rta[limit10][2] = {0};
double calc_peltier[limit10][2]={0}, calc_thermal_conductivity[limit10][2]={0} ;

double calc_mobility_pe[limit10][2] = {0}, calc_mobility_de[limit10][5] = {0}, calc_mobility_dis[limit10][2] = {0}, calc_mobility_ii[limit10][2] = {0};
double calc_mobility_po[limit10][2] = {0}, calc_mobility_to[limit10][2] = {0}, calc_mobility_alloy[limit10][2] = {0}, calc_mobility_iv[limit10][2] = {0}, calc_mobility_ir[limit10][2] = {0}, calc_mobility_skew[limit10][2] = {0};
double calc_mobility_neutral[limit10][2] = {0}, calc_mobility_npop[limit10][2] = {0}, calc_mobility_so_pop[limit10][2] = {0};

double mobility_hall_all[12]={0}, calc_mobility_hall[limit10][2] = {0}, calc_mobility_hall_rta[limit10][2] = {0};
double calc_sigma_hall[limit10][2] = {0}, calc_sigma_hall_rta[limit10][2] = {0};
double hall_factor[limit10][2] = {0}, hall_factor_rta[limit10][2] = {0};

double calc_mobility_hall_pe[limit10][2] = {0}, calc_mobility_hall_de[limit10][5] = {0}, calc_mobility_hall_dis[limit10][2] = {0};
double calc_mobility_hall_ii[limit10][2] = {0}, calc_mobility_hall_iv[limit10][2] = {0}, calc_mobility_hall_neutral[limit10][2] = {0}, calc_mobility_hall_npop[limit10][2] = {0};
double calc_mobility_hall_po[limit10][2] = {0}, calc_mobility_hall_so_po[limit10][2] = {0} ,calc_mobility_hall_to[limit10][2] = {0}, calc_mobility_hall_alloy[limit10][2] = {0}, calc_mobility_hall_ir[limit10][2] = {0}, calc_mobility_hall_skew[limit10][2] = {0}; 
double kcbm[3],kvbm[3];

double cond_band[limit2][2],val_band[limit2][2];

double k_grid[limit2]={0}, v_n[limit2]={0}, v_p[limit2]={0};	
double energy_n[limit2]={0}, energy_p[limit2]={0}, a_n[limit2]={0}, c_n[limit2]={0};

double a_p[limit2]={0}, c_p[limit2]={0}, Ds_n[limit2]={0}, Ds_p[limit2]={0};

double B_ii = 0, D_ii = 0, A_ii = 0;

double N_im_de, Nd_plus,N_im_modified;

int kk, ispin=1, TBS_material=0;
double lattice_constant;  

vector <double> e_pp;   // used for saving DOS for DOS for TBS
vector <double> e_nn;   // used for saving DOS for DOS for TBS
vector <double> xx1;     // used for saving energy for DOS


double T_array[limit10],epsilon_s[limit10],epsilon_inf[limit10],Bgap[limit10],P_piezo[limit10],C_piezo_h14[limit10],n_array[limit10],Nd[limit10]={0},Na[limit10]={0},N_im[limit10]={0};

double fermi;
double we[limit4],De[limit4];
int nfv[limit4];
double we_npop[limit5],De_npop[limit5];


int degree1;
double fraction[4];
int length_fraction;

int De_ionization,N_cb, N_vb, iterations=10, variation, scattering_mechanisms[15]={0}, iv_number, de_number=0, cel_number=0;
int fitting_1, fitting_2, fitting_3;

double rho=0, k_max, N_dis, omega_TO, E_deformation[3]={0}, C_long, C_trans, C_za=0, c_bar, C_11, C_12, C_44;
double C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

double Bfield=0;

int free_e=0;

//using namespace std;
 

void read_input_file()
{
	
    //char s[180];
    cout.precision(6);                        //set precision
	
	FILE *fid;	
	fid = fopen("input.dat", "r");
	if (fid==NULL)
	{
		cout<<"input.dat  file is not present. Exit from program";
		exit(EXIT_FAILURE);
	}
	fclose(fid);

        cout.setf(ios::scientific);
        cout<<"Data from input.dat file"<<endl;
        cout<< " -------------- "<<endl<<endl;
	
        //-------------------------------------------------------	
	int len_T,len_nn=0, len_na=0, len_nd=0;
	
        C_11=0;
        C_12=0;
        C_44=0;
      //-----------------------------------------------------
        for (int i=0;i<30;i++)
        {
		T_array[i]=0; epsilon_s[i]=0; epsilon_inf[i]=0; Bgap[i]=0; P_piezo[i] = 0;
		C_piezo_h14[i]=0; n_array[i]=0;  Nd[i]=0; Na[i]=0; N_im[i]=0;
        }

        for (int i=0;i<limit4;i++)
        {
	        we[i]=0; nfv[i]=0; De[i]=0;
        }

        De_ionization = 0;
        double Ed = 0.0092; // in eV   ionization energy  For ag 1 paper

        if (De_ionization==0)
	        Ed = 0; // in eV   ionization energy
        
        string str;
        string ss;
        
        ifstream in("input.dat");
        int count=0;

        while(!in.eof())
        {
		in>> str;
		
		// to check comments(if a line is not part of input then start tat line with #)
		if(str=="#")
		{
		  getline(in,ss);
		  cout<< "--COMMENT--" <<ss <<endl<<endl;
		}
		
		if(str=="TYPE")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>ss;
		  if(ss=="N")
		  	type="n";
		  else if(ss=="P")
		  	type="p";
		  else
		  {
		  	cout<<"Proper type is not given. Exit from program"<<endl;
		  	exit(EXIT_FAILURE);	
		  }	
		}
		
		if(str=="INPUT")     // to specify whether input is from VASP, tight binding band structure and from table form
		{			// In VASP five input files need to be given EIGENVAL, DOSCAR, PROCAR, OUTCAR and input.dat 
					// In tight binding only input.dat file
			getline(in,ss);
			stringstream tmp(ss);


			tmp>>ss;
			if(ss=="TABLE_FORM")      // // In table form ---> 4 files ---> EK_CB.dat  , EK_VB.dat --->    kx, ky, kz (unit is 1/cm )  and energy (eV) 
				VASP = 0;        // DOS.dat---> energy and DOS ---> to be mentioned// unit of energy is eV and unit of density of states is per eV per cm^3
			else if(ss=="VASP")
				VASP = 1;
			else if(ss=="TIGHT_BINDING_BAND_STRUCTURE")
				VASP = 2;				
			else
			{
				cout<<"Proper input are not given for VASP-INPUT. Exit from program"<<endl;
				exit(EXIT_FAILURE);	
			}			  
			//cout<< "VASP-INPUT  " <<VASP<<endl;
		}


		if(str=="GEOMETRY")
		{
			getline(in,ss);
			stringstream tmp(ss);

			tmp>>ss;
			if(ss=="3D")
				geometry = 1;
			else if(ss=="2D")
				geometry = 2;
			else
			{
				cout<<"Proper input are not given for GEOMETRY. Exit from program"<<endl;
				exit(EXIT_FAILURE);	
			}			  
			//cout<< "VASP-INPUT  " <<VASP<<endl;
		}


		if(str=="TEMPERATURE")   // unit K    line no. 1
		{
		  flag[0] = 1;  	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> T_array[i]) 
		  { count++; i++; } 
		  
		  len_T = count;                    
		  // To display temperature array
		  cout<<"Temp array"<<endl;
		  for(int i=0;i<len_T;i++)
		  {
		    cout << T_array[i]<<"K    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "TEMPERATURE " <<ss<<endl;
		}
		
		if(str=="DONOR_DOPING")  // unit per cm^3   or per cm^2    line no. 2
		{
		  flag[2] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> Nd[i]) 
		  { count++; i++; } 
		  
		  len_nd = count;                    
		}
		
		if(str=="CARRIER_CONC")  // unit per cm^3   or per cm^2    line no. 2
		{
		  flag[2] = 1;
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0;
		  while (tmp >> Nd[i])
		  { count++; i++; }

		  len_nd = count;
		}

		if(str=="ACCEPTOR_DOPING")   // unit per cm^3    line no. 3
		{
		  flag[2] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> Na[i]) 
		  { count++; i++; } 
		  
		  len_na = count;                    
		}
		

		if(str=="NEUTRAL_IMPURITY")    // unit per cm^3    line no. 4
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> N_im[i]) 
		  { count++; i++; } 
		  
		  len_nn = count;                    
		  // To display neutral impurity array 
		  cout<<"neutral impurity array "<<endl;
		  for(int i=0;i<len_nn;i++)
		  {
		    cout << N_im[i]<<" per cm^3    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "check saved"<< N_im[2]<<endl;
		  //cout<< "NEUTRAL IMPURITY CONCENTRATION " <<ss<<endl;
		}
		
		if(str=="DIELECTRIC_CONST_LF")    // unitless     line no. 5
		{
		  flag[2] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> epsilon_s[0]; 
			  
		  // To display low frequency dielectric constant
		  cout<<"low frequency dielectric constant "<<endl;
		  cout <<epsilon_s[0];

		  cout<<endl<<endl;
		  //cout<< "popW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="DIELECTRIC_CONST_HF")     // unitless     line no. 6
		{
		  flag[3] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> epsilon_inf[0]; 
		  
		  // To display high frequency dielectric constant
		  cout<<"high frequency dielectric constant "<<endl;
		  cout<<epsilon_inf[0];
		  
		  cout<<endl<<endl;
		 // cout<< "HIGH FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="REMOTE_IMPURITY_CONCENTRATION")   // unit per cm^2  used for 2D denotes remote impurity concnetration in substrate layer
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> rimp;
		  cout<< " REMOTE_IMPURITY_CONCENTRATION =  " <<rimp<<"  cm^-2 "<<endl;
		}

		if(str=="SUBSTRATE_DIELECTRIC_CONST_LF")    // unitless   used for 2D  // unitless  used for 2D    remote impurity or remote interfacial phonon
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> eps_sub_low; 
			  
		  // To display low frequency dielectric constant
		  cout<<"low frequency dielectric constant of substrate = "<<endl;
		  cout <<eps_sub_low;

		  cout<<endl<<endl;
		  //cout<< "popW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="SUBSTRATE_DIELECTRIC_CONST_HF")    // unitless  used for 2D    remote impurity or remote interfacial phonon
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> eps_sub_low; 
			  
		  // To display low frequency dielectric constant
		  cout<<"low frequency dielectric constant of substrate = "<<endl;
		  cout <<eps_sub_low;

		  cout<<endl<<endl;
		  //cout<< "popW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}

		
		if(str=="TOP_OXIDE_DIELECTRIC_CONST_LF")     // unitless  used for 2D // unitless  used for 2D    remote impurity or remote interfacial phonon
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> eps_up_low; 
			  
		  // To display low frequency dielectric constant
		  cout<<"low frequency dielectric constant of top oxide = "<<endl;
		  cout <<eps_up_low;

		  cout<<endl<<endl;
		  //cout<< "popW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="TOP_OXIDE_DIELECTRIC_CONST_HF")     // unitless  used for 2D  // unitless  used for 2D    remote impurity or remote interfacial phonon
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> eps_up_high; 
			  
		  // To display low frequency dielectric constant
		  cout<<"High frequency dielectric constant of top oxide = "<<endl;
		  cout <<eps_up_high;

		  cout<<endl<<endl;
		  //cout<< "High FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}


		if(str=="BAND_GAP")      // unit eV    line no. 7
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>Bgap[0];

		  // To display band-gap
		  cout<<"BAND-GAP "<<endl;
		  cout<<Bgap[0];
				
		  cout<<endl<<endl;
		  
		  //cout<< "ENERGY BANDGAP " <<ss<<endl;
		}
		
		 
		if(str=="TBS_MATERIAL")       // to fix tight binding material 1 means graphene, 2 means siliciene, 3 means germanene, 4 means Stanene
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  //cout<< "Tight binding bandstructure calculation are to be done for 2D material "<<endl;
		  tmp>>TB_material;
			
			// lattice_constant is not distance between nearest neighbour
			//  lattice_constant unit is in Angestron
			if(TB_material=="GRAPHENE") // graphene
			{
				TBS_material = 1;
				cout<<"Graphene is selected for tight binding bandstructure calculation"<<endl;
				lattice_constant = 2.46;     
				geometry = 2;   // for 2D materials
			}
			else if(TB_material == "SILICENE") // Silicene
			{
				TBS_material = 2;
				cout<<"Silicene is selected for tight binding bandstructure calculation"<<endl;
				lattice_constant = 3.88;
				geometry = 2;   // for 2D materials
			}
			else if(TB_material == "GERMANENE") // germanene
			{
				TBS_material = 3;
				cout<<"Germanene is selected for tight binding bandstructure calculation"<<endl;
				lattice_constant = 4.10;
				geometry = 2;   // for 2D materials			
			}
			else if(TB_material == "STANENE") //  Stanene
			{
				TBS_material = 4;
				cout<<"Stanene is selected for tight binding bandstructure calculation"<<endl;
				lattice_constant = 5.22;
				geometry = 2;   // for 2D materials				
			}
			else
			{
				cout<<"Wrong input is given for for tight binding bandstructure calculation. Exit from program"<<endl;
				exit(EXIT_FAILURE);

			}
		}
		
		
		if(str=="DENSITY_OF_SEMICONDUCTOR")    // density of semiconductor unit is g/cm^3 for 3D or gm/cm^2 for 2D    line no. 10
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> rho;	
		  //cout<<"Rho = "<<rho<<endl;	
		  //getchar();  
		}

		if(str=="DISLOCATION_DENSITY")     // unit cm^2 used dislocation scattering   line no. 12
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> N_dis;
		  cout<< "DISLOCATION DENSITY " <<N_dis<< " /cm^2 "<<endl;
		}

		if(str=="LATTICE_CONSTANT")     // // used for dislocation scattering 3D-ntype  Its unit is nm
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>c_lattice;
		  printf("\n Lattice constant   =   %e  nm \n", c_lattice);
		}

		if(str=="POP_FREQUENCY")   // unit is THZ    line no. 13
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we_pop[i]) 
			{ count++; i++; } 
		
			pop_number = count;

			cout<< "PHONON-FREQUENCY for POP  "<<endl;
			for(int i=0;i<pop_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we_pop[i]<<"THz  "<<endl;
				we_pop[i] = we_pop[i]*2*pi*1e12;
			}
			cout<<endl;
		}


		if(str=="POP_FREQUENCY_SUBSTRATE")      // unit is THz ; OPTICAL PHONON-FREQUENCY for pop frequency of substrate POP for only 2D
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we_to[i]) 
			{ count++; i++; } 
		
			so_pop_number = count;

			cout<< "OPTICAL PHONON-FREQUENCY for POP  "<<endl;
			for(int i=0;i<so_pop_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we_to[i]<<"THz  "<<endl;
				we_to[i] = we_to[i]*2*pi*1e12;
			}
			cout<<endl;
		}


		if(str=="ACOUSTIC")    //  DP unit is eV  and //  CEL unit is dyne/cm2  line no. 14
		{
			getline(in,ss);
			stringstream tmp(ss);

			int j=0,i=0, count=0;
			double dummy[6]={0}, CEL[3] = {0};			
			while (tmp >> dummy[i]) 
			{ 
				if((i+1)%2==1)
				{
					E_deformation[j] = dummy[i];
					//cout<<"E_deformation[j] = "<<E_deformation[j]<<endl;
					//cout<<"j = "<<j<<endl;
				}
				else
				{
					CEL[j] = dummy[i];
					//cout<<"CEL[j] = "<<CEL[j]<<endl;
					//cout<<"j = "<<j<<endl;
					//cout<<"i = "<<i<<endl;
					j++;	
				}
				count++; 
				i++; 
				//cout<<"count = "<<count<<endl;
				//cout<<"j = "<<j<<endl<<endl;
				//getchar();
			} 	
			
			
			de_number = count/2;   // 
			//cout<<"de_number = "<<de_number<<endl;
			//cout<<"count = "<<count<<endl;
			if(count%2!=0)
			{
				cout<<"Proper inputs are not given for ACOUSTIC scattering. Exit from program"<<endl;
				exit(EXIT_FAILURE);							
			}	
				
			cout<< "Sr. No.     ACOUSTIC DEFORMATION POTENTIAL (eV)           Elastic constants (dyne/cm^2)"<<endl;
			for(int i=0;i<de_number;i++)
			{
				cout<<i+1<<"          "<<E_deformation[i]<< " eV                                  "<<CEL[i]<<" dyne/cm^2 "<<endl;
			}
			
			if(de_number==1)
			{
				C_long = CEL[0];
			}
			else if(de_number==2)
			{
				C_long = CEL[0];
				C_trans = CEL[1];
			}
			else
			{
				C_long = CEL[0];
				C_trans = CEL[1];
				C_za = CEL[2];
			}

			cout<<endl<<endl;
		}



		if(str=="INTERVALLEY")    //  //  De_iv unit is 10^8 eV/cm  and //  we_iv unit is THz
		{
			getline(in,ss);
			stringstream tmp(ss);

			int j=0,i=0, count=0;
			double dummy[15]={0};			
			while (tmp >> dummy[i]) 
			{ 
				if((i+1)%3==1)
				{
					we[j] = dummy[i];
					//cout<<"we[j] = "<<we[j]<<"   THz"<<endl;
					//cout<<"j = "<<j<<endl;
					we[j] = we[j]*2*pi*1e12;
				}
				else if((i+1)%3==2)
				{
					De[j] = dummy[i];
					//cout<<"De[j] = "<<De[j]<<"  10^8 eV/cm "<<endl;
					//cout<<"j = "<<j<<endl;	
				}
				else
				{
					nfv[j] = dummy[i];
					//cout<<"nfv[j] = "<<nfv[j]<<endl;
					//cout<<"j = "<<j<<endl;
					j++;	
				}
				count++; 
				i++; 
				//cout<<"count = "<<count<<endl;
				//cout<<"j = "<<j<<endl;
			} 	
			
			
			iv_number = count/3;   // 
			//cout<<"iv_number = "<<iv_number<<endl;
			//cout<<"count = "<<count<<endl;
			if(count%3!=0)
			{
				cout<<"Proper inputs are not given for INTERVALLEY scattering. Exit from program"<<endl;
				exit(EXIT_FAILURE);							
			}	
				
			cout<< "Sr. No.     Frequency for intervally scattering (Hz)         Coupling constants (1e8 eV/cm)		Number of final valley for scattering "<<endl;
			for(int i=0;i<iv_number;i++)
			{
				cout<<i+1<<"          "<<we[i]<< " Hz                                  "<<De[i]<<"    e8 eV/cm              "<<nfv[i]<<endl;
			}
			
			cout<<endl<<endl;
		}

		if(str=="NPOP_3D")    //  De_npop unit is 10^8 eV/cm  and //  we_npop unit is THz
		{
			getline(in,ss);
			stringstream tmp(ss);

			int j=0,i=0, count=0;
			double dummy[10]={0};			
			while (tmp >> dummy[i]) 
			{ 
				if((i+1)%2==1)
				{
					we_npop[j] = dummy[i];
					//cout<<"we_npop[j] = "<<we_npop[j]<<"   THz"<<endl;
					//cout<<"j = "<<j<<endl;
					we_npop[j] = we_npop[j]*2*pi*1e12;
				}
				else
				{
					De_npop[j] = dummy[i];
					//cout<<"De_npop[j] = "<<De_npop[j]<<"  10^8 eV/cm "<<endl;
					//cout<<"j = "<<j<<endl;
					j++;	
				}

				count++; 
				i++; 
				//cout<<"count = "<<count<<endl;
				//cout<<"j = "<<j<<endl;
			} 	
			
			
			npop_number = count/2;   // 
			//cout<<"npop_number = "<<npop_number<<endl;
			//cout<<"count = "<<count<<endl;
			if(count%2!=0)
			{
				cout<<"Proper inputs are not given for NPOP scattering. Exit from program"<<endl;
				exit(EXIT_FAILURE);							
			}	
				
			cout<< "Sr. No.     Frequency for NPOP scattering (Hz)         Coupling constants (1e8 eV/cm)    "<<endl;
			for(int i=0;i<npop_number;i++)
			{
				cout<<i+1<<"          "<<we_npop[i]<< " Hz                                  "<<De_npop[i]<<"    e8 eV/cm              "<<endl;
			}
			
			cout<<endl<<endl;
		}


		if(str=="NPOP_2D")    //  De_npop unit is 10^8 eV/cm  and //  we_npop unit is THz   
		{
			getline(in,ss);
			stringstream tmp(ss);

			int j=0,i=0, count=0;
			double dummy[15]={0};			
			while (tmp >> dummy[i]) 
			{ 
				if((i+1)%3==1)
				{
					we_npop[j] = dummy[i];
					//cout<<"we_npop[j] = "<<we_npop[j]<<"   THz"<<endl;
					//cout<<"j = "<<j<<endl;
					we_npop[j] = we_npop[j]*2*pi*1e12;
				}
				else if((i+1)%3==2)
				{
					De_npop[j] = dummy[i];
					//cout<<"De_npop[j] = "<<De_npop[j]<<"  10^8 eV/cm "<<endl;
					//cout<<"j = "<<j<<endl;	
				}
				else
				{
					gd[j] = dummy[i];
					//cout<<"gd[j] = "<<gd[j]<<endl;
					//cout<<"j = "<<j<<endl;	
					j++;
				}

				count++; 
				i++; 
				//cout<<"count = "<<count<<endl;
				//cout<<"j = "<<j<<endl;
			} 	
			
			
			npop_number = count/3;   // 
			//cout<<"npop_number = "<<npop_number<<endl;
			//cout<<"count = "<<count<<endl;
			if(count%3!=0)
			{
				cout<<"Proper inputs are not given for NPOP scattering. Exit from program"<<endl;
				exit(EXIT_FAILURE);							
			}	
				
			cout<< "Sr. No.     Frequency for NPOP scattering (Hz)         Coupling constants (1e8 eV/cm)		Number of final valley for scattering "<<endl;
			for(int i=0;i<npop_number;i++)
			{
				cout<<i+1<<"          "<<we_npop[i]<< " Hz                                  "<<De_npop[i]<<"    e8 eV/cm              "<<gd[i]<<endl;
			}
			
			cout<<endl<<endl;
		}
			
			
		if(str=="PZ_COEFFICIENT")    // unitless if 3D or C/cm if 2D    line no. 15
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>P_piezo[0];		  
		  cout<< "PIEZOELECTRIC COEFFICIENT " <<P_piezo[0]<<endl;
		  // unitless if 3D or C/cm if 2D
		  
		}


		if(str=="CBAR")     // // used for npop 3D-ptype  Its unit is dyne/cm2   
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>c_bar;
		  printf("\n c_bar   =   %e dyne/cm^2 \n", c_bar);
		}
		


		
		if(str=="ALLOY")    // unit eV     for alloy scattering      line no. 18
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> Uall;
		  cout<< "ALLOY POTENTIAL " <<Uall<< "eV"<<endl;
		  tmp>> V0;
		  cout<< "VOLUME OF PRIMITIVE CELL OF ALLOY  " <<V0<< " (nm)^3 "<<endl;
		  tmp>>xx;
		  cout<< "FRACTION OF ATOM FOR ALLOY " <<xx<<endl;
		  cout<<endl<<endl;
		  //getchar();
		}
		
		
		if(str=="INTERFACE_ROUGHNESS")    // RMS height unit is nm  ;  CORRELATION_LENGTH unit is nm;  
		{				    // SHEET_CARRIER_DENSITY  unit is per cm^2 ; DOPING_CARRIER_DENSITY  unit is per cm^2
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> fl_L;
		  cout<< " Correlation length of the fluctuations =   " <<fl_L<<" nm"<<endl;	
		  tmp>> fl_delta;
		  cout<< " RMS Height of the fluctuations =   " <<fl_delta<<" nm"<<endl;
		  tmp>>ND;
    		  cout<< " Doping Carrier Density =   " <<ND<<" per cm^2    "<<endl;		  			
		  tmp>>ns;
		  cout<< " Sheet Carrier Density =   " <<ns<<" per cm^2    "<<endl;	
		  cout<<endl<<endl;
		  
		  //getchar();
		}
	
		if(str=="MAGNETIC_FIELD")     // for setting magnetic field unit is Tesla for fall factor calculation    line no. 28
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>> Bfield;
			cout<< "MAGNETIC-FIELD " <<Bfield<<"  Tesla "<<endl;
		}		
		
		/*
		if(str=="DOS")    // DOS 0    ---->  0 means DOSCAR is used otherwise ;  DOS 1  ---> means free electron desnity is used    line no. 24
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>free_e;
		  cout<< "-DENSITY OF STATES " <<free_e<<endl;
		  // free_e = 0, means DOSCAR is used, otherwise free electron desnity is used 
                 if (free_e == 0 )
			cout<<"free_e = false, DOSCAR is used for density of states"<<endl;
		  else
			cout<<"free_e = true, free electron density is used for calculation"<<endl;
		  
		}
		*/

		if(str=="SCATTERING_MECHANISM_CONTROL")       // line no. 
		{
		  flag[4] = 1;	
		  getline(in,ss);
		  stringstream smc(ss);
		  string a1[12];
		  for(int i=0;i<12;i++)
		    smc >> a1[i];


		  if(a1[0]=="true")
			  scattering_mechanisms[0] = 1;       // Ionized impurity
		  else
			  scattering_mechanisms[0] = 0;       // Ionized impurity

		  if(a1[1]=="true")
		      scattering_mechanisms[1] = 1;     // Polar Optical phonon scattering due to longitudinal phonon
		  else
		      scattering_mechanisms[1] = 0;     // Polar Optical phonon scattering due to longitudinal phonon

		  if(a1[2]=="true")
			  scattering_mechanisms[2] = 1;     // npop phonon
		  else
			  scattering_mechanisms[2] = 0;     // npop phonon

		  if(a1[3]=="true")
			  scattering_mechanisms[3] = 1;	// Acoustic deformation scattering
		  else
			  scattering_mechanisms[3] = 0;	// Acoustic deformation scattering

		  if(a1[4]=="true")
			  scattering_mechanisms[4] = 1;	// Piezoelectric scattering
		  else
			  scattering_mechanisms[4] = 0;	// Piezoelectric scattering

		  if(a1[5]=="true")
			  scattering_mechanisms[6] = 1;	// Dislocation scattering
		  else
			  scattering_mechanisms[6] = 0;	// Dislocation scattering

		  if(a1[6]=="true")
			  scattering_mechanisms[7] = 1;	// Alloy scattering
		  else
			  scattering_mechanisms[7] = 0;	// Alloy scattering

		  if(a1[7]=="true")
			  scattering_mechanisms[8] = 1;	// Inter-valley scattering
		  else
			  scattering_mechanisms[8] = 0;	// Inter-valley scattering

		  if(a1[8]=="true")
			  scattering_mechanisms[9] = 1;	// Neutral impurity scattering
		  else
			  scattering_mechanisms[9] = 0;	// Neutral impurity scattering

		  if(a1[9]=="true")
			  scattering_mechanisms[10] = 1;	// so pop scattering
		  else
			  scattering_mechanisms[10] = 0;	// so pop scattering

		  if(a1[10]=="true")
			  scattering_mechanisms[11] = 1;	// IR scattering
		  else
			  scattering_mechanisms[11] = 0;	// IR scattering

		  if(a1[11]=="true")
			  scattering_mechanisms[12] = 1;	// skew scattering
		  else
			  scattering_mechanisms[12] = 0;	// skew scattering


		  //scattering_mechanisms[5] = 0;  // Transverse optical POP scattering be default is zero
		  
		  cout<< " SCATTERING MECHANISM CONTROL " <<ss<<endl;
		  
		  //cout<<"scattering_mechanisms[1] = "<<scattering_mechanisms[1]<<"   a1[1] = "<<a1[1]<<endl;	
		//getchar();
		
		  //cout<< "alloy" << scattering_mechanisms[7]<<endl;
		}

		if(str=="NUMBER_OF_ITERATIONS")   // line no. 25
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> iterations;
		  cout<< "NUMBER OF ITERARTIONS " <<iterations<<endl;
		}

		if(str=="DISTANCE_FROM_SHEET")     // for remote impurity scattering for 2D   unit is nm
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> dist;
		  cout<< " DISTANCE FROM SHEET =  " <<dist<<"  nm "<<endl;
		  dist = dist*1e-9;  // converted from nm to m
		}

		if(str=="OMEGA")   // for frequency variation
		{
			freq_variation = 1;
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> freq[i]) 
			{ count++; i++; } 
			
			len_freq = count;

			// To display omega array
			cout<<"Frequency array"<<endl;
			for(int i=0;i<len_freq;i++)
			{
			  cout << freq[i]<<" Hz   ";
			}		
			cout<<endl<<endl;		
		}

		
		
		if(str=="IMPURITY_CARRIER_DENSITY")     // for skew scattering unit is per cm^2  Ns in expression
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>ni;
			cout<< " Impurity Carrier Density =   " <<ni<<" per cm^2    "<<endl;		  			
		}
		
		//----------------------------------------------------------------------------------------------------------------
		
		/*
		if(str=="PHONON_FREQUENCY_NPOP")     // unit is THz 
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we_npop[i]) 
			{ count++; i++; } 
		
			if(npop_number!=0)
			{
				if(npop_number!=count)
				{
					cout<<"Error same number of npop constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				npop_number = count;

			cout<< "PHONON FREQUENCY for NPOP  "<<endl;
			for(int i=0;i<npop_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we_npop[i]<<"THz  "<<endl;
				we_npop[i] = we_npop[i]*2*pi*1e12;
			}
			cout<<endl;
		}

		
		if(str=="COUPLING_CONSTANT_NPOP")    // unit is 10^8 eV/cm
		{
			getline(in,ss);
			stringstream tmp(ss);
			int count = 0;
			int i=0; 
			while (tmp >> De_npop[i]) 
			{ count++; i++; } 
		
			if(npop_number!=0)
			{
				if(npop_number!=count)
				{
					cout<<"Error same number of npop constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				npop_number = count;
				
			cout<< "COUPLING CONSTANT for NPOP  "<<endl;
			for(int i=0;i<npop_number;i++)
				cout<<"Coupling constant["<<i<<"] = "<<De_npop[i]<<"e8 eV/cm  "<<endl;			
		}

		if(str=="NUMBER_OF_VALLEY_NPOP")      // for npop scatering 2D
		{
			getline(in,ss);
			stringstream tmp(ss);
			int count = 0; 
			int i=0;
			while (tmp >> gd[i]) 
			{  count++;   i++;   }

			if(npop_number!=0)
			{
				if(npop_number!=count)
				{
					cout<<"Error same number of npop constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				npop_number = count;

			for(int i=0;i<npop_number;i++)
			{
				cout<<"number of final valleys for npop scatttering ["<<i<<"] = "<<gd[i]<<endl;
			}
			
		}


		*/	
		
		/*
		if(str=="RMS_HEIGHT")     // for interface roughness scattering unit is nm
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>fl_delta;
			cout<< " RMS Height of the fluctuations =   " <<fl_delta<<" nm"<<endl;		  			
		}
		
		if(str=="CORRELATION_LENGTH")      // for interface roughness scattering unit is nm
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>fl_L;
			cout<< " Correlation length of the fluctuations =   " <<fl_L<<" nm"<<endl;		  			
		}
		
		
		if(str=="SHEET_CARRIER_DENSITY")     // for interface roughness scattering unit is per cm^2  Ns in expression
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>ns;
			cout<< " Sheet Carrier Density =   " <<ns<<" per cm^2    "<<endl;		  			
		}
		
		if(str=="DOPING_CARRIER_DENSITY")    // for interface roughness scattering unit is per cm^2  Nd in expression
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>ND;
			cout<< " Doping Carrier Density =   " <<ND<<" per cm^2    "<<endl;		  			
		}
		
		*/


		//*/

		/*
		// used for npop 3D ptype
		c_bar = (1.0/3.0)*C_long + (2.0/3)*C_trans;   // in dyne/cm^2
		printf("\n c_bar   =   %e dyne/cm^2 \n", c_bar);
		*/
		
		/*
		if(str=="CL")     // spherically averaged elastic constant for longitudinal modes. Its unit is dyne/cm2   For piezoelectric scattering  line no. 16
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> C_long;
		}

		if(str=="CT")    // spherically averaged elastic constant for transverse modes. Its unit is dyne/cm2     For piezoelectric scattering    line no. 17
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>C_trans;
		}
		

		if(str=="CZA")	     // spherically averaged elastic constant for ZO modes. Its unit is dyne/cm2     For piezoelectric scattering   
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>C_za;
		}
		
		
		
		if(str=="ALLOY_POTENTIAL")    // unit eV     for alloy scattering      line no. 18
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> Uall;
		  cout<< "ALLOY POTENTIAL " <<Uall<< "eV"<<endl;
		}

		if(str=="VOLUME_PRIMITIVE_CELL")    // unit nm^3       for alloy scattering    line no. 19
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> V0;
		  cout<< "VOLUME OF PRIMITIVE CELL OF ALLOY  " <<V0<< " (nm)^3 "<<endl;
		}

		if(str=="FRACTION_OF_ATOM")        // for alloy scattering    line no. 20
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>xx;
		  cout<< "FRACTION OF ATOM FOR ALLOY " <<xx<<endl;
		}
		
		
		if(str=="INTERVALLEY_PHONON_FREQUENCY")    // for INTER VALLEY scattering      // phonon frequency   unit is THZ   line no. 21
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we[i]) 
			{ count++; i++; } 
		
			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of inter valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			cout<< "PHONON FREQUENCY FOR INTERVALLEY SCATTERING "<<endl;
			for(int i=0;i<iv_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we[i]<<"THz  "<<endl;
				we[i] = we[i]*2*pi*1e12;
			}
			cout<<endl<<endl;

		}

		if(str=="COUPLING_CONSTANT_INTERVALLEY")    // for INTER VALLEY scattering  3D  unit is 10^8 eV/cm     line no. 22
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> De[i]) 
			{ count++; i++; } 
		
			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of inter valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			cout<< "COUPLING CONSTANT FOR INTERVALLEY SCATTERING  "<<endl;
			for(int i=0;i<iv_number;i++)
			{
				cout<<"Coupling constant["<<i<<"] = "<<De[i]<<"e8 eV/cm  "<<endl;
			}
			cout<<endl;


		}

		if(str=="NUMBER_OF_VALLEY")      // // for INTER VALLEY scattering   NUMBER OF final VALLEY    line no. 23
		{
			getline(in,ss);
			stringstream tmp(ss);
			int count = 0; 
			int i=0;
			while (tmp >> nfv[i]) 
			{  count++;   i++;   }

			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of inter valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			for(int i=0;i<iv_number;i++)
			{
				cout<<"number of final valleys ["<<i<<"] = "<<nfv[i]<<endl;
			}
			
		}
		*/

		//------------------------------------------- --------------------------------------------------------------------
		if(str=="FITTING_DOP")       //  line no. 26
		{
		  getline(in,ss);
		  cout<< "DEGREE OF POLYNOMIAL USED FOR FITTING " <<ss<<endl;
		}
		
		if(str=="K_SEGMENT")	     //  line no. 27
		{
			getline(in,ss);
			stringstream tmp(ss);
			fraction[0] = 0;
			fraction[2] = 0;
			fraction[2] = 0;
			fraction[3] = 0; 

			int count = 0;
			int i=0; 
			while (tmp >> fraction[i]) 
			{ count++; i++; } 

			cout<< "FRACTION array AT WHICH K-POINT IS DIVIDED "<<endl;

			for(int i=0;i<count;i++)
			{
			  cout << fraction[i]<<"    "<<endl;
			}

			length_fraction = count;
			//cout<<"Length of fraction array given by user is  "<<length_fraction<<endl;
			if (count==1 && fraction[0]==0)
			{
				//cout<<"Code will automatically calculate points for division of wavevector for minimum discontinuity"<<endl;
			}	
			else if ((count) > 4 || (count) < 3)
			{
				cout<<"Three or four points should be given for division of k segment. Exit from program "<<endl;
				exit (EXIT_FAILURE);		
			}		
			cout<<endl;
		}
	
		if(str=="N_VALENCE_BAND_VALLEYS")      // NUMBER OF VALENCE-BAND VALLIES IN BZ refere earlier manual  line no. 8
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>N_vb;
		  cout<< "NUMBER OF VALENCE-BAND VALLIES IN BZ " <<N_vb<<endl;
		}
		
		if(str=="N_CONDUCTION_BAND_VALLEYS")    // NUMBER OF CONDUCTION-BAND-BAND VALLIES IN BZ refere earlier manual   line no. 9
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>N_cb;
		  cout<< "NUMBER OF CONDUCTION-BAND VALLIES IN BZ " <<N_cb<<endl;
		}

  
		if(str=="SCREENING")      // leave it
		{
			screening = 1;
		}
		
		if(str=="SAVE_DATA")     // leave it
		{
			save_data = 1;
		}

		if(str=="LINEAR_FITTING")      // leave it
		{	
		  linear_fit = 1;	// linear fitting without intercept or intercept = 0	  
		  cout<< "Linear fitting is selected for band structure. Maximum value of energy given is "<<Emax<<endl;
		}

		if(str=="LINEAR_FITTING_2")
		{	
		  linear_fit = 2;	// linear fitting with intercept	  
		  cout<< "Linear fitting with intercept is selected for band structure"<<endl;
		}

		if(str=="EMAX")     // leave it
		{	
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> Emax;
		  cout<< "Maximum value of energy given for linear fitting is  "<<Emax<<"  eV "<<endl;
		}
	  
		if(str=="SORTING") // sorting is done for only distance from reference point      leave it
		{
			SORT = 1;		
			cout<< "SORTING is selected   "<<endl;
		}
		
		if(str=="Effective_spin_orbit_coupling")     // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>lambda_so;
			cout<< " Effective spin orbit coupling =   " <<lambda_so<<"  A^2 "<<endl;	
			
			lambda_so = lambda_so*1e-20;  //converted from A^2 to m^2
			  
		}

		//------------------------------------------- --------------------------------------------------------------------
		
		//---------------------time variation input read part ---------------------------------------------------		
		if(str=="TIME_VARIATION")    // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>time_variation;
			cout<< " Time_variation  =   " <<time_variation<<endl;		  
		}

		if(str=="OMEGA_S")     // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>omega_s;
			cout<< "OMEGA-S " <<omega_s<<endl;
		}
		
		if(str=="TIME_STEPS")   // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>time_limit;
			cout<< "STEPS  for time variations   =  " <<time_limit<<endl;
		}


		if(str=="INITIAL")    // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>initial;
			cout<< " INITIAL time  =   " <<initial<<endl;		  
		}
		
		if(str=="TAU")     // leave it
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>tau;
			cout<< " TAU =   " <<tau<<endl;		  
		}
		
		//------------------------time variation input completed--------------------------------------------------------------		

        } // end of while loop 

	//eps_sub_low = eps_sub_low*epsilon_0;
	//eps_sub_high = eps_sub_high*epsilon_0; 
	//eps_up_low = eps_up_low*epsilon_0;
	//eps_up_high = eps_up_high*epsilon_0;
	eps_avg_low = (eps_sub_low + eps_up_low)/2.0;

	  if(geometry==1)  // 3D
	  {
	  	cout<< "DENSITY OF SEMICONDUCTOR " <<rho<< " gm/cm3 "<<endl;
	  	rho = rho*1000;
		  //converted from g/cm^3 to kg/m^3
	  }
	  else if(geometry==2) // 2D
	  {
	  	cout<< "DENSITY OF SEMICONDUCTOR " <<rho<< " gm/cm2 "<<endl;
	  	rho = rho*10000;
		  //converted from g/cm^2 to kg/m^2
	  }

	if(geometry==1)  // 3D
	{
		// To display donor doping array
		cout<<"donor doping array "<<endl;
		for(int i=0;i<len_nd;i++)
		{
			cout << Nd[i]<<" per cm^3    ";
		}		
		cout<<endl<<endl;
		//cout<< "DONOR DOPING CONCENTRATION " <<ss<<endl;

		// To display acceptor doping array 
		cout<<"acceptor doping array "<<endl;
		for(int i=0;i<len_na;i++)
		{
			cout << Na[i]<<" per cm^3    ";
		}		
		cout<<endl<<endl;
		//cout<< "ACCEPTOR DOPING CONCENTRATION " <<ss<<endl;

	}
	else if(geometry==2) // 2D
	{
		// To display donor doping array
		cout<<"donor doping array "<<endl;
		for(int i=0;i<len_nd;i++)
		{
			cout << Nd[i]<<" per cm^2    ";
		}		
		cout<<endl<<endl;
		//cout<< "DONOR DOPING CONCENTRATION " <<ss<<endl;

		// To display acceptor doping array 
		cout<<"acceptor doping array "<<endl;
		for(int i=0;i<len_na;i++)
		{
			cout << Na[i]<<" per cm^2    ";
		}		
		cout<<endl<<endl;
		//cout<< "ACCEPTOR DOPING CONCENTRATION " <<ss<<endl;

	}
	
	/*
	if(geometry==1)	
	{
		cout<< "ELASTIC CONSTANT FOR LONGITUDINAL MODE " <<C_long<<" dyne/cm^2 "<<endl;
		cout<< "ELASTIC CONSTANT FOR TRANSVERSE MODE " <<C_trans<<" dyne/cm^2 "<<endl;

		if(de_number==3)		
			cout<< "ELASTIC CONSTANT FOR ZA MODE " <<C_za<<" dyne/cm^2 "<<endl;
	}
	else
	{
		cout<< "ELASTIC CONSTANT FOR LONGITUDINAL MODE " <<C_long<<" dyne/cm^2 "<<endl;
		
		if(de_number==2)
			cout<< "ELASTIC CONSTANT FOR TRANSVERSE MODE " <<C_trans<<" dyne/cm^2 "<<endl;

		if(de_number==3)		
		{
			cout<< "ELASTIC CONSTANT FOR TRANSVERSE MODE " <<C_trans<<" dyne/cm^2 "<<endl;
			cout<< "ELASTIC CONSTANT FOR ZA MODE " <<C_za<<" dyne/cm^2 "<<endl;
		}	
	}
	*/
	
	if(type=="n" && geometry==1)  // 3D and n type
	{
		  scattering_mechanisms[10] = 0;	// so pop scattering	
	}
	else if(type=="p" && geometry==1)	// 3D and p-type
	{
		  scattering_mechanisms[4] = 0;	// Piezoelectric scattering          
		  //scattering_mechanisms[5] = 0;  // Transverse optical POP scattering be default is zero
		  scattering_mechanisms[6] = 0;	// Dislocation scattering 
		  scattering_mechanisms[7] = 0;	// Alloy scattering
		  scattering_mechanisms[8] = 0;	// Inter-valley scattering
		  scattering_mechanisms[9] = 0;	// Neutral impurity scattering
		  scattering_mechanisms[10] = 0;	// so pop scattering
		  scattering_mechanisms[11] = 0;	// IR scattering
		  scattering_mechanisms[12] = 0;	// skew scattering
	}
	else // 2D 
	{
		  //scattering_mechanisms[5] = 0;  // Transverse optical POP scattering be default is zero
		  scattering_mechanisms[6] = 0;	// Dislocation scattering 
		  scattering_mechanisms[7] = 0;	// Alloy scattering
		  scattering_mechanisms[8] = 0;	// Inter-valley scattering
		  scattering_mechanisms[9] = 0;	// Neutral impurity scattering
	}	
		
	if(flag[0]==0)
	{
		cout<<"Temeprature is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	
	if(flag[2]==0)
	{
		cout<<"Doping is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	
	}
	
	/*
	else if(flag[2]==0)
	{
		cout<<"Low frequency dielectric constant is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	else if(flag[3]==0)
	{
		cout<<"High frequency dielectric constant is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	*/
	
	else if(flag[4]==0)
	{
		cout<<"Scattering Mechanism control is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	
	
	if(len_na ==0)
		len_na = len_nd;
	else if(len_nd ==0)
		len_nd = len_na;

	
	// To check array size of donor, neutral impurity and acceptor concentration
	if(len_nd!=len_na)
	{
		cout<<"Error Donor array length and acceptor array length array length must be same. Exit from program "<<endl;
		exit (EXIT_FAILURE);
	}

	// to fix ntype or ptype
	if(type =="n")
	{
		cout<<"Net doping array"<<endl;
		for(int j=0;j<len_nd;j++)
		{
			
			n_array[j] = Nd[j] - Na[j];
			cout<<n_array[j]<<"   "<<"cm^(-3)";

		}
		cout<<endl;		
	}
	else
	{
		cout<<"Net doping array"<<endl;
		for(int j=0;j<len_nd;j++)
		{
			
			n_array[j] = Na[j] - Nd[j];
			cout<<n_array[j]<<"   "<<"cm^(-3)";

		}
		cout<<endl;		
	}		
        
	// to fix value of variation
	if (len_nd==1 & len_T==1)
		variation = 0;    // Any variation can be taken temperature variation
	else if (len_nd==1)
		variation = 0;   // temp variation
	else if (len_T==1)
		variation=1;   // Doping variation
	else
	{	
		cout<<"Error, out of doping and temperature array one should be one length exit from the program"<<endl;
		exit (EXIT_FAILURE);
	}	


	count_d = len_nd;
	count_t = len_T;


	N_cb = 1;
	N_vb = 1;


	// copy bandgap dielectric constant in all elements of array	
	if (len_T > 1)
 	{
		for (int i=1;i<len_T;i++)
		{
			epsilon_s[i]=epsilon_s[0]; 
			epsilon_inf[i]=epsilon_inf[0]; 
			Bgap[i]=Bgap[0]; 
		}
	}

	if(VASP!=1)
	{
		if(rho==0)
		{
			cout<<"Density of material is required.  "<<endl;
			cout<<"Exit from program"<<endl;
			exit(EXIT_FAILURE);
		}
	}	

	if(TBS_material == 0 && VASP==2)
	{
		cout<<"Material input is not given for tight binding band strcture calculation. Programming is running by assuming graphene as input material"<<endl;
		TBS_material = 1;   // now graphene is selected as input material
	}
	for(int i=0;i<len_T;i++)
	C_piezo_h14[i] = C_piezo_c14/(epsilon_s[i]*epsilon_0);   // unit - N/C
	
	/*
	if (C_trans == 0)
	{
		C_trans =  (C_11 - C_12 + 3 * C_44)/5;  // in dyne/cm2   Equation 99 from rode book
		printf("\n C_trans =   %e dyne/cm^2 \n", C_trans);
	}
	
	if (C_long == 0)
	{
		C_long =  (3*C_11 + 2 * C_12 + 4 * C_44)/5;         // in dyne/cm2   Equation 100 from rode book
		printf("\n C_long  =   %e dyne/cm^2 \n", C_long);
	}
	
	*/
	
	
	/*
	if (P_piezo[0] == 0 && scattering_mechanisms[4] ==1)
	{
		for(int i=0;i<len_T;i++)
		{
		    P_piezo[i] = (pow(C_piezo_h14[i],2)*epsilon_0*epsilon_s[i]*(12/C_long+16/C_trans)/35*1e1 );
		    P_piezo[i] = pow(P_piezo[i],0.5);
		// P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
		}
		printf("\n P_piezo =   %e \n",P_piezo[0]);
		
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
	*/
	
	if (len_T>1)
	{		
		for(int i=1;i<len_T;i++)
		{
		    P_piezo[i] = P_piezo[0];
		// P_piezo is unitless for 3D and its unit C/m for 2D
		}
	}
	
	k_max = 6;
	T_trans = 40;
	
	//cout<<"End of Reading input.dat file"<<endl;

	//printf("\n P_piezo =   %e \n",P_piezo);
	
	if(geometry==1)
	{
		if (scattering_mechanisms[1] == 0)
		{
			if (iterations != 1 )
			{
			    iterations = 1;
			    cout<<"Since polar optical phonon scattering is not used, so iterations is set to 1"<<endl;
			}
			T_trans = 0;
		}
	}
	else
	{
		if (scattering_mechanisms[1] == 0 && scattering_mechanisms[10] == 0)
		{
			if (iterations != 1 )
			{
			    iterations = 1;
			    cout<<"Since polar optical phonon scattering and remote polar optical phonon scattering are not used, so iterations is set to 1"<<endl;
			}
			T_trans = 0;
		}
	}
	
	if (De_ionization == 1)
	{	// neutral impurity scattering
		if (scattering_mechanisms[9] == 0)
		{
		    scattering_mechanisms[9] = 1;
		    cout<<"De-ionization is to be included so neutral impurity scattering is to included in simulation";
		}
	}

	h_bar =  h_planck / (2 * pi);    // Planck constant divided by 2pi[eV.s]

	fitting_1 = 0;  // FOR BAND
	fitting_2 = 0;  // FOR PROCAR
	fitting_3 = 1;  // FOR DOSCAR
  
}
