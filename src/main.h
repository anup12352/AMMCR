
#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <cstring>
#include<iostream>
using namespace std;

#include<fstream>

#include<stdio.h>
#include <math.h>
#include <time.h>
#include<stdlib.h>
#include <vector>
#include <sys/stat.h>   // this for FileExists function
#include<bits/stdc++.h>
#include <algorithm>


//#include <opencv/cv.h>

//#include <Python.h>
//#include "pyhelper.hpp"



//using namespace std;
//using namespace cv;

//------------ - *** Some Global constants ***----------------------------------------------------
#define m_e 9.10938291e-31   // Electron mass[kg]
#define e 1.60217657e-19     // Electron charge[C]
#define pi 3.14159265359      // Pi number
#define h_planck 4.135667516e-15    // Planck constant[eV.s]
#define k_B 8.6173324e-5       // Boltzmann constant[eV / K]
#define epsilon_0 8.854187817e-12    // Absolute value of dielectric constant in vacuum[C ^ 2 / m ^ 2N]
//#define h_bar

#define k_min0 0.001
#define k_trans0 0.1
#define k_step_fine0  0.001
#define k_step0 0.01
#define E 1
#define dTdz 10
#define LORBIT 11

//------------ - ------------------------------------------------------------------------------

void read_OUTCAR();
void read_input_file();
void copyright(); 
void generate_output_files();

void find_cbm_vbm(int , int);
void calculate_mobility(double T, int ii);

void make_band(int type);
int* analytical_fitting(double band[2000][2],int counta, int aa);
double * polyfit1(double x[],double y[], int n, int N);   // n-degree   N - length of x[] and y[] array
void polyval(double p[], double x[], int nn, int mm);
double * linspace(double a, double b, int numbers);
double rsquare(double data[], double fitted[], int length);
int decompose_wave();
int decompose_wave_p();
void n_DOS();
int FindMinInd(double arr[],int length);

double conduction_dispersion(double k,double coefficients[5][7],double kindex[4],int aa[2]);
double dedk(double k, double coeff[5][7], double kindex[4], int aa[2]);
double admixture_value(double k,int coloumn);
double admixture_value_p(double k,int coloumn);
double DOS_value(double energy, int a);
void find_fermi(double n, double T, int ii);
double f0(double E1, double E_f, double T);
double beta(double T,int T_loop);
double N_poph(double omega, double T);
double kplus(int counter, double omega, int points);
double kminus(int counter, double omega, int points);
double nu_ii(double k, int counter, double beta_constant, double v, double epsilon_s);

double Aplus(int counter, double omega, int points);
double Aminus(int counter, double omega, int points);

double betaplus(int counter, double omega, double epsilon_s, double epsilon_inf , int points);
double betaminus(int counter, double omega, double epsilon_s, double epsilon_inf, int points);

double lambda_i_plus(int counter, double omega,double A_plus,double epsilon_s, double epsilon_inf, int points);
double lambda_i_minus(int counter,double omega,double A_minus,double epsilon_s, double epsilon_inf,int points);

double lambda_o_plus(int counter,double omega,double A_plus, double epsilon_s,double epsilon_inf,int points);
double lambda_o_minus(int counter,double omega,double A_minus, double epsilon_s,double epsilon_inf,int points);

double lambda_e_minus(int counter,double omega,double rho,double De,int nfv,int points);
double lambda_e_plus(int counter,double omega,double rho,double De,int nfv,int points);

double nu_de(double k,int counter,double T,double v);                    
double nu_pe(double k,int counter,double T,double P_piezo,double epsilon_s,double v);
double nu_dis(double k, int counter, double T, double beta_constant, double epsilon_s, double v);
double nu_alloy1(double k,double v);
double nu_im(double k_dum, int counter, double epsilon_s, double N_im, double v);

void generate_required_data(double T);

double df0dk(double k,double T,double E_f, double coefficients[5][7],double kindex[], int aa[]);
double df0dz(double k, double E_f, double T, double df0dz_integral, double coefficients[5][7], double kindex[], int aa[]);


double f(double k,double E_f,double T,double  coefficients[5][7],double kindex[],double g[], int points, int aa[]);
double fH(double k,double E_f,double T,double  coefficients[5][7],double kindex[],double g[], double h[],int points, int aa[]);

double mu_elastic(double E_f,double T,double coefficients[5][7],double kindex[],double nu_elastic[],double g[],int points, int aa[]);
double mu_elasticH(double E_f,double T,double coefficients[5][7],double kindex[],
        double nu_elastic[], int points, int aa[]);

double mu_overall(double E_f,double T,double coefficients[5][7],double kindex[],
        double g[],double nu_el[],int points,int aa[]);
double mu_overallH(double E_f,double T,double coefficients[5][7],double kindex[],
        double g[], double h[], double nu_el[],int points,int aa[]);
void save_scattering_rate();
void save_perturbation();
double mu_po(double E_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double g[],double nu_el[],int points, int aa[]);
double mu_poH(double E_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double h_LO[],double nu_el[],int points, int aa[]);
void components_BTE(double T, int T_loop, double efef_n, int ii);

double J(double T,double m,double g_th[],  int points);
double DOS_value1(double energy, int a);
void read_ispin();
void fitting_band();
void solve_g(double T);
void save_results();
void initialize_array();

/*
void find_cbm(int spinorbit);
void bubbleSort(double arr[], int n);
void swap(double *xp, double *yp);
void analytical_fitting(int type);
int FindMinInd(double *array, int size);
void polyFit(double *x, double *y, double *z, int n, int l1, int l2);
void polyval(double *p, double *x, double *y, int p1, int x1);
double rsquare(double *data, double *fitted, int n);
double *linspace(double initial_value, double final_value, int samples);
void write_files(char type);
double conduction_dispersion(double k);
double dedk(double k);
double find_fermi();
*/
//------------ - *** Some Global variables ***----------------------------------------------------
extern double h_bar;
//--------------------------------------------------------------------------------------------
extern int cc;
extern double lm[10][10], volume1, ion_mass1[5];
extern int ion_numbers1[5], spin_orbit_coupling;

extern double CBM, kcbm[3];
extern double kvbm[3], VBM;
//extern double conduction_band[10000][2];
//------------------------------------ make_band ----------------------------------------
extern double cond_band[2000][2], val_band[2000][2],poly[4000];
extern int count1,count2,countx;    //count1 for conduction band count2 for valence band
extern int count_orbital, count_orbital_p;
extern int count_DOS_n,count_DOS_p;
extern double E_F,n_e,n_h;
//extern double short_range_band_num[1000][3];
//extern int length_dd, maximum_degree, points ;
//extern double coefficients[10][10], kindex[10];
extern int number;
extern int fitting_1;
extern void read_procar();
//--------------------------------------- analytical_fitting----------------------------------
extern double coefficients[5][7],kindex[4],coefficients_cond[5][7],kindex_cond[4],coefficients_val[5][7],kindex_val[4];
extern int a11[2],b11[2];


extern int degree1,length_fraction;
extern double fraction[4];

//------------------------------------- inside read_input----------------------------------------------
extern double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];
extern double we[10],De[10];
extern int nfv[10];
extern int variation;
extern double kcbm[3],kvbm[3];
extern int kk, count_d ,count_t;
extern double k_min, k_trans, k_step_fine, k_step;
extern int points, points1, points2;

//extern double k_min0,k_trans0,k_step_fine0,k_step0;

extern int De_ionization,N_cb, N_vb, iterations, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

extern double Ed, c_lattice, rho, k_max, N_dis, omega_LO, omega_TO, E_deformation_n, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

extern double Bfield;
extern double kcbm[3],kvbm[3];

extern double beta1[2000], gH[2000], hH[2000], gH_rta[2000], hH_rta[2000];
extern double gH_LO[2000], hH_LO[2000], S_i_grid_g[2000], S_i_grid_h[2000];
extern double S_iLO_grid_g[2000], S_iLO_grid_h[2000], S_o_gridH[2000], S_o_grid_totalH[2000];
            
extern string  type;
extern int free_e;
extern int len_T,len_n;
//----------------------------find_cbm_vbm-------------------------------------------------------------------------
extern double ecbm,evbm,kcbm1[3],kvbm1[3];
extern double kcbm[3],kvbm[3],VBM,CBM;
extern int NKPTS,NBVAL,NBTOT;
extern double orbital_decomposedd[1000][4], orbital_decomposedd_p[1000][4];
extern int ispin;
extern double DOS_n[2000][2],DOS_p[2000][2];
extern int ispin;

extern double k_grid[2000], v_n[2000], v_p[2000], energy_n[2000], energy_p[2000];
extern double a_n[2000], c_n[2000], a_p[2000], c_p[2000], Ds_n[2000], Ds_p[2000];
extern double electric_driving_force[2000], thermal_driving_force[2000];
extern double B_ii, D_ii, A_ii;

extern double df0dk_grid[2000], f0x1_f0[2000], electric_driving_force[2000], thermal_driving_force[2000], f_dist[2000];

extern double kplus_grid[2000], kminus_grid[2000], betaplus_grid[2000], betaminus_grid[2000];

extern double  Aminus_grid[2000], Aplus_grid[2000], lambda_i_plus_grid[2000], lambda_o_plus_grid[2000];
extern double lambda_i_minus_grid[2000], lambda_o_minus_grid[2000], lambda_e_plus_grid[2000][5], lambda_e_minus_grid[2000][5];

extern double nu_deformation[2000], nu_piezoelectric[2000], nu_ionizedimpurity[2000], nu_dislocation[2000], nu_alloy[2000];
extern double nu_neutralimpurity[2000], nu_iv[2000][5], nu_iv_total[2000], nu_el[2000];

extern double beta1[2000], gH[2000], hH[2000], gH_rta[2000], hH_rta[2000];
extern double gH_LO[2000], hH_LO[2000], S_i_grid_g[2000], S_i_grid_h[2000];
extern double S_iLO_grid_g[2000], S_iLO_grid_h[2000], S_o_gridH[2000], S_o_grid_totalH[2000];

extern double N_poph_atT, df0dz_integral_n, N_e[10], beta_constant;
extern double kplus_grid[2000], kminus_grid[2000]; 

extern double g[2000], g_rta[2000], g_old[2000], g_LO[2000], g_iv[2000], g_th[2000], g_th_old[2000], g_LO_th[2000];
extern double S_o_grid[2000], S_o_grid_total[2000], S_i_grid[2000], S_iLO_grid[2000], S_i_th_grid[2000], S_iLO_th_grid[2000];
extern double result_g[2000][15+1], result_g_LO[2000][15+1], result_g_th[2000][15+1], result_f[2000][15+1];

extern double mobility_ii, mobility_po, mobility_to, mobility_npop, mobility_de, mobility_pe, mobility_dis;
extern double mobility_alloy, mobility_iv, mobility_neutral, mobility_avg, mobility, mobility_rta;
extern double mobility_hall_ii, mobility_hall_po, mobility_hall_to, mobility_hall_npop, mobility_hall_de;
extern double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv;
extern double mobility_hall_neutral, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
extern double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;

extern double mobility_all[10] , calc_mobility[30][2] , calc_mobility_rta[30][2] ;
extern double calc_thermopower[30][2] , calc_sigma[30][2] , calc_sigma_rta[30][2] ;

extern double calc_mobility_pe[30][1] , calc_mobility_de[30][1] , calc_mobility_dis[30][1] , calc_mobility_ii[30][1] ;
extern double calc_mobility_po[30][1] , calc_mobility_to[30][1] , calc_mobility_alloy[30][1] , calc_mobility_iv[30][1] ;
extern double calc_mobility_neutral[30][1] ;

extern double mobility_hall_all[10], calc_mobility_hall[30][2] , calc_mobility_hall_rta[30][2];
extern double calc_sigma_hall[30][2] , calc_sigma_hall_rta[30][2] , hall_factor[30][2] , hall_factor_rta[30][2] ;

extern double calc_mobility_hall_pe[30][1] , calc_mobility_hall_de[30][1] , calc_mobility_hall_dis[30][1] ;
extern double calc_mobility_hall_ii[30][1] , calc_mobility_hall_iv[30][1] , calc_mobility_hall_neutral[30][1] ;
extern double calc_mobility_hall_po[30][1], calc_mobility_hall_to[30][1], calc_mobility_hall_alloy[30][1];
extern double n0,Nd1,Na1, efef_n, efef_p, N_ii;
//-------------------------------------------------------------------------------------------------------------
#endif // MAIN_H_INCLUDED
