
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
void read_input();
void read_input_file();
void copyright(); 

void find_cbm_vbm(int , int);

void make_band(int type);
int* analytical_fitting(double band[2000][2],int counta, int aa);
double * polyfit1(double x[],double y[], int n, int N);   // n-degree   N - length of x[] and y[] array
void polyval(double p[], double x[], int nn, int mm);
double * linspace(double a, double b, int numbers);
double rsquare(double data[], double fitted[], int length);
int decompose_wave();
void n_DOS();
int FindMinInd(double arr[],int length);

double conduction_dispersion(double k,double coefficients[5][7],double kindex[4],int a[2]);
double dedk(double k, double coeff[5][7], double kindex[4], int a[2]);
double admixture_value(double k,int coloumn);
double DOS_value(double energy, int a);
void find_fermi(double k_grid[], double n, double T, int ii, double energy_n[], double energy_p[], double Ds_n[], double Ds_p[], int points);
double f0(double E1, double E_f, double T);
double beta(double T,double Ds[], double energy[], double k_grid[], double v[],  int points, int T_loop);
double N_poph(double omega, double T);
double kplus(int counter, double k_grid[], double omega, double energy[], int points);
double kminus(int counter, double k_grid[], double omega, double energy[], int points);
double nu_ii(double k, int counter, double B_ii, double D_ii, double beta_constant, double N_ii, double v, double epsilon_s);

double Aplus(int counter,double k_grid[],double omega, double a[],double c[],double energy[], int points);
double Aminus(int counter,double k_grid[],double omega, double a[],double c[],double energy[], int points);

double betaplus(int counter,double k_grid[],double omega,double epsilon_s, double epsilon_inf , double energy[], double v[],int points);
double betaminus(int counter,double k_grid[],double omega,double epsilon_s, double epsilon_inf, double energy[], double v[], int points);

double lambda_i_plus(int counter,double k_grid[],double omega,double A_plus,double c[],double epsilon_s,
                     double epsilon_inf,double energy[],double v[],int points);
double lambda_i_minus(int counter,double k_grid[],double omega,double A_minus,double c[],double epsilon_s,
                     double epsilon_inf,double energy[],double v[],int points);

double lambda_o_plus(int counter,double k_grid[],double omega,double A_plus,double a[],double c[],
        double epsilon_s,double epsilon_inf,double energy[],double v[],int points);
double lambda_o_minus(int counter,double k_grid[],double omega,double A_minus,double a[],double c[],
        double epsilon_s,double epsilon_inf,double energy[],double v[],int points);

double nu_de(double k,int counter,double T,double c[],double v);
double nu_pe(double k,int counter,double T,double c[],double P_piezo,double epsilon_s,double v);
double nu_dis(double k, int counter, double T, double beta_constant, double epsilon_s, double v);
double nu_alloy1(double k,double v);

double lambda_e_minus(int counter,double k_grid[],double omega,double rho,double De,int nfv,double energy[],double v[],int points);
double lambda_e_plus(int counter,double k_grid[],double omega,double rho,double De,int nfv,double energy[],double v[],int points);

double nu_im(double k_dum,int counter, double epsilon_s,double N_im,double v);

double df0dk(double k_grid[],double k,double T,double E_f, double coefficients[5][7],double kindex[], int a[]);
double df0dz(double k,double E_f,double T,double df0dz_integral,double coefficients[5][7],double kindex[], int a[]);


double f(double k,double k_grid[],double E_f,double T,double  coefficients[5][7],double kindex[],double g[], int points, int a[]);
double fH(double k,double k_grid[],double E_f,double T,double  coefficients[5][7],double kindex[],double g[], double h[],int points, int a[]);

double mu_elastic(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],
                  double Ds[],double v[],double nu_elastic[],double g[],int points, int a[]);
double mu_elasticH(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],
                  double Ds[],double v[],double nu_elastic[], int points, int a[]);

double mu_overall(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],double Ds[],
                  double v[],double g[],double nu_el[],int points,int a[]);
double mu_overallH(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],double Ds[],
                  double v[],double g[], double h[], double nu_el[],int points,int a[]);

double mu_po(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],double Ds[],double v[],
             double g_LO[],double g[],double nu_el[],int points, int a[]);
double mu_poH(double k_grid[],double energy[],double E_f,double T,double coefficients[5][7],double kindex[],double Ds[],double v[],
             double g_LO[],double h_LO[],double nu_el[],int points, int a[]);


double J(double k_grid[],double T,double m,double g_th[], double Ds[],double energy[],double v[], int points);
double DOS_value1(double energy, int a);
void read_ispin();


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

extern double lm[10][10], volume1, ion_mass1[5];
extern int ion_numbers1[5], spin_orbit_coupling;

extern double CBM, kcbm[3];
extern double kvbm[3], VBM;
//extern double conduction_band[10000][2];
//------------------------------------ make_band ----------------------------------------
extern double cond_band[2000][2], val_band[2000][2],poly[4000];
extern int count1,count2,countx;    //count1 for conduction band count2 for valence band
extern int count_orbital;
extern int count_DOS_n,count_DOS_p;
extern double E_F,n_e,n_h;
//extern double short_range_band_num[1000][3];
//extern int length_dd, maximum_degree, points ;
//extern double coefficients[10][10], kindex[10];
extern int number;
extern int fitting_1;

//--------------------------------------- analytical_fitting----------------------------------
extern double coefficients[5][7],kindex[4],coefficients_cond[5][7],kindex_cond[4],coefficients_val[5][7],kindex_val[4];
extern int a[2],b[2];


extern int degree1,length_fraction;
extern double fraction[4];

//------------------------------------- inside read_input----------------------------------------------
extern double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];
extern double we[10],De[10];
extern int nfv[10];
extern int variation;
extern double kcbm[3],kvbm[3];
extern int kk;

//extern double k_min0,k_trans0,k_step_fine0,k_step0;

extern int De_ionization,N_cb, N_vb, iterations, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

extern double Ed, c_lattice, rho, k_max, N_dis, omega_LO, omega_TO, E_deformation_n, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

extern double Bfield;

extern string  type;
extern int free_e;
extern int len_T,len_n;
//----------------------------find_cbm_vbm-------------------------------------------------------------------------
extern double ecbm,evbm,kcbm1[3],kvbm1[3];
extern double kcbm[3],kvbm[3],VBM,CBM;
extern int NKPTS,NBVAL,NBTOT;
extern double orbital_decomposedd[1000][4];
extern int ispin;
extern double DOS_n[2000][2],DOS_p[2000][2];
extern int ispin;


//-------------------------------------------------------------------------------------------------------------
#endif // MAIN_H_INCLUDED
