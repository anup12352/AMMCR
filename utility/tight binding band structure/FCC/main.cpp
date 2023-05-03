#include <iostream>
#include<bits/stdc++.h>
#include <cmath>
#include<complex>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include<vector>
using namespace std;
using namespace Eigen;

/*
double rand_num() 
{
	double num;
	num = double(rand()) / double(RAND_MAX);	
	//printf("%lf \n",num);
	return num;
}
*/

void sort(double a[], int n) 
{
   int i, j, min;
   double temp;
   
   for (i = 0; i < n - 1; i++) 
   {
      min = i;
      for (j = i + 1; j < n; j++)
      {
	      if (a[j] < a[min])
	      {
		      min = j;	
		      temp = a[i];
		      a[i] = a[min];
		      a[min] = temp;
		}
	}
   }
}
/*
vector<double> linspace(double min, double max, int n)
{
	vector<double> result;
	int iterator = 0;
	for (int i = 0; i <= n-2; i++)	
	{
		double temp = min + i*(max-min)/(floor((double)n) - 1);
		result.insert(result.begin() + iterator, temp);
		iterator += 1;
	}

	result.insert(result.begin() + iterator, max);
	return result;
}			/// function to create a linear dimention subspace grid
*/

int main()
{   
	FILE *fid, *fid1, *fid2, *fid3;
	
	const complex<double> ii(0.0,1.0);
	//int datapoints = 150, steps = 5000;

	
	
	char line[1000];
	double m_e, hbar, pi, d, e, l, constant;
	double Eta_sss = 0, Eta_sps = 0, Eta_pps = 0, Eta_ppp = 0, Eta_sds = 0, Eta_pds = 0, Eta_pdp = 0, Eta_dds = 0, Eta_ddp = 0, Eta_ddd = 0, lattice_constant = 0;
	double Vsss = 0, Vsps = 0, Vpps = 0, Vppp = 0, Ess = 0, Esp = 0, Epp = 0, Epxpy = 0, Epxpz = 0;
	double Es0, Ep0, Ed0;
	
	pi = 3.14159265359;
	hbar = 6.62607015e-34;
	m_e = 9.109e-31;
	e = 1.60217657e-19;
	
	fid = fopen("input.dat","r");
	if (fid == NULL)
	{
        	cout<<"input.dat is not present. Exit from program";
        	exit(EXIT_FAILURE);
	}

	
	
        cout.setf(ios::scientific);
        cout<<"Data from input.dat file"<<endl;
        cout<< " -------------- "<<endl<<endl;

        string str;
        string ss;
        
        ifstream in("input.dat");

        while(!in.eof())
        {
		in>> str;
		
		// to check comments(if a line is not part of input then start tat line with #)
		if(str=="#")
		{
		  getline(in,ss);
		  cout<< "--COMMENT--" <<ss <<endl<<endl;
		}
		
		if(str=="Eta_sss")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_sss; 
			  
		  cout<< "Eta_sss =  " <<Eta_sss<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_sps")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_sps; 
			  
		  cout<< "Eta_sps =  " <<Eta_sps<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_pps")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_pps; 
			  
		  cout<< "Eta_pps =  " <<Eta_pps<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_ppp")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_ppp; 

		  cout<< "Eta_ppp =  " <<Eta_ppp<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_sds")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_sds; 
			  
		  cout<< "Eta_sds =  " <<Eta_sds<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_pds")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_pds; 
			  
		  cout<< "Eta_pds =  " <<Eta_pds<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_pdp")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_pdp; 
			  
		  cout<< "Eta_pdp =  " <<Eta_pdp<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_dds")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_dds; 
			  
		  cout<< "Eta_dds =  " <<Eta_dds<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_ddp")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_ddp; 
			  
		  cout<< "Eta_ddp =  " <<Eta_ddp<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Eta_ddd")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Eta_ddd; 
			  
		  cout<< "Eta_ddd =  " <<Eta_ddd<<endl;
		  cout<<endl<<endl;
		}

		if(str=="lattice_constant")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> lattice_constant; 
			  
		  cout<<"lattice_constant =  " <<lattice_constant<<"  A "<<endl;
		  lattice_constant = lattice_constant *1e-10;  // converted from angestron to meter
		  
		  cout<<endl<<endl;
		}
		
		if(str=="Es0")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Es0; 
			  
		  cout<< "Es0 =  " <<Es0<<" eV"<<endl;
		  cout<<endl<<endl;
		}
		
		if(str=="Ep0")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Ep0; 
			  
		  cout<< "Ep0 =  " <<Ep0<<" eV"<<endl;
		  cout<<endl<<endl;
		}
		
		if(str=="Ed0")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Ed0; 
			  
		  cout<< "Ed0 =  " <<Ed0<<" eV"<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Vsss")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Vsss; 
			  
		  cout<< "Vsss =  " <<Vsss<<" eV"<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Vsps")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Vsps; 
			  
		  cout<< "Vsps =  " <<Vsps<<" eV"<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Vpps")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Vpps; 
			  
		  cout<< "Vpps =  " <<Vpps<<" eV"<<endl;
		  cout<<endl<<endl;
		}

		if(str=="Vppp")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> Vppp; 
			  
		  cout<< "Vppp =  " <<Vppp<<" eV"<<endl;
		  cout<<endl<<endl;
		}
	}
	
	fclose(fid);

	d = lattice_constant/sqrt(2);
	
	cout<<"lattice_constant =  " <<lattice_constant<<"  meter "<<endl;
	cout<<"d = "<<d<<"  meter "<<endl;
	
	
	constant = 7.62e-20;
	
	if(Vsss==0)	
		Vsss = Eta_sss*constant/(d*d);   
	
	if(Vsps==0)
		Vsps = Eta_sps*constant/(d*d);
	
	if(Vpps==0)
		Vpps = Eta_pps*constant/(d*d);
	
	if(Vppp==0)	
		Vppp = Eta_ppp*constant/(d*d);

	
	//constant = hbar*hbar/(m_e*e);
	
	//cout<<"constant = "<<constant<<endl;
	// divided by e to convert joule to eV
	//Vsss = Eta_sss*constant/(e*d*d);   
	//Vsps = Eta_sps*constant/(e*d*d);
	//Vpps = Eta_pps*constant/(e*d*d);
	//Vppp = Eta_ppp*constant/(e*d*d);

	cout<<"Vsss = "<<Vsss<<endl;
	cout<<"Vsps = "<<Vsps<<endl;
	cout<<"Vpps = "<<Vpps<<endl;
	cout<<"Vppp = "<<Vppp<<endl;

	//cout<<"Press key to continue"<<endl;
	//getchar();


	//l = 1/sqrt(2);
	
	//cout<<"l = "<<l<<endl;
	//cout<<"Press key to continue"<<endl;
	//getchar();

	/*
	Ess = Vsss;
	Esp = l*Vsps;
	Epp = l*l*Vpps + (1-l*l) *Vppp;
	Epxpy = l*l*Vpps - l*l*Vppp;
	Epxpz = l*l*Vpps - l*l*Vppp;
	
	
	cout<<"Ess    =   "<<Ess<<"  eV"<<endl;
	cout<<"Esp   =   "<<Esp<<"  eV"<<endl;
	cout<<"Epp   =   "<<Epp<<"  eV"<<endl;
	cout<<"Epxpy   =   "<<Epxpy<<"  eV"<<endl;
	cout<<"Epxpz   =   "<<Epxpz<<"  eV"<<endl;
	*/
	
	cout<<"lattice_constant    =   "<<lattice_constant<<" meter "<<endl<<endl;
	//getchar();
	
	double kpoints[5000][3]={0},kx,ky,kz;
	int number;
	
	MatrixXcd h(4,4);	
	
	fid = fopen("kpoints_file","r");
	if (fid == NULL)
	{
        	cout<<"kpoints_file is not present. Exit from program";
        	exit(EXIT_FAILURE);
	}

	fid1 = fopen("energy_file","w");
	fid2 = fopen("EK_VB.dat","w");
	fid3 = fopen("EK_CB.dat","w");
	
	
	fprintf(fid1,"# Sr. No.   kx            ky            kz               e[1]           e[2]           e[3]           e[4] \n");  
	
	
	// valence band
	fprintf(fid2,"# kx     	   ky     	kz (1/cm)   	 Energy (eV) \n");  
	
	// conduction band
	fprintf(fid3,"# kx     	   ky     	kz (1/cm)   	 Energy (eV) \n");

	// read kpoint_file, pass three starting lines
	fgets(line, 1000, (FILE*)fid);
	fgets(line, 1000, (FILE*)fid);
	sscanf(line, "%d ",&number);
	cout<<"Total kpoints = "<<number<<endl;
	fgets(line, 1000, (FILE*)fid);
	
	double ev[3];
	//const complex<double> ii(0.0,1.0);
	
	for(int i=0;i<number;i++)
	{
		fgets(line, 1000, (FILE*)fid);
		sscanf(line, "%lf %lf  %lf ", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2]);
		
		//cout<<"i   =   "<<i<<"    kpoints[i][0] = "<<kpoints[i][0]<<"   kpoints[i][1] = "<<kpoints[i][1]<<"  kpoints[i][2] = "<<kpoints[i][2]<<endl;
		//getchar();
		
		kx = kpoints[i][0]*2*pi;
		ky = kpoints[i][1]*2*pi;
		kz = kpoints[i][2]*2*pi;
		
		
		h(0,0) = Es0 + 4*Vsss*(cos(kx/2.0)*cos(ky/2.0)+cos(kx/2.0)*cos(kz/2.0)+cos(ky/2.0)*cos(kz/2.0));
		h(1,1) = Ep0 + 2*Vpps*(cos(kx/2.0)*cos(ky/2.0)+cos(kx/2.0)*cos(kz/2.0));
		h(2,2) = Ep0 + 2*Vpps*(cos(ky/2.0)*cos(kx/2.0)+cos(ky/2.0)*cos(kz/2.0));
		h(3,3) = Ep0 + 2*Vpps*(cos(kz/2.0)*cos(kx/2.0)+cos(kz/2.0)*cos(ky/2.0));
		
		
		h(0,1) = 2*sqrt(2)*ii*(sin(kx/2.0)*cos(ky/2.0)+sin(kx/2.0)*cos(kz/2.0))*Vsps;
		h(1,0) = conj(h(0,1));
		//h(1,0) = -1*2*sqrt(2)*ii*(sin(kx/2.0)*cos(ky/2.0)+sin(kx/2.0)*cos(kz/2.0))*Vsps;

		h(0,2) = 2*sqrt(2)*ii*(sin(ky/2.0)*cos(kx/2.0)+sin(ky/2.0)*cos(kz/2.0))*Vsps;
		h(2,0) = conj(h(0,2));

		h(0,3) = 2*sqrt(2)*ii*(sin(kz/2.0)*cos(kx/2.0)+sin(kz/2.0)*cos(ky/2.0))*Vsps;
		h(3,0) = conj(h(0,3));;
		
		h(1,2) = -2*(Vpps-Vppp)*sin(kx/2.0)*sin(ky/2.0);
		h(2,1) = h(1,2);

		h(1,3) = -2*(Vpps-Vppp)*sin(kx/2.0)*sin(kz/2.0);
		h(3,1) = h(1,3);
		
		h(2,3) = -2*(Vpps-Vppp)*sin(ky/2.0)*sin(kz/2.0);
		h(3,2) = h(2,3);
		//*/
		
		Eigen::ComplexEigenSolver<MatrixXcd> ces;
		ces.compute(h);    // computes complex eigenvalues and eigenvectors
		
		//cout<<ces<<endl;
		//getchar();
		
		for(int j=0;j<=3;j++)
		{
			ev[j] = (ces.eigenvalues()[j]).real();
			//cout<<"j   =     "<<j<<"  ev[]   =    "<<ev[j]<<endl;
		}
		
		// 4 is size of matrix or number of bands
		sort(ev,4);

		fprintf(fid1,"  %d    %e   %e  %e    %e   %e    %e    %e \n", i, kpoints[i][0], kpoints[i][1], kpoints[i][2], ev[0], ev[1], ev[2], ev[3]);  
		
		//Valence Band
		fprintf(fid2,"  %e   %e  %e    %e   \n", 
		kpoints[i][0]*(2*pi)/(lattice_constant*100), kpoints[i][1]*(2*pi)/(lattice_constant*100), kpoints[i][2]*(2*pi)/(lattice_constant*100), ev[0]); 
		
		// Conduction band
		fprintf(fid3,"  %e   %e  %e    %e   \n", 
		kpoints[i][0]*(2*pi)/(lattice_constant*100), kpoints[i][1]*(2*pi)/(lattice_constant*100), kpoints[i][2]*(2*pi)/(lattice_constant*100), ev[1]); 

	}		
	
	fclose(fid);
	fclose(fid1);
	fclose(fid2);
	fclose(fid3);
	//cout<<"Reached here"<<endl;	
	
	cout<<"Band structure saved"<<endl;
		
	exit(0);
	
	return 0;
	//cout<<"Reached here"<<endl;
}

