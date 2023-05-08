#include <iostream>
#include<bits/stdc++.h>
#include <cmath>
#include<complex>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include<vector>
using namespace std;
using namespace Eigen;

double rand_num() 
{
	double num;
	num = double(rand()) / double(RAND_MAX);	
	//printf("%lf \n",num);
	return num;
}

void sort(double a[]) 
{
   int i, j, min, n = 10;
   double temp;
   
   for (i = 0; i < n - 1; i++) 
   {
      min = i;
      for (j = i + 1; j < n; j++)
      if (a[j] < a[min])
      min = j;
      temp = a[i];
      a[i] = a[min];
      a[min] = temp;
   }
}

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


int main()
{   
	FILE *fid, *fid1, *fid2, *fid3;
	
	const complex<double> ii(0.0,1.0);
	int datapoints = 150, steps = 5000;
	vector <double> energy, DOS, xx1;
	
	
	char line[1000],dummy;
	double pi = 3.14159265359;
	
	fid = fopen("material_DB_TB_Vogl.dat","r");
	if (fid == NULL)
	{
        	cout<<"material_DB_TB_Vogl.dat is not present. Exit from program";
        	exit(EXIT_FAILURE);
	}

	
	cout<<"Select Material: Enter the number corrresponding to it"<<endl;
	cout<<"3C-SiC :    1"<<endl;
	cout<<"AlAs :      2"<<endl;
	cout<<"AlP:        3"<<endl;
	cout<<"AlSb:       4"<<endl;
	cout<<"C:          5"<<endl;
	cout<<"CdS:        6"<<endl;
	cout<<"CdSe:       7"<<endl;
	cout<<"GaAs:       8"<<endl;
	cout<<"GaP:        9"<<endl;
	cout<<"GaSb:       10"<<endl;
	cout<<"Ge:         11"<<endl;
	cout<<"InP:        12"<<endl;
	cout<<"InAs:       13"<<endl;
	cout<<"InSb:       14"<<endl;
	cout<<"Si:         15"<<endl;
	cout<<"Sn:	    16"<<endl;
	cout<<"ZnSe:	    17"<<endl;
	cout<<"ZnTe:	    18"<<endl;
	cout<<"ZnS:	    19"<<endl;
	cout<<endl;
	
	int j;
	scanf("%d", &j);
	
	double Esa, Epa, Esc, Epc, Esea, Esec, Vss, Vxx, Vxy, Vsapc, Vpasc, Vseapc, Vpasec, b_length;

	for(int i=1;i<=j;i++)
	{
		fgets(line, 1000, (FILE*)fid);
		//cout<<line<<endl;
	}
	
	//cout<<line[0]<<"   "<<line[1]<<"   "<<line[2]<<"   "<<line[3]<<"   "<<line[4]<<"   "<<line[5]<<"   "<<line[6]<<"   "<<line[7]<<"   "<<line[8]<<"   "<<line[9]<<"   "
	//<<line[10]<<"   "<<line[11]<<"   "<<line[12]<<"   "<<line[13]<<"   "<<line[14]<<"   "<<endl;
	fgets(line, 1000, (FILE*)fid);	

	//cout<<line1[0]<<"   "<<line1[1]<<"   "<<line1[2]<<"   "<<line1[3]<<"   "<<line1[4]<<"   "<<line1[5]<<"   "<<line1[6]<<"   "<<line1[7]<<"   "<<line1[8]<<"   "<<line1[9]<<"   "
	//<<line1[10]<<"   "<<line1[11]<<"   "<<line1[12]<<"   "<<line1[13]<<"   "<<line1[14]<<"   "<<endl;

	sscanf(line, "%s  %lf %lf  %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf  %lf %lf", &dummy, &Esa, &Epa, &Esc, &Epc, &Esea, 
	&Esec, &Vss, &Vxx, &Vxy, &Vsapc, &Vpasc, &Vseapc, &Vpasec, &b_length);

	
	cout<<"Esa    =   "<<Esa<<endl;
	cout<<"Epa   =   "<<Epa<<endl;
	cout<<"Esc    =   "<<Esc<<endl;
	cout<<"Epc    =   "<<Epc<<endl;
	cout<<"Esea    =   "<<Esea<<endl;
	cout<<"Esec    =   "<<Esec<<endl;
	cout<<"Vss    =   "<<Vss<<endl;
	cout<<"Vxx    =   "<<Vxx<<endl;
	cout<<"Vxy    =   "<<Vxy<<endl;
	cout<<"Vsapc    =   "<<Vsapc<<endl;
	cout<<"Vpasc    =   "<<Vpasc<<endl;
	cout<<"Vseapc    =   "<<Vseapc<<endl;
	cout<<"Vpasec    =   "<<Vpasec<<endl;
	cout<<"b_length    =   "<<b_length<<endl<<endl;
	getchar();
	
	fclose(fid);
	
	int number;
	double kpoints[10000][3], lattice_constant;
	
	lattice_constant = b_length*4.0/sqrt(3); // unit is Angestron
	lattice_constant = lattice_constant * 1e-10;   // converted from Angestron to m 
	
	fid = fopen("kpoints_file","r");
	if (fid == NULL)
	{
        	cout<<"kpoints_file is not present. Exit from program";
        	exit(EXIT_FAILURE);
	}

	fid1 = fopen("energy_file","w");
	fid2 = fopen("EK_VB.dat","w");
	fid3 = fopen("EK_CB.dat","w");
	
	
	fprintf(fid1,"# kx            ky            kz               e[1]           e[2]           e[3]           e[4]           e[5]           e[6]           e[7]           e[8]           e[9]           e[10] \n");  
	
	
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
	
	double d1[] = {0.25, 0.25, 0.25}, d2[] =  {0.25, -0.25, -0.25}, d3[] = {-0.25, 0.25, -0.25}, d4[] = {-0.25, -0.25, 0.25};
	
	double k1,k2,k3,k4;
	double ev[10];
	complex<double> p1,p2,p3,p4,g0,g1,g2,g3, a(4.0,0);
	
	MatrixXcd h(10,10), ht(10,10),h1(10,10), H(10,10);

	//cout<<"i =  "<<ii<<endl
	for(int i=0;i<number;i++)
	{
		fgets(line, 1000, (FILE*)fid);
		sscanf(line, "%lf %lf  %lf ", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2]);
		
		//cout<<"kpoints = "<<kpoints[i][0]<<"     "<<kpoints[i][1]<<"     "<<kpoints[i][2]<<endl;
		
		
		//kpoints[i][0]  = kpoints[i][0]*(2*pi);   
		//kpoints[i][1]  = kpoints[i][1]*(2*pi);
		//kpoints[i][2]  = kpoints[i][2]*(2*pi); 
		
		k1 = (d1[0]*kpoints[i][0] + d1[1]*kpoints[i][1] + d1[2]*kpoints[i][2])*(2*pi);
		k2 = (d2[0]*kpoints[i][0] + d2[1]*kpoints[i][1] + d2[2]*kpoints[i][2])*(2*pi);
		k3 = (d3[0]*kpoints[i][0] + d3[1]*kpoints[i][1] + d3[2]*kpoints[i][2])*(2*pi);
		k4 = (d4[0]*kpoints[i][0] + d4[1]*kpoints[i][1] + d4[2]*kpoints[i][2])*(2*pi);
		
		
		//getchar();
		
		p1 = exp(ii*k1);
		p2 = exp(ii*k2);
		p3 = exp(ii*k3);
		p4 = exp(ii*k4);
		
				
		g0 = (p1+p2+p3+p4)/a;
		g1 = (p1+p2-p3-p4)/a;
		g2 = (p1-p2+p3-p4)/a;
		g3 = (p1-p2-p3+p4)/a;
				
		
		//MatrixXcd H = MatrixXcd::Zero(10,10);
		
		h << Esa/2, Vss*g0, 0, 0, 0, Vsapc*g1, Vsapc*g2, Vsapc*g3, 0, 0,
		     0, Esc/2, -Vpasc*conj(g1), -Vpasc*conj(g2), -Vpasc*conj(g3), 0 , 0, 0, 0, 0,
		     0, 0, Epa/2, 0, 0, Vxx*g0, Vxy*g3, Vxy*g2, 0, -Vpasec*g1,  
		     0, 0, 0, Epa/2, 0, Vxy*g3, Vxx*g0, Vxy*g1, 0, -Vpasec*g2, 
		     0, 0, 0, 0, Epa/2, Vxy*g2, Vxy*g1, Vxx*g0, 0, -Vpasec*g3,	
		     0, 0, 0, 0, 0, Epc/2, 0, 0, Vseapc*g1, 0, 
		     0, 0, 0, 0, 0, 0, Epc/2, 0, Vseapc*g2, 0, 
		     0, 0, 0, 0, 0, 0, 0, Epc/2, Vseapc*g3, 0, 
		     0, 0, 0, 0, 0, 0, 0, 0, Esea/2, 0,
		     0, 0, 0, 0, 0, 0, 0, 0, 0, Esec/2;
		
		ht = h.transpose();
		h1 = ht.conjugate();
		
		H = h + h1;
		
		  
		//cout<<"k1   =  "<<k1<<"    k2   =  "<<k2<<"    k3   =  "<<k3<<"   k4   =  "<<k4<<endl<<endl;
		
		/*
		cout<<"p1   =   "<<p1<<endl;
		cout<<"p2   =   "<<p2<<endl;
		cout<<"p3   =   "<<p3<<endl;
		cout<<"p4   =   "<<p4<<endl<<endl;

		cout<<"g0   =   "<<g0<<endl;
		cout<<"g1   =   "<<g1<<endl;
		cout<<"g2   =   "<<g2<<endl;
		cout<<"g3   =   "<<g3<<endl<<endl;

		cout<<"h   =   "<<endl<<h<<endl;  
		cout<<"h1   =   "<<endl<<h1<<endl;  
		cout<<"No. of rows in H = "<<H.rows()<<endl;
		cout<<"No. of Coloumns in H = "<<H.cols()<<endl;		
		cout<<"H   =   "<<endl<<H<<endl;     
		//getchar();
		//*/
		
		Eigen::ComplexEigenSolver<MatrixXcd> ces;
		ces.compute(H);    // computes complex eigenvalues and eigenvectors
		
		//cout<<ces<<endl;
		//getchar();
		
		for(int i=0;i<=9;i++)
		{
			ev[i] = (ces.eigenvalues()[i]).real();
			//cout<<"i   =     "<<i<<"  ev[]   =    "<<ev[i]<<endl;
		}
		
		/*
		cout<<"Before sorting "<<endl;
		for(int i=0;i<=9;i++)
		{
			cout<<ev[i]<<"       ";
		}
		cout<<endl;
		*/
		
		sort(ev);
		
		/*
		cout<<"After sorting "<<endl;
		for(int i=0;i<=9;i++)
		{
			cout<<ev[i]<<"       ";
		}
		cout<<endl;
		getchar();
		//*/
		
		
		// save results
		fprintf(fid1,"  %e   %e  %e    %e   %e  %e    %e   %e  %e    %e   %e  %e  %e  \n",
		kpoints[i][0], kpoints[i][1], kpoints[i][2], ev[0], ev[1], ev[2], ev[3], ev[4], ev[5], ev[6], ev[7], ev[8], ev[9]);  
		
		//Valence Band
		fprintf(fid2,"  %e   %e  %e    %e   \n",
		kpoints[i][0]*(2*pi)/(lattice_constant*100), kpoints[i][1]*(2*pi)/(lattice_constant*100), kpoints[i][2]*(2*pi)/(lattice_constant*100), ev[3]); 
		
		// Conduction band
		fprintf(fid3,"  %e   %e  %e    %e   \n",
		kpoints[i][0]*(2*pi)/(lattice_constant*100), kpoints[i][1]*(2*pi)/(lattice_constant*100), kpoints[i][2]*(2*pi)/(lattice_constant*100), ev[4]); 
		
		
		 
		
	}
	
	fclose(fid);
	fclose(fid1);
	fclose(fid2);
	fclose(fid3);
	

	/*		
	for(int l=0;l<N;l++)
	{    
		fprintf(fid1,"  %e               %e               %e               %e\n",ep1[l], en1[l], ep2[l], en2[l]);        
	}
	fclose(fid1);
	*/
	
	//-------------------------------------------------------------------------------------------------------------------------------------------------------------
	/*
	// Calculation of DOS
	int N = 0;
	double kx,ky,kz;
	
	vector <double> xl =linspace(-1, 1,  datapoints );
	vector <double> yl =linspace(-1, 1,  datapoints );
	vector <double> zl =linspace(-1, 1,  datapoints );

	cout<<"Loop for DOS running "<<endl;
	//for (int i=1;i<=datapoints;i++)
	int iteration=1;
	
	for(int i1=0;i1<datapoints;i1++)
	{	
		for(int i2=0;i2<datapoints;i2++)
		{
			for(int i3=0;i3<datapoints;i3++)
			{
				/*
				cout<<"i = "<<i<<endl;
					
				kx = rand_num();
				ky = rand_num();
				kz = rand_num();
				
				if(rand_num() < 0.5)
					kx = -kx;
					
				if(rand_num()<0.5)
					ky = -ky;

				if(rand_num()<0.5)
					kz = -kz;
				*/
				
				/*
				cout<<"i1   =    "<<i1<<"   i2   =    "<<i2<<"   i3   =    "<<i3<<endl;
				//cout<<"x[]   =    "<<xl[i1]<<"   yl   =    "<<yl[i2]<<"   zl[]   =    "<<zl[i3]<<endl;
				cout<<"Iteration     =   "<<iteration<<endl;
				//getchar();
				iteration = iteration+1;
				
				kx = xl[i1];
				ky = yl[i2];
				kz = zl[i3];
				
				
				if( abs(kx)+abs(ky)+abs(kz) < 1.5 )
				{	
					N = N + 1;
					
					k1 = (d1[0]*kx + d1[1]*ky + d1[2]*kz)*(2*pi);
					k2 = (d2[0]*kx + d2[1]*ky + d2[2]*kz)*(2*pi);
					k3 = (d3[0]*kx + d3[1]*ky + d3[2]*kz)*(2*pi);
					k4 = (d4[0]*kx + d4[1]*ky + d4[2]*kz)*(2*pi);
					
					p1 = exp(ii*k1);
					p2 = exp(ii*k2);
					p3 = exp(ii*k3);
					p4 = exp(ii*k4);
							
					g0 = (p1+p2+p3+p4)/a;
					g1 = (p1+p2-p3-p4)/a;
					g2 = (p1-p2+p3-p4)/a;
					g3 = (p1-p2-p3+p4)/a;
							
							
					h << Esa/2, Vss*g0, 0, 0, 0, Vsapc*g1, Vsapc*g2, Vsapc*g3, 0, 0,
					     0, Esc/2, -Vpasc*conj(g1), -Vpasc*conj(g2), -Vpasc*conj(g3), 0 , 0, 0, 0, 0,
					     0, 0, Epa/2, 0, 0, Vxx*g0, Vxy*g3, Vxy*g2, 0, -Vpasec*g1,  
					     0, 0, 0, Epa/2, 0, Vxy*g3, Vxx*g0, Vxy*g1, 0, -Vpasec*g2, 
					     0, 0, 0, 0, Epa/2, Vxy*g2, Vxy*g1, Vxx*g0, 0, -Vpasec*g3,	
					     0, 0, 0, 0, 0, Epc/2, 0, 0, Vseapc*g1, 0, 

					     0, 0, 0, 0, 0, 0, Epc/2, 0, Vseapc*g2, 0, 
					     0, 0, 0, 0, 0, 0, 0, Epc/2, Vseapc*g3, 0, 
					     0, 0, 0, 0, 0, 0, 0, 0, Esea/2, 0,
					     0, 0, 0, 0, 0, 0, 0, 0, 0, Esec/2;
					
					ht = h.transpose();
					h1 = ht.conjugate();
					
					H = h + h1;
					
					  
					//cout<<"k1   =  "<<k1<<"    k2   =  "<<k2<<"    k3   =  "<<k3<<"   k4   =  "<<k4<<endl<<endl;
							
					Eigen::ComplexEigenSolver<MatrixXcd> ces;
					ces.compute(H);    // computes complex eigenvalues and eigenvectors

					for(int i=0;i<=9;i++)
					{
						energy.push_back((ces.eigenvalues()[i]).real());
						//cout<<"i   =     "<<i<<"  ev[]   =    "<<ev[i]<<endl;
					}
				}
			}
		}
	}
		
	cout<<"Loop for DOS completed "<<endl;
	cout<<"No. of cell N = "<<N<<endl;
	
	sort(energy.begin(),energy.end());
	
	double energy_min = energy[0]; 
	
	double energy_max = energy[energy.size()-1];		
	
	double de = (energy_max-energy_min)/steps;

	int kk=0;
	
	// volume of primitive unit cell   unit is m^3
	double volume = lattice_constant*lattice_constant*lattice_constant/4;
	
	cout<<"Second Loop for DOS started "<<endl;
	for(int i=0; i<energy.size();i++) 
	{	
		if(energy[i] <  energy_min)
		{
			kk=kk+1;  
		}
		else 
		{
			DOS.push_back(kk);	
			kk=0;
			xx1.push_back(energy_min);
			energy_min = de + energy_min;	
		}
	}
	cout<<"Second Loop for DOS completed "<<endl;
			
	for(int i=0; i < DOS.size(); i++)
	{  

		// DOS (per eV per cm^2) multiplied with 1e-6 to convert from per m^3 to per cm^3
		DOS[i] = (2*DOS[i]*1e-6)/(de*N*volume);    
		// multiplied with 2 to account for spin
		// unit is per eV per cm^3
	}
	
	fid = fopen("DOS.dat","w");
	fprintf(fid,"# Energy (eV)		DOS  (per eV per cm^3 )\n");  
	
	for(int i=0; i < DOS.size(); i++)
	{  
		fprintf(fid,"%e	%e \n", xx1[i],DOS[i]);  
	}
	fclose(fid);	
	*/
	
	
	return 0;
}
