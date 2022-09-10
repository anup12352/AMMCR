#include"main.h"


void tight_binding_band_structure()
{	


	// bond_length = distance between two nearest atoms
	// relationship between lattic_constant and bond_length => 
	
	
//------------------------------------ generating kpoints ---------------------------------------------- 

	double center_kpoint[3] = {0,0,0};
	
	double factor[3]; 
	
	/*
	factor[0] = 9;
	factor[1] = 9; 
	factor[2] = 10;

	double step[3]={0.001, 0.01, 0.02};
	*/
	
	factor[0] = 9;
	factor[1] = 9; 
	factor[2] = 10;
	factor[3] = 10;

	double step[4]={0.0001, 0.001, 0.01, 0.05};
	
	// center point lies at 
	// K point  =    (2*pi/(3*a), 2*pi/(3*sqrt(3)*a)) (here a= bond length, distance between two atoms) 
	// or K point =  (2*pi/(a1*sqrt(3)), 2*pi/(3*a1)) (here a1= lattice constant)
	// relationship between lattice constant (a1) and bond length (a) is a1 = a*sqrt(3)) 
	
	// wrt 2*pi/a1   -- a1 = lattice constant

	center_kpoint[0] = 0.577350269;   //
	center_kpoint[1] = 0.333333334;
	
	double kx,ky,kz;

	kz = 0;
	
	double t = 2*pi*10/(lattice_constant);  // unit is now 1/nm
	//cout<<"lattice_constant   =  "<<lattice_constant<<"   unit Angestron"<<endl;
	//cout<<"t  =  "<<t<<"  unit 1/nm  "<<endl;
	//getchar();
	
	kpoints[0][0] = center_kpoint[0]*t;
	kpoints[0][1] = center_kpoint[1]*t;
	kpoints[0][2] = kz;

	int j=1;
	
	
	for (int i=0;i<4;i++)
	{
	       	kx=center_kpoint[0]-step[i]*factor[i];

		while (kx<=center_kpoint[0]+step[i]*factor[i]+0.00001)   
		// 0.00001 is added because in place of zero it taking 1.38778e-17 in c++
		{
		    ky=center_kpoint[1]-step[i]*factor[i];
		    while(ky<=center_kpoint[1]+step[i]*factor[i]+0.00001)
		    {
			    if (pow((kx-center_kpoint[0]),2)+pow((ky-center_kpoint[1]),2)>0.0000000000001)
			    {
			    


				kpoints[j][0] = kx*t;
				kpoints[j][1] = ky*t;
				kpoints[j][2] = kz*t;
				// unit is 1/nm

				/*
		    		cout<<"i   =   "<<i<<" kx =  "<<kx<<"  normalized"<<endl;
				cout<<"j   =   "<<j<<" ky =  "<<ky<<"  "<<endl;
				cout<<"j   =   "<<j<<" kz =  "<<kz<<"  "<<endl;
		    		cout<<"j   =   "<<i<<"   kpoints[j][0] =  "<<kpoints[j][0]<<"  unit 1/nm"<<endl;
				cout<<"j   =   "<<j<<"   kpoints[j][1] =  "<<kpoints[j][1]<<"  "<<endl;
				cout<<"j   =   "<<j<<"   kz =  "<<kz<<"  "<<endl;
				getchar();
				*/
					
				j = j+1;
			    }
			ky = ky +step[i];
		    }
		    kx = kx + step[i];
		}
	}
	
	NKPTS= j;
	/*
	if(NKPTS!=1161)
	{
		cout<<"Error in NKPTS calculation"<<endl;
		exit(EXIT_FAILURE);
	}
	*/
	//cout<<"NKPTS =  "<<NKPTS<<endl;
	
//------------------------------------ generating kpoints completed ----------------------------------------------



//-------------------------- calculating energy at different kpoints----------------------------------------------
//      Relationship used 
//	E = t*sqrt(1 + 4 cos((kx*sqrt(3)*a0)/2) cos(a0*ky/2) + 4 * cos(a0*ky/2)*cos(a0*ky/2))
	
	// nearest neighbour hopping parameter unit eV
	double t0 = 3.033 * 2.46 * 2.46 / (lattice_constant * lattice_constant);  // unit is eV
				//  lattice_constant unit is in Angestron
				
	double a0 = lattice_constant*1e-10;   // unit is m
	// a0 is lattice_constant in meter
	
	for (int i=0;i<NKPTS;i++)
	{
		//cout<<"i   =   "<<i<<" kpoints[i][0] =  "<<kpoints[i][0]<<"  unit 1/nm"<<endl;
		//cout<<"i   =   "<<i<<" kpoints[i][1] =  "<<kpoints[i][1]<<"  unit 1/nm"<<endl;
		
		kx = kpoints[i][0]*1e9;  // unit is 1/m
		ky = kpoints[i][1]*1e9;  // unit is 1/m
		
		energy_n[i] = t0*sqrt(1 + 4 * cos((kx*sqrt(3)*a0)/2.0) * cos(a0*ky/2) + 4 * cos(a0*ky/2)* cos(a0*ky/2) );	
		energy_p[i] = -1*energy_n[i];
		
		//cout<<"i =   "<<i<<"  energy_n[] =  "<<energy_n[i]<<endl;
		//cout<<"i =   "<<i<<"  energy_p[] =  "<<energy_p[i]<<endl;
		//getchar();
		
	}

//--------------------------------------------------- save band data ------------------------------------------------------------	
	/*
	FILE *fid1, *fid2;
	fid1 = fopen("EK_CB.dat","w");
	fid2 = fopen("EK_VB.dat","w");
	
	fprintf(fid1,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");
	fprintf(fid2,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");

	for (int i = 0; i < NKPTS; i++)
	{
		fprintf(fid1,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energy_n[i]);
		fprintf(fid2,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energy_p[i]);
	}
	
	fclose(fid1);
	fclose(fid2);
	//*/
//--------------------------------------------------- saved band data ------------------------------------------------------------


}

