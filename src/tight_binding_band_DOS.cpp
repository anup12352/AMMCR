#include"main.h"


vector<double> linspace_vector(double min, double max, int n)
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

void tight_binding_band_DOS()
{	


	//cout<<"Inside tight binding dos "<<endl;
	// nearest neighbour hopping parameter unit eV
	double t0 = 3.033 * 2.46 * 2.46 / (lattice_constant * lattice_constant);  // unit is eV
				//  lattice_constant unit is in Angestron

	double a0 = lattice_constant*1e-10;   // unit is m
	// a0 is lattice_constant in meter

	//cout<<"lattice_constant    =   "<<lattice_constant<<"   unit is angestron "<<endl;
	//cout<<"a0    =   "<<a0<<"   unit is meter "<<endl;
	
	// bond_length = distance between two nearest atoms
	// relationship between lattic_constant and bond_length => 

	vector <double> e_p;
	vector <double> e_n;
	double kx,ky,kz;
	
// ------------------------------------------------- DOS calculation -----------------------------------------------------------				

	FILE *fid1;
	int dataPoints=2000, steps=500;
	
	// a0 unit is meter  // a0 is lattice_constant in meter
	double gM=(2*pi/(sqrt(3)*a0));
	double mK=(2*pi/(3*a0));
	double kG=(4*pi/(3*a0));

	double temp_S[2];
	
	double area;
	
	area = 3*sqrt(3)*(a0*a0)/(3*2);
	
	//cout<<"Area of unit cell = "<<area<<"   m^2"<<endl;

	//double k;
	vector <double> xl =linspace_vector(-gM, gM,  dataPoints*2 );
	vector <double> yl =linspace_vector(-kG, kG,  dataPoints*2 );
	
	//cout << "Size of xl " << xl.size()<<endl;
	//cout << "Size of yl " << yl.size()<<endl;
	//getchar();

	/*
	for(int i=0;i<=dataPoints*2-1;i++)
	{
	    cout << "i  =  "<<i<<"  xl[i]  =  " << xl[i]<<endl;
	    cout << "i  =  "<<i<<"  yl[i]  =  " << yl[i]<<endl<<endl;
	
	}
	getchar();
	*/
		
	double cell_no = 0;
	
	for (int i=0; i <2*dataPoints;i++)
	{

	    for (int j=0; j<2*dataPoints;j++)
	    { 
	    			
	    //cout << "i  =  "<<i<<"  xl[i]  =  " << xl[i]<<endl;
	    //cout << "j  =  "<<i<<"  yl[j]  =  " << yl[j]<<endl<<endl;

		if( (xl[i]>0 && yl[j] < xl[i]/abs(sqrt(3))) && (yl[j]>0.00))
		{
				cell_no = cell_no + 1;
				//n=xl[i];   p=yl[j];  
				
				kx = xl[i];  // unit is 1/m
				ky = yl[j];  // unit is 1/m

				temp_S[0] = t0*sqrt(1 + 4 * cos((kx*sqrt(3)*a0)/2.0) * cos(a0*ky/2) + 4 * cos(a0*ky/2)* cos(a0*ky/2) );

				temp_S[1] = -1*temp_S[0];
				
			//cout<<" i    =   "<<i<<"  j    =   "<<j<<"   temp_S[0]  = "<<temp_S[0]<<"   temp_S[1] = "<<temp_S[1]<<endl;
		    	//cout<<"kx   = "<<kx<<"   ky   =  "<<ky<<"  unit 1/m "<<endl;
		    	//getchar();
		    	
		    	
		    	e_p.push_back((temp_S[1]));
		    	e_n.push_back((temp_S[0]));		
		    	
		    	
		}
	    }
	}	

	//cout<<"  cell_no = "<<cell_no<<endl;
	//getchar();

	sort(e_p.begin(),e_p.end());
	sort(e_n.begin(),e_n.end());
	
	/*
	fid1 = fopen("e_p_sort.dat","w");
	fprintf(fid1,"#  index       Energy (eV)  \n");

	for(int i=0;i < e_p.size()  ;i++)
		fprintf(fid1,"  %d       %e  \n ",i, e_p[i]);
		
	fclose(fid1);
	
	
	fid1 = fopen("e_n_sort.dat","w");
	fprintf(fid1,"#  index       Energy (eV)  \n");

	for(int i=0;i < e_n.size()  ;i++)
		fprintf(fid1,"  %d       %e  \n ",i, e_n[i]);
		
	fclose(fid1);
	//*/

	/*
	cout<<"e_n.size() =  "<<e_n.size()<<endl;
	cout<<"e_p.size() =  "<<e_p.size()<<endl;
	//*/
		
	double P_min = e_p[0]; 
	
	double P_max = e_p[e_p.size()-1];		
	
	double N_min =e_n[0];	
	double N_max = e_n[e_n.size()-1];	
	
	/*
	cout<<"P_max   =   "<<P_max<<endl;
	cout<<"P_min   =   "<<P_min<<endl;
	cout<<"N_max   =   "<<N_max<<endl;
	cout<<"N_min   =   "<<N_min<<endl;
	getchar();
	//*/
	
	/*
	for(int i=0;i<e_n.size();i++)
	{
		cout<<" i    =   "<<i<<"   e_n[i]  = "<<e_n[i]<<"   e_p[i]  = "<<e_p[i]<<endl;
		getchar();
	}
	//*/
	
	/*
	double sim= 0;
	sim = accumulate(e_n.begin(),   e_n.end(), sim );
	cout<<"sim   = "<<sim<<endl;

	//double sem = accumulate(e_p.rbegin() ,  e_p.rend() , sem );  
	//cout<<"sem   = "<<sem<<endl;
	
	//getchar();
	//*/
	
	double dp = (P_max-P_min)/steps;	// delta E for valence band 	
	double de = (N_max-N_min)/steps;  	// delta E for conduction band 
	
	/*
	cout<<"de =  "<<de<<"   eV  "<<endl;
	cout<<"dp =  "<<dp<<"   eV  "<<endl;
	getchar();
	//*/

	int k=1;
	for(int i=0; i<e_p.size();i++) 
	{	
		if(e_p[i]<  P_min)
		{
			k=k+1;  
		}
		else 
		{
			e_pp.push_back(k);	
			k=1;
			xx1.push_back(P_min);
			P_min = dp+P_min;	
		}
	}
	k=1;



	for(int i=0; i<e_n.size();i++)
	{		
		if(e_n[i]<N_min)
		{
			k=k+1; 
		}
		else 
		{
			e_nn.push_back(k);
			k=1;
			N_min = de+ N_min;
			xx1.push_back(N_min);
		}
	}
	
	
	
//------------------------------------------------------------ save DOS ------------------------------------------------------------------------------------------------------	 
	/*	
	fid1 = fopen("valence_band_DOS.dat","w");
	fprintf(fid1,"#  Energy (eV)  DOS (per eV per atom)  DOS (per eV per unit cell)     DOS (per eV per m^2)  \n");

	for(int i=0;i < e_pp.size()  ;i++)
	{  
		fprintf(fid1,"  %e        %e        %e           %e \n ",xx1[i], 2*e_pp[i]/(2*dp*cell_no), 2*e_pp[i]/(dp*cell_no),   2*e_pp[i]/(dp*cell_no*area));
		// multiplied with 2 to account for spin and divide with 2 to calculate DOS in units of per eV per atom (since there are two atoms in unit cell)
		//fprintf(fid2,"  %e   %e  \n ",-xx1[i], (steps*e_pp[i])/(2*sim));  
		
	}   
	
	// e_nn[] and e_pp[] contains DOS with energy for CB and VB
	
	fclose(fid1);
	
	int nn = e_pp.size();
	
	fid1 = fopen("conduction_band_DOS.dat","w");
	fprintf(fid1,"#  Energy (eV)  DOS (per eV per atom)   DOS (per eV per unit cell)       DOS (per eV per m^2)\n ");
	  
	for(int i=0;i < e_nn.size()  ;i++)
	{
		fprintf(fid1,"  %e           %e         %e        %e \n ", xx1[i+nn], 2*e_nn[i]/(2*de*cell_no), 2*e_nn[i]/(de*cell_no),   2*e_nn[i]/(de*cell_no*area));
		// multiplied with 2 to account for spin and divide with 2 to calculate DOS in units of per eV per atom (since there are two atoms in unit cell)
		//fprintf(fid1,"  %e   %e  \n ",xx1[i], steps*(e_nn[i])/(2*sim));   
		
	}
	
	fclose(fid1);
	//*/

//------------------------------------------------------------ save DOS ------------------------------------------------------------------------------------------------------



	for(int i=0;i < e_pp.size()  ;i++)
	{  

		// DOS (per eV per cm^2) multiplied with 1e-4 to convert from per m^2 to per cm^2
		e_pp[i] = 2*e_pp[i]*1e-4/(dp*cell_no*area);    
	}
	
	for(int i=0;i < e_nn.size()  ;i++)
	{  
		// DOS (per eV per cm^2) multiplied with 1e-4 to convert from per m^2 to per cm^2
		e_nn[i] = 2*e_nn[i]*1e-4/(de*cell_no*area); 
	}
	

}

