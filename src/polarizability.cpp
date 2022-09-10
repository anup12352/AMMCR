#include"main.h"

void polarizability(double T, int ii)
{
	double de, dk2, En[points]={0}, k[points]={0}, v[points]={0}, mu;
	double const1[points]={0};
	double Ef;  // fermi level at T=0
	
	de = 0.001;
	int gv=1;
	
	//double Ef = find_fermi(n_array[ii], 0, 0);   // first zero for 0 temp and second for T_loop
	// fermi level at T= 0 K
	
	for(int i=0;i<points;i++)
	{	
		En[i] = energy_n[i];
	        k[i] = k_grid[i]*1e9;   // converted from 1/nm to 1/m					
    		v[i] = v_n[i]*1e-2;	      // converted from cm/s to m/s	
    		//const1[i] = Ds_n[i]/(e*volume1*1e-6);    
    		// DOS at energy E;  1e-6 is multiplied in denominator to convert volume from cm^3 to m^3  
    		// e is multiplied to convert from per eV to per J
    		const1[i] = gv*k[i]/(pi*h_bar*e*v[i]) ; 
	}
	
	//cout<<"At E = 0 const1 = "<<const1[0]<<endl;
	
	dk2 = 2*k[points-1]/limit6;
	
	// loop for q change in wave vector variation from 0 to 2k
	for(int i=0;i<=limit6;i++)
	{
		q[i] = dk2*i;   // change in wave vector
		pz[i] = 0;
	}

	double pze, term, denom, q1, aa;
	
	mu = efef_n;
	//cout<<endl<<"Temp = "<<T<<endl;
	//cout<<"mu = "<<mu<<endl;
		
		
	for (int j=1;j<=limit6;j++)     // loop for q change in wave vector variation from 0 to 2k
	{
		q1 = q[j];
		
		//cout<<"j = "<<j<<endl;
		//cout<<"q1 = "<<q1<<endl;
		//getchar();
		
		for(int i=1;i<points;i++)  // loop for integration wrt E
		{
							
			if(q1 > 2*k[i])
				pze = const1[i]*(1 - sqrt(1 - (2*k[i]/q1)*(2*k[i]/q1) ));
			else
				pze = const1[i];
				
			de = En[i] - En[i-1];
			aa = cosh(((mu-En[i])*e)/(2*k_B*e*T));
			denom = 4*k_B*e*T*aa*aa;    
			
			term = pze/denom;
			
			
			pz[j] = pz[j] + e*de*term; // polarizibility
			
			/*
			cout<<"i = "<<i<<endl;
			cout<<"k = "<<k<<endl;
			cout<<"numer - pze = "<<pze<<endl;
			cout<<"aa  =   "<<aa<<endl;
			cout<<"denom = "<<denom<<endl;
			cout<<"term = "<<term<<endl;
			cout<<"pz[j]  =   "<<pz[j]<<endl;  
			getchar();
			//*/
			
		}  
		//cout<<endl<<endl<<"final pz for j = "<<j<<endl;
		//cout<<"q1 = "<<q1<<endl;
		//cout<<"pz[j] =   "<<pz[j]<<endl;  
		//cout<<"Press key to continue "<<endl;
		//getchar();
	}
	
	pz[0] = pz[1];   
	
	//kf = sqrt(2*m_eff*e*abs(mu))/(q*h_bar);

	//pz_q0 = abs(-gs*gv*m_eff/(2*pi*e*h_bar*e*h_bar));
	
	double pz_q0 = const1[0];
	
	cout<<" pz_q0 = "<<pz_q0<<endl;					
	
	
	FILE *fid1;
	fid1 = fopen("pz.dat","w");
	
	fprintf(fid1,"# Change in wave vector \n");
	    		
	for (int j=0;j<=limit6;j++)     // loop for q change in wave vector variation from 0 to 2k
		fprintf(fid1," %e   %e \n",q[j], pz[j]);		

	//*/
	
}


