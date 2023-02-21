# include "main.h"

int main()
{
		
    cout<<"Inside main function"<<endl;	
    //cout<<"iterations = "<<iterations<<endl;

    cout.precision(15);
    FILE *fid1;
    FILE *fid2;
	
	
	if(iterations > 20)
	{
		cout<<"No. of iterations are higher than maximul limit of 20. No. of itertions set to maximum limit of 20"<<endl;
		iterations = 20;
	}	
	
//------------------// read input file  ------------------------------------------------------------------------------
    	
    copyright();	
    read_input_file();	
    cout<<endl;
	    
   		
    //VASP=1;	// read VASP files
    //VASP=0;	// read from table
    //VASP=2;	// tight binding

    int jj=0;
	
			    	
    if(VASP==1)    	
    {
	//--------------------------------// reading outcar file ----------------------------------  	
	    read_OUTCAR();
	//--------------------------------------// reading ispin ----------------------------------
	    jj = read_ispin();	    	
	//--------------------------------Reading ispin completed  --------------------------------
    }
    
    else if(VASP==2)	
	tight_binding_band_structure();

    //cout<<"thickness  =   "<<thickness<<endl;
    //getchar();
      	    
//---------------------------------------------output file generated---------------------------------------------------------

    // running code for different ispin 		
    for (kk=0; kk<= jj ; kk++)
    {
	    if (ispin == 2 && kk==0)    	
		cout<<"-------------------- For up spin ----------------------------------------- "<<endl;

	    if (ispin == 2 && kk==1)
	    {
	    cout<<"-------------------- For down spin ---------------------------------------"<<endl;
	    cc  = -1;
	    }	
	    //--------------------------------------------- finding cbm and vbm --------------------------------------------
	    	    		    		
	    find_cbm_vbm(spin_orbit_coupling);   
	    	
	    kcbm[0]=kcbm1[0];
	    kcbm[1]=kcbm1[1];
	    kcbm[2]=kcbm1[2];
			
	    kvbm[0]=kvbm1[0];
	    kvbm[1]=kvbm1[1];
	    kvbm[2]=kvbm1[2];
	    
	    
	    // unnormalized form of kcbm and kvbm unit is 1/nm
	    cout<<endl<<"kcbm = "<<"   "<<kcbm[0]<<"   "<<kcbm[1]<<"   "<<kcbm[2]<<"   1/nm "<<endl;
	    cout<<endl<<"kvbm = "<<"   "<<kvbm[0]<<"   "<<kvbm[1]<<"   "<<kvbm[2]<<"   1/nm "<<endl;
	    
	    if(VASP==1)
	    {
		    cout<<endl<<"In normalized form kcbm = "<<"   "<<temp1[cbm_index][0]<<"   "<<temp1[cbm_index][1]<<"   "<<temp1[cbm_index][2]<<endl;
		    cout<<endl<<"In normalized form  kvbm = "<<"   "<<temp1[vbm_index][0]<<"   "<<temp1[vbm_index][1]<<"   "<<temp1[vbm_index][2]<<endl;
	    }


	    //------------------------ reading procar file for wave function admixture ---------------------------------------
            read_procar();

	    //------------------------ reading EIGENVAL file for band structure ---------------------------------------		    		    
	    // reading conduction band	
	    make_band(1); // for conduction band
	    count1 = countx;   // length of conduction band
	    //cout<<"count1 = "<<count1<<endl;

	    // reading valence band	
	    make_band(2); // for valence band
	    count2 = countx;   // valence band elements number
	    //cout<<"count2 = "<<count2<<endl;

	    //------------------------------------------------ reading EIGENVAL file is completed -------------------------------------
	    		    
	    cout<<"In main.cpp"<<endl;	
	    cout<<endl<<"CBM = "<<CBM<<" eV"<<endl;
	    cout<<endl<<"VBM = "<<VBM<<" eV"<<endl;
	    

	    if (Bgap[0]==0)
	    {
		int l = (sizeof(Bgap)/sizeof(*Bgap));
		for (int i=0;i<l;i++)
		    Bgap[i] = abs(CBM - VBM);
	    }
            cout<<endl<<"Bgap = "<<Bgap[0]<<" eV"<<endl;

// ------------------------  analytical fitting of CB and VB bands -------------------------------------------------------- 	
	    fitting_band();
		
//---------------------------------------------------Reading DOSCAR -----------------------------------------------
	    if (free_e==0)
		n_DOS();
//--------------------------------------------------- reading doscar completed --------------------------------------

//------------------------------------------   effective mass calculation ---------------------------
		
	    // Electron effective mass calculation 
	    double qq = coefficients_cond[0][a11[0]-2]*2*1e-18/(pow(h_bar,2)*e);
	    m = pow(qq,-1)/m_e;
	    
	    // hole effective mass calculation
	    qq = coefficients_val[0][b11[0]-2]*2*1e-18/(pow(h_bar,2)*e);
	    m_h = pow(qq,-1)/m_e;
	    
	    cout<<"Average electron effective mass, m*_e/m =  "<<m<<endl;
	    cout<<"Average hole effective mass, m*_h/m =  "<<m_h<<endl;
//--------------------------------------------------------------------------------------------------------------------		    
    //cout<<"thickness    =  "<<thickness<<endl;
    //getchar();
	    
	initialize_array();	
//---------------------------------------------------- Runing loop for different doping variation ------------------
		
	    double T;
	    for(int ii=0;ii<count_d;ii++)  // For doping variation
	    {
		n0 = n_array[ii];
		Nd1 = Nd[ii];
		Na1 = Na[ii];
		cout<<"Nd1 = "<<Nd1<<" per cm^3    "<<endl;
		cout<<"Na1 = "<<Na1<<" per cm^3    "<<endl;

//---------------------------------------------------- Runing loop for different temperature variation ------------------
		for(int T_loop=0;T_loop<count_t;T_loop++)
		{
		     T = T_array[T_loop];


		     printf("\nFor Temperature = %f ",T);
		     printf("\nFor Doping = %e \n",n_array[ii]);
	 		
//----------------------------------------energy band important data ----------------------------------------------------------
		    // generate k-grid, energy_n, v_n, a_n, c_n, Ds_n, energy_p, v_p, a_p, c_p, Ds_p 		
		    generate_required_data(T);
			
//---------------------------------------------------------------------------------------------------------------

		    if (energy_n[points-1]<6*k_B*T)
		        cout<<"The selected value of k_max is too low"<<endl;
		    //cout<<"energy_n[points-1] = "<<energy_n[points-1]<<endl;
		    
		    		
//------------------------------------finding fermi level ------------------------------------------------------------
			
		    if (De_ionization==0)
			find_fermi(n_array[ii], T, T_loop);
		    
		    	                		
		    efef_n = 0; // effective ef (fermi energy)
		    efef_p = 0; // effective ef (fermi energy)
		    
		    if(type=="n")	
		    {
			    efef_n = E_F;
			    efef_p = -(E_F+Bgap[T_loop]);
		    }
		    else
		    {
			    efef_n = -(E_F+Bgap[T_loop]);
			    efef_p = E_F;
		    }	
//-------------------------------------------------------------------------------------------------------------------

		     components_BTE(T, T_loop, efef_n, efef_p, ii);
			

		     solve_g(T, T_loop, ii);	


		     save_scattering_rate();


		     save_perturbation();	


		     calculate_mobility(T, T_loop, ii);
				

		}   // temperature variation loop
	    }   //  doping variation loop
	    
	    save_results();
    
    }	// end of loop ue to ispin

    return 0;
}
