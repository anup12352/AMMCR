#include"main.h"

void find_cbm_vbm(int spin_orbit_coupling)
{
	evbm = -1000;
	ecbm = 1000;

	char line[1000];
	int a[10],i,j ;
	FILE *fid;

	double temp[3];
	
	if(VASP==1)
	{
		fid = fopen("EIGENVAL","r");
		if (fid==NULL)
		{
			cout<<"EIGENVAL is not present. Exit from program";
			exit(EXIT_FAILURE);
		}

		fgets(line, 1000, fid);
		sscanf(line, "%d %d %d %d", &a[0], &a[1], &a[2], &a[3]);

		for(i=1;i<=5;i++)
		{
			fgets(line,1000,fid);   // passing next 4 unimportant line and reading next 5th line
		}

		sscanf(line, "%d %d %d", &a[0], &a[1], &a[2] );

		NKPTS = a[1];

		if (spin_orbit_coupling==0)  //spin-orbit_coupling == false
			NBVAL = int(a[0]/2.0) ;
		else
			NBVAL = int(a[0]) ;    //spin-orbit_coupling == true

		NBTOT = a[2];

		//cout<<" NKPTS = "<<NKPTS<<"    NBVAL = "<<NBVAL<<"  NBTOT   ="<<NBTOT<<endl;
		//getchar();
		    
		for (i=0; i<NKPTS; i++)
		{
			fgets(line,1000,fid);   //passing one empty line
			fgets(line,1000,fid);    // reading kpoint line
			sscanf(line, "%lf %lf %lf", &temp1[i][0], &temp1[i][1], &temp1[i][2]);

		kpoints[i][0] = (temp1[i][0] * lm[0][3] + temp1[i][1] * lm[1][3] + temp1[i][2] * lm[2][3]) * 2. * 3.14159265359 * 10;
		kpoints[i][1] = (temp1[i][0] * lm[0][4] + temp1[i][1] * lm[1][4] + temp1[i][2] * lm[2][4]) * 2. * 3.14159265359 * 10;
		kpoints[i][2] = (temp1[i][0] * lm[0][5] + temp1[i][1] * lm[1][5] + temp1[i][2] * lm[2][5]) * 2. * 3.14159265359 * 10;
			// unit is 1/nm
			
		//cout<<"line = "<<line<<endl;
		//cout<<"kpoints = "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"   "<<kpoints[i][3]<<endl;
		//getchar();
		
			for (j = 0; j<NBTOT; j++)
			{
			    //cout<<"j = "<<j<<endl;
			    fgets(line,1000,fid);
			    //cout<<"line = "<<line<<endl;
			    sscanf(line, "%lf %lf %lf ", &temp[0], &temp[1], &temp[2]); 
			    
			    if(j==NBVAL-1)    // VB           
			    	energies[i][0] = temp[1];     
			    else if(j==NBVAL)  // CB
			    	energies[i][1] = temp[1];            
			    
			    //energies[i][0]  contains valence band
			    // energies[i][1]  contains CB
			    					
			    if (ispin == 2 && kk==1)    // reading for spin down
			    {	
				    if(j==NBVAL-1)    // VB           
				    	energies[i][0] = temp[2];     
				    else if(j==NBVAL)  // CB
				    	energies[i][1] = temp[2];
			    }            
			    /*
			     if(j == NBVAL-1 )
			     	cout<<"Valence Energy = "<<energies[i][0]<<endl;		

			     if(j == NBVAL )
			     	cout<<"Conduction Energy = "<<energies[i][1]<<endl;		
			    //getchar();
			    
			    */
			} // loop for reading energy for bands is completd here
		//getchar();
		} // loop for all kpoints is completed here
	}   // if condiction for VASP==1 completed
	else if(VASP==0)    // reading from table
	{
		//cout<<"Reading CB for CBM "<<endl;
		fid = fopen("EK_CB.dat","r");
		if (fid==NULL)
		{
			cout<<"EK_CB.dat file is not present";
			exit(EXIT_FAILURE);
		}

		fgets(line, 1000, fid);   // pass first line
		i=0;
		while (fgets(line, 1000, fid)!= NULL)
		{
			//fgets(line, 1000, fid);
		       sscanf(line, "%lf %lf %lf %lf ", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2], &energies[i][1]);
			
			// converted from 1/cm to 1/nm; all calculation are done by assuming k unit is 1/nm
			kpoints[i][0] = kpoints[i][0]/1e7;
			kpoints[i][1] = kpoints[i][1]/1e7;
			kpoints[i][2] = kpoints[i][2]/1e7;
			//cout<<"i = "<<i<<"   "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"  "<<energies[i][1]<<endl;
			//getchar();			
			i++;
		}
		NKPTS = i;
		//cout<<"NKPTS =     "<<NKPTS<<endl;

		//cout<<"Reading VB for VBM "<<endl;
		fid = fopen("EK_VB.dat","r");
		if (fid==NULL)
		{
			cout<<"EK_VB.dat file is not present";
			exit(EXIT_FAILURE);
		}

		fgets(line, 1000, fid);   // pass first line
		i=0;
		while (fgets(line, 1000, fid)!= NULL)
		{
			//fgets(line, 1000, fid);
		       sscanf(line, "%lf %lf %lf %lf ", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2], &energies[i][0]);
			kpoints[i][0] = kpoints[i][0]/1e7;
			kpoints[i][1] = kpoints[i][1]/1e7;
			kpoints[i][2] = kpoints[i][2]/1e7;
			// converted from 1/cm to 1/nm; 
			// all calculation are done by assuming k unit is 1/nm
			//cout<<"i = "<<i<<"   "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"  "<<energies[i][0]<<endl;
			//getchar();			
			i++;
		}
		NKPTS = i;
		//cout<<"NKPTS =     "<<NKPTS<<endl;
	}   // reading from table else if ended
	else if(VASP==2)    // For tight_binding_band_structure
	{
		tight_binding_band_structure();	
		for (i=0; i<NKPTS; i++)
		{
			energies[i][0] = energy_p[i];
			energies[i][1] = energy_n[i];
			
		}		    	
	} // end of if else condiction for VASP =0, 1 and 2
		
	//cout<<"Comparing VB and CB for VBM and CBM "<<endl;		
	for (i=0; i<NKPTS; i++)
	{
		if(energies[i][0] > evbm)
		{
			evbm = energies[i][0];      
			kvbm1[0] = kpoints[i][0];    // unit is 1/nm
			kvbm1[1] = kpoints[i][1];
			kvbm1[2] = kpoints[i][2];

			if(VASP==1)
				vbm_index = i;

		}

		if(energies[i][1] < ecbm )
		{
			ecbm = energies[i][1];
			kcbm1[0] = kpoints[i][0];
			kcbm1[1] = kpoints[i][1];
			kcbm1[2] = kpoints[i][2];

			if(VASP==1)
				cbm_index = i;

		}

	}
	
	CBM = ecbm;
	VBM = evbm;
	
	
	if(save_data==1)
	{
		fid = fopen("EK_CB.dat","w");
		fprintf(fid,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");
	
		for (int i = 0; i < NKPTS; i++)
			fprintf(fid,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energies[i][1]);
		//kpoints are save in 1/cm
	
		fclose(fid);

		fid = fopen("EK_VB.dat","w");
		fprintf(fid,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");
	
		for (int i = 0; i < NKPTS; i++)
			fprintf(fid,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energies[i][0]);
		//kpoints are save in 1/cm
	
		fclose(fid);
	}

	//cout<<"evbm = "<<evbm<<"  kvbm1[0] = "<<kvbm1[0]<<"   kvbm1[1] = "<<kvbm1[1]<<"   kvbm1[2]   = "<<kvbm1[2]<<endl;
	//cout<<"ecbm = "<<ecbm<<"   kcbm1[0] = "<<kcbm1[0]<<"   kcbm1[1] = "<<kcbm1[1]<<"   kcbm1[2]   = "<<kcbm1[2]<<endl;

	//cout<<"Outside evbm = "<<evbm<<"  kvbm1[0] = "<<kvbm1[0]<<"   kvbm1[1] = "<<kvbm1[1]<<"   kvbm1[2]   = "<<kvbm1[2]<<endl;
	//cout<<"Outside ecbm = "<<ecbm<<"  kcbm1[0] = "<<kcbm1[0]<<"   kcbm1[1] = "<<kcbm1[1]<<"   kcbm1[2]   = "<<kcbm1[2]<<endl;
	    	
}



