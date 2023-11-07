#include "main.h"
int countx;


// typee==1 means CB
// typee==2 means VB

void make_band(int typee)
{
	
	char line[1000];
	FILE *fid,*fid1,*fid2;
	
	// reading is already done in find_cbm_vbm so removed from here
	/*
	if(VASP==1)  // reading from VASP files
	{
	
		int VBM_band_number, CBM_band_number,dummy;
		double temp, temp1[4],temp2;
		int num1;


		if (typee==1) // 1 represent conduction band
		{
			fid = fopen("EIGENVAL_n", "r");
			if (fid==NULL)
			{
			    fid = fopen("EIGENVAL","r");
			    if (fid==NULL)
			    {
				cout<<"EIGENVAL is not present. Exit from program";
				exit(EXIT_FAILURE);
			    }
			}
		}
		else   // 2 reprsent valence band
		{
			fid = fopen("EIGENVAL_p", "r");
			if (fid==NULL)
			{
			    fid = fopen("EIGENVAL","r");
			    if (fid==NULL)
			    {
				cout<<"EIGENVAL is not present. Exit from program";
				exit(EXIT_FAILURE);
			    }
			}
		}

		fgets(line, 1000, (FILE*)fid);
		sscanf(line, "%d %d %d %d", &dummy, &dummy, &dummy, &dummy);  
		//printf("%d ", ispin);
		//getchar();

		for (int i = 0; i < 5; i++)
			fgets(line, 1000, (FILE*)fid);   // passing 5 unwanted line

		//printf("\n%s", line);
		//getchar();

		int a[3];

		sscanf(line, "%d %d %d", &a[0], &a[1], &a[2] );

		NKPTS = a[1];

		if (spin_orbit_coupling==0)  //spin-orbit_coupling == false
			NBVAL = int(a[0]/2.0) ;
		else
			NBVAL = int(a[0]) ;    //spin-orbit_coupling == true

		NBTOT = a[2];

		VBM_band_number = NBVAL;

		CBM_band_number = VBM_band_number + 1; // +1

		//printf("\nCBM_band_number = %d", CBM_band_number);
		
		/*
		if (ispin == 1)
		{
			
			for (int i = 0; i < NKPTS; i++)
			{
			    fgets(line, 1000, (FILE*)fid);   //passing blank line
			    fgets(line, 1000, (FILE*)fid);   // reading next kpoint line
			    sscanf(line, "%lf %lf %lf %lf", &temp1[0], &temp1[1], &temp1[2], &temp1[3]);
			    
			    	
			    kpoints[i][0] = (temp1[0] * lm[0][3] + temp1[1] * lm[1][3] + temp1[2] * lm[2][3]) * 2. * 3.14159265359 * 10;
			    kpoints[i][1] = (temp1[0] * lm[0][4] + temp1[1] * lm[1][4] + temp1[2] * lm[2][4]) * 2. * 3.14159265359 * 10;
			    kpoints[i][2] = (temp1[0] * lm[0][5] + temp1[1] * lm[1][5] + temp1[2] * lm[2][5]) * 2. * 3.14159265359 * 10;


			    for (int j = 0; j < VBM_band_number; j++)
				fgets(line, 1000, (FILE*)fid);

			    sscanf(line, "%d %lf", &num1, &temp);   // save values of valence band only
			    energies[i][0] = temp;   // valence band energy

			    fgets(line, 1000, (FILE*)fid);  // read conduction band
			    sscanf(line, "%d %lf", &num1, &temp);

			    energies[i][1] = temp;     // conduction band energy

			    //cout<<"k_points = "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<endl;
			    //cout<<energies[i][0]<<"   "<<energies[i][1]<<endl;
			    //getchar();

			    for (int j = VBM_band_number+1; j < NBTOT; j++)
				fgets(line, 1000, (FILE*)fid);
			    //printf("\n %d %e %e %e %e", i, kpoints[i][0], kpoints[i][1], kpoints[i][2], energies[i]);
			}
			
		}

		if (ispin == 2)
		{
			for (int i = 0; i < NKPTS; i++)
			{
			    fgets(line, 1000, (FILE*)fid);   //passing blank line
			    fgets(line, 1000, (FILE*)fid);   // reading next kpoint line
			    sscanf(line, "%lf %lf %lf %lf", &temp1[0], &temp1[1], &temp1[2], &temp1[3]);
				
			//cout<<"k_points  =  "<<temp1[0]<<"      "<<temp1[1]<<"       "<<temp1[2]<<"      "<<temp1[3]<<endl;
			    
			    // unit is 1/nm		
			    kpoints[i][0] = (temp1[0] * lm[0][3] + temp1[1] * lm[1][3] + temp1[2] * lm[2][3]) * 2. * 3.14159265359 * 10;
			    kpoints[i][1] = (temp1[0] * lm[0][4] + temp1[1] * lm[1][4] + temp1[2] * lm[2][4]) * 2. * 3.14159265359 * 10;
			    kpoints[i][2] = (temp1[0] * lm[0][5] + temp1[1] * lm[1][5] + temp1[2] * lm[2][5]) * 2. * 3.14159265359 * 10;


			    for (int j = 0; j < VBM_band_number; j++)
				fgets(line, 1000, (FILE*)fid);

			    sscanf(line, "%d %lf %lf ", &num1, &temp, &temp2);   // save values of valence band only

				//cout<<"valence Band    :=  "<<num1<<"    "<<temp<<" "<<temp2<<endl;
				
			    if (kk == 0)	
			    	energies[i][0] = temp;   // valence band energy for up spin
			    else	
			    	energies[i][0] = temp2;   // valence band energy for down spin


			    fgets(line, 1000, (FILE*)fid);  // read conduction band
			    sscanf(line, "%d %lf %lf", &num1, &temp, &temp2);

			    if (kk == 0)	
				    energies[i][1] = temp;     // conduction band energy for up spin
			     else
				    energies[i][1] = temp2;     // conduction band energy for down spin

				//cout<<"Conduction band :=  "<<num1<<"    "<<temp<<"   "<<temp2<<endl;
			   	//getchar();
			   	
			    //cout<<"k_points = "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<endl;
			    //cout<<energies[i][0]<<"   "<<energies[i][1]<<endl;
			    //getchar();

			    for (int j = VBM_band_number+1; j < NBTOT; j++)
				fgets(line, 1000, (FILE*)fid);
			    //printf("\n %d %e %e %e %e", i, kpoints[i][0], kpoints[i][1], kpoints[i][2], energies[i]);
			}  	// end of for loop
		}  		// end of ispin==2
		*/
		
		//--------------------------------------- save band data ---------------------------------
		//cout<<"Saving band data EK  getchar(); two times"<<endl;
		//getchar();
		/*
		if(typee==1)	// for CB // save CB band data	
		{
			fid2 = fopen("EK_CB.dat","w");
			fprintf(fid2,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");
		
			for (int i = 0; i < NKPTS; i++)
			fprintf(fid2,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energies[i][1]);
		
			fclose(fid2);
		}
		else     // for VB  // save valence band data 		
		{
			fid2 = fopen("EK_VB.dat","w");
			fprintf(fid2,"# kx(1/cm)   ky(1/cm)    kz(1/cm)    energy  \n");
		
			for (int i = 0; i < NKPTS; i++)
			fprintf(fid2,"%e  %e  %e   %e \n", kpoints[i][0]*1e7, kpoints[i][1]*1e7, kpoints[i][2]*1e7, energies[i][0]);
		
			fclose(fid2);
		}
		//--------------------------------------------- band data saved ----------------------------------------------
		//
	}
	else if(TBS==1)
	{
				
	}
	else if(VASP==0)// for table form reading band structure 
	{
		//cout<<"Reading band data getchar(); three times"<<endl;
		//getchar(); getchar(); getchar();
		// reading CB
		int i;
		
		if(typee==1)    
		{
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
			       sscanf(line, "%lf %lf %lf %lf", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2], &energies[i][1]);

				// converted from 1/cm to 1/nm; all calculation are done by assuming k unit is 1/nm
				kpoints[i][0] = kpoints[i][0]/1e7;
				kpoints[i][1] = kpoints[i][1]/1e7;
				kpoints[i][2] = kpoints[i][2]/1e7;
			//cout<<"i = "<<i<<"   "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"  "<<energies[i][1]<<endl;
				//getchar();
				i++;
			}

		}
		// reading VB
		else    // reading VB
		{
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
			       sscanf(line, "%lf %lf %lf %lf", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2], &energies[i][0]);

				// converted from 1/cm to 1/nm; all calculation are done by assuming k unit is 1/nm
				kpoints[i][0] = kpoints[i][0]/1e7;
				kpoints[i][1] = kpoints[i][1]/1e7;
				kpoints[i][2] = kpoints[i][2]/1e7;
			//cout<<"i = "<<i<<"   "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"  "<<energies[i][0]<<endl;
				//getchar();	
				i++;
			}
		}
		NKPTS = i;		
	// for table form reading band structure completed	
	}    // if and else for reading VASP or from table completed

	fclose(fid);
	
	
	
	
	/* Already find in find_cbm_vbm.cpp
	if (typee==1)
	{
		double ECBM = 1000;
		
		for (int i = 0; i < NKPTS; i++)
		{
		    if (energies[i][1] < ECBM)
			ECBM = energies[i][1];    // energies[][1]  contains conduction band
						
		}
		CBM = ECBM;
		
	}
	else
	{
		double EVBM = -1000, EVBM1 = 1000;
		for (int i = 0; i < NKPTS; i++)
		{
		    if (energies[i][0] > EVBM)
			EVBM = energies[i][0];    // energies[][0]  contains valence band
		}
		VBM = EVBM;
	}
	*/
	
	/*
	cout<<"Inside make_band.cpp "<<endl;
	if (typee == 1)
	printf("\n CBM = %lf \n", CBM);
	else
	printf("\n VBM = %lf \n", VBM);
	//getchar();
	*/	


	// this below portion is common for all VASP, TBS and table form 
	// only kcbm kvbm, NKPTS  are required  
	double temp_reference[3], reference_point[3];
	
	// unit of kcbm and kvbm is 1/nm
	if (typee==1) // means CB
	{
		temp_reference[0] = kcbm[0];
		temp_reference[1] = kcbm[1];
		temp_reference[2] = kcbm[2];
	}
	else  // means VB
	{
		temp_reference[0] = kvbm[0];
		temp_reference[1] = kvbm[1];
		temp_reference[2] = kvbm[2];
	}
	
	
	double sum;
	//printf("\n Reference point = ");

	for (int i = 0; i < 3; i++)  // // i = 0 gives kx,   i = 1 gives ky and i = 2 gives kz
	{
		reference_point[i] = temp_reference[i];
		//cout<<"reference_point[i] = "<<reference_point[i]<<endl;
		//printf(" %e", reference_point[i]);
	}
	//getchar();
	//*/


	double conduction_band[NKPTS][2];
	double distance, kx, ky, kz;

	if(typee==1 || (typee==2 && flagg ==0))
	{
		for (int i = 0; i < NKPTS; i++)
		{
			kx = pow(kpoints[i][0] - reference_point[0], 2);
			ky = pow(kpoints[i][1] - reference_point[1], 2);
			kz = pow(kpoints[i][2] - reference_point[2], 2);
			// unit 1/nm

			//cout<<"i = "<<i<<endl;
			//cout<<"kz = "<<kz<<endl;
			//cout<<"kpoints[i][2] = "<<kpoints[i][2]<<endl;
			//getchar();

			distance = sqrt(kx + ky + kz);   // unit is 1/nm
			conduction_band[i][0] = distance;
			conduction_band[i][1] = energies[i][1];    // CB
			//printf("\n %e %lf", distance,energies[i]);
		}
	}
	else if(typee==2 && flagg ==1)
	{
		for (int i = 0; i < NKPTS; i++)
		{
			kx = pow(kpoints_p[i][0] - reference_point[0], 2);
			ky = pow(kpoints_p[i][1] - reference_point[1], 2);
			kz = pow(kpoints_p[i][2] - reference_point[2], 2);
			// unit 1/nm

			//cout<<"i = "<<i<<endl;
			//cout<<"kz = "<<kz<<endl;
			//cout<<"kpoints[i][2] = "<<kpoints[i][2]<<endl;
			//getchar();

			distance = sqrt(kx + ky + kz);   // unit is 1/nm
			conduction_band[i][0] = distance;
			conduction_band[i][1] = energies[i][1];    // VB
			//printf("\n %e %lf", distance,energies[i]);
		}
	}
		
	if (typee == 2)
	{
		for (int i = 0; i < NKPTS; i++)
		    conduction_band[i][1] = - energies[i][0];  // VB
	}

	//checked

	double offset = abs(conduction_band[0][1]);
	//printf("\nThe offset = %lf", offset);

	double conduction_band_num_dum[NKPTS][2];

	double conduction_band1[NKPTS], conduction_band2[NKPTS];

	//fid1 = fopen("conduction_band","w");

	for (int i = 0; i < NKPTS; i++)
	{
		conduction_band1[i] = conduction_band[i][0];  // distance
		conduction_band2[i] = conduction_band[i][1];  // energy
	}

	if(SORT==0)	// by default this part is slected
	{
		sort(conduction_band1, conduction_band1+NKPTS);
		sort(conduction_band2, conduction_band2+NKPTS);
	}
	else // if sorting is selected then only this part will run 
	{	
		double temp3;
		// bubble sort according to distance 
		for(int i=0;i<NKPTS-1;i++)
		{
			for(int j=0;j<NKPTS-i;j++)
			{
				if(conduction_band1[j+1] < conduction_band1[j])
				{
					// swap elements
					temp3 = conduction_band1[j+1]; 
					conduction_band1[j+1] = conduction_band1[j];
					conduction_band1[j] = temp3;
					
					temp3 = conduction_band2[j+1]; 
					conduction_band2[j+1] = conduction_band2[j];
					conduction_band2[j] = temp3;
				}   
			}
		}		 
	}   // sorting completed
	
	
	for (int i = 0; i < NKPTS; i++)
	{
		conduction_band_num_dum[i][0] = conduction_band1[i];    // distance
		conduction_band_num_dum[i][1] = conduction_band2[i];    // energy
	// checked //cout<< conduction_band_num_dum[i][0] <<"    "<<conduction_band_num_dum[i][1]<<endl;
	}

	double conduction_band_num[NKPTS][2];   // this 2d array is passed to main

	conduction_band_num[0][0] = conduction_band_num_dum[0][0];
	conduction_band_num[0][1] = conduction_band_num_dum[0][1];
	//printf("\n %e %e", conduction_band_num[0][0], conduction_band_num[0][1]);

//--------------------// doing average kpoints started----------------------------------------------------------
	double avg[1][2] = {0}, z = 0.0001;
	int start = 0, stop = 0;
	countx = 0;

	for (int i = 1; i < NKPTS - 1; i++)
	{
		if (abs(conduction_band_num_dum[i][0] - conduction_band_num_dum[i - 1][0]) < z)
		{
		    if (start == 0)
			start = i - 1;

		    if ((conduction_band_num_dum[i+1][0] - conduction_band_num_dum[i][0]) > z)
		    {
			stop = i;
			for (int j = start; j <= stop; j++)
			{
			    avg[0][0] = avg[0][0]+conduction_band_num_dum[j][0];
			    avg[0][1] = avg[0][1]+conduction_band_num_dum[j][1];
			}
			avg[0][0] = avg[0][0]/(stop - start + 1);
			avg[0][1] = avg[0][1]/(stop - start + 1);

			countx++;
			conduction_band_num[countx][0] = avg[0][0];
			conduction_band_num[countx][1] = avg[0][1];

			start = 0;
			stop = 0;
			avg[0][0] = 0;
			avg[0][1] = 0;
		    }
		}
		if( (abs(conduction_band_num_dum[i][0] - conduction_band_num_dum[i - 1][0]) > z) &
		   ((conduction_band_num_dum[i + 1][0] - conduction_band_num_dum[i][0]) > z) )
		{
		    countx++;
		    conduction_band_num[countx][0] = conduction_band_num_dum[i][0];
		    conduction_band_num[countx][1] = conduction_band_num_dum[i][1];

		}
	}

//--------------------// doing average of kpoints completed----------------------------------------------------------

	if (conduction_band_num_dum[NKPTS-1][0] > conduction_band_num_dum[NKPTS-2][0])
	{
		countx++;
		conduction_band_num[countx][0] = conduction_band_num_dum[NKPTS-1][0];
		conduction_band_num[countx][1] = conduction_band_num_dum[NKPTS-1][1];
	}

	if ((conduction_band_num[1][0] - conduction_band_num[0][0]) < z)
	{
		for (int i = 1; i < countx+1; i++)
		{
		    conduction_band_num[i - 1][0] = conduction_band_num[i][0];
		    conduction_band_num[i - 1][1] = conduction_band_num[i][1];
		}
		countx--;
	}


//------- -----------remove offset of energy 	 ------------------------------------------------------
	double normalization_factor;
	if(SORT==0)
	{
		normalization_factor = 0 - conduction_band_num_dum[0][1];
		for (int i = 0; i < countx+1; i++)
			conduction_band_num[i][1] = conduction_band_num[i][1] + normalization_factor;
	}
	else  // if buuble sort is done for only distance 
	{
		cout<<"SORTING is selected for bands"<<endl;
		// find minimum energy
		double  minimum=100;
		for(int i=0;i<NKPTS;i++)
		{
			if(minimum > conduction_band_num[i][1])
				minimum = conduction_band_num[i][1];	
		}
	
		normalization_factor = -1*minimum;
		for (int i = 0; i < countx+1; i++)
			conduction_band_num[i][1] = conduction_band_num[i][1] + normalization_factor;			
	}
//------- -----------remove offset of energy completed ------------------------------------------------------	

//---------------------------------------------- save band ----------------------------------------------------		
	/*
	if (typee==1)
		fid1 = fopen("conduction_band_num.txt","w");
	else
		fid1 = fopen("valence_band_num.txt","w");
	
	fprintf(fid1,"index    kpoint     energy (eV) \n");
	for (int i = 0; i < countx+1; i++)
		fprintf(fid1,"%d  %e  %e   \n", i+1, conduction_band_num[i][0], conduction_band_num[i][1]);
	

	fclose(fid1);	
	//*/
//--------------------------------------------------------------------------------------------------------------------------	
	countx = countx+1;

	if (typee==1)  //  Conduction Band 
	{
		for (int i=0;i<countx;i++)
		{
		    cond_band[i][0] = conduction_band_num[i][0];
		    cond_band[i][1] = conduction_band_num[i][1];
		}
	}
	else    // Valence Band
	{
		for (int i=0;i<countx;i++)
		{
		    val_band[i][0] = conduction_band_num[i][0];
		    val_band[i][1] = conduction_band_num[i][1];
		}
	}


//-------------------------------------------------------------------------------------------------------------------


		
//------------------------------------- linear fitting of band is done here for ----------------------------------

	if(linear_fit==1) // without intercept
	{
		//cout<<endl<<"Linear fitting without intercept is selected"<<endl;
		//cout<<"Emax = "<<Emax<<"  eV"<<endl;
		//cout<<"countx = "<<countx<<endl;
		
		double sum_xy=0,sum_x2=0, a;
		int n=0;
		for(int i=0;i<countx;i++)
		{	
			if(conduction_band_num[i][1] < Emax)
			{
				n++;
				sum_xy = sum_xy + conduction_band_num[i][1] * conduction_band_num[i][0]*1e9*e;
				sum_x2 = sum_x2 + conduction_band_num[i][0] * conduction_band_num[i][0]*1e9*1e9;  

				// multiplied with 1e9 to converted from 1/nm to 1/m
				/*
				if(typee!=1)
				{
					cout<<"i = "<<i<<endl;
					cout<<"sum_xy = "<<sum_xy<<endl;
					cout<<"sum_x2 = "<<sum_x2<<endl;
					cout<<"conduction_band_num[i][0] = "<<conduction_band_num[i][0]<<endl;
					cout<<"conduction_band_num[i][1] = "<<conduction_band_num[i][1]<<endl;
					getchar();
				}
				*/
			}			
		}
		
		//cout<<"sum_xy = "<<sum_xy<<endl;
		//cout<<"sum_x2 = "<<sum_x2<<endl;
		cout<<"n = "<<n<<endl;
		// slope multilped 
		a = (sum_xy/sum_x2);  // unit is J-m     
		
		if(typee==1)
		{
			vf_cb = a/(h_bar*e);  // unit is m/s, h_bar unit is eV-s
			cout<<endl<<"Calculated fermi velocity for conduction band = "<<vf_cb<<"   m/s"<<endl;		
		}
		else
		{
			vf_vb = a/(h_bar*e);  // unit is m/s, h_bar unit is eV-s
			vf = (vf_cb + vf_vb)/2.0;
			cout<<endl<<"Calculated fermi velocity for valence band = "<<vf_vb<<"   m/s"<<endl;		
			cout<<endl<<"Average fermi velocity = "<<vf<<"   m/s"<<endl;	
		}
		
		//getchar();	
		
	}
	else if(linear_fit==2) //with intercept
	{
		//cout<<endl<<"Linear fitting with intercept is selected"<<endl;
		//cout<<"Emax = "<<Emax<<"  eV"<<endl;
		//getchar();
		double sum_xy=0, sum_x2=0, sum_x=0, sum_y=0, delta=0, a, b;
		int n=0;
		for(int i=0;i<countx;i++)
		{
			if(conduction_band_num[i][1] < Emax)
			{
				n++;
				sum_xy = sum_xy + conduction_band_num[i][1] * conduction_band_num[i][0]*1e9*e;
				sum_x2 = sum_x2 + conduction_band_num[i][0] * conduction_band_num[i][0]*1e9*1e9;
				sum_x = sum_x + conduction_band_num[i][0] * 1e9;
				sum_y = sum_y + conduction_band_num[i][1] * e;
				// multiplied with q to converted from eV to joule  
				// multiplied with 1e9 to converted from 1/nm to 1/m
			}
			//cout<<"sum_xy = "<<sum_xy<<endl;
			//cout<<"sum_x2 = "<<sum_x2<<endl;
			//cout<<"sum_x = "<<sum_x<<endl;
			//cout<<"sum_y = "<<sum_y<<endl;
		}
		delta = n*sum_x2 - sum_x*sum_x;
		a = (n*sum_xy - sum_x * sum_y)/delta ;  // slope 
		b = (sum_x2 * sum_y - sum_x * sum_xy)/delta;
		if(typee==1)
		{
			vf_cb = a/(h_bar*e);  // unit is m/s, h_bar unit is eV-s
			cout<<endl<<"Calculated fermi velocity for conduction band = "<<vf_cb<<"   m/s"<<endl;		
		}
		else
		{
			vf_vb = a/(h_bar*e);  // unit is m/s, h_bar unit is eV-s
			vf = (vf_cb + vf_vb)/2.0;
			cout<<endl<<"Calculated fermi velocity for valence band = "<<vf_vb<<"   m/s"<<endl;		
			cout<<endl<<"Average fermi velocity = "<<vf<<"   m/s"<<endl;	
		}

		//cout<<"NKPTS= "<<NKPTS<<endl;
		//cout<<"countx  = "<<countx<<endl;
		//cout<<"n = "<<n<<endl;
		//cout<<"a = "<<a<<endl;
		//cout<<"b = "<<b<<endl;
		//cout<<"delta = "<<delta<<endl;	
		//cout<<"With intercept "<<endl;	
		//getchar();	

	}	
//--------------------------------------- linear fitting completed ---------------------------------------------	

}
