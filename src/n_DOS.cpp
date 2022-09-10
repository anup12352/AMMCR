# include "main.h"
double DOS_n[2000][2],DOS_p[2000][2];
int count_DOS_n,count_DOS_p;

int min_index(double band[][2],int county, int coloumn)
{
    double min1 = 1000;
    int min_indexx = -1;
	for (int i = 0; i < county; i++)
	{
		if (abs(band[i][coloumn-1]) < min1)
		{
		    min_indexx = i;
		    min1 = abs(band[i][coloumn-1]);
		}
	}
	return min_indexx;
}

int max_index(double band[][2],int county,int coloumn)
{
    double max1 = -1000;
    int max_indexx = -1;
	for (int i = 0; i < county; i++)
	{
		if (abs(band[i][coloumn-1]) > max1)
		{
		    max_indexx = i;
		    max1 = abs(band[i][coloumn-1]);
		}
	}
	return max_indexx;
}

double max_value(double band[][2],int county, int coloumn)
{
    double max1 = -100000000;
    int max_indexx = -1;
	for (int i = 0; i < county; i++)
	{
		if ((band[i][coloumn-1]) > max1)
		{
		    max1 = band[i][coloumn-1];
		    max_indexx = i;
		}
	}
	return max1;
}

void n_DOS()
{
	//cout<<CBM<<"  "<<VBM<<"  "<<ispin<<"   "<<volume1<<endl;
	char line[1000];
	FILE *fid;

	double dummy;

	int num_points;		
	double energy[limit3],states[limit3], states_spin[limit3]={0};

	if(VASP==1)
	{
		//------------------reading DOS from DOSCAR -----------------------------
		fid = fopen("DOSCAR","r");
		if (fid==NULL)
		{
			cout<<"DOSCAR is not present. Program is running with free_e = true"<<endl;
			free_e = 1;
			return ;
		}

		for(int j=1;j<=6;j++)
		fgets(line, 1000, (FILE*)fid);

		sscanf(line, "%lf %lf %d %lf", &dummy, &dummy, &num_points,&fermi);  //ispin is useful

		cout<<"dos fermi = "<<fermi<<endl;
		//cout<<"num_points = "<<num_points<<endl;
		//cout<<"energy[0] = "<<energy[0]<<endl;

		for(int i=0;i<num_points;i++)
		{
			fgets(line, 1000, (FILE*)fid);
			sscanf(line, "%lf %lf %lf", &energy[i], &states[i], &states_spin[i]);

			if (ispin==2 && kk == 1)    // reading for downspin 
				states[i] = states_spin[i];

			if (states[i]<0)
			    states[i]=0;
		    
		}
		fclose(fid);		
		
	}
	else  if(VASP==0)//VASP=0 table form is used
	{
		//cout<<"fermi is not give till now getchar(); four times"<<endl;
		//getchar();  getchar(); getchar(); getchar();
		//fermi = 2.952756400000000e+00;
		fermi = VBM;
		
		cout<<"fermi = "<<fermi<<endl;
		//getchar();
		//fermi = 1.301916870000000e+00;  // for ZnSe
		//volume1 = 1/(1e24);  // converted from 1 cm^3 to 1 A^3 volume is just now 
		volume1 = 1;
				
		fid = fopen("DOS.dat","r");
		if (fid==NULL)
		{
			cout<<"DOS.dat file is not present. ";
			cout<<"Program is running by assuming free desnity of states"<<endl;
			free_e = 1;
			return ;
		}

		fgets(line, 1000, fid);   // pass first line
		int i=0;
		while ((fgets(line, 1000, fid)!= NULL))
		{
		       //fgets(line, 1000, fid);
			sscanf(line, "%lf %lf ", &energy[i], &states[i]);   
			// unit of energy is eV and unit of density of states is per eV per cm^3
			
			//states[i] = states[i]/(1e24);  // converted from A^-3 to cm^-3 
			i++;
		}
		num_points = i;
	}
	else if(VASP==2)  // tight binding band structure
	{
		fermi = VBM;
		volume1 = 1;

		cout<<"fermi = "<<fermi<<endl;
		//getchar();
		
		tight_binding_band_DOS();
		
		num_points = xx1.size();
		
		for(int i=0;i<xx1.size();i++)
		{
			energy[i] = xx1[i];
			
			// unit of e_pp and e_nn is per eV per cm^2
			if(i < e_pp.size() )
				states[i] = e_pp[i];
			else
				states[i] = e_nn[i - e_pp.size()];
		}
			
	}

	/*
	fid = fopen("DOS_new.dat","w");
	
	fprintf(fid,"# energy (eV)  	DOS (per eV per cm^2) \n");			
	for (int i = 0; i < xx1.size(); i++)
		fprintf(fid,"%d     %e     %e \n", i+1, energy[i], states[i]);
	
	fclose(fid);
	//*/
	
	// this portion is common for VASP=0,1,2 
	// energy[i], states[] and fermi, Bgap[0], num_points are required
	
	double DOS[num_points][2];   		// used for CB
	double DOS_temp[num_points][2];       // used for VB
	
	for (int i=0;i<num_points;i++)
	{
		DOS[i][0] = energy[i]-(fermi+Bgap[0]);   // this is used for DOS_n Now CBM is at 0 eV
		DOS[i][1] = states[i];
		
		DOS_temp[i][0] = energy[i]-fermi;  // this is used for DOS_p Now VBM is at 0 eV
		DOS_temp[i][1] = states[i];
	}
	//------------------ printing DOS -----------------------------
	/*
	for(int i=0;i<count_DOS;i++)
	{
	cout<<DOS[i][0]<<"   "<<DOS[i][1]<<endl;
	if (i==50||i==200||i==300||i==2000||i==3000||i==4000||i==5000)
	getchar();

	}
	*/
	//------------------------------------------------------------------

	//----------------------- Separating DOS_n and DOS_p ------------------------

	int end_index,start_index;


	double max1;

//----------------------------------------------------------------------------------------------------------------------------------	
	// for VB
	max1 = max_value(val_band,count2,2);
	//cout<<"max value = valence band ="<<max1<<endl;
	
	double dum[num_points][2]={0};
	
	for(int i=0;i<num_points;i++)
	{ 
		dum[i][0] = DOS_temp[i][0]+ 1.1 * max1;
	}  
	
	//cout<<"dum[][] = "<<dum[i][0]<<endl; }
	//getchar();

	start_index = min_index(dum,num_points,1);    //  Find a maximum for energy range in DOS for valence band
	end_index = min_index(DOS_temp,num_points,1);   // Find the closest point to VBM Since VBM is at 0 energy

	//cout<<"start_index = "<<start_index<<"    end_index = "<<end_index<<endl;
	//getchar();
	
	if(VASP==1)
	{
		while (DOS_temp[end_index][1] < 1e-6)
		{
			end_index = end_index - 1;
			start_index = start_index - 1;
		}
	}
	
	for(int i=0;i<num_points;i++)
		DOS_temp[i][0] = -1*DOS_temp[i][0];

	for(int i=0;i<num_points;i++)
	    DOS_temp[i][0] = DOS_temp[i][0] - DOS_temp[end_index][0];

	int j=0;
	for(int i=end_index;i>=start_index;i--)
	{
		DOS_p[j][0] = DOS_temp[i][0];
		DOS_p[j][1] = DOS_temp[i][1];
		j++;
	}
	count_DOS_p = j;
//----------------------------------------------------------------------------------------------------------------------------------
	
	
//----------------------------------------------------------------------------------------------------------------------------------
	// for CB
	max1 = max_value(cond_band,count1,2);
	//cout<<"max1 = conduction band = "<<max1<<endl;
	//cout<<"next";
	for(int i=0;i<num_points;i++)
	{ 
		dum[i][0] = DOS[i][0] - 1.1 * max1;      //  Find a maximum for energy range in DOS for conduction band
		//cout<<"dum[][] = "<<dum[i][0]<<endl; 
	}

	end_index = min_index(dum,num_points,1);    //  Find a maximum for energy range in DOS for valence band
	start_index = min_index(DOS,num_points,1);   // Find the closest point to CBM, since CBM is at 0 energy
	//cout<<"start_index = "<<start_index<<"    end_index = "<<end_index<<endl;

	if(VASP==1)
	{
		while (DOS[start_index][1] < 1e-6 )
		{
		    end_index = end_index + 1;
			start_index = start_index + 1;
		}
	}
	
	//cout<<"DOS[start_index][0] = "<<DOS[start_index][0]<<endl;
	//cout<<"DOS[start_index+1][0] = "<<DOS[start_index+1][0]<<endl;

	double dummy1 = DOS[start_index][0];

	for(int i=0;i<num_points;i++)
	DOS[i][0] = DOS[i][0] - dummy1;

	//cout<<"DOS[start_index][0] = "<<DOS[start_index][0]<<endl;
	//cout<<"DOS[start_index+1][0] = "<<DOS[start_index+1][0]<<endl;



	
	//-------------------------------------------------------------------------	
	j=0;
	
	for(int i=start_index;i<=end_index;i++)
	{
		DOS_n[j][0] = DOS[i][0];
		DOS_n[j][1] = DOS[i][1];
		//cout<<"DOS[i][0 and 1] = "<<DOS[i][0]<<"   "<<DOS[i][1]<<endl;
		//getchar();
		j++;
		//cout<<"DOS[i][0] = "
	}
	count_DOS_n = j;


//----------------------------------------------------------------------------------------------------------------------------------


	if(VASP==1)
	{
		fid = fopen("DOS_n.dat","w");
		
		//volume1 unit is cm^3
		
		fprintf(fid,"# energy (eV)  	DOS (per eV per cm^3) \n");			
		for (int i = 0; i < count_DOS_n; i++)
		fprintf(fid,"%d     %e     %e \n", i+1, DOS_n[i][0], DOS_n[i][1]/volume1);

		fid = fopen("DOS_p.dat","w");
		fprintf(fid,"# energy (eV)  	DOS  (per eV per cm^3) \n");			
		for (int i = 0; i < count_DOS_p; i++)
		fprintf(fid,"%d     %e     %e \n", i+1, DOS_p[i][0], DOS_p[i][1]/volume1);
	}
	
	else if(VASP==2)  // FOR tight binding material 2D
	{
		fid = fopen("DOS_n.dat","w");
		
		
		fprintf(fid,"# energy (eV)  	DOS (per eV per cm^2) \n");			
		for (int i = 0; i < count_DOS_n; i++)
		fprintf(fid,"%d     %e     %e \n", i+1, DOS_n[i][0], DOS_n[i][1]);

		fid = fopen("DOS_p.dat","w");
		fprintf(fid,"# energy (eV)  	DOS  (per eV per cm^2) \n");			
		for (int i = 0; i < count_DOS_p; i++)
		fprintf(fid,"%d     %e     %e \n", i+1, DOS_p[i][0], DOS_p[i][1]);
	}

	
	if(save_data==1)
	{	
		fid = fopen("DOS.dat","w");

		fprintf(fid,"# energy (eV)  	DOS (per eV per cm^3)\n");			
		for (int i = 0; i < num_points; i++)
			fprintf(fid," %e     %e \n", energy[i], states[i]/volume1);
	}

	//*/
	
	/*
	for(int i=0;i<count_DOS_n;i++)
	cout<<DOS_n[i][0]<<"   "<<DOS_n[i][1]<<endl;

	//cout<<"DOS_n =  "<<"    "<<"count_n = "<<count_DOS_n<<endl;

	for(int i=0;i<count_DOS_p;i++)
	cout<<DOS_p[i][0]<<"   "<<DOS_p[i][1]<<endl;
	cout<<"DOS_p =  "<<"    "<<"count_p = "<<count_DOS_p<<endl;
	*/
    
    
    
}
