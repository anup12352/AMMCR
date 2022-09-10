
#include "main.h"

#include<iostream>
using namespace std;

double orbital_decomposedd_p[1000][4];

int decompose_wave_p()
{
	FILE *fid, *fid2;
	int num_kpoints;
	char line[1000],data[100];
    
	double kpoints_BZ[limit3][3], s_total[limit3]={0}, p_total[limit3]={0}, d_total[limit3]={0};

	if(VASP==1)		
	{
		int num_bands1, num_ions,band_number;
		band_number = NBVAL;


		fid = fopen("PROCAR_p", "r");
		if (fid==NULL)
		fid = fopen("PROCAR","r");

		if (fid==NULL)
		{
		cout<<"PROCAR is not present. Program is running by assuming valence band to be p-like"<<endl;
		return 0;
		}
		fgets(line, 1000, fid);
		fgets(line, 1000, fid);
		//cout<<"line = "<<line<<endl;

		int i=0,j,l;
		l = strlen(line);
		while(line[i]!=':')
		{
		i++;
		}
		i++;
		j=i;
		while(line[i]!='#')
		{
		data[i-j] = line[i];
		i++;
		}

		sscanf(data, "%d", &num_kpoints);

		while(line[i]!=':')
		i++;
		i++;
		j=i;
		while(line[i]!='#')
		{
		data[i-j] = line[i];
		i++;
		}

		sscanf(data, "%d", &num_bands1);

		while(line[i]!=':')
		i++;
		i++;
		j=i;
		while(i!=l)
		{
		data[i-j] = line[i];
		i++;
		}

		sscanf(data, "%d", &num_ions);
		//cout<<"num_kpoints = "<<num_kpoints<<endl;
		//cout<<"num_bands1 =  "<<num_bands1<<endl;
		//cout<<"num_ions = "<<num_ions<<endl;
		//getchar();


		double **kpoints = new double*[num_kpoints];
		for (int i = 0; i < num_kpoints; i++)
		kpoints[i] = new double[4];


		double **band_energies = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		band_energies[i] = new double[num_kpoints];

		double ***s = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		s[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    s[i][j] = new double[num_kpoints];

		double **s_total1 = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		s_total1[i] = new double[num_kpoints];

		double ***py = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		py[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    py[i][j] = new double[num_kpoints];

		double **py_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		py_total[i] = new double[num_kpoints];

		double ***pz = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		pz[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    pz[i][j] = new double[num_kpoints];

		double **pz_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		pz_total[i] = new double[num_kpoints];

		double ***px = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		px[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    px[i][j] = new double[num_kpoints];

		double **px_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		px_total[i] = new double[num_kpoints];

		double **p = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		p[i] = new double[num_kpoints];

		double **p_total1 = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		p_total1[i] = new double[num_kpoints];

		double ***dxy = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		dxy[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    dxy[i][j] = new double[num_kpoints];

		double **dxy_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		dxy_total[i] = new double[num_kpoints];

		double ***dyz = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		dyz[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    dyz[i][j] = new double[num_kpoints];

		double **dyz_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		dyz_total[i] = new double[num_kpoints];

		double ***dz2 = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		dz2[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    dz2[i][j] = new double[num_kpoints];

		double **dz2_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		dz2_total[i] = new double[num_kpoints];

		double ***dxz = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		dxz[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    dxz[i][j] = new double[num_kpoints];

		double **dxz_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		dxz_total[i] = new double[num_kpoints];


		double ***dx2 = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		dx2[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    dx2[i][j] = new double[num_kpoints];

		double **dx2_total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		dx2_total[i] = new double[num_kpoints];

		double **d = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		d[i] = new double[num_kpoints];

		double **d_total1 = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		d_total1[i] = new double[num_kpoints];

		double ***ion_total = new double**[num_ions];
		for (int i = 0; i < num_ions; i++)
		ion_total[i] = new double*[num_bands1];
		for (int i = 0; i < num_ions; i++)
		for(int j=0; j< num_bands1; j++)
		    ion_total[i][j] = new double[num_kpoints];

		double **total = new double*[num_bands1];
		for (int i = 0; i < num_bands1; i++)
		total[i] = new double[num_kpoints];

		int nb=0,dummy;

		for(int k=0;k<num_kpoints;k++)
		{
		//cout<<"kpoint no. = "<<k<<endl;
		fgets(line, 1000, fid);
		fgets(line, 1000, fid);
		//cout<<"line = "<<line<<endl;
		//getchar();
		i=0;
		l = strlen(line);
		while(line[i]!=':')
		    i++;
		i++;

		j=i;

		while(line[i]!='w')
		{
		    data[i-j] = line[i];
		    i++;
		}

		sscanf(data, "%lf %lf %lf", &kpoints[k][0], &kpoints[k][1], &kpoints[k][2]);

		while(line[i]!='=')
		    i++;
		i++;
		j=i;
		while(i!=l)
		{
		    data[i-j] = line[i];
		    i++;
		}

		sscanf(data, "%lf ", &kpoints[k][3]);

		//cout<<"kpoints = "<<kpoints[k][0]<<"   "<<kpoints[k][0]<<"   "<<kpoints[k][1]
		//<<"   "<<kpoints[k][2]<<"   "<<kpoints[k][2]<<endl;
		//getchar();

		for (nb=0;nb<num_bands1;nb++)
		{
		    //cout<<"band number =  "<<nb<<endl;
		    fgets(line, 1000, fid);
		    fgets(line, 1000, fid);
		    i =19;
		    j=i;
		    while(line[i]!='#')
		    {
			data[i-j]=line[i];
			i++;
		    }
		    data[i-j]='\0';
		    sscanf(data, "%lf ", &band_energies[nb][k]);
		    //cout<<"energy = "<<band_energies[nb][k]<<"     end"<<endl;
		    //getchar();
		    fgets(line, 1000, fid);
		    fgets(line, 1000, fid);
		    fgets(line, 1000, fid);
		    //cout<<line<<endl<<"ssss";
		    //getchar();
		    for(i=0;i<num_ions;i++)
		    {
			sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy, &s[i][nb][k],
			   &py[i][nb][k], &pz[i][nb][k], &px[i][nb][k], &dxy[i][nb][k], &dyz[i][nb][k], &dz2[i][nb][k],
			   &dxz[i][nb][k], &dx2[i][nb][k], &ion_total[i][nb][k]);

			/*
			cout<<dummy<<"  "<<s[i][nb][k]<<"  "<<py[i][nb][k]<<"  "<<pz[i][nb][k]<<"  "<<px[i][nb][k]<<"   "<<dxy[i][nb][k]
			 <<"  "<<dyz[i][nb][k]<<"  "<<dz2[i][nb][k]<<"  "<<dxz[i][nb][k]<<"  "<<dx2[i][nb][k]<<"  "<<ion_total[i][nb][k]<<endl;
			getchar();
			*/
			//cout<<"check  "<<endl;
			fgets(line, 1000, fid);
			//cout<<line<<"xxxx"<<endl;
			//getchar();
		    }
		    l = strlen(line);
		    for (i=3;i<l;i++)
			data[i-3] = line[i];
		    sscanf(data, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &s_total1[nb][k], &py_total[nb][k], &pz_total[nb][k],
			   &px_total[nb][k], &dxy_total[nb][k], &dyz_total[nb][k], &dz2_total[nb][k],
			   &dxz_total[nb][k], &dx2_total[nb][k], &total[nb][k]);

		    p_total1[nb][k] = py_total[nb][k]+ px_total[nb][k] + pz_total[nb][k];
		    d_total1[nb][k] = dxy_total[nb][k]+ dyz_total[nb][k] + dz2_total[nb][k] + dxz_total[nb][k] + dx2_total[nb][k];

		    //cout<<"Here total"<<endl;
		    //cout<<s_total1[nb][k]<<"    "<<p_total1[nb][k]<<"    "<<d_total1[nb][k]<<"    "<<endl;
		    //getchar();

		    if (spin_orbit_coupling == 1)
		    {
			//cout<<"ssshoww"<<endl;
			for (int skp=1;skp<=3*(num_ions+1);skp++)
			    fgets(line, 1000, fid);
		    }
		}
		fgets(line, 1000, fid);
		}


		//cout<<"test";
		//getchar();

		//cout<<"check";

		for(int i=0;i<num_kpoints;i++)
		{
		kpoints_BZ[i][0] = (kpoints[i][0] * lm[0][3] + kpoints[i][1] * lm[1][3] + kpoints[i][2] * lm[2][3]) * 2. * 3.14159265359 * 10;
		kpoints_BZ[i][1] = (kpoints[i][0] * lm[0][4] + kpoints[i][1] * lm[1][4] + kpoints[i][2] * lm[2][4]) * 2. * 3.14159265359 * 10;
		kpoints_BZ[i][2] = (kpoints[i][0] * lm[0][5] + kpoints[i][1] * lm[1][5] + kpoints[i][2] * lm[2][5]) * 2. * 3.14159265359 * 10;
		//cout<<kpoints_BZ[i][0]<<"    "<<kpoints_BZ[i][1]<<"    "<<kpoints_BZ[i][2]<<endl;
		//getchar();
		}
		
		for(int i=0;i<num_kpoints;i++)
		{
			s_total[i] = s_total1[band_number-1][i];
			p_total[i] = p_total1[band_number-1][i];
			d_total[i] = d_total1[band_number-1][i];
		}
		/*
		//--------------------------------------- save wave function admixture ---------------------------------
		fid2 = fopen("WAVE_ADMIXTURE_VB.dat","w");
		fprintf(fid2,"# kx   ky    kz    s	p  \n");
		for (int i = 0; i < num_kpoints; i++)
			fprintf(fid2,"%e  %e  %e   %e  %e \n", kpoints_BZ[i][0], kpoints_BZ[i][1], kpoints_BZ[i][2], s_total[i], p_total[i]);
		
		fclose(fid2);				
		//--------------------------------------------- wave function admixture saved ------------------------------	
		*/
		
	}	
	else if(VASP==0)    // reading wave function admixture from table form
	{
		
		fid = fopen("WAVE_ADMIXTURE_VB.dat","r");
		if (fid==NULL)
		{
			cout<<"WAVE_ADMIXTURE_VB.dat file is not present";
			cout<<"Program is running by assuming valnece band to be p-like"<<endl;
			return 0;
		}

		fgets(line, 1000, fid);   // pass first line
		int i=0;
		while (fgets(line, 1000, fid)!= NULL)
		{
		       //fgets(line, 1000, fid);
		       sscanf(line, "%lf %lf %lf %lf %lf ", &kpoints_BZ[i][0], &kpoints_BZ[i][1], &kpoints_BZ[i][2],                      					&s_total[i], &p_total[i]);
			i++;
		}
		num_kpoints = i;
				
	}    
	else if(VASP==2)   // tight binding band strructure
	{
		//cout<<"For tight binding band structure program is running by assuming conduction band to be s-like and valence band to be p-like"<<endl;
		return 0;
		// In generate_required_data,.cpp if  count_orbital_n and count_orbital_p are zero then it will make conduction band to be s-like and valence band to be p-like
		// By default count_orbital_p = 0 and count_orbital_n = 0
		
	}
	// if and else condition for VASP==1 completed
		
	// if and else condition for VASP==1 completed
        //--------------------------------------------------------------------------------------------------------------------------


	double temp_reference[3], reference_point[3];

	temp_reference[0] = kvbm[0];
	temp_reference[1] = kvbm[1];
	temp_reference[2] = kvbm[2];

	//cout<<"reached here";
	double sum;
	//printf("\n Reference point = ");
	for (int i = 0; i < 3; i++)
	{
		reference_point[i] = temp_reference[i];
	//cout<<"reference_point[i] = "<<reference_point[i]<<endl;
	//printf(" %e", reference_point[i]);
	}
	//getchar();


	double **orbital_decomposed = new double*[num_kpoints];
	for (int i = 0; i < num_kpoints; i++)
	orbital_decomposed[i] = new double[4];

	//cout<<"num_kpoints = "<<num_kpoints<<endl;
	//getchar();
	double kx,ky,kz,distance;
	for (int i = 0; i < num_kpoints; i++)
	{
	kx = pow(kpoints_BZ[i][0] - reference_point[0], 2);
	ky = pow(kpoints_BZ[i][1] - reference_point[1], 2);
	kz = pow(kpoints_BZ[i][2] - reference_point[2], 2);
	distance = sqrt(kx + ky + kz);
	orbital_decomposed[i][0] = distance;
	orbital_decomposed[i][1] = s_total[i];
	orbital_decomposed[i][2] = p_total[i];
	orbital_decomposed[i][3] = d_total[i];
	//printf("\n %lf %lf %lf %lf", orbital_decomposed[i][0], orbital_decomposed[i][1], orbital_decomposed[i][2], orbital_decomposed[i][3]);
	//getchar();
	}
	//getchar();

	/*
	//------------------------------Readin  ----------------------------------------------------
	fid = fopen("orbital_decomposed1.txt", "r");
	if (fid==NULL)
	cout<<"File is not present "<<endl;
	for(int i=0;i<num_kpoints;i++)
	{
	    fgets(line, 1000, fid);
	    //cout<<line<<endl;
	    //getchar();
	sscanf(line, "%lf %lf %lf %lf", &orbital_decomposed[i][0], &orbital_decomposed[i][1],
	       &orbital_decomposed[i][2], &orbital_decomposed[i][3]);
	}

	cout<<"After Reading"<<endl;

	for(int i=0;i<num_kpoints;i++)
	{
	cout<<"i+1 = "<<i+1<<endl;
	cout.precision(8);
	cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]<<"    "<<
	orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<endl;
	getchar();
	}
	cout<<"After Reading"<<endl;
	*/

	//-----------------------------------------------------------------------------------
	// Sorting array orbital_decomposed according to distance from reference point
	int dum = num_kpoints;
	double dummy1;

	for(int i=0;i<dum-1;i++)
	{
	for (int j=0; j<dum-1; j++)
	{
	    if(orbital_decomposed[j][0] > orbital_decomposed[j+1][0])
	    {
		dummy1 = orbital_decomposed[j][0];
		orbital_decomposed[j][0] = orbital_decomposed[j+1][0];
		orbital_decomposed[j+1][0] = dummy1;

		dummy1 = orbital_decomposed[j][1];
		orbital_decomposed[j][1] = orbital_decomposed[j+1][1];
		orbital_decomposed[j+1][1] = dummy1;

		dummy1 = orbital_decomposed[j][2];
		orbital_decomposed[j][2] = orbital_decomposed[j+1][2];
		orbital_decomposed[j+1][2] = dummy1;

		dummy1 = orbital_decomposed[j][3];
		orbital_decomposed[j][3] = orbital_decomposed[j+1][3];
		orbital_decomposed[j+1][3] = dummy1;
	    }
	}
	}

	/*
	cout<<"After sorting"<<endl;

	for(int i=0;i<num_kpoints;i++)
	{
	cout<<"i+1 = "<<i+1<<endl;
	cout.precision(8);
	cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]<<"    "<<
	orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<endl;
	getchar();
	}
	cout<<"After sorting"<<endl;

	getchar();
	*/
	//------------------------- AFter sorting save orbital_decomposed ------------------------------
	/*
	fid = fopen("orbital_decomposed_sort.txt", "w");
	for(int i=0;i<num_kpoints;i++)
	{
	fprintf(fid,"%e    %e    %e    %e \n", orbital_decomposed[i][0], orbital_decomposed[i][1],
	orbital_decomposed[i][2], orbital_decomposed[i][3]);
	}
	fclose(fid);
	*/


	//------------------------Read orbital_decomposed from txt file after sorting--------------------------------------------------------
	/*
	fid = fopen("orbital_decomposed2.txt", "r");
	if (fid==NULL)
	cout<<"File is not present "<<endl;
	for(int i=0;i<num_kpoints;i++)
	{
	    fgets(line, 1000, fid);
	    //cout<<line<<endl;
	    //getchar();
	sscanf(line, "%lf %lf %lf %lf", &orbital_decomposed[i][0], &orbital_decomposed[i][1],
	       &orbital_decomposed[i][2], &orbital_decomposed[i][3]);
	}

	cout<<"After rreading sorting array from matlab "<<endl;
	for(int i=0;i<num_kpoints;i++)
	{
	cout<<"i+1 = "<<i+1<<endl;
	cout.precision(8);
	cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]<<"    "<<
	orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<endl;
	getchar();
	}
	*/
	//-----------Normalizing s and p coefficients
	double factor;
	//cout<<"After normalization = "<<endl;

	for(int i=0;i<num_kpoints;i++)
	{
		factor=1/(pow((orbital_decomposed[i][1]),2)+pow((orbital_decomposed[i][2]),2)+0.0000000001);
		orbital_decomposed[i][1]=orbital_decomposed[i][1]*pow(factor,0.5);
		orbital_decomposed[i][2]=orbital_decomposed[i][2]*pow(factor,0.5);

		//cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]<<"    "
		//<<orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<"    "<<endl;
		//getchar();

	}
	//getchar();

	double **orbital_decomposed_dum = new double*[num_kpoints];
	for (int i = 0; i < num_kpoints; i++)
	orbital_decomposed_dum[i] = new double[4];

	for(int i = 0; i < num_kpoints; i++)
	{
	orbital_decomposed_dum[i][0] = orbital_decomposed[i][0];    // distance
	orbital_decomposed_dum[i][1] = orbital_decomposed[i][1];
	orbital_decomposed_dum[i][2] = orbital_decomposed[i][2];
	orbital_decomposed_dum[i][3] = orbital_decomposed[i][3];
	}
	/*
	cout<<"After normalization"<<endl;

	for(int i=0;i<num_kpoints;i++)
	{
	cout<<"i+1 = "<<i+1<<endl;
	cout.precision(8);
	cout<<orbital_decomposed_dum[i][0]<<"    "<<orbital_decomposed_dum[i][1]<<"    "<<
	orbital_decomposed_dum[i][2]<<"    "<<orbital_decomposed_dum[i][3]<<endl;
	getchar();
	}
	cout<<"After norm"<<endl;

	getchar();
	*/

	orbital_decomposed[0][0] = orbital_decomposed_dum[0][0];    // distance
	orbital_decomposed[0][1] = orbital_decomposed_dum[0][1];
	orbital_decomposed[0][2] = orbital_decomposed_dum[0][2];
	orbital_decomposed[0][3] = orbital_decomposed_dum[0][3];


	double avg[1][4] = {0}, z = 0.0001;
	int start = 0, stop = 0;
	countx = 0;
	//cout<<"Before loop"<<endl;
	double a1,a2;
	for (int i = 1; i < num_kpoints - 1; i++)
	{
	//cout<<"i = "<<i<<endl;

	//cout<<"orbital-decopmose[i][0] = "<<orbital_decomposed_dum[i][0]<<endl;
	//cout<<"orbital-decopmose[i+1][0] = "<<orbital_decomposed_dum[i+1][0]<<endl;
	//cout<<"orbital-decopmose[i-1][0] = "<<orbital_decomposed_dum[i-1][0]<<endl;
	//a1 = orbital_decomposed_dum[i+1][0] - orbital_decomposed_dum[i][0];
	//a2 = orbital_decomposed_dum[i][0] - orbital_decomposed_dum[i-1][0];

	//cout<< "i+1 - i  =       "<< a1<<endl;
	//cout<< "i - i-1 =        "<< a2<<endl;
	//if(i==50||i==100||i==150||i==200||i==1000||i==2000||i==3000)
	 //   { cout<<"testing "<<endl; getchar(); }

	if (abs(orbital_decomposed_dum[i][0] - orbital_decomposed_dum[i - 1][0]) < z)
	{
	    if (start == 0)
		start = i - 1;

	    //cout<<"start = "<<start<<endl;

	    if ((orbital_decomposed_dum[i+1][0] - orbital_decomposed_dum[i][0]) > z)
	    {
		//cout<<"Before stop "<<endl;
		stop = i;
		for (int j = start; j <= stop; j++)
		{
		    avg[0][0] = avg[0][0]+orbital_decomposed_dum[j][0];
		    avg[0][1] = avg[0][1]+orbital_decomposed_dum[j][1];
		    avg[0][2] = avg[0][2]+orbital_decomposed_dum[j][2];
		    avg[0][3] = avg[0][3]+orbital_decomposed_dum[j][3];
		}
		avg[0][0] = avg[0][0]/(stop - start + 1);
		avg[0][1] = avg[0][1]/(stop - start + 1);
		avg[0][2] = avg[0][2]/(stop - start + 1);
		avg[0][3] = avg[0][3]/(stop - start + 1);

		countx++;
		orbital_decomposed[countx][0] = avg[0][0];
		orbital_decomposed[countx][1] = avg[0][1];
		orbital_decomposed[countx][2] = avg[0][2];
		orbital_decomposed[countx][3] = avg[0][3];

		//cout<<"In betweenn"<<endl;
		//cout<<orbital_decomposed[countx][0]<<"    "<<orbital_decomposed[countx][1]
		//<<"     "<<orbital_decomposed[countx][2]<<"    "<<orbital_decomposed[countx][3]<<endl;
		//cout<<"countx = "<<countx<<endl;

		start = 0;
		stop = 0;
		avg[0][0] = 0;
		avg[0][1] = 0;
		avg[0][2] = 0;
		avg[0][3] = 0;
	    }
	}
	if( (abs(orbital_decomposed_dum[i][0] - orbital_decomposed_dum[i - 1][0]) > z) &
	   ((orbital_decomposed_dum[i + 1][0] - orbital_decomposed_dum[i][0]) > z) )
	{
	    countx++;
	    orbital_decomposed[countx][0] = orbital_decomposed_dum[i][0];
	    orbital_decomposed[countx][1] = orbital_decomposed_dum[i][1];
	    orbital_decomposed[countx][2] = orbital_decomposed_dum[i][2];
	    orbital_decomposed[countx][3] = orbital_decomposed_dum[i][3];
	    //cout<<"countx++ inside next if = "<<countx++<<endl;
	}
	//getchar();

	}

	//cout<<"Reached outside   "<<endl;
	//getchar();

	if (orbital_decomposed_dum[num_kpoints-1][0] > orbital_decomposed_dum[num_kpoints-2][0])
	{
	countx++;
	orbital_decomposed[countx][0] = orbital_decomposed_dum[num_kpoints-1][0];
	orbital_decomposed[countx][1] = orbital_decomposed_dum[num_kpoints-1][1];
	orbital_decomposed[countx][2] = orbital_decomposed_dum[num_kpoints-1][2];
	orbital_decomposed[countx][3] = orbital_decomposed_dum[num_kpoints-1][3];
	}

	if ((orbital_decomposed[1][0] - orbital_decomposed[0][0]) < z)
	{
	for (int i = 1; i < countx+1; i++)
	{
	    orbital_decomposed[i - 1][0] = orbital_decomposed[i][0];
	    orbital_decomposed[i - 1][1] = orbital_decomposed[i][1];
	    orbital_decomposed[i - 1][2] = orbital_decomposed[i][2];
	    orbital_decomposed[i - 1][3] = orbital_decomposed[i][3];
	}
	countx--;
	}

	//FILE *fid1;
	//fid1 = fopen("orbital_decomposedd_p.txt","w");

	for (int i = 0; i < countx+1; i++)
	{
		orbital_decomposedd_p[i][0] = orbital_decomposed[i][0];
		orbital_decomposedd_p[i][1] = orbital_decomposed[i][1];
		orbital_decomposedd_p[i][2] = orbital_decomposed[i][2];
		orbital_decomposedd_p[i][3] = orbital_decomposed[i][3];


		/*fprintf(fid1,"%e  %e  %e  %e \n", orbital_decomposed[i][0], orbital_decomposed[i][1],
			orbital_decomposed[i][2], orbital_decomposed[i][3]);
			/*
			cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]<<"    "
			<<orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<"    "<<endl;
			getchar();
			*/
	}
	countx = countx+1;
	fclose(fid);
	//fclose(fid1);
	/*
	cout<<"At the end"<<endl;
	for (int i=0;i<countx;i++)
	{
	cout<<"i+1 = "<<i+1<<endl;
	cout<<orbital_decomposed[i][0]<<"    "<<orbital_decomposed[i][1]
	<<"     "<<orbital_decomposed[i][2]<<"    "<<orbital_decomposed[i][3]<<endl;
	getchar();
	}
	//*/

	return countx;
}

