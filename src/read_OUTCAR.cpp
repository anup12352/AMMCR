#include"main.h"

void read_OUTCAR()
{
	char line[1000], temp[100], temp1[100], temp2[100], temp3[200], dummy;
	char ion_mass[100], ion_numbers[100], volume[100], lattice_constant[100], c_a_ratio[100],str[100];
	double lattice_constant1, c_a_ratio1;
	int ions;	
	FILE *fid;
	fid = fopen("OUTCAR", "r");
	int flag = 1, data_num=0, i;
    number = 0;
    cout<<"\n Inside read_OUTCAR\n";

	if (fid == NULL)
	{
        	cout<<"OUTCAR file is not present. Exit from program";
        	exit(EXIT_FAILURE);
	}

			
		
	int count = -1;

	while (flag)
	{
		while (fgets(line, 1000, (FILE*)fid))
		{
			//printf("%s \n", line);
			//printf("\n %lu \n", strlen(line));
			if (strlen(line) > 18)
            		{	
				//cout<<"SSSSSSSSSSS"<<endl;
				strncpy(temp, line, 18);
				temp[18] = '\0';
				//printf("%s\n", temp); getchar();
				if (strcmp(temp, "   LSORBIT =      ")==0)
                		{
					if (line[18] == 'T')
                    			{
						spin_orbit_coupling = 1;
					}
					else
                    			{
						spin_orbit_coupling = 0;
					}
				}
				
				if (strcmp(temp, "   number of dos  ")==0)
				{
					//cout<<"found"<<endl;
					sscanf(line, "%s %s %s %s %s %s %s %s %s %s %s %d", &dummy, &dummy, &dummy, &dummy, &dummy,
					&dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &ions);
					cout<<"Total no. of ions = "<<ions<<endl;
					//getchar();
				}
				
				strncpy(temp, line, 14);
				temp[14] = '\0';

				//printf("hereee    %s\n", temp); //getchar();
				//cout<<"ttt"<<temp[0]<<"ttt"<<temp[1]<<"ttt"<<temp[2]<<"ttt"<<endl;
				//cout<<"ttt"<<temp[3]<<"ttt"<<temp[4]<<"ttt"<<temp[5]<<"ttt"<<endl;
				//cout<<"ttt"<<temp[6]<<"ttt"<<temp[7]<<"ttt"<<temp[8]<<"ttt"<<endl;
				//cout<<"ttt"<<temp[9]<<"ttt"<<temp[10]<<"ttt"<<temp[11]<<"ttt"<<endl;
				//cout<<"ttt"<<temp[12]<<"ttt"<<temp[13]<<"ttt"<<temp[14]<<"ttt"<<endl;
				//cout<<"aaaaaaaaaaa"<<endl;
				if (scattering_mechanisms[6]==1)   // dislocation scattering 
				{	

					if (strcmp(temp, " ALAT       = ") == 0)   // 13 elememts upto equal to including =
					{	
						//cout<<"Inside dislocation "<<endl;
						//cout<<"Enter some value for getchar"<<endl; getchar();
						//cout<<"temp = "<<temp<<endl;
						
						for (i = 14; i < strlen(line); i++)
						{
							lattice_constant[i - 14] = line[i];
						}
						lattice_constant[strlen(line) - 14] = '\0';

						sscanf(lattice_constant, "%lf", &lattice_constant1);
						double aa = lattice_constant1/10;
						cout<<"lattice_constant a = "<<aa<<" nm"<<endl;
					
						fgets(line, 1000, (FILE*)fid);    // reading next line for c-a ratio
						strncpy(temp, line, 14);
						temp[14] = '\0';
						//printf("%s\n", temp); //getchar();

						for (i = 14; i < strlen(line); i++)
						{
							c_a_ratio[i - 14] = line[i];
						}
						c_a_ratio[strlen(line) - 14] = '\0';					
																				
						sscanf(c_a_ratio, "%lf", &c_a_ratio1);
						cout<<"c_by_a_ratio = "<<c_a_ratio1<<endl;
						c_lattice = (c_a_ratio1 * lattice_constant1)/10;  // divided with 10 to convert A into nm  
						cout<<"c_lattice = "<<c_lattice<<"nm"<<endl; 
						//cout<<"Enter some value for getchar"<<endl; getchar();
					}
				}


				//strncpy(temp1, line, 12);
				//temp1[12] = '\0';
				//printf("%s", temp1); getchar();
				//if (strcmp(temp1, "   POMASS = ")== 0)
				
				strncpy(temp1, line, 20);
				temp1[20] = '\0';
				
				if (strcmp(temp1,"  Mass of Ions in am")==0)
                		{
                			//printf("%s", temp1); getchar();
                			fgets(line, 1000, (FILE*)fid);  // read next line  	
					count++;
					for (i = 12; i < strlen(line); i++)
					{
						ion_mass[i - 12] = line[i];
					}
					ion_mass[i-12] = '\0';
				}


				strncpy(temp2, line, 18);
				temp2[18] = '\0';

				if (strcmp(temp2, "   ions per type =")==0)
                		{
					for (i = 18; i < strlen(line); i++)
					{
						float d = strlen(line);
						ion_numbers[i - 18] = line[i];
					}
					ion_numbers[strlen(line) - 18] = '\0';
				}

				if (strcmp(temp2, "  volume of cell :")==0)
                		{
					data_num = data_num + 1;
					for (i = 18; i < strlen(line); i++)
					{
						volume[i - 18] = line[i];
					}
					volume[i - 18] = '\0';

					fgets(line, 1000, (FILE*)fid);
					//fgets(line, 1000, (FILE*)fid);
					//fgets(line, 1000, (FILE*)fid);

					//printf("\n lattice vectors:");
					//printf("\n %s", line);
					fgets(line, 1000, (FILE*)fid);
					sscanf(line, "%lf %lf %lf  %lf %lf %lf", &lm[0][0], &lm[0][1], &lm[0][2],
						&lm[0][3], &lm[0][4], &lm[0][5]);
					//printf("\n %s", line);
					fgets(line, 1000, (FILE*)fid);
					sscanf(line, "%lf %lf %lf  %lf %lf %lf", &lm[1][0], &lm[1][1], &lm[1][2],
						&lm[1][3], &lm[1][4], &lm[1][5]);

					//printf("\n %s", line);
					fgets(line, 1000, (FILE*)fid);
					sscanf(line, "%lf %lf %lf  %lf %lf %lf", &lm[2][0], &lm[2][1], &lm[2][2],
						&lm[2][3], &lm[2][4], &lm[2][5]);

					//printf("\n %s", line);
					flag = 0;
				}
				if(geometry==2)   // 2D material
				{
					strncpy(temp3, line, 42);
					temp3[42] = '\0';
					if (strcmp(temp3, " position of ions in cartesian coordinates")==0)
					{
						// reading ions position z ccordinate
						double z[ions], max=-1e6,min = 1e6, dummy1; 
						for(int i=0;i<ions;i++)
						{
							fgets(line, 1000, (FILE*)fid);
							sscanf(line, "%lf %lf %lf", &dummy1, &dummy1, &z[i]);
							//cout<<"z[i] = "<<z[i]<<endl;
						}
						// finding max and min z coordinate
						for(int i=0;i<ions;i++)
						{
							if(z[i]>max)					
								max = z[i];
							if(z[i]<min)
								min = z[i];
						}
						thickness = max - min;
						thickness = thickness*1e-10; // converted from A to m	
						//cout<<"max = "<<max<<"   min = "<<min<<"  thickness = "<<thickness<<"  m "<<endl;
						cout<<"  Thickness along z direction = "<<thickness<<"  m "<<endl;
						//getchar();
					}
				}			
			}
		}

	}
	fclose(fid);

	//cout<<"Outside looop"<<endl;
	//printf("%lf", ion_mass1[0]); getchar();
	int c=0;
	i=0;
	while(ion_numbers[i]!='\0' && ion_numbers[i]!=char(10) )
    	{
		if(ion_numbers[i] != char(32) )
		    number++;
		i++;
    	}
	
	
	
	if (number==1)
	{	
		sscanf(ion_numbers, "%d", &ion_numbers1[0]);
	 	sscanf(ion_mass, "%lf", &ion_mass1[0]);
	}
	else if (number==2)
	{
		sscanf(ion_numbers, "%d  %d", &ion_numbers1[0], &ion_numbers1[1]);
		sscanf(ion_mass, "%lf %lf", &ion_mass1[0], &ion_mass1[1]);
	}        
	else if (number==3)
	{
		sscanf(ion_numbers, "%d %d %d ", &ion_numbers1[0], &ion_numbers1[1], &ion_numbers1[2]);
		sscanf(ion_mass, "%lf %lf %lf", &ion_mass1[0], &ion_mass1[1], &ion_mass1[2]);	
	}
	else if (number==4)
	{
		    sscanf(ion_numbers, "%d %d %d %d", &ion_numbers1[0], &ion_numbers1[1], &ion_numbers1[2], &ion_numbers1[3]);
		sscanf(ion_mass, "%lf %lf %lf %lf", &ion_mass1[0], &ion_mass1[1], &ion_mass1[2],&ion_mass1[3]);			
	}        
	else if(number==5)
	{
		sscanf(ion_numbers, "%d %d %d %d %d", &ion_numbers1[0], &ion_numbers1[1], &ion_numbers1[2], &ion_numbers1[3], &ion_numbers1[4]);
    		sscanf(ion_mass, "%lf %lf %lf %lf %lf",&ion_mass1[0], &ion_mass1[1],&ion_mass1[2],&ion_mass1[3],&ion_mass1[4]);			
	}
	else
	{
    		exit(EXIT_FAILURE);
	}
	//*/
	
	/*
	for(int i=0;i<numbers-1;i++)
	{
		sscanf(ion_numbers, "%d", &ion_numbers1[i]);
		sscanf(ion_mass, "%lf", &ion_mass1[i]);
	}
	*/
	
	sscanf(volume, "%lf", &volume1);

//------------------------// showing outcar data ---------------------------------------------------------------------------------
//------------------------// showing outcar data ---------------------------------------------------------------------------------
    			
    //cout<<"number = "<<number<<endl;

    if (number==1)
    {
        printf("\nion_mass =  %lf \n", ion_mass1[0] );   // unit amu
        printf("\nion_numbers =  %d \n", ion_numbers1[0] );
    }
    else if (number==2)
    {
        printf("\nion_mass =  %lf %lf \n", ion_mass1[0], ion_mass1[1]);   // unit amu
        printf("\nion_numbers =  %d %d\n", ion_numbers1[0], ion_numbers1[1]);
    }
    else if (number==3)
    {
        printf("\nion_mass =  %lf %lf %lf \n", ion_mass1[0], ion_mass1[1],ion_mass1[2]);   // unit amu
        printf("\nion_numbers =  %d %d %d \n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2]);
    }
    else if (number==4)
    {
        printf("\nion_mass =  %lf %lf %lf %lf\n", ion_mass1[0], ion_mass1[1],ion_mass1[2], ion_mass1[3]);   // unit amu
        printf("\nion_numbers =  %d %d %d %d\n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2], ion_numbers1[3]);
    }
    else if (number==5)
    {
        printf("\nion_mass =  %lf %lf %lf %lf %lf \n", ion_mass1[0], ion_mass1[1],ion_mass1[2], ion_mass1[3], ion_mass1[4]);     // unit amu
        printf("\nion_numbers =  %d %d %d %d %d \n", ion_numbers1[0], ion_numbers1[1], ion_numbers1[2], ion_numbers1[3], ion_numbers1[4]);
    }
    	
    	volume1 = volume1/1e24;   // converted from A^3 to cm^3
	printf("\nvolume =  %e cm^(-3)\n", volume1);
	//getchar();
	// unit 1/cm^3  volume = a^3/4   (a -> lattice constant)

    //cout<<"spin_orbit_coupling = "<<spin_orbit_coupling<<endl;
    if(spin_orbit_coupling == 0)
        cout<<endl<<"spin_orbit_coupling = false"<<endl;
    if(spin_orbit_coupling == 1)
        cout<<endl<<"spin_orbit_coupling = true"<<endl;

	cout<<"Lattice matrix : "<<endl;

	for (int i = 0; i < 3; i++)
        {
		for (int j = 3; j < 6; j++)
        	{
			printf(" %lf", 10 * lm[i][j]);
		}
		printf("\n");
	}

	
	if (rho==0)
    	{
		double sum=0;
		for (int i = 0; i < number; i++)
		{
		    sum += (ion_mass1[i] * ion_numbers1[i]);
			//cout<<"sum = "<<sum<<endl;
		}
		rho = (sum*1.67377e-27) / (volume1*1e-6);    // unit kg/m^3
		if(geometry==1)   // for 3D bulk
	    		cout<<endl<<"Calculated density = "<<rho/1000<<" g/cm^3"<<endl;
	        else  // for 2D bulk
	        {
	        	rho = rho * thickness;  // unit kg/m^2
	        	cout<<endl<<"Calculated density = "<<rho/10<<" g/cm^2"<<endl;
	        }	   
	}

	if (scattering_mechanisms[6]==1)   // dislocation scattering 
	{
		cout<<endl<<"c_lattice constant = "<<c_lattice<<" nm"<<endl;       		
	}	
	//	
	//return 0;

	//getchar();
	
    
//----------------// showing outcar data completed  -----------------------------------------------------------------------
	
}
