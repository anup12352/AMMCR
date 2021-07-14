#include"main.h"

int spin_orbit_coupling;
double ion_mass1[5];
double lm[10][10];
int ion_numbers1[5];
double volume1;
double c_lattice;
int number;

void read_OUTCAR()
{
	char line[1000], temp[100], temp1[100], temp2[100];
	char ion_mass[100], ion_numbers[100], volume[100], lattice_constant[100], c_a_ratio[100],str[100];
	double lattice_constant1, c_a_ratio1;	
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


				strncpy(temp1, line, 12);
				temp1[12] = '\0';

				//printf("%s", temp1); getchar();
				if (strcmp(temp1, "   POMASS = ")== 0)
                		{
					if (line[20] != ';')
					{
						count++;
						for (i = 12; i < strlen(line); i++)
						{
							ion_mass[i - 12] = line[i];
						}
						ion_mass[i-12] = '\0';
					}
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
	{	sscanf(ion_numbers, "%d", &ion_numbers1[0]);
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
	sscanf(volume, "%lf", &volume1);
}
