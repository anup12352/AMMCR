#include "main.h"
int countx;
double cond_band[2000][2],val_band[2000][2];
void make_band(int type)
{
	char line[1000];
	int VBM_band_number, CBM_band_number,dummy;
	double temp, temp1[4],temp2;
	int num1;

	FILE *fid,*fid1;

	if (type==1) // 1 represent conduction band
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

    double energies[NKPTS][2], k_points[NKPTS][3];

    if (ispin == 1)
    {
        for (int i = 0; i < NKPTS; i++)
        {
            fgets(line, 1000, (FILE*)fid);   //passing blank line
            fgets(line, 1000, (FILE*)fid);   // reading next kpoint line
            sscanf(line, "%lf %lf %lf %lf", &temp1[0], &temp1[1], &temp1[2], &temp1[3]);

            k_points[i][0] = (temp1[0] * lm[0][3] + temp1[1] * lm[1][3] + temp1[2] * lm[2][3]) * 2. * 3.14159265359 * 10;
            k_points[i][1] = (temp1[0] * lm[0][4] + temp1[1] * lm[1][4] + temp1[2] * lm[2][4]) * 2. * 3.14159265359 * 10;
            k_points[i][2] = (temp1[0] * lm[0][5] + temp1[1] * lm[1][5] + temp1[2] * lm[2][5]) * 2. * 3.14159265359 * 10;



            for (int j = 0; j < VBM_band_number; j++)
                fgets(line, 1000, (FILE*)fid);

            sscanf(line, "%d %lf", &num1, &temp);   // save values of valence band only
            energies[i][0] = temp;   // valence band energy

            fgets(line, 1000, (FILE*)fid);  // read conduction band
            sscanf(line, "%d %lf", &num1, &temp);

            energies[i][1] = temp;     // conduction band energy

            //cout<<"k_points = "<<k_points[i][0]<<"   "<<k_points[i][1]<<"   "<<k_points[i][2]<<endl;
            //cout<<energies[i][0]<<"   "<<energies[i][1]<<endl;
            //getchar();

            for (int j = VBM_band_number+1; j < NBTOT; j++)
                fgets(line, 1000, (FILE*)fid);
            //printf("\n %d %e %e %e %e", i, k_points[i][0], k_points[i][1], k_points[i][2], energies[i]);
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
			
            k_points[i][0] = (temp1[0] * lm[0][3] + temp1[1] * lm[1][3] + temp1[2] * lm[2][3]) * 2. * 3.14159265359 * 10;
            k_points[i][1] = (temp1[0] * lm[0][4] + temp1[1] * lm[1][4] + temp1[2] * lm[2][4]) * 2. * 3.14159265359 * 10;
            k_points[i][2] = (temp1[0] * lm[0][5] + temp1[1] * lm[1][5] + temp1[2] * lm[2][5]) * 2. * 3.14159265359 * 10;


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
	   	
            //cout<<"k_points = "<<k_points[i][0]<<"   "<<k_points[i][1]<<"   "<<k_points[i][2]<<endl;
            //cout<<energies[i][0]<<"   "<<energies[i][1]<<endl;
            //getchar();

            for (int j = VBM_band_number+1; j < NBTOT; j++)
                fgets(line, 1000, (FILE*)fid);
            //printf("\n %d %e %e %e %e", i, k_points[i][0], k_points[i][1], k_points[i][2], energies[i]);
        }
    }

    fclose(fid);

    double temp_reference[3], reference_point[3];

    if (type==1)
    {
        temp_reference[0] = kcbm[0];
        temp_reference[1] = kcbm[1];
        temp_reference[2] = kcbm[2];
    }
    else
    {
        temp_reference[0] = kvbm[0];
        temp_reference[1] = kvbm[1];
        temp_reference[2] = kvbm[2];
    }

    double sum;
    //printf("\n Reference point = ");
    for (int i = 0; i < 3; i++)
    {
        sum = 0;
        for (int j = 0; j < 3; j++)
        {
            sum = sum + (temp_reference[j] * lm[j][i+3]*10 );
        }
        reference_point[i] = sum* 2.0 * 3.14159265359;
        //cout<<"reference_point[i] = "<<reference_point[i]<<endl;
        //printf(" %e", reference_point[i]);
    }
    //getchar();

    if (type==1)
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
        double EVBM = -1000;
        for (int i = 0; i < NKPTS; i++)
        {
            if (energies[i][0] > EVBM)
                EVBM = energies[i][0];    // energies[][1]  contains conduction band
        }
        VBM = EVBM;
    }
    /*
    cout<<"Inside make_band.cpp "<<endl;
    if (type == 1)
    	printf("\n CBM = %lf \n", CBM);
    else
    	printf("\n VBM = %lf \n", VBM);
    //getchar();
    */	
	
    double conduction_band[NKPTS][2];
    double distance, kx, ky, kz;
    for (int i = 0; i < NKPTS; i++)
    {
        kx = pow(k_points[i][0] - reference_point[0], 2);
        ky = pow(k_points[i][1] - reference_point[1], 2);
        kz = pow(k_points[i][2] - reference_point[2], 2);
        distance = sqrt(kx + ky + kz);
        conduction_band[i][0] = distance;
        conduction_band[i][1] = energies[i][1];
        //printf("\n %e %lf", distance,energies[i]);
    }
    if (type == 2)
    {
        for (int i = 0; i < NKPTS; i++)
            conduction_band[i][1] = - energies[i][0];
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


    sort(conduction_band1, conduction_band1+NKPTS);
    sort(conduction_band2, conduction_band2+NKPTS);

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

    double normalization_factor = 0 - conduction_band_num_dum[0][1];
    for (int i = 0; i < countx+1; i++)
    {
        conduction_band_num[i][1] = conduction_band_num[i][1] + normalization_factor;
    }

    //FILE *fid1;
    /*
    if (type==1)
        fid1 = fopen("conduction_band_num.txt","w");
    else
        fid1 = fopen("valence_band_num.txt","w");

    for (int i = 0; i < countx+1; i++)
    {
        fprintf(fid1,"\n%d  %e  %e", i+1, conduction_band_num[i][0], conduction_band_num[i][1]);
    }
    fclose(fid);
    fclose(fid1);
    */
    countx = countx+1;

    if (type==1)
    {
        for (int i=0;i<countx;i++)
        {
            cond_band[i][0] = conduction_band_num[i][0];
            cond_band[i][1] = conduction_band_num[i][1];
        }
    }
    else
    {
        for (int i=0;i<countx;i++)
        {
            val_band[i][0] = conduction_band_num[i][0];
            val_band[i][1] = conduction_band_num[i][1];
        }
    }

}
