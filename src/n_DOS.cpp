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
    int num_points;
    double dummy;
    double fermi;
    FILE *fid;
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

    //cout<<"num_points = "<<num_points<<endl;
    //cout<<"dos fermi = "<<fermi<<endl;

    double energy[num_points],states[num_points],states_spin[num_points]={0};
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
    
    double DOS[num_points][2],DOS_temp[num_points][2];
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
    max1 = max_value(val_band,count2,2);
    //cout<<"max value = valence band ="<<max1<<endl;
    double dum[num_points][2]={0};
    for(int i=0;i<num_points;i++)
        { dum[i][0] = DOS_temp[i][0]+ 1.1 * max1;}  // cout<<"dum[][] = "<<dum[i][0]<<endl; }
    //getchar();

    start_index = min_index(dum,num_points,1);    //  Find a maximum for energy range in DOS for valence band
    end_index = min_index(DOS_temp,num_points,1);   // Find the closest point to VBM Since VBM is at 0 energy

    //cout<<"start_index = "<<start_index<<"    end_index = "<<end_index<<endl;

    while (DOS_temp[end_index][1] < 1e-6)
    {
        end_index = end_index - 1;
        start_index = start_index - 1;
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



    max1 = max_value(cond_band,count1,2);
    //cout<<"max1 = conduction band = "<<max1<<endl;
    //cout<<"next";
    for(int i=0;i<num_points;i++)
        { dum[i][0] = DOS[i][0] - 1.1 * max1;  }    //  Find a maximum for energy range in DOS for conduction band
        //cout<<"dum[][] = "<<dum[i][0]<<endl; }

    end_index = min_index(dum,num_points,1);    //  Find a maximum for energy range in DOS for valence band
    start_index = min_index(DOS,num_points,1);   // Find the closest point to CBM, since CBM is at 0 energy
    //cout<<"start_index = "<<start_index<<"    end_index = "<<end_index<<endl;

    while (DOS[start_index][1] < 1e-6 )
	{
	    end_index = end_index + 1;
		start_index = start_index + 1;
	}

    //cout<<"DOS[start_index][0] = "<<DOS[start_index][0]<<endl;
    //cout<<"DOS[start_index+1][0] = "<<DOS[start_index+1][0]<<endl;

    double dummy1 = DOS[start_index][0];

    for(int i=0;i<num_points;i++)
        DOS[i][0] = DOS[i][0] - dummy1;

    //cout<<"DOS[start_index][0] = "<<DOS[start_index][0]<<endl;
    //cout<<"DOS[start_index+1][0] = "<<DOS[start_index+1][0]<<endl;

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

    /*
    FILE *fid1;
    fid1 = fopen("DOS_n.txt","w");

    for (int i = 0; i < count_DOS_n; i++)
        fprintf(fid1,"%d     %e     %e \n", i+1, DOS_n[i][0], DOS_n[i][1]);

    fid1 = fopen("DOS_p.txt","w");
    for (int i = 0; i < count_DOS_p; i++)
        fprintf(fid1,"%d     %e     %e \n", i+1, DOS_p[i][0], DOS_p[i][1]);


    /*
    for(int i=0;i<count_DOS_n;i++)
        cout<<DOS_n[i][0]<<"   "<<DOS_n[i][1]<<endl;

    //cout<<"DOS_n =  "<<"    "<<"count_n = "<<count_DOS_n<<endl;

    for(int i=0;i<count_DOS_p;i++)
        cout<<DOS_p[i][0]<<"   "<<DOS_p[i][1]<<endl;
    cout<<"DOS_p =  "<<"    "<<"count_p = "<<count_DOS_p<<endl;
    */
}
