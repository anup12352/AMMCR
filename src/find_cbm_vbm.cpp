
#include"main.h"

double ecbm,evbm,kcbm1[3],kvbm1[3];
int NKPTS,NBVAL,NBTOT;


void find_cbm_vbm(int aa,int spin_orbit_coupling)
{
	char line[1000];
	int a[10],i,j ;
	FILE *fid;
	if (aa==1)
    {
        fid = fopen("EIGENVAL_n","r");
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
    else
    {
        fid = fopen("EIGENVAL_p","r");
        if (fid==NULL)
        {
            fid = fopen("EIGENVAL","r");
            if (fid==NULL)
            {   cout<<"EIGENVAL is not present. Exit from program";
                exit(EXIT_FAILURE);
            }
        }
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
    
    double energies[NKPTS][NBTOT], kpoints[NKPTS][4], temp[3];

    for (i=0; i<NKPTS; i++)
    {
        fgets(line,1000,fid);   //passing one empty line
        fgets(line,1000,fid);    // reading kpoint line
        sscanf(line, "%lf %lf %lf %lf ", &kpoints[i][0], &kpoints[i][1], &kpoints[i][2], &kpoints[i][3] );

        //cout<<"line = "<<line<<endl;
        //cout<<"kpoints = "<<kpoints[i][0]<<"   "<<kpoints[i][1]<<"   "<<kpoints[i][2]<<"   "<<kpoints[i][3]<<endl;

        for (j = 0; j<NBTOT; j++)
        {
            //cout<<"j = "<<j<<endl;
            fgets(line,1000,fid);
            //cout<<"line = "<<line<<endl;
            sscanf(line, "%lf %lf %lf ", &temp[0], &temp[1], &temp[2]);            
            energies[i][j] = temp[1];     
                        
		
	    if (ispin == 2 && kk==1)    // reading for spin down
	            energies[i][j] = temp[2];
	    /*
	     if(int(temp[0]) == NBVAL-1 )
	     	cout<<"Valence Energy = "<<energies[i][j]<<endl;		

	     if(int(temp[0]) == NBVAL )
	     	cout<<"Conduction Energy = "<<energies[i][j]<<endl;		
	    //getchar();
	    
            //cout<<"energy = "<<energies[i][j]<<endl;
            //getchar();
            */
        }
        //getchar();
    }
    
    evbm = -1000;
    ecbm = 1000;

    for (i=0; i<NKPTS; i++)
    {
            if(energies[i][NBVAL-1] > evbm)
            {
                evbm = energies[i][NBVAL-1];
                kvbm1[0] = kpoints[i][0];
                kvbm1[1] = kpoints[i][1];
                kvbm1[2] = kpoints[i][2];
            }
    }
    for (i = 0; i< NKPTS; i++)
    {
            if(energies[i][NBVAL] < ecbm )
            {
                ecbm = energies[i][NBVAL];
                kcbm1[0] = kpoints[i][0];
                kcbm1[1] = kpoints[i][1];
                kcbm1[2] = kpoints[i][2];
            }
    }

}



