#include"main.h"

int read_ispin()
{
	char line[1000];
	int a[10],i,j,jj;
	FILE *fid;

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

    fgets(line, 1000, fid);
    sscanf(line, "%d %d %d %d", &a[0], &a[1], &a[2], &a[3]);
    ispin=a[3]; // %4th number is ispin (ispin =1: non-magnetic calculation, ispin=2: magnetic)

    cout<<"ispin = "<<ispin<<endl;
    //cout<<"Next "<<endl;
    if (ispin == 1)    	
    {
	cout<<"Material is non magnetic"<<endl;
	jj = 0;
    }
    else
    {
    	cout<<"Material is magnetic"<<endl;
	jj = 1;
    }

    //cout<<"ispin = "<<ispin<<endl;
    //cout<<"Next Line"<<endl;
    //getchar();
    return jj;	
}



