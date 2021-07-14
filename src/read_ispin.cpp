#include"main.h"

void read_ispin()
{
	char line[1000];
	int a[10],i,j ;
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

    //cout<<"ispin = "<<ispin<<endl;
    //cout<<"Next Line"<<endl;
    //getchar();	
}



