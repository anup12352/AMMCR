#include<iostream>
using namespace std;

#include<fstream>
#include <math.h>

#include<stdio.h>

int main()
{
	

    FILE *fid1;

    double k[5000][3] = {0};
    int steps = 1000;
    double delta = 0.001;
	
//---------------------------------------------------------------------------------	
    // gamma to X path		
    k[0][0] = 0;
    k[0][1] = 0;
    k[0][2] = 0;

    
    // gamma to X path
    for(int i=1;i<=steps;i++)
    {
	k[i][0] =  k[i-1][0] + delta;
	k[i][1] =  0;
	k[i][2] =  0;
    }
//---------------------------------------------------------------------------------
    
//---------------------------------------------------------------------------------    
    // For gamma to L path
    k[steps+1][0] = 0;
    k[steps+1][1] = 0;
    k[steps+1][2] = 0;

    // gamma to L path
    for(int i=steps+2;i<2*(steps+1);i++)
    {
	k[i][0] =  k[i-1][0] + delta/2;
	k[i][1] =  k[i-1][1] + delta/2;
	k[i][2] =  k[i-1][2] + delta/2;
    }
//--------------------------------------------------------------------------------- 

    	
	cout<<"Total number of kpoints = "<<(steps+1)*2<<endl;

	fid1 = fopen("kpoints_file","w");
	fprintf(fid1,"Total k points  =  \n 2002  \n");
	fprintf(fid1,"Reciprocal Lattice \n");
	for(int i = 0; i<2*(steps+1);i++)
	{
		fprintf(fid1,"%lf    %lf    %lf \n",k[i][0],k[i][1],k[i][2]);
		//cout<<k[i][0]<<"   "<<k[i][1]<<"   "<<k[i][2]<<"   "<<k[i][3]<<"   "<<endl;
	}
	
	
	fclose(fid1);
	return 0;

}
