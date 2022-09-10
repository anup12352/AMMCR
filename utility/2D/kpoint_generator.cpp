#include<iostream>
using namespace std;

#include<fstream>
#include <math.h>

#include<stdio.h>

int main()
{
    double center_kpoint[4] = {0,0,0,1};
    int a;

    cout<<"Enter center point for kx and ky for generating k-point file"<<endl;
    cin>>center_kpoint[0]>>center_kpoint[1];

    cout<<"Enter a fixed value of kz for generating k-point file"<<endl;
    cin>>center_kpoint[2];

    cout<<"Select no. of k_points"<<endl<<"1161"<<endl<<"1025"<<endl<<"897"<<endl;
    cin>>a;
		
	double factor[3]; 
	if(a ==1161){
    		factor[0] = 9;
		factor[1] = 9; 
		factor[2] = 10;}
	else if(a==1025){
    		factor[0] = 9;
		factor[1] = 7; 
		factor[2] = 10;}
	else if(a==897){
    		factor[0] = 8;
		factor[1] = 6; 
		factor[2] = 10;}
			
	

    FILE *fid1;
	char line[1000];

    double kx,ky,kz;
    double k[10000][4] = {0};
	
    kz = center_kpoint[2];	
	
    k[0][0] = center_kpoint[0];
    k[0][1] = center_kpoint[1];
    k[0][2] = center_kpoint[2];
    k[0][3] = center_kpoint[3];
    
    int j=1;
    
    if (a==1161 || a==1025 || a==897)
    {
	double step[3]={0.001, 0.01, 0.02};

	for (int i=0;i<3;i++)
	{
	       	kx=center_kpoint[0]-step[i]*factor[i];

		while (kx<=center_kpoint[0]+step[i]*factor[i]+0.00001)   // 0.00001 is added because in place of zero it taking 1.38778e-17 in  									// c++
		{
		    ky=center_kpoint[1]-step[i]*factor[i];
		    while(ky<=center_kpoint[1]+step[i]*factor[i]+0.00001)
		    {
			    if (pow((kx-center_kpoint[0]),2)+pow((ky-center_kpoint[1]),2)>0.0000000000001)
			    {
				k[j][0] = kx;
				k[j][1] = ky;
				k[j][2] = kz;
				k[j][3] = 1.000000;
				j = j+1;
			    }
			ky = ky +step[i];
		    }
		    kx = kx + step[i];
		}
	}
    }

	    cout<<"Total number of kpoints = "<<j<<endl;

	    fid1 = fopen("kpoints_file","w");
	    fprintf(fid1,"Total k points \n");
	    fprintf(fid1,"%d \n",j);
	    fprintf(fid1,"Reciprocal Lattice \n");
	    for(int i = 0; i<j;i++)
	    {
		fprintf(fid1,"%lf    %lf    %lf   %lf \n",k[i][0],k[i][1],k[i][2], k[i][3]);
		//cout<<k[i][0]<<"   "<<k[i][1]<<"   "<<k[i][2]<<"   "<<k[i][3]<<"   "<<endl;
	    }
	    //cout<<j<<endl;
	    fclose(fid1);
	    return 0;

}
