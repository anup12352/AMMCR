#include<iostream>
using namespace std;

#include<fstream>
#include <math.h>

#include<stdio.h>

int main()
{
    double center_kpoint[4] = {0,0,0,1};
    int a;

    cout<<"Enter center point for generating k-point file"<<endl;
    cin>>center_kpoint[0]>>center_kpoint[1]>>center_kpoint[2];

    cout<<"Select no. of k_points"<<endl<<"8531"<<endl<<"6585"<<endl<<"5047"<<endl<<"3869"<<endl<<"3267"<<endl<<"833"
	<<endl<<"391"<<endl;
    cin>>a;
		
	double factor[3]; 
	if(a ==8531){
    		factor[0] = 9;
		factor[1] = 3; 
		factor[2] = 5;}
	else if(a==6585){
    		factor[0] = 8;
		factor[1] = 3; 
		factor[2] = 5;}
	else if(a==5047){
    		factor[0] = 7;
		factor[1] = 3; 
		factor[2] = 5;}
	else if(a==3869){
    		factor[0] = 6;
		factor[1] = 3; 
		factor[2] = 5;}
	else if(a==3267){
    		factor[0] = 6;
		factor[1] = 3; 
		factor[2] = 4;}
			
	

    FILE *fid1;
	char line[1000];

    double kx,ky,kz;
    double k[10000][4] = {0};
	
	
    k[0][0] = center_kpoint[0];
    k[0][1] = center_kpoint[1];
    k[0][2] = center_kpoint[2];
    k[0][3] = center_kpoint[3];
    int j=1;
    if (a==8531 || a==6585 || a==5047 || a==3869 || a==3267)
    {
	double step[3]={0.001, 0.01, 0.04};

	for (int i=0;i<3;i++)
	{
	       	kx=center_kpoint[0]-step[i]*factor[i];

		while (kx<=center_kpoint[0]+step[i]*factor[i]+0.00001)   // 0.00001 is added because in place of zero it taking 1.38778e-17 in  									// c++
		{
		    ky=center_kpoint[1]-step[i]*factor[i];
		    while(ky<=center_kpoint[1]+step[i]*factor[i]+0.00001)
		    {
			kz=center_kpoint[2]-step[i]*factor[i];
			while(kz<=center_kpoint[2]+step[i]*factor[i]+0.00001)
			{
			    if (pow((kx-center_kpoint[0]),2)+pow((ky-center_kpoint[1]),2)+pow((kz-center_kpoint[2]),2)>0.0000000000001)
			    {
				k[j][0] = kx;
				k[j][1] = ky;
				k[j][2] = kz;
				k[j][3] = 1.000000;
				j = j+1;
			    }
			    kz = kz+step[i];
			}
			ky = ky +step[i];
		    }
		    kx = kx + step[i];
		}
	}
    }
	else
	{

 	    double step[40]={0}; int length;
	    if (a==391)
	    {
		double step1[15]={0.001, 0.002, 0.003, 0.004, 0.008, 0.012, 0.016, 0.024, 0.032, 0.048, 0.064, 0.096, 0.128, 0.192, 0.256};
		length = 15;
	  	for (int i=0;i<length;i++)
		    step[i] = step1[i];
		
	    }
	    else if (a==833)
	    {
	    	double step1[32]= {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.012, 0.016, 0.018, 0.020, 0.022, 0.024, 					   0.026, 0.028, 0.030, 0.032, 0.036, 0.040, 0.048, 0.052, 0.064, 0.080, 0.096, 0.112, 0.128, 0.150, 0.192, 					   0.216, 0.256};
		length = 32;
	  	for (int i=0;i<length;i++)
		    step[i] = step1[i];
	    }
	    else
	    {
		cout<<"Exit from program. Run program again and enter right number of kpoints from above given choices"<<endl;
		return 0;
	    }	    
				
	    for (int i=0;i<length;i++)
	    {

		kx=center_kpoint[0]-step[i];

		while (kx<=center_kpoint[0]+step[i]+0.00001)   // 0.00001 is added because in place of zero it taking 1.38778e-17 in c++
		{
		    ky=center_kpoint[1]-step[i];
		    while(ky<=center_kpoint[1]+step[i]+0.00001)
		    {
		        kz=center_kpoint[2]-step[i];

		        while(kz<=center_kpoint[2]+step[i]+0.00001)
		        {
		            if (pow((kx-center_kpoint[0]),2)+pow((ky-center_kpoint[1]),2)+pow((kz-center_kpoint[2]),2)>0.0000000000001)
		            {
		                k[j][0] = kx;
		                k[j][1] = ky;
		                k[j][2] = kz;
		                k[j][3] = 1.000000;
		                j = j+1;
		            }
		            kz = kz+step[i];
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
		fprintf(fid1,"%lf    %lf    %lf    %lf  \n",k[i][0],k[i][1],k[i][2],k[i][3]);
		//cout<<k[i][0]<<"   "<<k[i][1]<<"   "<<k[i][2]<<"   "<<k[i][3]<<"   "<<endl;
	    }
	    //cout<<j<<endl;
	    fclose(fid1);
	    return 0;

}
