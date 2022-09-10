#include"main.h"

void save_perturbation()
{
		FILE *fid1,*fid2;
//----------------------------------save perturbation g(E)--------------------------------------------------------
	    if (ispin == 1 )			    
		    fid1 = fopen("g.dat","w");

	    if (ispin == 2 && kk == 0)			    
		    fid1 = fopen("g_up_spin.dat","w");

	    if (ispin == 2 && kk == 1)			    
		    fid1 = fopen("g_down_spin.dat","w");

       	    fprintf(fid1,"# S.No.    Energy (eV)       g(E)\n");

		if(type=="n")
		{
		    for (int i = 0; i < points; i++)
			{
				//cout<<"i+1 = "<<i+1<<"    g[i] = "<<g[i]<<endl;			
				fprintf(fid1,"  %d        %e      %e  \n", i+1, energy_n[i], g[i]);
				//getchar();
			}
			fclose(fid1);
		}		    
		else
		{
		    for (int i = 0; i < points; i++)
			{
				//cout<<"i+1 = "<<i+1<<"    g[i] = "<<g[i]<<endl;			
				fprintf(fid1,"  %d        %e      %e  \n", i+1, energy_p[i], g[i]);
				//getchar();
			}
			fclose(fid1);
		}		    


	//cout<<"g saved"<<endl;
	//getchar();    

//------------------------------------Save gH(E) and hH(E)-----------------------------------------------------
	if(Bfield!=0)
	{

	   if (ispin == 1 )
	   {			    
		    fid1 = fopen("gH.dat","w");
		    fid2 = fopen("hH.dat","w");
	   }	
	   if (ispin == 2 && kk == 0)			    
	   {
		    fid1 = fopen("gH_up_spin.dat","w");
		    fid2 = fopen("hH_up_spin.dat","w");
	   }	
	    if (ispin == 2 && kk == 1)			    
	     {
	     	    fid1 = fopen("gH_down_spin.dat","w");
		    fid2 = fopen("hH_down_spin.dat","w");
	     }	
       	    fprintf(fid1,"# S.No.    Energy (eV)       gH(E)\n");
       	    fprintf(fid2,"# S.No.    Energy (eV)       hH(E)\n");


	    for (int i = 0; i < points; i++)
	        {
			//cout<<"i+1 = "<<i+1<<"    gH[i] = "<<gH[i]<<endl;			
			//cout<<"i+1 = "<<i+1<<"    hH[i] = "<<hH[i]<<endl;			
			fprintf(fid1,"  %d        %e      %e  \n", i+1, energy_n[i], gH[i]);
			fprintf(fid2,"  %d        %e      %e  \n", i+1, energy_n[i], hH[i]);
			//getchar();
		}
		fclose(fid1);
	    	fclose(fid2);
	}
	
}
