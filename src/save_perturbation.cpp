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
				fprintf(fid1,"  %d        %e    ", i+1, energy_n[i]);

				for(int j=0;j<iterations;++j)
					fprintf(fid1,"  %e     ", result_g[i][j+1]);

				fprintf(fid1,"\n");

				//getchar();
			}
			fclose(fid1);
		}		    
		else
		{
		    for (int i = 0; i < points; i++)
			{
				//cout<<"i+1 = "<<i+1<<"    g[i] = "<<g[i]<<endl;			
				fprintf(fid1,"  %d        %e      ", i+1, energy_p[i]);


				for(int j=0;j<iterations;++j)
					fprintf(fid1,"  %e     ", result_g[i][j+1]);

				fprintf(fid1," \n ");
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
			fprintf(fid1,"  %d        %e     ", i+1, energy_n[i]);
			fprintf(fid2,"  %d        %e     ", i+1, energy_n[i]);


			for(int j=0;j<iterations;++j)
				fprintf(fid1,"  %e     ", result_gH[i][j+1]);

			fprintf(fid1," \n ");

			for(int j=0;j<iterations;++j)
				fprintf(fid2,"  %e     ", result_hH[i][j+1]);

			fprintf(fid2," \n ");

			//getchar();
		}
		fclose(fid1);
	    fclose(fid2);
	}
	
}
