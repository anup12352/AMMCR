#include"main.h"

void nu_ii(double epsilon_s)
{
	double k_dum, B_ii,D_ii, xx, yy, zz, k, v;
	
	for (int counter = 0;counter < points;counter++)
	{
		k_dum = k_grid[counter];
		k = k_dum;
		v = v_n[counter];
	    //cout<<"k_dum = "<<k_dum<<endl;
	    //cout<<"beta_constant =  "<<beta_constant<<endl;
	    //cout<<"c_n[counter] =  "<<c_n[counter]<<endl;
	    //cout<<"counter =  "<<counter+1<<endl;

	    B_ii = (4*(k_dum*k_dum)/(beta_constant*beta_constant))/(1+4*k_dum*k_dum/(beta_constant*beta_constant))
	    +8*(beta_constant*beta_constant+2*k_dum*k_dum)/(beta_constant*beta_constant+4*k_dum*k_dum)*(pow(c_n[counter],2))
	    +(3*pow(beta_constant,4)+
	      6*pow(beta_constant,2)*k_dum*k_dum-8*k_dum*k_dum*k_dum*k_dum)/((beta_constant*beta_constant+4*k_dum*k_dum)*k_dum*k_dum)*pow(c_n[counter],4);
	    // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)


	    D_ii = 1+(2*pow(beta_constant,2)*(pow(c_n[counter],2)/pow(k_dum,2))+
	              (3*pow(beta_constant,4)*(pow(c_n[counter],4))/(4*pow(k_dum,4))));
	              // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)

		xx = (pow(e,4)*abs(N_ii));
		yy = (8*pi*v*epsilon_s*epsilon_s*epsilon_0*epsilon_0*h_bar*h_bar*k*k);
		zz = (D_ii*log(1+4*k*k/(beta_constant*beta_constant))-B_ii);

		//cout<<"xx = "<<xx<<endl;
		//cout<<"yy = "<<yy<<endl;
		//cout<<"zz = "<<zz<<endl;
		nu_ionizedimpurity[counter] = abs( xx/yy*zz*3.89564386e27 );
		// unit 1/second	
		//cout<<"B_ii = "<<B_ii<<endl;
		//cout<<"D_ii = "<<D_ii<<endl;
		//cout<<"N_ii = "<<N_ii<<endl;
		//cout<<"nu_ionizedimpurity[counter] = "<<nu_ionizedimpurity[counter]<<endl;

	}

	/*
	fid1 = fopen("nu_ionizedimpurity.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, nu_ionizedimpurity[i]);
	fclose(fid1);
	*/
	
	/*
	fid1 = fopen("nu_ionizedimpurity.txt","r");
	for (int i = 0; i < points; i++)
	{
	fgets(line, 1000, fid1);
	sscanf(line, "%lf", &nu_ionizedimpurity[i]);
	}
	fclose(fid1);
	
	*/
}
