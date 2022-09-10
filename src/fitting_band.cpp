#include"main.h"

void fitting_band()
{
    FILE *fid1;    	
    char line[1000];
    if (fitting_1==0)
    {
	//cout<<"Reached here"<<endl; getchar();
        cout<<endl<<"Analytical fitting of conduction band"<<endl;
        int *c = analytical_fitting(cond_band,count1,1);    // coloumn - a11[0]+1   row - a11[1]+1    of coefficients
        a11[0] = c[0];
        a11[1] = c[1];

        //cout<<"a11[0] = "<<a11[0]<<endl;
        //cout<<"a11[1] = "<<a11[1]<<endl;
        //getchar();
        // a11[0] --contain mamimum degree   coloumn - a11[0]+1
        // a11[1] = contain length_dd        row - a11[1]+1

//--------------------------- save coefficient and kindex -----------------------------
	/*
        fid1 = fopen("coefficients_conduction_band.txt","w");
        for (int i = 0; i <a11[1]+1 ; i++)
        {
            for (int j = 0; j<a11[0]+1 ; j++)
                fprintf(fid1,"%e    ", coefficients_cond[i][j]);
            fprintf(fid1,"\n");
        }
	fclose(fid1);
	
        fid1 = fopen("kindex_conduction_band.txt","w");
        for (int i = 0; i <a11[1] ; i++)
        {
            fprintf(fid1,"%e \n", kindex_cond[i]);

        }
	fclose(fid1);
	*/
//--------------------------- ---------------------------- -----------------------------
    }
    else    // this part is for debugging code only
    {
        fid1 = fopen("coefficients_conduction_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"coefficients_conduction_band.txt `f\ECle is not present exit from program"<<endl;
            exit(EXIT_FAILURE);        
        }
        //a11[2] and b11[2] contains size of coefficient of conduction and valence band
        a11[0] = 6;  // a11[0]+1 coloumn of conduction band
        a11[1] = 4;  // a11[1]+1 row of conduction band   and a11[1] number of elements in k_index

        for (int i = 0; i < 5; i++)
        {
            fgets(line, 1000, fid1);

                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &coefficients_cond[i][0], &coefficients_cond[i][1], &coefficients_cond[i][2],
                       &coefficients_cond[i][3],&coefficients_cond[i][4],&coefficients_cond[i][5], &coefficients_cond[i][6]);
        }
	fclose(fid1);

        fid1 = fopen("kindex_conduction_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"kindex_conduction_band.txt `f\ECle is not present exit from program"<<endl;
            exit(EXIT_FAILURE);        
        }
    
	for (int i = 0; i < 4; i++)
	{
	    fgets(line, 1000, fid1);

		sscanf(line, "%lf", &kindex_cond[i]);
	}
	fclose(fid1);
    }

    if (fitting_1==0)
    {
        cout<<endl<<"Analytical fitting of valence band"<<endl;
        int *c = analytical_fitting(val_band,count2,2);
        cout<<endl;
        b11[0] = c[0];
        b11[1] = c[1];

        //cout<<"b11[0] = "<<b11[0]<<endl;
        //cout<<"b11[1] = "<<b11[1]<<endl;
        //getchar();
        // b11[0] --contain mamimum degree   coloumn - b11[0]+1
        // b11[1] = contain length_dd         row - b11[1]+1

//--------------------------- save coefficient and kindex -----------------------------
	/*
        fid1 = fopen("coefficients_valence_band.txt","w");
        for (int i = 0; i <b11[1]+1 ; i++)
        {
            for (int j = 0; j<b11[0]+1 ; j++)
                fprintf(fid1,"%e    ", coefficients_val[i][j]);
            fprintf(fid1,"\n");
        }
	fclose(fid1);

        fid1 = fopen("kindex_valence_band.txt","w");
        for (int i = 0; i <b11[1] ; i++)
        {
            fprintf(fid1,"%e \n", kindex_val[i]);
        }
	fclose(fid1);
	*/
//----------------------------------------------------------- -----------------------------
    }
    else                 // this part is for debugging code only
    {
        b11[0] = 6;  // b11[0]+1 coloumn of conduction band
        b11[1] = 4;  // b11[1]+1 row of conduction band  and b11[1] number of elements in k_index
        fid1 = fopen("coefficients_valence_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"coefficients_valence_band.txt `f\ECle is not present exit from program"<<endl;
            exit(EXIT_FAILURE);        
        }
        
        for (int i = 0; i < 5; i++)
        {
            fgets(line, 1000, fid1);

                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &coefficients_val[i][0], &coefficients_val[i][1], &coefficients_val[i][2],
                       &coefficients_val[i][3],&coefficients_val[i][4],&coefficients_val[i][5], &coefficients_val[i][6]);
        }
	fclose(fid1);

        fid1 = fopen("kindex_valence_band.txt", "r");
        if (fid1==NULL)
        {
            cout<<"kindex_valence__band.txt `f\ECle is not present exit from program"<<endl;
            exit(EXIT_FAILURE);        
        }

        for (int i = 0; i < 4; i++)
        {
            fgets(line, 1000, fid1);
            sscanf(line, "%lf", &kindex_val[i]);
        }
	fclose(fid1);
    }
//-------------------------------------------Save coefficient result in txt file -----------------------

//------------------------------- print coefficient -and kindex --------------------------
    cout<<"coefficient of conduction band"<<endl;

    for (int i = 0; i < a11[1]+1; i++)
	{
		for (int j = 0; j < a11[0]+1; j++)
			cout<<coefficients_cond[i][j]<<"     " ;
        cout<<endl;
	}
	cout<<endl;
    cout<<"kindex of conduction band"<<endl;
	for (int i = 0; i < a11[1]; i++)
    {
			cout<<kindex_cond[i]<<endl;
    }
    cout<<endl;


    cout<<"coefficient of valence band"<<endl;
    for (int i = 0; i < b11[1]+1; i++)
	{
		for (int j = 0; j < b11[0]+1; j++)
			cout<<coefficients_val[i][j]<<"     ";
        cout<<endl;
    }
    cout<<endl;

    cout<<"kindex of valence band"<<endl;
	for (int i = 0; i < b11[1]; i++)
        cout<<kindex_val[i]<<endl;
    cout<<endl;

}



