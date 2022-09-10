#include"main.h"

int FindMinInd(double arr[],int length);
double * linspace(double a, double b, int numbers);
double rsquare(double data[], double fitted[], int length);
void polyval(double p[], double x[], int nn, int mm);

//double coefficients[5][7] = {0};
//double kindex[4]={1,1,1,1};

double coefficients[5][7]={0},kindex[4]={1,1,1,1},coefficients_cond[5][7]={0},kindex_cond[4]={1,1,1,1},
coefficients_val[5][7]={0},kindex_val[4]={1,1,1,1};

int max(int a, int b)
{
	if (a >= b)
		return a;
	else
		return b;
}

int min(int a, int b)
{
	if (a <= b)
		return a;
	else
		return b;
}

double max_band(double band[][2],int county)
{
    double max1 = -1000;
	for (int i = 0; i < county; i++)
	{
		if (band[i][1] > max1)
			max1 = band[i][1];
	}
	return max1;
}

int* analytical_fitting(double band[2000][2],int county, int aa)   // if aa = 1 then its conduction band   if aa=2 then its for valence band
{
	//double suggested_division_bs[10] = { 0.0000766, 0.034909090535536, 0.213199787119379, 0.849825392674487 };
    /*
    for(int i=0;i<county;i++){
        cout<<"band = "<<band[i][0]<<"    "<<band[i][1]<<endl;
        if(i==50||i==100||i==150||i==200||i==250||i==300)
            getchar();
        }
        */
	double max_range, growth_base=1.0, division_default[10], division_start[10],
		coefficients_selected[10][10], kindex_selected[10],
		division[10], division_selected[10], r2_penalty = 10 ;

    double short_range_band_num[2000][2];

	int merge_curves, accuracy_factor = 1,maximum_degree;
	int degree[10], degree_selected[10];
	int start_index, end_index, counter, i_start = -1, start_merge_curves;
	double growth, temp[10000][5];
	double r2[10], r2_selected[10], tolerance = 5, test_negative = -10;
	double step, maxitr = 4, average_discontinuity_selected, average_velocity_discontinuity_selected;
	double total_fitting_time, coeffs[1000];
	double array_data[10000];
	double polycoeffs_short_range[1000], src1[2000], src2[2000], y[1000];
	int length1, length2, number_of_curves_fitted, count_iter = 0, count_i = 0, count_r2 = 0;
	double average_velocity_discontinuity, average_discontinuity;
	double *powers, energy_left, energy_right, derivative_left, derivative_right;
	//std::vector<double> v1, v2;

	clock_t start, end1;
	start = clock();

	max_range = max_band(band,county);
	//printf("\n max_range = %lf ", max_range);
    //cout<<endl;

	merge_curves = max(40, floor( (county) / 4) - 1);
	//printf("merge curves = %d \n", merge_curves);
    //cout<<endl;



	if (degree1 == 0)
	{		
		maximum_degree = 6;
		//cout<<"default degree"<<endl;
	}
	else
		maximum_degree = degree1;


	if ((county) < 41)
		maximum_degree = max(floor((county)/8), 2);

	growth_base = 1.0;

	int length_dd = 4;

	if (fraction[0]==0 && fraction[1]==0 && fraction[2]==0 && fraction[3]==0)
	{
		//cout<<"default division_default is running. Press Enter to continue"<<endl;
		//getchar();
		division_default[0] = 0.1*max_range;
		division_default[1] = 0.2*max_range;
		division_default[2] = 0.5*max_range;
		division_default[3] = 1.2*max_range;
	}
	else
	{
		if (length_fraction==4)		
		{
			division_default[0] = fraction[0]*max_range;
			division_default[1] = fraction[1]*max_range;
			division_default[2] = fraction[2]*max_range;
			division_default[3] = fraction[3]*max_range;	
		}
		else if (length_fraction == 3)
		{

			division_default[0] = fraction[0]*max_range;
			division_default[1] = fraction[1]*max_range;
			division_default[2] = fraction[2]*max_range;
			division_default[3] = 0.0;	
			length_dd = 3;
						
		}
	}



    //cout<<"division_default = "<<division_default[0]<<"    "<<division_default[1]<<"    "
    //<<division_default[2]<<"    "<<division_default[3]<<endl;

	if ((county) < (3 * maximum_degree))
    {
		printf("\n ERROR!!!:more data points are required for band structure fitting");
		end1 = clock();
		total_fitting_time = ((double)(end1 - start) / (CLOCKS_PER_SEC * 60));

		printf("\n Curve fitting time = %lf s \n", total_fitting_time);
		return 0;
	}
	else
	{
		if ((county ) <  length_dd * (maximum_degree + 4))
        {
			division_default[0] = 0.1*max_range;
			division_default[1] = 0.5*max_range;
			division_default[2] = 1.2*max_range;
			division_default[3] = 0.0;
			length_dd = 3;
		}
	}

    //cout<<"division_default= "<<division_default[0]<<"    "
    //<<division_default[1]<<"   "<<division_default[2]<<"   "<<division_default[3];
    //cout<<"length_dd = "<<length_dd<<endl;

	merge_curves = max(0, min(merge_curves, floor((county) / length_dd) - 1));
	//printf("\n merge curves = %d", merge_curves); getchar();

	for (int i = 0; i < length_dd; i++)
    {
		division_start[i] = division_default[i];
		division[i] = division_start[i];
		division_selected[i] = division_start[i];
	}

	for (int i = 0; i < 5; i++)
    {
		degree[i] = maximum_degree;
		degree_selected[i] = maximum_degree;
	}

	int d1_initial = maximum_degree;
	int d2_initial = maximum_degree;
	int d3_initial = maximum_degree;

	for (int i = 0; i < 5; i++)
    {
		r2[i] = 0.5;
		r2_selected[i] = 0.5;
	}
	step = 0.02*max_range;
    //cout<<"step0 = "<<step<<endl;   getchar();

	for (int iteration = 0; iteration <= maxitr; iteration++)
    {
		step = step/pow(2,iteration);
        //cout<<"step1 = "<<step<<endl;  getchar();
		while (test_negative < 0)
        {
			for (int i = 0; i < length_dd; i++)
				division_start[i] = division_selected[i];
			average_discontinuity_selected = 1e30;
			average_velocity_discontinuity_selected = 1e30;
			int div1,div2,div3;
			for (int d1 = d1_initial; d1 <= maximum_degree; d1++)   // there is no need of this loop
			{
				for (int d2 = d2_initial; d2 <= maximum_degree; d2++)   // there is no need of this loop
				{
					for (int d3 = d3_initial; d3 <= maximum_degree; d3++)   // there is no need of this loop
					{
						for (div1 = -tolerance; div1 <= tolerance * 2 + floor(sqrt(1 / accuracy_factor)); div1++)
						{
                        //cout<<"running "<<endl;

							for (div2 = -tolerance; div2 <= tolerance * 2; div2++)
							{
								for (div3 = -tolerance * 2; div3 <= tolerance; div3++)
								{

                                    //cout<<endl<<"div3   =    "<<div3<<endl;
                                    //cout<<endl;
									degree[0] = d1;
									degree[1] = d2;
									degree[2] = d3;
									for (int i = 0; i < length_dd+1; i++)
									{
										for (int j = 0; j < maximum_degree+1; j++)
										{
											coefficients[i][j] = 0;
										}
									}

									//cout <<"wwwwwwwwwwwwwwww"<<endl;

									for (int i = 0; i < length_dd; i++) {
										kindex[i] = 1;
									}
									for (int i = 0; i < length_dd; i++)
                                    {
										division[i] = division_start[i];
                                        //cout<<"i =  "<<i<<" division[i] = "<<division[i]<<endl;
									}

                                    /*
                                    cout<<"step = "<<step<<endl;
                                    cout<<"accuracy factor = "<<accuracy_factor<<endl;
                                    cout<<"div2 =  "<<div2<<endl;
                                    cout<<"div1 =  "<<div1<<endl;
                                    */

									division[0] = division_start[0] + div1*step*accuracy_factor; //accuracy factor = 1 no meaning
									division[1] = division_start[1] + div2*step;
									division[2] = division_start[2] + div3*step;


                                    /*
                                    cout<<"division[0] = "<<division[0]<<endl;
                                    cout<<"division[1] = "<<division[1]<<endl;
                                    cout<<"division[2] = "<<division[2]<<endl;
                                    */

									for (int i = 0; i < county; i++)
										array_data[i] = abs(band[i][1] - division[0]);
									start_index = FindMinInd(array_data, county);

                                    //cout<<"start_index = "<<start_index<<endl;
                                    //getchar();

									counter = 0;
                                    //cout<<"Before first while loop line no. 205"<<endl;
									while (start_index + 1 < (maximum_degree + 1 + 3))
                                    {
										growth = pow(growth_base, counter);   // it is always 1
										division[0] = division[0] + step*accuracy_factor*growth;
										for (int i = 0; i < county; i++)
											array_data[i] = abs(band[i][1] - division[0]);
										start_index = FindMinInd(array_data, county);
										counter++;

										if (counter > 50)
											break;
									}
									//cout <<"wwwwwwwwwwwwwwww"<<endl;
                                    /*
                                    cout<<"After first while loop line no. 220"<<endl;
                                    cout<<"start_index = "<<start_index<<endl;
                                    cout<<"division[0] = "<<division[0]<<endl;
                                    //cout<<"press key to continue "<<endl;  getchar();
                                    */
									for (int i = 1; i <= start_index; i++)
									{
										for (int j = 0; j < 2; j++)
											temp[i - 1][j] = band[i][j];
									}

                                    /*
                                    cout<<"start_index = "<<start_index<<endl;
                                    cout<<endl<<"Tempppppppppppp"<<endl;
                                    for (int i = 0; i < start_index; i++)
                                        {
                                            cout<<"i+1 = "<<i+1<<"  temp[i][0] =  "<<temp[i][0]<<"  temp[i][1] =  "<<temp[i][1]<<endl;
                                            if (i==100 || i==200 || i==300 || i==400 || i==500)
                                            getchar();
                                        }
                                    cout<<"Press any key to continue"<<endl; getchar();
                                    */
									for (int i = 0; i < start_index; i++)
									{
										short_range_band_num[i][0] = -temp[i][0];
										short_range_band_num[i][1] = temp[i][1];
										//cout<<"i+1 = "<<i+1<<"   srbn[i][0] = "<<short_range_band_num[i][0]
										//<<"srbn[i][1] = "<<short_range_band_num[i][1]<<endl;
                                    }

									short_range_band_num[start_index][0] = band[0][0];
									short_range_band_num[start_index][1] = band[0][1];
									for (int i = start_index + 1; i <= 2 * start_index; i++)
									{
										short_range_band_num[i][0] = temp[i - start_index - 1][0];
										short_range_band_num[i][1] = temp[i - start_index - 1][1];
									}

									for (int i = 0; i <= 2 * start_index; i++)
                                    {
										src1[i] = short_range_band_num[i][0];
										src2[i] = short_range_band_num[i][1];
										//printf("\n %lf %lf %d", src1[i], src2[i], degree[0]);
									}
                                    int ss = 2 * start_index + 1;
                                    //cout<<"ss = "<<ss<<endl;
                                    //cout<<"degree[0] = "<<degree[0]<<endl;
                                    //getchar();
									double *coeff = polyfit1(src1, src2, degree[0], ss);

                                    /*
                                    cout<<"ss = "<<ss<<endl;
                                    cout<<"degree[0] = "<<degree[0]<<endl;
                                    cout<<"src1 and src2 of polyfit "<<endl;
                                    for(int i=0;i<ss;i++)
                                        { cout<<"i+1 = "<<i+1<<"    src1[i] =  "<<src1[i]<<"    src2[i] = "<<src2[i]<<endl;
                                        if (i==200 || i==300 || i==400 || i==600 || i==800)
                                            getchar();
                                        }
                                    cout<<"Press any key to continue"<<endl; getchar();
                                    */

									for (int i = 0; i < maximum_degree + 1; i++)
										coefficients[0][i] = coeff[i];

                                    //count_i = 1;

                                    //cout<<" src1 = "<<endl;

									for (int i = 0; i < degree[0] + 1; i++)
									{
									    src1[i] = coefficients[0][i];
									    //cout<<src1[i]<< "    ";
									}
                                    //getchar();

                                    //cout<<" src2 = "<<endl;

                                    //cout<<"short_range_band_num"<<endl;
                                    //cout<<"length of short range band num = "<<2 * start_index+1<<endl;
									for (int i = 0; i <= 2 * start_index; i++)
                                    {
										src2[i] = short_range_band_num[i][0];
                                        //cout<<"i+1 = "<<i+1<<"    "<<src2[i]<<endl;
                                        //if (i==100 || i==200 || i==300 || i==400 || i==500 || i==600 || i==700 || i==800 || i==900)
                                        //        getchar();
                                    }
                                    //getchar();

									polyval(src1, src2, degree[0] + 1 , 2 * start_index + 1);   // poly array is changed

									/*
                                    cout<<"poly = "<<endl;
                                    for(int i=0; i < 2 * start_index + 1;i++)
                                    {
                                        cout<<poly[i]<<endl;
                                        if (i==50 || i==100 || i==150 || i==200 || i==250 || i==300)
                                                getchar();
                                    }
                                    cout<<"Press any key to continue"<<endl;  getchar();
									*/

									for (int i = 0; i <= 2 * start_index; i++)
										src1[i] = short_range_band_num[i][1];

									r2[0] = rsquare(src1, poly, 2 * start_index+1);
									kindex[0] = band[start_index][0];

                                    /*
                                    cout <<"Coefficient after first while loop"<<endl;
                                    cout<<"coeff = "<<coeff[0]<<endl<<coeff[1]<<endl<<coeff[2]<<endl
                                    <<coeff[3]<<endl<<coeff[4]<<endl<<coeff[5]<<endl<<coeff[6]<<endl;
                                    cout<<" r2[0] = "<<r2[0]<<endl;
                                    cout<<"kindex[0] = "<<kindex[0];
                                    cout<<"Press any key to continue "<<endl;
                                    getchar();
                                    */
									for (int i =0; i < length_dd-1; i++)
                                    {
										//count_i++;
										for (int j = 0; j <= county; j++)
											array_data[j] = abs(band[j][1] - division[i + 1]);
										end_index = FindMinInd(array_data, county);
										kindex[i + 1] = band[end_index][0];

                                        //cout<<"I m inside for loop of i_start line no. 328"<<endl;
                                        //cout<<"i = "<<i<<endl;
                                        //cout<<"end_index = "<<end_index<<endl;
                                        //cout<<"kindex[i+1] = "<<kindex[i+1]<<endl;
                                        //cout<<"press key to continue "<<endl; getchar();

										counter = -1;
                                        //cout<<"Before second while loop line no. 338"<<endl;
										while ((end_index - start_index) < (maximum_degree + merge_curves + 1))
										{
											growth = pow(growth_base, counter); // it is always 1
											division[i + 1] = division[i + 1] + step*growth;
											for (int j = 0; j < county; j++)
												array_data[j] = abs(band[j][1] - division[i + 1]);
											end_index = FindMinInd(array_data, county);
											kindex[i + 1] = band[end_index][0];

											counter++;
											if (counter > 50)
												break;

										}///end second while
										/*
                                        cout<<"After second while loop line no. 353"<<endl;
                                        cout<<"i    =  "<<i<<endl;
                                        cout<<"end_index = "<<end_index<<endl;
                                        cout<<"kindex[i+1] = "<<kindex[i+1]<<endl;
                                        cout<<"press key to continue "<<endl; getchar();
                                        */
										start_merge_curves = merge_curves;
										counter = 0;
								        //cout<<"Before third while loop line no. 332"<<endl;
                                        //cout<<"start_merge_curves = "<<start_merge_curves<<endl;

										while ((start_index+1) < (maximum_degree + start_merge_curves))
                                        {
											start_merge_curves = floor(start_merge_curves / 2.0);
                                            //cout<<"Inside third while loop line no. 406"<<endl;
                                            //printf("start_merge_curves  = %d", start_merge_curves);

											counter++;
											if (counter > 50)
												break;

										}/// end third while
                                        /*
                                        cout<<"After third while loop line no. 370"<<endl;
                                        cout<<"start_merge_curves =  "<<start_merge_curves<<endl;
                                        cout<<"Press any key to contunie"<<endl;  getchar();
                                        cout<<"before if condiction after third while loop"<<endl;
                                        */
										if (end_index < (county  - (merge_curves + maximum_degree)))
                                        {
                                            //cout<<"Inside if loop line no. 376"<<endl;
											//getchar();
											if ((end_index + 1  + merge_curves - (start_index + 1 - start_merge_curves)) < maximum_degree)
												continue;

											/*
                                            printf("\n i = %d", i);
                                            //cout<<"Press any key to continue"<<endl;getchar();
                                            if (i + 1 == 1) printf("I am after  if condition line no. 431 %d", i+1);
                                            printf("\nif start_index = %d %d %d %d\n", start_index, end_index, start_merge_curves, merge_curves); getchar();
                                            */
											for (int j = (start_index - start_merge_curves); j <= (end_index + merge_curves); j++)
                                            {
												src1[j - start_index + start_merge_curves] = band[j][0];
											}

											for (int j = start_index - start_merge_curves ; j <= end_index + merge_curves; j++)
                                            {
												src2[j - start_index + start_merge_curves] = band[j][1];
											}

											length1 = end_index + merge_curves - start_index + start_merge_curves + 1;
											//length2 = end_index + merge_curves - start_index + start_merge_curves + 1  ;

											//printf("\n if now length = %d %d", length1);getchar();
                                            /*
                                            cout<<"src1[]     src2[]"<<endl;
											for (int j = 0; j < length1; j++)
												printf("\n %lf     %lf", src1[j], src2[j]);
											*/
											//getchar();

											coeff = polyfit1(src1, src2, degree[i + 1], length1);

                                            //cout.precision(15);

                                        	for (int j = 0; j < maximum_degree + 1; j++)
												coefficients[i + 1][j] = coeffs[j];
											//count_i = i + 2 ;
											for (int j = 0; j < degree[i + 1] + 1; j++)
                                            {
												src1[j] = coefficients[i + 1][j];
												//printf("\n src1 %lf ", src1[j]);
											}
											//getchar();

											for (int j = start_index - start_merge_curves; j <= end_index + merge_curves; j++)
												src2[j - start_index + start_merge_curves] = band[j][0];

											length1 = degree[i + 1] + 1 ;
											length2 = end_index + merge_curves - start_index + start_merge_curves + 1;
											//printf(" length1 and length2   %d   %d", length1, length2);
                                            //getchar();

											polyval(src1, src2, length1, length2);

                                    		//polyval always return in poly function and it is of length2(same as of length src2)
											/*printf("y[k] = ");
											for (int k = 0; k < length2; k++) {
												printf(" %lf", y[k]);
											}*/

											for (int j = start_index - start_merge_curves; j <= end_index + merge_curves; j++)
												src1[j - start_index + start_merge_curves] = band[j][1];
											length1 = end_index + merge_curves- start_index + start_merge_curves + 1;

											r2[i + 1] = rsquare(src1, poly, length1);
										}///end if
										else
                                        {
                                            //cout<<"Inside else line no. 476"<<endl;
											if ((end_index - (start_index  - start_merge_curves) ) < maximum_degree)
												continue;
											for (int j = start_index - start_merge_curves; j <= end_index; j++)
                                            {
												src1[j - start_index + start_merge_curves] = band[j][0];
                                                //cout<<"src1 = "<<src1[j - start_index + start_merge_curves]<<endl;

                                                src2[j - start_index + start_merge_curves] = band[j][1];
                                                //cout<<"src2 = "<<src2[j - start_index + start_merge_curves]<<endl;
											}

											length1 = end_index - start_index + start_merge_curves + 1;
											//length2 = end_index - start_index + start_merge_curves + 1;
											//printf("\n length =  %d %d %d", end_index,start_index, start_merge_curves);
											//printf("\n coefficients = ");
											coeff = polyfit1(src1, src2, degree[i + 1], length1);

											for (int j = 0 ; j < maximum_degree + 1; j++)
                                            {
												coefficients[i + 1][j] = coeff[j];
												//printf(" %e", coefficients[i+1][j]);
                                            }
											//count_i = i + 2;

											for (int j = 0; j < degree[i + 1] + 1; j++)
                                            {
                                                src1[j] = coefficients[i + 1][j];
                                                //cout<<"src1 = "<<src1[j]<<"     "<<endl;
                                            }
                                            //cout<<endl; getchar();


                                            //cout<<"start_index   =  "<<start_index<<endl;
                                            //cout<<"start merge curves  = "<<start_merge_curves<<endl;
                                            //cout<<"end_index    =  "<<end_index<<endl;
                                            //cout<<"press key to continue "<<endl; getchar();

											for (int j = start_index - start_merge_curves; j <= end_index; j++)
                                            {
												src2[j - start_index + start_merge_curves] = band[j][0];
                                                //cout<<"src2 = "<<src2[j - start_index + start_merge_curves]<<endl;
                                                //cout<<"band  "<<band[j][0]<<endl;
                                            }
                                            // getchar();
                                            //getchar();

											length1 = degree[i + 1]  + 1;
											length2 = end_index - start_index + start_merge_curves + 1;
											//src2[0] =0.171600053278014;
											//length2 =1;
											//cout<<"length ===="<<endl;
											//printf(" %d      %d", length1, length2);
											//getchar();

                                            //cout<<"Before poly"<<endl;
                                            //length2 = 1;
											//src2[0] = 0.171600053278014;
											polyval(src1, src2, length1, length2);

                                            /*cout<<"poly  =="<<endl;
                                            for(int i=0; i < length2; i++)
                                            {
                                                cout<<poly[i]<<endl;
                                                if (i==50 || i==100 || i==150 || i==200 || i==250 || i==300)
                                                        getchar();
                                            } getchar();
                                            */

											for (int j = start_index - start_merge_curves; j <= end_index; j++)
                                            {
											    src1[j - start_index + start_merge_curves] = band[j][1];
                                                //cout<<"src1 = "<<src1[j - start_index + start_merge_curves]<<endl;
                                            }

											length1 = end_index - start_index + start_merge_curves + 1;

                                            /*
                                            FILE *fid;
                                            fid = fopen("poly.txt", "r");

                                            char line[1000];
                                            for(int i=0; i < length1; i++)
                                            {
                                                fgets(line, 1000, (FILE*)fid);
                                                sscanf(line, "%lf", &poly[i]);
                                                cout<<"poly new =   "<<poly[i]<<endl;
                                            }
                                            getchar();
                                            */

											r2[i + 1] = rsquare(src1, poly, length1);
											//cout<<"r2 check inside else "<<endl;
											//count_r2++;
											//printf("\n %d %lf ", count_r2, r2[i + 1]);

											/*for (int j = 0; j < length1; j++)
												printf(" %lf", src1[j]);*/
												/*printf("\n else start end = %d %d", start_index, end_index);
												getchar();*/

										}  ///end if else
                                        /*
                                        cout<<"After if else condiction after third whhile loop"<<endl;

                                        cout<<"coefficient "<<endl;
                                        for (int j = 0; j <= degree[i + 1]; j++)
                                            printf(" %lf    ", coeff[j]);

                                        cout<<"r2[i + 1] =  "<<r2[i + 1]<<endl;

                                    //count_iter++;
                                        //cout<<"At end of for loop:    start index: "<<endl;
                                        //printf("\n %d  \n", start_index);
                                        //cout<<"press key to continue "<<endl; getchar();
                                        */
                                        start_index = end_index;
									}///end for i loop


                                    //cout<<endl<<"After for loop Press key to continue"<<endl;
                                    //cout<<"Press any keey to continue"<<endl;
                                    //getchar();

									number_of_curves_fitted = length_dd+1;
									//printf("\n  number of curves fitted = %d\n", number_of_curves_fitted); getchar();

									average_discontinuity = 0;
									average_velocity_discontinuity = 0;
									/*
									for (int j = 0; j < number_of_curves_fitted - 1; j++)
                                    {
										printf("\n r2 = %lf %d %lf", r2[j], r2_penalty, max_range);
									}
									*/

									for (int j = 0; j < number_of_curves_fitted-1; j++)
                                    {
										powers = linspace(maximum_degree, 0, maximum_degree + 1);
                                        /*
                                        cout<<"powersssss"<<endl;
                                        getchar();
										//printf("\nr2[j] = %e", r2[j]);
										printf("\n powers= ");
										for(int k = 0; k < maximum_degree + 1; k++)
											cout<<powers[k]<<endl;
										//printf(" size coefficients = %d %d", count_i , maximum_degree+1);
										getchar();
                                        cout<<"aaaaaaaaaaaaa";
                                        */
										energy_left = 0;
										energy_right = 0;
										for (int k = 0; k < maximum_degree + 1; k++)
										{
										    energy_left = energy_left + coefficients[j][k] * pow(kindex[j], powers[k]);
											energy_right = energy_right + coefficients[j + 1][k] * pow(kindex[j], powers[k]);
                                        }
                                        //printf("\n Energy_left = %e", energy_left); getchar();
                                        //printf("\n Energy_right= %e", energy_right); getchar();

                                        //cout<<"reached in between  "<<endl;
                                        //getchar();

										derivative_left = 0;
										derivative_right = 0;
                                        for (int k = 0; k < maximum_degree + 1; k++)
                                        {
                                            derivative_left = derivative_left + coefficients[j][k] * powers[k] * pow(kindex[j], powers[k] - 1);
                                            derivative_right = derivative_right + coefficients[j + 1][k] * powers[k] * pow(kindex[j], powers[k] - 1);
                                        }
                                        //printf("\n derivative_left= %e", derivative_left); getchar();
                                        //printf("\n derivative_left= %e", derivative_right); getchar();

                                        //cout<<"reached in between  222222"<<endl;
                                        //getchar();

										average_discontinuity += 1 / (pow(r2[j], r2_penalty))*abs((energy_left - energy_right) / (energy_left + energy_right) / 2)*exp(-energy_left / max_range);
										//printf(" check = %lf %d", r2[j], j);
										//printf("\n average_discontinuity= %e", average_discontinuity); getchar();
										average_velocity_discontinuity += 1 / pow(r2[j], r2_penalty)*abs((derivative_left - derivative_right) / (derivative_left + derivative_right) / 2)*exp(-energy_left / max_range);
										//printf("\n average_velocity_discontinuity= %e", average_velocity_discontinuity); getchar();

									}//end for j loop

									//cout<<"average_discontinuity = "<<average_discontinuity<<endl;
									//cout<<"average_velocity_discontinuity  = "<<average_velocity_discontinuity<<endl;
                                    //cout<<"press key to continue "<<endl; getchar();
                                    //cout<<"average_velocity_discontinuity = "<<average_velocity_discontinuity<<endl;
                                    //cout<<"average_velocity_discontinuity_selected = "<<average_velocity_discontinuity_selected<<endl;
                                    //cout<<"average_discontinuity = "<<average_discontinuity<<endl;
                                    //cout<<"average_discontinuity_selected = "<<average_discontinuity_selected<<endl;

									if (average_velocity_discontinuity < average_velocity_discontinuity_selected)
                                    {
										//printf("\n I am here %e %e", average_velocity_discontinuity, average_velocity_discontinuity_selected); getchar();

										average_velocity_discontinuity_selected = average_velocity_discontinuity;
										average_discontinuity_selected = average_discontinuity;
										//printf("\n \n Division_selected = ");
										for (int k = 0; k < length_dd; k++)
                                        {
											division_selected[k] = division[k];
											kindex_selected[k] = kindex[k];
											r2_selected[k] = r2[k];
											degree_selected[k] = degree[k];
											//printf("\n kindex = %lf", kindex[k]);
										} // end if
										//getchar();
										for (int i = 0; i < length_dd+1; i++)
										{
											for (int k = 0; k < maximum_degree + 1; k++) {
												coefficients_selected[i][k] = coefficients[i][k];
											} // end for i loop
										}	//end for k loop
									}//end if


                                /*
                                for (int i = 0; i < length_dd+1; i++)
                                {
                                    for (int k = 0; k < maximum_degree + 1; k++)
                                                cout<<"coeff= "<<coefficients_selected[i][k]<<endl;
                                    cout<<endl;
                                }

                                for (int k = 0; k < length_dd; k++)
                                    cout<<"kindex = "<<kindex_selected[k]<<endl;
                                cout<<endl;
                                cout<<"Press any key to continue...."<<endl;
                                cout<<"div3 = "<<div3<<endl;
                                cout<<"At the end of div3 loop after one iteration ..  Press any key to continue"<<endl<<endl<<endl;
                                getchar(); getchar();
                                */

								}/// end for div3
                                /*
                                cout<<"after completing div3 loop";

                                for (int i = 0; i < length_dd+1; i++)
                                {
                                    for (int k = 0; k < maximum_degree + 1; k++)
                                        cout<<"coefficient selected = "<<coefficients_selected[i][k]<<endl;
                                    cout<<endl;
                                }

                                for (int k = 0; k < length_dd; k++)
                                    cout<<"kindex = "<<kindex_selected[k]<<endl;
                                cout<<endl;
                                cout<<"Press any key to continue...."<<endl;
                                cout<<"div3 = "<<div3<<endl;
                                cout<<"div2 = "<<div2<<endl;
                                cout<<"At the end of div2 loop after one iteration "<<endl<<endl;

                                getchar();
                                */
							}/// end for div2
                            /*
                            cout<<"after completing div2 loop"<<endl;
                            for (int i = 0; i < length_dd+1; i++)
                            {
                                for (int k = 0; k < maximum_degree + 1; k++)
                                    cout<<"coefficient = "<<coefficients_selected[i][k]<<endl;
                                cout<<endl;
                            }

                            for (int k = 0; k < length_dd; k++)
                                cout<<"kindex = "<<kindex_selected[k]<<endl;
                            cout<<endl;
                            cout<<"Press any key to continue...."<<endl;
                            cout<<"div3 = "<<div3<<endl;
                            cout<<"div2 = "<<div2<<endl;
                            cout<<"div1 = "<<div1<<endl;
                            cout<<"At the end of div1 loop after one iteration "<<endl;
                            //getchar();
                            */

                        }/// end for div1
                        /*
                        cout<<"after completing div1 loop";

                        for (int i = 0; i < length_dd+1; i++)
                        {
                            for (int k = 0; k < maximum_degree + 1; k++)
                                cout<<"coefficient selected = "<<coefficients_selected[i][k]<<endl;
                            cout<<endl;
                        }

                        for (int k = 0; k < length_dd; k++)
                            cout<<"kindex_selected = "<<kindex_selected[k]<<endl;
                        cout<<endl;
                        cout<<"div3 = "<<div3<<endl;
                        cout<<"div2 = "<<div2<<endl;
                        cout<<"div1 = "<<div1<<endl;
                        cout<<"At the ned of d3 loop"<<endl;
                        cout<<"Press any key to continue...."<<endl;getchar();
                        */
					} ///end for d3
				}/// end for d2
			}/// end for d1
			test_negative = 0;
            /*
            cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<endl;
            for (int k = 0; k < maximum_degree + 1; k++)
                cout<<"coefficient = "<<coefficients[1][k]<<endl;
			cout<<"coefficient after 6 loops"<<endl; getchar();
            */

			powers = linspace(maximum_degree, 0, maximum_degree + 1);
			for (int k = 0; k < maximum_degree + 1; k++)
				test_negative = test_negative + coefficients[1][k] * pow(kindex[0], powers[k]);
			//printf("\n test_negative = %lf", test_negative);
		}// end while loop
	}// end for iteration loop
	//cout<<"ssssssssssssssaaaaaaaaaaaaaaaa"<<endl;
	end1 = clock();

	printf("\n");
	/*
	cout<<"coefficients final   =   "<<endl;

	for (int i = 0; i < length_dd+1; i++)
	{
		for (int j = 0; j < maximum_degree+1; j++) {
			coefficients[i][j] = coefficients_selected[i][j] ;
			cout<<coefficients[i][j]<<endl;
		}
		cout<<endl;
	}

    cout<<" kindex   =  "<<endl;
	for (int i = 0; i < length_dd; i++)
    {
			kindex[i] = kindex_selected[i];
			cout<<kindex[i]<<"    ";
    }
    */

    if (aa==1)
    {
        for (int i = 0; i < length_dd+1; i++)
        {
            for (int j = 0; j < maximum_degree+1; j++)
                coefficients_cond[i][j] = coefficients_selected[i][j] ;
        }

        for (int i = 0; i < length_dd; i++)
            kindex_cond[i] = kindex_selected[i];
    }
    else
    {
        for (int i = 0; i < length_dd+1; i++)
        {
            for (int j = 0; j < maximum_degree+1; j++)
                coefficients_val[i][j] = coefficients_selected[i][j] ;
        }

        for (int i = 0; i < length_dd; i++)
            kindex_val[i] = kindex_selected[i];
    }

	total_fitting_time = ((double)(end1 - start) / (CLOCKS_PER_SEC * 60));
	printf("Curve fitting time = %lf s \n", total_fitting_time);

    int *c = new int[2];
    c[0] = maximum_degree;
    c[1] = length_dd;

    return c;
}

