
#include "main.h"

double dedk(double k, double coeff[5][7], double kindex[4], int a[2])
{
    int segment_of_the_curve=0;

    for (int j=0;j<a[0];j++)
    {
        if (k<=kindex[j])
        {
            segment_of_the_curve=j;
            break;
        }
        else
            segment_of_the_curve = a[0];
    }
    //cout<<"segment_of_the_curve = "<<segment_of_the_curve<<endl;
    double *powers;
    powers = linspace(a[0],0,a[0]+1);

    //for(int i=0;i<a[0]+1;i++)
    //cout<<"powers = "<<powers[i]<<endl;

    double derivative=0;
    for(int i=0;i<a[0];i++)
        derivative = derivative +  powers[i] * (pow(k,powers[i]-1))*coeff[segment_of_the_curve][i];
    //cout<<"derivative = "<<derivative<<endl;

    return derivative;
}



