
#include "main.h"
double rsquare(double data[], double fitted[], int length)
{
    double sum1=0,sum2=0,sum3=0,r2,num;
    for(int i=0;i<length;i++)
    {
        sum1 = sum1 + data[i];
    }
    double average;
    average= sum1/length;

    for(int i=0;i<length;i++)
    {
         sum2 = sum2 + pow(data[i]-average,2);
         sum3 = sum3 + pow(fitted[i]-data[i],2);
    }
    double denom = sum2;

    num = sum3;

    r2 = 1 - num/denom;
    return r2;
}
