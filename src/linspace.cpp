
#include "main.h"

double * linspace(double a, double b, int numbers)
{
    double step;
    double *arr = new double[numbers];

    step = (b-a)/(numbers-1);
    arr[0] = a;
    arr[numbers-1] = b;
    for(int i=1;i<=numbers-1;i++)
    {
        arr[i] = arr[i-1] + step;
    }
    return arr;
}

