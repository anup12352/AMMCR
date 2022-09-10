#include "main.h"

double admixture_value_p(double k,int coloumn)
{
    double arr[count_orbital_p];
    for(int i=0;i<count_orbital_p;i++)
    {
        arr[i] = abs(orbital_decomposedd_p[i][0]- k);
    }

    int index;

    index = FindMinInd(arr, count_orbital_p);

    //cout<<"index = "<<index<<endl;

    double c = orbital_decomposedd_p[index][coloumn-1];

    return c;
}
