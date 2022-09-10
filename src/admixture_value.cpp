#include "main.h"

double admixture_value(double k,int coloumn)
{
    double arr[count_orbital];
    for(int i=0;i<count_orbital;i++)
    {
        arr[i] = abs(orbital_decomposedd[i][0]- k);
    }

    int index;

    index = FindMinInd(arr, count_orbital);

    //cout<<"index = "<<index<<endl;

    double c = orbital_decomposedd[index][coloumn-1];

    return c;
}
