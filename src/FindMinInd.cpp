
#include "main.h"
int FindMinInd(double arr[],int length)
{
    double min1 = arr[0];
    int index = 0;
    for(int i=0;i<length;i++)
    {
        if (arr[i] < min1)
        {
            min1 = arr[i];
            index = i;
        }
    }
    return index;
}
