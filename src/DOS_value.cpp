

#include "main.h"

double DOS_value(double energy, int a)
{
    double DOS[2000][2];
    int countx;
    if (a==1)
    {
        for(int i=0;i<count_DOS_n;i++)
        {
            DOS[i][0] = DOS_n[i][0];
            DOS[i][1] = DOS_n[i][1];
        }
        countx = count_DOS_n;
    }
    else
    {
        for(int i=0;i<count_DOS_p;i++)
        {
            DOS[i][0] = DOS_p[i][0];
            DOS[i][1] = DOS_p[i][1];
        }
        countx = count_DOS_p;
    }

    double arr[countx],ds,slope,energy_earlier;
    for(int i=0;i<countx;i++)
        arr[i] = abs(DOS[i][0]- energy) ;

    int index;
    index = FindMinInd(arr, countx);

    double temp;
    temp = DOS[index][0];

    int index_earlier,index_after;

    if (temp > energy)
    {
        index_earlier = index-1;
        index_after = index;
    }
    else
    {
        index_earlier = index;
        index_after = index+1;
    }

    if (index_earlier == -1)
    {
        slope = (DOS[index_after][1])/(DOS[index_after][0]);
        // (y2 - y1)/(x2 - x1)   delta DOS/ delta energy
        energy_earlier = 0;
        ds = slope * (energy - energy_earlier) + 0 ;
    }
    else
    {
        slope = (DOS[index_after][1] - DOS[index_earlier][1])/(DOS[index_after][0] - DOS[index_earlier][0]);
        // (y2 - y1)/(x2 - x1)   delta DOS/ delta energy
        energy_earlier = DOS[index_earlier][0];
        ds = slope * (energy - energy_earlier) + DOS[index_earlier][1];
    }
    
    return ds;
}
