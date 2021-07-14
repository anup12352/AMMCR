#include"main.h"

double max_array(double band[],int county)
{
    double max1 = -3e20;
	for (int i = 0; i < county; i++)
	{
		if (band[i] > max1)
			max1 = band[i];
	}
	return max1;
}


double Aplus(int counter,double k_grid[],double omega, double a[],double c[],double energy[], int points)
{
    double k_plus = kplus(counter,k_grid,omega,energy,points);
    double AA;
    double max1 = max_array(energy,points);
    if (max1 < energy[counter] + h_bar*omega)
        AA =0;
    else
    {
        double arr[points];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k_plus);
        int plus_index =FindMinInd(arr,points);
        double k = k_grid[counter];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k);
        int index =FindMinInd(arr,points);

        AA = a[index]*a[plus_index]+(k_plus*k_plus+k*k)/(2*k_plus*k)*c[index]*c[plus_index];
    }

    return AA;
}

