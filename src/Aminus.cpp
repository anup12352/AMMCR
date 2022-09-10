#include"main.h"

double Aminus(int counter, double omega, int points, double energy[])
{
    double k_minus =kminus(counter,omega,points,energy);
    double AA;
    if (energy[counter] < h_bar*omega)
        AA =0;
    else
    {
        double arr[points];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k_minus);
        int minus_index =FindMinInd(arr,points);
        double k = k_grid[counter];
        for (int i=0;i<points;i++)
            arr[i] = abs(k_grid[i] - k);
        int index =FindMinInd(arr,points);

        AA = a_n[index]*a_n[minus_index]+(k_minus*k_minus+k*k)/(2*k_minus*k)*c_n[index]*c_n[minus_index];
        
        /*
        cout<<"minus_index = "<<minus_index<<endl;
        cout<<"k = "<<k<<endl;
        cout<<" index = "<<index<<endl;
        cout<<"a[index] = "<<a[index]<<endl;
        cout<<"a[minus_index] = "<<a[minus_index]<<endl;
        cout<<" = "<<<<endl;
        */
    }
    
    //cout<<"AA = "<<AA<<endl;
    return AA;
}
