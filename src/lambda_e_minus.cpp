#include"main.h"

double lambda_e_minus(int counter,double omega,double rho,double De,int nfv,int points)
{
    //cout<<endl<<"Inside lambda_e_minus"<<endl;

    double k_minus = kminus(counter,omega,points,energy_n);
    //cout<<"k_minus = "<<k_minus<<endl;

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_minus);
    int minus_index =FindMinInd(arr,points);
    //cout<<"minus_index = "<<minus_index<<endl;

    double k = k_grid[counter];
    //cout<<" k = "<<k<<endl;

    double l;

    if ((energy_n[counter]<h_bar*omega)||(k_minus==k))
        l = 0;
         //If the energy of electron is lower than the phonon, there will be no
         //elmission so lambda_minus terms will be zero
         //(page 40 of Semiconductors and Semimetals, volume 10)
    else
            l = pow((e*De*1e10),2)*nfv*(k_minus*(1e9))*(k*(1e9))/(2*pi*rho*omega*(h_bar*e)*v_n[counter]*1e-2);

    //cout<<"l = "<<l<<endl;
    //cout<<"End of  lambda_e_minus"<<endl;
    //getchar();
    return l;
}

