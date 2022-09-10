
#include"main.h"

double lambda_e_plus(int counter, double omega,double rho,double De,int nfv,int points)
{
    //cout<<endl<<"Inside lambda_e_plus"<<endl;

    double k_plus = kplus(counter,omega,points,energy_n);
    //cout<<"k_plus = "<<k_plus<<endl;

    double arr[points];
    for (int i=0;i<points;i++)
        arr[i] = abs(k_grid[i] - k_plus);
    int plus_index =FindMinInd(arr,points);
    //cout<<"plus_index = "<<plus_index<<endl;

    double k = k_grid[counter];
    //cout<<" k = "<<k<<endl;
    double l;

    if(k_plus == k)
        l = 0;
    else
        l = pow((e*De*1e10),2)*nfv*(k_plus*(1e9))*(k*(1e9))/(2*pi*rho*omega*(h_bar*e)*v_n[counter]*1e-2);

   //cout<<"rho = "<<rho<<endl;
   //getchar();	

   //cout<<"pow((e*De*1e10),2)  =  "<<pow((e*De*1e10),2)<<endl;
   //cout<<"nfv  =  "<<nfv<<endl;
   //cout<<"(k_plus*(1e9))*(k*(1e9))  =  "<<(k_plus*(1e9))*(k*(1e9))<<endl;
   //cout<<"2*pi*rho*omega*(h_bar*e)*v_n[counter]*1e-2  =  "<<2*pi*rho*omega*(h_bar*e)*v_n[counter]*1e-2<<endl;	

    //cout<<"l = "<<l<<endl;
    //cout<<"End of  lambda_e_plus"<<endl;
    //getchar();
    return l;
}

//l = pow((e*De*1e10),2)*nfv*(k_plus*(1e9))*(k*(1e9))/(2*pi*rho*omega*(h_bar*e)*v_n[counter]*1e-2);  

