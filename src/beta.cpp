#include"main.h"

double beta(double T,int T_loop)
// it gives inverse screening length (beta) in units of (1/nm) in certain Fermi energy (efef_n),
// and at Temperature (T)
{
    // According to equation (70) of Rode's book (book8):

    double k , de, k_step;
    double integral = 0;

    for(int counter = 0;counter <= points-2; counter++)
    {
        k = k_grid[counter];
        de = energy_n[counter+1]-energy_n[counter];
        if(de<0)
        	de = 0;
        k_step = k_grid[counter+1]-k_grid[counter];

        if (free_e ==0)  // DOSCAR is used
        {
           integral = integral+ de*((1e6*Ds_n[counter])/volume1)*f0(energy_n[counter],efef_n,T)*(1-f0(energy_n[counter],efef_n,T));
            
            // unit is (1/m)^3; 1e6 is multiplied to convert volume from (1/cm)^3 to (1/m)^3
        }
        else
        {
            integral = integral+k_step*(k/pi)*(k/pi)*f0(energy_n[counter],efef_n,T)*(1-f0(energy_n[counter],efef_n,T))*1e27;
            // unit is (1/m)^3
            // Part of equation (70) or Rode's book (book8)
        }
        //cout<<"counter = "<<counter+1<<endl;
        //cout<<"k = "<<k<<endl;
        //cout<<"de = "<<de<<endl;
        //cout<<"k_step = "<<k_step<<endl;
        //cout<<"integral = "<<integral<<endl;
        //getchar();
    }
    //cout<<"integral = "<<integral<<endl;
    //getchar();
    
    //cout<<"epsilon_s[T_loop] = "<<epsilon_s[T_loop]<<endl;

    double bet = (e*e/(epsilon_s[T_loop]*epsilon_0*(k_B*e)*T)*integral);   // unit (1/m)^2
    											// 	
    bet = pow(bet,0.5);   // unit (1/m) 
    
    bet = bet*1e-9;
    //converted from 1/m to 1/nm
    return bet;
    // Equation (70) of Rode's book (book8), conversion constant is to get beta in (1/nm) unit
}
