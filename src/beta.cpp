#include"main.h"

double beta(double T,double Ds[], double energy[], double k_grid[], double v[],  int points, int T_loop)
// it gives inverse screening length (beta) in units of (1/nm) in certain Fermi energy (E_F),
// and at Temperature (T)
{
    // According to equation (70) of Rode's book (book8):

    double k , de, k_step;
    double integral = 0;

    for(int counter = 0;counter <= points-2; counter++)
    {
        k = k_grid[counter];
        de = energy[counter+1]-energy[counter];
        k_step = k_grid[counter+1]-k_grid[counter];

        if (free_e ==0)
            integral = integral+ de*((Ds[counter]*1000)/volume1)*f0(energy[counter],E_F,T)*(1-f0(energy[counter],E_F,T));
            // unit is (1/nm)^3, 1000 is multiplied to convert volume from (1/Angstron)^3 to (1/nm)^3
        else
        {
            integral = integral+k_step*(k/pi)*(k/pi)*f0(energy[counter],E_F,T)*(1-f0(energy[counter],E_F,T));
            // unit is (1/nm)^3
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

    double bet = (e*e/(epsilon_s[T_loop]*epsilon_0*k_B*T)*integral*6.241509324e27);   // unit (1/nm)^2
    											// 	
    bet = pow(bet,0.5);   // unit (1/nm) 
    return bet;
    // Equation (70) of Rode's book (book8), conversion constant is to get beta in (1/nm) unit
}
