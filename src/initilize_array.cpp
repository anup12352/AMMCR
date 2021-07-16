#include"main.h"

void initialize_array()
{

//------------------------------------- generating varaibles for doping or temperature variation -------------------------    
	    if (variation==1)
	    {
		for (int i=0;i<count_d;i++)
		{
		    calc_mobility[i][0] = n_array[i];
		    calc_mobility_rta[i][0] = n_array[i];
		    calc_thermopower[i][0] = n_array[i];
		    calc_sigma[i][0] = n_array[i];
		    calc_sigma_rta[i][0] = n_array[i];

		    calc_mobility_ii[i][0] = n_array[i];
		    calc_mobility_po[i][0] = n_array[i];
		    calc_mobility_de[i][0] = n_array[i];
		    calc_mobility_pe[i][0] = n_array[i];
		    calc_mobility_dis[i][0] = n_array[i];
		    calc_mobility_to[i][0] = n_array[i];
		    calc_mobility_alloy[i][0] = n_array[i];
		    calc_mobility_iv[i][0] = n_array[i];
		    calc_mobility_neutral[i][0] = n_array[i];


		    calc_mobility_hall[i][0] = n_array[i];
		    calc_mobility_hall_rta[i][0] = n_array[i];
		    calc_sigma_hall[i][0] = n_array[i];
		    calc_sigma_hall_rta[i][0] = n_array[i];

		    calc_mobility_hall_ii[i][0] = n_array[i];
		    calc_mobility_hall_po[i][0] = n_array[i];
		    calc_mobility_hall_de[i][0] = n_array[i];
		    calc_mobility_hall_pe[i][0] = n_array[i];
		    calc_mobility_hall_dis[i][0] = n_array[i];
		    calc_mobility_hall_to[i][0] = n_array[i];
		    calc_mobility_hall_alloy[i][0] = n_array[i];
		    calc_mobility_hall_iv[i][0] = n_array[i];
		    calc_mobility_hall_neutral[i][0] = n_array[i];

		    hall_factor[i][0] = n_array[i];
		    hall_factor_rta[i][0] = n_array[i];
		}
	    }
	    else
	    {
		for (int i=0;i<count_t;i++)
		{
		    calc_mobility[i][0] = T_array[i];
		    calc_mobility_rta[i][0] = T_array[i];
		    calc_thermopower[i][0] = T_array[i];
		    calc_sigma[i][0] = T_array[i];
		    calc_sigma_rta[i][0] = T_array[i];

		    calc_mobility_ii[i][0] = T_array[i];
		    calc_mobility_po[i][0] = T_array[i];
		    calc_mobility_de[i][0] = T_array[i];
		    calc_mobility_pe[i][0] = T_array[i];
		    calc_mobility_dis[i][0] = T_array[i];
		    calc_mobility_to[i][0] = T_array[i];
		    calc_mobility_alloy[i][0] = T_array[i];
		    calc_mobility_iv[i][0] = T_array[i];
		    calc_mobility_neutral[i][0] = T_array[i];


		    calc_mobility_hall[i][0] = T_array[i];
		    calc_mobility_hall_rta[i][0] = T_array[i];
		    calc_sigma_hall[i][0] = T_array[i];
		    calc_sigma_hall_rta[i][0] = T_array[i];

		    calc_mobility_hall_ii[i][0] = T_array[i];
		    calc_mobility_hall_po[i][0] = T_array[i];
		    calc_mobility_hall_de[i][0] = T_array[i];
		    calc_mobility_hall_pe[i][0] = T_array[i];
		    calc_mobility_hall_dis[i][0] = T_array[i];
		    calc_mobility_hall_to[i][0] = T_array[i];
		    calc_mobility_hall_alloy[i][0] = T_array[i];
		    calc_mobility_hall_iv[i][0] = T_array[i];
		    calc_mobility_hall_neutral[i][0] = T_array[i];
		    hall_factor[i][0] = T_array[i];
		    hall_factor_rta[i][0] = T_array[i];

		}
	    }
//----------------------------------------------- varaibles for doping or temperature variation generated-------------------------
	
}
