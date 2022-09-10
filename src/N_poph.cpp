
#include"main.h"

double N_poph(double omega, double T)  // Polar Optical Phonon Bose-Einstein Distribution; equation (114) of Rode's book
{
    double N = 1/(exp(h_bar*omega/(k_B*T))-1);
    return N;

}




