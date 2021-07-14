#include"main.h"

double nu_ii(double k, int counter, double B_ii, double D_ii, double beta_constant, double N_ii, double v, double epsilon_s)
{
    /*
    cout<<"k = "<<k<<endl;
    cout<<"counter = "<<counter<<endl;
    cout<<"B_ii = "<<B_ii<<endl;
    cout<<"D_ii = "<<D_ii<<endl;
    cout<<"beta_constant = "<<beta_constant<<endl;
    cout<<"N_ii = "<<N_ii<<endl;
    cout<<"v = "<<v<<endl;
    cout<<"epsilon_s = "<<epsilon_s<<endl;
    */

    //B_ii = 1.436011541671097e+04;
    //D_ii =  7.097643979710515e+09;
    //beta_constant =  0.140607315390570;
    //N_ii = 2.999793089963898e+16;
    //v =   2.181315148394955e+04;


    double xx = (pow(e,4)*abs(N_ii));
    double yy = (8*pi*v*epsilon_s*epsilon_s*epsilon_0*epsilon_0*h_bar*h_bar*k*k);
    double zz = (D_ii*log(1+4*k*k/(beta_constant*beta_constant))-B_ii);

    //cout<<"xx = "<<xx<<endl;
    //cout<<"yy = "<<yy<<endl;
    //cout<<"zz = "<<zz<<endl;

    double ii = abs( xx/yy*zz*3.89564386e27 );


    return ii;
}
