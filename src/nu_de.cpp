
#include"main.h"

double nu_de(double k,int counter,double T,double c[],double v)
//deformation potential scattering rate according to equation (112) or Rode's book
{

// From equation (112) of Rode's book (book8):
    double de =(k_B*T*pow(E_deformation_n,2)*k*k)/(3*pi*h_bar*h_bar*v*C_long)*(3-8*pow(c[counter],2)
                +6*pow(c[counter],4))*1e10*1.60217657/1e8;
    /*
    if (kk==1)
    {
        cout<<endl;
        cout<<"E_deformation_n = "<<E_deformation_n<<endl;
        cout<<"v = "<<v<<endl;
        cout<<"C_long = "<<C_long<<endl;
        cout<<"c[counter] = "<<c[counter]<<endl;
        cout<<"k =  "<<k<<endl;
        getchar();
    }
    */
    return de;
// 1e10 is coming from unit conversion (take a look at OneNote notes in Deformation potential section) and *1.60217657/1e8 is to get from cm/s (v) to hk/(md)
}

