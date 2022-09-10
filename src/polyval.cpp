
#include "main.h"

double poly[4000]={0};

void polyval(double p[], double x[],  int n, int m)  // n -length of coefficient array p   m--- length of x array
{
    //cout<<"Inside polyval"<<endl;
    //cout<<"m =   " <<m<<endl;
    //cout<<"n =   "<<n<<endl;

   for(int j=0;j<m;j++)
    {
        //cout<<"x[j] = "<<x[j]<<endl;
        poly[j]=0;
        for (int i=0; i<n; i++)
            {
                poly[j] = poly[j] + p[i]*pow(x[j],(n-i-1));
                //cout<<"poly in between = "<<poly[j]<<endl;
                //cout<<"(n-i-1) = "<<(n-i-1)<<endl;
                //cout<<"p[i] = "<<p[i]<<endl;
            }
            //cout<<"poly -----  "<<poly[j]<<endl;
            //getchar();
    }
}

