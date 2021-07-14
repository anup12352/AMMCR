
#include"main.h"

//double k_min0,k_trans0,k_step_fine0,k_step0;
double Ed;

void read_input()
{
    read_input_file();
    fitting_1 = 0;  // FOR BAND
    fitting_2 = 0;  // FOR PROCAR
    fitting_3 = 1;  // FOR DOSCAR
	
    /*
    len_T =1;
    len_n = 3;

    variation = 1;  // 0 - means temperature and 1 means doping variation
    iterations = 10;


    for (int i=0;i<30;i++)
    {
        T_array[i]=0; epsilon_s[i]=0; epsilon_inf[i]=0; Bgap[i]=0; P_piezo[i] = 0;
        C_piezo_h14[i]=0; n_array[i]=0;  Nd[i]=0; Na[i]=0; N_im[i]=0;
    }

    for (i=0;i<10;i++)
    {
        we[i]=0; nfv[i]=0; De[i]=0;
    }

    T_array[0] = 200;
    n_array[0] = 3e15;   n_array[1] = 3e16;    n_array[2] = 3e17;
    Nd[0] = 3e15;        Nd[1] = 3e16;         Nd[2] = 3e17;
    Na[0] = 0;           Na[1] = 0;            Na[2] = 0;
    N_im[0] = 3e15;      N_im[1] = 3e15;       N_im[2] = 3e15;

    De_ionization = 0;
    Ed = 0.0092; // in eV   ionization energy  For ag 1 paper

    if (De_ionization==0)
        Ed = 0; // in eV   ionization energy

    /*
    // General Constant
        epsilon_s[0] = 7.45;
        epsilon_inf[0] = 3.45;
        m =0;   //
        m_h = 0;  //
        Bgap[0]= 0;       // 2.78;  // At 300
        N_vb = 1;  // No. of valley at valence band maxima
        N_cb = 1;  // No. of valley at conduction band minima along L diection [111]
        c_lattice = 0.562;  // in nm
        spin_orbit_coupling = '0';
        k_max = 6;
        rho = 0;

        scattering_mechanisms[0] = 1;  // Ionized impurity
        scattering_mechanisms[1] = 1;  // POP
        scattering_mechanisms[2] = 0;  // NPOP
        scattering_mechanisms[3] = 1;   // deformation
        scattering_mechanisms[4] = 1;   // piezoelectric
        scattering_mechanisms[5] = 1;   // dislocation
        scattering_mechanisms[6] = 0;   // To phonon
        scattering_mechanisms[7] = 1;   // alloy
        scattering_mechanisms[8] = 1;   // intervalley
        scattering_mechanisms[9] = 1;   // neutral impurity

        // [ionized impurity(1,1); polar optical phonon (POP)(2,1); NPOP(3,1); deformation(4,1); piezoelectric(5,1);
        // dislocation(6,1); TO phonon(7,1); alloy (8,1); intervalley(9,1); neutral impurity(10,1)];

        // dislocation scattering
            N_dis =100;

        // POP scattering constant
            omega_LO =  5.88*2*pi*1e12; //7.50*pi*1e12;     // from paper    // 5.88 - abinitio calculated
            omega_TO = 0.0;

        // Acoustic deformation potential constant
            E_deformation_n =12 ;

        // Piezoelectric scattering constant
            P_piezo[0] =  0.0392;

        // Spherically averaged elastic constants for longitudinal and transverse
        //modes
            C_long = 0;     // in dyne/cm^2  //10.34e11;
            C_trans = 0;    // in dyne/cm^2  //3.29e11;


            C_11 = 7.99e11; // abinitio value // elastic constant in dyne/cm2
            C_12 = 4.54e11; // abinitio value  // elastic constant in dyne/cm2
            C_44 = 3.71e11; // abinitio value  // elastic consstant in dyne/cm2
            C_piezo_c14 = -0.08162 ;  // in C/m2 ab initio calculated = ( -0.72919 + 0.64757)

            for(i=0;i<len_T;i++)
                C_piezo_h14[i] = C_piezo_c14/(epsilon_s[i]*epsilon_0);   // unit - N/C

            // Spherically averaged elastic constants for longitudinal and transverse
            // modes
            if (C_trans == 0)
                C_trans =  (C_11 - C_12 + 3 * C_44)/5;  // in dyne/cm2   Equation 99 from rode book

            if (C_long == 0)
                C_long =  (3*C_11 + 2 * C_12 + 4 * C_44)/5;         // in dyne/cm2   Equation 100 from rode book

            c_bar = (1.0/3.0)*C_long + (2.0/3)*C_trans;   // in dyne/cm^2

            if (P_piezo[0] == 0)
                for(i=0;i<len_T;i++)
                {
                    P_piezo[i] = (pow(C_piezo_h14[i],2)*epsilon_0*epsilon_s[i]*(12/C_long+16/C_trans)/35*1e1 );
                    P_piezo[i] = pow(P_piezo[i],0.5);
                // P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
                }

        // Alloy scattering constants
            Uall = 0.53;  // Alloy Potential in eV
            V0 = 0.562*0.562*0.562; // volume of the primitive cell unit (nm)^3
            xx = 0.47; //

        // Intervalley scattering constants
            iv_number = 1;     // No. of intervalley scattering to include   // 3 for f type and 3 for g type

            // conductivity effective mass and density of states effective mass
            if (scattering_mechanisms[9]==1)
            {
                we[0] = 2*pi*7.96*1e12;      // unit in 1/s
                De[0] = 3;      // unit (*10^8 eV/cm )
                nfv[0] = 3;    // number of final valley's
            }

            free_e = 0;  //  0 means - false   1 means - true
            type = "n";
            T_trans = 40;

            kcbm[0]=0; kcbm[1]=0; kcbm[2]=0;
            kvbm[0]=0; kvbm[1]=0; kvbm[2]=0;

            fitting_1 = 1;  // FOR BAND
            fitting_2 = 0;  // FORPROCAR
            fitting_3 = 1;  // FOR DOSCAR
        */
}
