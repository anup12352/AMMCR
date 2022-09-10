#include"main.h"

void read_procar()
{

	//------------------------------------  Reading Procar file not for conduction band dependent on ispin ----------
	//cout<<"fitting_2 = "<<fitting_2<<endl;
	if (fitting_2 == 0)
	{
		count_orbital = decompose_wave();
		//cout<<"count_orbital = "<<count_orbital<<endl;
	}
	else     // this part is for debugging code only
	{
		int i;
		ifstream cin("orbital_decomposedd.txt");

		i=0;
		while(!cin.eof())
		{
		    cin>>orbital_decomposedd[i][0]>>orbital_decomposedd[i][1]>>orbital_decomposedd[i][2]>>orbital_decomposedd[i][3];
		    /*

		    cout<<orbital_decomposedd[i][0]<<"    "<<orbital_decomposedd[i][1]<<"    "<<
		    orbital_decomposedd[i][2]<<"    "<<orbital_decomposedd[i][3]<<endl;
		    getchar();
		    */
		    i++;
		}
		count_orbital = i;
	}

	/*
	if(fitting_2 == 0)
	{
	    fid1 = fopen("orbital_decomposedd.txt","w");

	    for (int i = 0; i < count_orbital; i++)
		fprintf(fid1,"%d    %e    %e    %e    %e\n", i+1, orbital_decomposedd[i][0],
		orbital_decomposedd[i][1], orbital_decomposedd[i][2],orbital_decomposedd[i][3]);
	    fclose(fid1);

	}
	*/

	//---------------------------- Reading Procar file not for valence band dependent on ispin  -------------------

	//cout<<"fitting_2 = "<<fitting_2<<endl;
	if (fitting_2 == 0)
	{
		count_orbital_p = decompose_wave_p();
		//cout<<"count_orbital_p = "<<count_orbital_p<<endl;
	}
	else     // this part is for debugging code only
	{
		int i;
		ifstream cin("orbital_decomposedd_p.txt");

		i=0;
		while(!cin.eof())
		{
		    cin>>orbital_decomposedd_p[i][0]>>orbital_decomposedd_p[i][1]>>orbital_decomposedd_p[i][2]>>orbital_decomposedd_p[i][3];
		    /*
		    cout<<orbital_decomposedd_p[i][0]<<"    "<<orbital_decomposedd_p[i][1]<<"    "<<
		    orbital_decomposedd_p[i][2]<<"    "<<orbital_decomposedd_p[i][3]<<endl;
		    getchar();
		    */
		    i++;
		}
		count_orbital_p = i;
	}

	//-----------------------------------------------------------------------------------------------------------

}



