#include"main.h"

double read_fermi()
{
	FILE *fid;
	char line[1000];
    double dummy,fermi;
    int dumm;

    fid = fopen("DOSCAR", "r");
    if (fid==NULL)
        cout<<"DOSCAR is not present";
    else
    {
        fgets(line, 1000, (FILE*)fid);
        fgets(line, 1000, (FILE*)fid);
        fgets(line, 1000, (FILE*)fid);
        fgets(line, 1000, (FILE*)fid);
        fgets(line, 1000, (FILE*)fid);
        fgets(line, 1000, (FILE*)fid);
        //cout<<line;
        sscanf(line, "%lf %lf %d %lf", &dummy, &dummy, &dumm, &fermi);  //ispin[3] is useful
        //cout<<"fermi = "<<fermi;
    }
    return fermi;
}
