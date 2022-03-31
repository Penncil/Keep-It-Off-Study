#include <R.h>

void mylikC(int *p, int *n, double *yvec, int *index, int *count, double *xall, int *neach, double *mypar, double *sigma2, double *result)
{
	int i, ii, j, jj, ini;
	result[0] =0;
	double temp;
	for (i=0; i<(n[0]-1); i++)
	for (ii=i+1; ii<n[0]; ii++)
	{
		for (j=0; j<count[i]; j++)
		for (jj=0; jj<count[ii]; jj++)
		{
			temp=0;
			for (ini=0; ini<p[0]; ini++)
			{
				temp +=(xall[ini*neach[0]+index[i]+j]-xall[ini*neach[0]+index[ii]+jj])*mypar[ini];
			}
			result[0] += -log(1+exp(-(yvec[index[i]+j]-yvec[index[ii]+jj])*temp/sigma2[0]));
		}
	}
}