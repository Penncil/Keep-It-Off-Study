// same with cov3C.c

// make a minor change compared with covC.c code in normal folder to adjust the generated large number in real data

#include <R.h>

void covcalC(int *p, int *n, double *yvec, int *index, int *count, double *xall, int *neach, double *sigma2, double *beta, double *mat1, double *mat2, double *phi)
{
	int i, ii, iii, j, jj, count1, ini, count2;
	count1=0;
	count2=0;
	double phis[p[0]];
	int istart[n[0]+1];
	double veca[p[0]];
	double phiall[n[0]*(n[0]-1)*p[0]];
	double temp;
	for (ini=0; ini<p[0]; ini++)
	{
		phi[ini]=0;
	}
	for (i=0; i<n[0]; i++)
	{
		istart[i]=count1;
		for (ii=0; ii<n[0]; ii++)
		{
			if (i !=ii)
			{
				for (ini=0; ini<p[0]; ini++)
				{
					phis[ini]=0;
				}
				for (j=0; j<count[i]; j++)
				for (jj=0; jj<count[ii]; jj++)
				{
					for (ini=0; ini<p[0]; ini++)
					{
						veca[ini]=(yvec[index[i]+j]-yvec[index[ii]+jj])*(xall[ini*neach[0]+index[i]+j]-xall[ini*neach[0]+index[ii]+jj])/sigma2[0];
					}
					temp=0;
					for (ini=0; ini<p[0]; ini++)
					{
						temp=temp+veca[ini]*beta[ini];
					}
					temp=exp(-temp);
					for (ini=0; ini<p[0]; ini++)
					{
						phi[ini] +=(1-1/(1+temp))*veca[ini];
						phis[ini] +=(1-1/(1+temp))*veca[ini];
					}
					for (ini=0; ini<(p[0]*p[0]); ini++)
					{
						mat1[ini]+=-veca[ini/p[0]]*veca[ini%p[0]]*temp/((1+temp)*(1+temp));
					}
				}
				for (ini=0; ini<p[0]; ini++)
				{
					phiall[count1*p[0]+ini]=phis[ini];
				}
				count1+=1;
			}
		}

	}
	for (ini=0; ini<p[0]*p[0]; ini++)
	{
		mat1[ini]=mat1[ini]/count1;
	}
	for (ini=0; ini<p[0]; ini++)
	{
		phi[ini]=phi[ini]/count1;
	}

	istart[n[0]]=count1;

	for (i=0; i<n[0]; i++)
	{
		for (ii=istart[i]; ii<istart[i+1]-1; ii++)
		{
			for (iii=ii+1; iii<istart[i+1]; iii++)
			{
				for (ini=0; ini<p[0]*p[0]; ini++)
				{
					count2 +=2;
					mat2[ini]+=4*(phiall[ii*p[0]+ini/p[0]]*phiall[iii*p[0]+ini%p[0]]-phi[ini/p[0]]*phi[ini%p[0]])/n[0];
					mat2[ini]+=4*(phiall[iii*p[0]+ini/p[0]]*phiall[ii*p[0]+ini%p[0]]-phi[ini/p[0]]*phi[ini%p[0]])/n[0];
				}
			}
		}
	}
	for (ini=0; ini<p[0]*p[0]; ini++)
	{
		mat2[ini]=mat2[ini]/(n[0]-1)/(n[0]-2);
	}
	p[0]=count2;
}





