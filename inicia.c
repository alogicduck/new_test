#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "nrutil.h"
#include "nrutil.c"
#define EPS 3.0e-11
#define TINY 1.0e-30
#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#define cons_t (0.5*1/(pow(PI,2)))
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif
#define B 4.0
#define alf -0.216162
#define NR_END 1
/****************************************************************/

void ludcmp(double **a, int n, int *indx, double *d){
   int i,imax,j,k;
   double   big,dum,sum,temp;
   double   *vv;

   vv = (double*)malloc((unsigned)(n*sizeof(double)));
   if (!vv) 
     {
      fprintf(stderr,"Error Allocating Vector Memory\n");
      exit(1);
     }
   *d = 1.0;
   for (i=0;i<n;i++)      {
       big = 0.0;
       for (j=0;j<n;j++){
        if ((temp=fabs(a[i][j])) > big) big = temp;
       }
       if (big == 0.0){
       fprintf(stderr,"Singular Matrix in Routine LUDCMP\n");
       for (j=0;j<n;j++){ printf(" %f ",a[i][j]); printf("/n");

	 }
         }
       vv[i] = 1.0/big;
      }
   for (j=0;j<n;j++)
      {
       for (i=0;i<j;i++)
       {
        sum = a[i][j];
        for (k=0;k<i;k++) sum -= a[i][k] * a[k][j];
        a[i][j] = sum;
       }
       big = 0.0;
       for (i=j;i<n;i++)
          {
        sum = a[i][j];
        for (k=0;k<j;k++) sum -= a[i][k] * a[k][j];
        a[i][j] = sum;
        if ((dum=vv[i]*fabs(sum)) >= big)
          {
           big = dum;
           imax = i;
          }
       }
       if (j != imax)
         {
       for (k=0;k<n;k++)
          {
           dum = a[imax][k];
           a[imax][k] = a[j][k];
           a[j][k] = dum;
          }
          *d = -(*d);
       vv[imax] = vv[j];
         }
       indx[j] = imax;
       if (a[j][j] == 0.0) a[j][j] = TINY;
       if (j != n-1)
         {
       dum = 1.0 / a[j][j];
       for (i=j+1;i<n;i++) a[i][j] *= dum;
         }
      }
   free(vv);
  }

void gauleg(double x1,double x2,double *x,double *w,int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
	
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++){
		z=cos(PI*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++){
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
		
	}	
}
#undef EPS


double transf(double x){
double q;
q = x*pow(1-x,-1);
return q;	
	}

double jacob(double x){
double a;
a = pow(pow(x-1,2),-1);
return a;
}
double fp(double p){
double g,u=0.5;

g = pow((pow(transf(p),2)+pow(u,2)),-1);
return g;	
}

double g_q(double x){
double g,v=0.3;
g = pow((pow(x,2)+pow(v,2)),-1);
return g;	
}


double v(double x,double p){
double r;
r = fp(p)*g_q(transf(x));
return r;
}

double G_q(double x){
double a,mu = 1.3;

a = -1/(pow(2*mu,-1)*pow(transf(x),2)+B);
return a;
}

double fun(double x,double p)
{
double result;
result = cons_t*pow(transf(x),2)*v(x,p)*G_q(x);	
return result;
}

int main(int argc, char* argv[])
{
FILE *saida;
int i,j;
int n;
double **M,**K;
int *indx;
double *x,*w;
double d,var,a,y;

saida = fopen("saida0.txt","w");
printf("dimensao da matriz :");
scanf("%d",&n);
printf("\n");

K = dmatrix(0,n,0,n);
M = dmatrix(0,n,0,n);
indx = ivector(0,n);
w = dvector(0,n);
x = dvector(0,n);

gauleg(0,1,x,w,n);

/*montando a matriz K */

for(i=1;i<=n;i++){
	for(j=1;j<=n;j++){
		a = x[j];
		y = x[i];
		K[i][j] = alf*fun(transf(a),transf(y))*jacob(transf(y))*w[j];
		printf("fun : %.10e \n",fun(transf(a),transf(y)));
		printf("x[%d] : %lf \n",j,x[j]);
		printf("w[%d] : %lf \n",j,w[j]);
		printf("K[%d][%d] : %.10e \n", i,j,K[i][j]);
		printf("\n");
		fprintf(saida,"%lf  %.10e \n",x[j],K[i][j]);
		}
	}
fclose(saida);
	
	
for(i=1;i<=n;i++){
	for(j=1;j<=n;j++){
		if (i!= j) M[i][j] = -K[i][j];
		else M[i][j] = 1 - K[i][j];
		}
	}

//~ for(i=1;i<=n;i++){
	//~ for(j=1;j<n;j++){
		//~ printf("M[%d][%d]: %.10e \n",i,j,M[i][j]);
		//~ }
	//~ }

ludcmp(M,n,indx,&d);


//~ for(i=0;i<n;i++){
	//~ for(j=0;j<n;j++){
		//~ printf("M[%d][%d]: %.10e \n",i,j,M[i][j]);
		//~ }
	//~ }


/*Calcula o determinante da matriz M*/

var = d;
for(j=1;j<=n;j++){
 var *= M[j][j]; 
}
printf("\n");
printf("determinante : %lf ",var);


free_dmatrix(K,0,n,0,n);
free_dmatrix(M,0,n,0,n);
free_ivector(indx,0,n);
free_dvector(x,0,n);
free_dvector(w,0,n);

return 0;
}
