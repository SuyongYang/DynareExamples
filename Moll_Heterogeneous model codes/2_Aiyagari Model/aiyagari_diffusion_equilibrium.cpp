#include <stdio.h>
#include <math.h>
#include "umfpack.h"

//using namespace std;

const double ga = 2;
const double alpha = 0.35;
const double delta = 0.1;
const double zmean = 1.0;
const double sig2 = 0.01;
const double Corr = exp(-0.3);
const double rho = 0.05;

double K = 3.8;             //Initial Guess of Capital Level
double r = alpha*pow(K,alpha-1)-delta;
double w = (1-alpha)*pow(K,alpha);
const double relax = 0.99;
const int J = 40;
const double zmin = 0.5;
const double zmax = 1.5;
const double amin = -1;
const double amax = 30;
const int I = 100;
const double da = (amax-amin)/(I-1);
const double dz = (zmax-zmin)/(J-1);
const double dz2 = dz*dz;
double a[I];
double z[J];
const int n = I*J;

const int maxit = 100;
const int maxitK = 100;
const double crit = 1e-6; 
const double critK = 1e-5;
const double Delta = 1000;

const double the = -log(Corr);
const double Var = sig2/(2*the);

//*** PREALLOCATION ***//
double v[n];
double vold[n];
double gg[n];
double b[n];
double S;

double mu[J];
double s2[J];

int mat(int i, int j){
  return j*I+i;
};

int main(){
  int Ap[n+1];
  int Ai[5*n-2*I-2*J];
  double Ax[5*n-2*I-2*J];
  double x[n];
  double centdiag[n];
  double tmpmax;
  double g_sum;
  double Axt[5*n-2*I-2*J];
  double Info [UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  double vaf, vab, caf, cab, ca0, va0, va_upwind=0, u, c, sf, sb;

  // a=linspace(amin,amax,I);
  for (int i=0;i<I;i++){
    a[i]=i*(amax-amin)/(I-1)+amin;
  }
  
  // z=linspace(zmin,zmax,J);
  // mu=the*(zmean-z);
  for (int j=0;j<J;j++){
    z[j]=j*(zmax-zmin)/(J-1)+zmin;
    mu[j]=the*(zmean-z[j]);
  }

  for (int i=0;i<I;i++) {
    for (int j=0;j<J;j++){
      v[j*I+i]=pow(w*z[j]+r*a[i],1-ga)/(1-ga)/rho;
    }
  }

  // Getting the Diagonal Shape of Our matrix //
  // Afterward fill in offdiagnonal elements that depend only on External Skill Level //
  Ap[0]=0;
  Ai[0]=0;
  Ai[1]=1;
  Ai[2]=I;
  Ax[2]=-mu[0]/dz-sig2/2/dz2;
  centdiag[0]=-mu[0]/dz-sig2/2/dz2;
  Ap[1]=3;
  for (int i=1;i<I-1;i++){
    Ai[4*i-1]=i-1;
    Ai[4*i]=i;
    Ai[4*i+1]=i+1;
    Ai[4*i+2]=I+i;
    Ax[4*i+2]=-mu[0]/dz-sig2/2/dz2;
    centdiag[i]=-mu[0]/dz-sig2/2/dz2;    
    Ap[i+1]=4*(i+1)-1;
  }
  Ai[4*(I-1)-1]=I-2;
  Ai[4*(I-1)]=I-1;
  Ai[4*(I-1)+1]=2*I-1;
  Ax[4*(I-1)+1]=-mu[0]/dz-sig2/2/dz2;
  centdiag[I-1]=-mu[0]/dz-sig2/2/dz2;
  Ap[I]=4*I-2;
  for (int j=1;j<J-1;j++){
    Ai[j*(5*I-2)-I]=I*(j-1);
    Ai[j*(5*I-2)-I+1]=I*j;
    Ai[j*(5*I-2)-I+2]=I*j+1;
    Ai[j*(5*I-2)-I+3]=I*(j+1);
    Ax[j*(5*I-2)-I]=-sig2/2/dz2;
    Ax[j*(5*I-2)-I+3]=-mu[j]/dz-sig2/2/dz2;
    centdiag[j*I]=-mu[j]/dz-sig2/dz2;
    Ap[I*j+1]=j*(5*I-2)-I+4;
    for (int i=1;i<I-1;i++){
      Ai[j*(5*I-2)-I+5*i-1]=I*(j-1)+i;
      Ai[j*(5*I-2)-I+5*i]=I*j+i-1;
      Ai[j*(5*I-2)-I+5*i+1]=I*j+i;
      Ai[j*(5*I-2)-I+5*i+2]=I*j+i+1;
      Ai[j*(5*I-2)-I+5*i+3]=I*(j+1)+i;
      Ax[j*(5*I-2)-I+5*i-1]=-sig2/2/dz2;
      Ax[j*(5*I-2)-I+5*i+3]=-mu[j]/dz-sig2/2/dz2;
      centdiag[j*I+i]=-mu[j]/dz-sig2/dz2;
      Ap[I*j+i+1]=j*(5*I-2)-I+5*i+4;
    }
    Ai[j*(5*I-2)+4*I-6]=I*j-1;
    Ai[j*(5*I-2)+4*I-5]=I*(j+1)-2;
    Ai[j*(5*I-2)+4*I-4]=I*(j+1)-1;
    Ai[j*(5*I-2)+4*I-3]=I*(j+2)-1;
    Ax[j*(5*I-2)+4*I-6]=-sig2/2/dz2;
    Ax[j*(5*I-2)+4*I-3]=-mu[j]/dz-sig2/2/dz2;
    centdiag[(j+1)*I-1]=-mu[j]/dz-sig2/dz2;
    Ap[I*j+I]=j*(5*I-2)+5*I-2-I;
  }
  Ai[(J-1)*(5*I-2)-I]=(J-2)*I;
  Ai[(J-1)*(5*I-2)-I+1]=(J-1)*I;
  Ai[(J-1)*(5*I-2)-I+2]=(J-1)*I+1;
  Ax[(J-1)*(5*I-2)-I]=-sig2/2/dz2;
  centdiag[(J-1)*I]=-sig2/2/dz2;
  Ap[I*(J-1)+1]=(J-1)*(5*I-2)-I+3;
  for (int i=1;i<I-1;i++){
    Ai[(J-1)*(5*I-2)-I+4*i-1]=(J-2)*I+i;
    Ai[(J-1)*(5*I-2)-I+4*i]=(J-1)*I+i-1;
    Ai[(J-1)*(5*I-2)-I+4*i+1]=(J-1)*I+i;
    Ai[(J-1)*(5*I-2)-I+4*i+2]=(J-1)*I+i+1;
    Ax[(J-1)*(5*I-2)-I+4*i-1]=-sig2/2/dz2;
    centdiag[(J-1)*I+i]=-sig2/2/dz2;
    Ap[I*(J-1)+i+1]=(J-1)*(5*I-2)-I+4*i+3;
  }
  Ai[(J-1)*(5*I-2)-I+4*(I-1)-1]=(J-1)*I-1;
  Ai[(J-1)*(5*I-2)-I+4*(I-1)]=J*I-2;
  Ai[(J-1)*(5*I-2)-I+4*(I-1)+1]=J*I-1;
  Ax[(J-1)*(5*I-2)-I+4*(I-1)-1]=-sig2/2/dz2;
  centdiag[J*I-1]=-sig2/2/dz2;
  Ap[I*J]=J*(5*I-2)-2*I;

  // Prepare Matrix for Fokker-Plank Equation //
  for (int i=0;i<5*J*I-2*J-2*I;i++){
    Axt[i]=-Ax[i];
  }

  (void) umfpack_di_symbolic(n,n,Ap,Ai,Ax,&Symbolic,null,Info);

  // Capital Iteration //
  for (int iter=0; iter<maxitK; iter++) {
    
    // Solve HBJ Equation //
    for (int it=0; it<maxit; it++) {
      // vold=v;
      for (int i=0;i<n;i++){
	vold[i]=v[i];
      }

      // Set Solution Equation //
      va_upwind=0;
      vaf=(v[1]-v[0])/da;
      caf=pow(vaf,-1/ga);
      sf=w*z[0]+r*a[0]-caf;
      ca0=w*z[0]+r*a[0];
      va0=pow(ca0,-ga);
      if (sf>0) {
	va_upwind+=vaf;
	Ax[0]=1/Delta+rho-centdiag[0]+sf/da;
	Ax[1]=-sf/da;
	Axt[0]=centdiag[0]-sf/da;
	Axt[1]=sf/da;
      }
      else {
	va_upwind+=va0;
	Ax[0]=1/Delta+rho-centdiag[0];
	Ax[1]=0;
	Axt[0]=centdiag[0];
	Axt[1]=0;
      }
      c=pow(va_upwind,-1/ga);
      u=pow(c,(1-ga))/(1-ga);
      b[0]=u+v[0]/Delta;
      for (int i=1;i<I-1;i++) {
	va_upwind=0;
	vaf=(v[i+1]-v[i])/da;
	caf=pow(vaf,-1/ga);
	sf=w*z[0]+r*a[i]-caf;
	vab=(v[i]-v[i-1])/da;
	cab=pow(vab,-1/ga);
	ca0=w*z[0]+r*a[i];
	va0=pow(ca0,-ga);
	sb=w*z[0]+r*a[i]-cab;
	Ax[4*i-1]=0;
	Ax[4*i]=1/Delta+rho-centdiag[i];
	Ax[4*i+1]=0;
	Axt[4*i-1]=0;
	Axt[4*i]=centdiag[i];
	Axt[4*i+1]=0;
	if (sf>0) {
	  va_upwind+=vaf;
	  Ax[4*i]+=sf/da;
	  Ax[4*i+1]-=sf/da;
	  Axt[4*i]-=sf/da;
	  Axt[4*i+1]+=sf/da;
	}
	if (sb<0) {
	  va_upwind+=vab;
	  Ax[4*i-1]=sb/da;
	  Ax[4*i]-=sb/da;
	  Axt[4*i-1]=-sb/da;
	  Axt[4*i]+=sb/da;
	}
	if ((sf<=0)&&(sb>=0)) va_upwind+=va0;
	c=pow(va_upwind,-1/ga);
	u=pow(c,(1-ga))/(1-ga);
	b[i]=u+v[i]/Delta;
      } // End of i Iteration //
      va_upwind=0;
      vab=(v[I-1]-v[I-2])/da;
      cab=pow(vab,-1/ga);
      sb=w*z[0]+r*a[I-1]-cab;
      ca0=w*z[0]+r*a[I-1];
      va0=pow(ca0,-ga);
      if (sb<0) {
	Ax[4*(I-1)-1]=sb/da;
	Ax[4*(I-1)]=1/Delta+rho-centdiag[I-1]-sb/da;
	Axt[4*(I-1)-1]=-sb/da;
	Axt[4*(I-1)]=centdiag[I-1]+sb/da;
	va_upwind+=vab;
      }
      else {
	va_upwind+=va0;
	Ax[4*(I-1)-1]=0;
	Ax[4*(I-1)]=1/Delta+rho-centdiag[I-1];
	Axt[4*(I-1)-1]=0;
	Axt[4*(I-1)]=centdiag[I-1];
      }
      c=pow(va_upwind,-1/ga);
      u=pow(c,(1-ga))/(1-ga);
      b[I-1]=u+v[I-1]/Delta;

      for (int j=1;j<J-1;j++) {
	va_upwind=0;
	vaf=(v[j*I+1]-v[j*I])/da;
	caf=pow(vaf,-1/ga);
	sf=w*z[j]+r*a[0]-caf;
	ca0=w*z[j]+r*a[0];
	va0=pow(ca0,-ga);
	if (sf>0) {
	  va_upwind+=vaf;
	  Ax[j*(5*I-2)-I+1]=1/Delta+rho-centdiag[j*I]+sf/da;
	  Ax[j*(5*I-2)-I+2]=-sf/da;
	  Axt[j*(5*I-2)-I+1]=centdiag[j*I]-sf/da;
	  Axt[j*(5*I-2)-I+2]=sf/da;
	}
	else {
	  va_upwind+=va0;
	  Ax[j*(5*I-2)-I+1]=1/Delta+rho-centdiag[j*I];
	  Ax[j*(5*I-2)-I+2]=0;
	  Axt[j*(5*I-2)-I+1]=centdiag[j*I];
	  Axt[j*(5*I-2)-I+2]=0;
	}
	c=pow(va_upwind,-1/ga);
	u=pow(c,(1-ga))/(1-ga);
	b[j*I]=u+v[j*I]/Delta;
	for (int i=1;i<I-1;i++) {
	  va_upwind=0;
	  int tmp=mat(i,j);
	  vaf=(v[tmp+1]-v[tmp])/da;
	  caf=pow(vaf,-1/ga);
	  sf=w*z[j]+r*a[i]-caf;
	  vab=(v[tmp]-v[tmp-1])/da;
	  cab=pow(vab,-1/ga);
	  ca0=w*z[j]+r*a[i];
	  va0=pow(ca0,-ga);
	  sb=w*z[j]+r*a[i]-cab;
	  Ax[j*(5*I-2)-I+5*i]=0;
	  Ax[j*(5*I-2)-I+5*i+1]=1/Delta+rho-centdiag[j*I+i];
	  Ax[j*(5*I-2)-I+5*i+2]=0;
	  Axt[j*(5*I-2)-I+5*i]=0;
	  Axt[j*(5*I-2)-I+5*i+1]=centdiag[j*I+i];
	  Axt[j*(5*I-2)-I+5*i+2]=0;
	  if (sf>0) {
	    va_upwind+=vaf;
	    Ax[j*(5*I-2)-I+5*i+1]+=sf/da;
	    Ax[j*(5*I-2)-I+5*i+2]-=sf/da;
	    Axt[j*(5*I-2)-I+5*i+1]-=sf/da;
	    Axt[j*(5*I-2)-I+5*i+2]+=sf/da;
	  }
	  if (sb<0) {
	    va_upwind+=vab;
	    Ax[j*(5*I-2)-I+5*i]=sb/da;
	    Ax[j*(5*I-2)-I+5*i+1]-=sb/da;
	    Axt[j*(5*I-2)-I+5*i]=-sb/da;
	    Axt[j*(5*I-2)-I+5*i+1]+=sb/da;
	  }
	  if ((sf<=0)&&(sb>=0)) va_upwind+=va0;
	  c=pow(va_upwind,-1/ga);
	  u=pow(c,(1-ga))/(1-ga);
	  b[tmp] = u+v[tmp]/Delta;
	} // End of i Iteration //
	va_upwind=0;
	vab=(v[j*I+I-1]-v[j*I+I-2])/da;
	cab=pow(vab,-1/ga);
	sb=w*z[j]+r*a[I-1]-cab;
	ca0=w*z[j]+r*a[I-1];
	va0=pow(ca0,-ga);
	if (sb<0) {
	  Ax[j*(5*I-2)+4*I-5]=sb/da;
	  Ax[j*(5*I-2)+4*I-4]=1/Delta+rho-centdiag[j*I+I-1]-sb/da;
	  Axt[j*(5*I-2)+4*I-5]=-sb/da;
	  Axt[j*(5*I-2)+4*I-4]=centdiag[j*I+I-1]+sb/da;
	  va_upwind+=vab;
	}
	else {
	  va_upwind+=va0;
	  Ax[j*(5*I-2)+4*I-5]=0;
	  Ax[j*(5*I-2)+4*I-4]=1/Delta+rho-centdiag[j*I+I-1];
	  Axt[j*(5*I-2)+4*I-5]=0;
	  Axt[j*(5*I-2)+4*I-4]=centdiag[j*I+I-1];
	}
	c=pow(va_upwind,-1/ga);
	u=pow(c,(1-ga))/(1-ga);
	b[j*I+I-1] = u+v[j*I+I-1]/Delta;
      } // End of j Iteration //

      va_upwind=0;
      vaf=(v[(J-1)*I+1]-v[(J-1)*I])/da;
      caf=pow(vaf,-1/ga);
      sf=w*z[J-1]+r*a[0]-caf;
      ca0=w*z[J-1]+r*a[0];
      va0=pow(ca0,-ga);
      if (sf>0) {
	va_upwind+=vaf;
	Ax[(J-1)*(5*I-2)-I+1]=1/Delta+rho-centdiag[(J-1)*I]+sf/da;
	Ax[(J-1)*(5*I-2)-I+2]=-sf/da;
	Axt[(J-1)*(5*I-2)-I+1]=centdiag[(J-1)*I]-sf/da;
	Axt[(J-1)*(5*I-2)-I+2]=sf/da;
      }
      else {
	va_upwind+=va0;
	Ax[(J-1)*(5*I-2)-I+1]=1/Delta+rho-centdiag[(J-1)*I];
	Ax[(J-1)*(5*I-2)-I+2]=0;
	Axt[(J-1)*(5*I-2)-I+1]=centdiag[(J-1)*I];
	Axt[(J-1)*(5*I-2)-I+2]=0;
      }
      c=pow(va_upwind,-1/ga);
      u=pow(c,(1-ga))/(1-ga);
      b[(J-1)*I] = u + v[(J-1)*I]/Delta;
      for (int i=1;i<I-1;i++) {
	va_upwind=0;
	vaf=(v[(J-1)*I+i+1]-v[(J-1)*I+i])/da;
	caf=pow(vaf,-1/ga);
	sf=w*z[J-1]+r*a[i]-caf;
	vab=(v[(J-1)*I+i]-v[(J-1)*I+i-1])/da;
	cab=pow(vab,-1/ga);
	ca0=w*z[J-1]+r*a[i];
	va0=pow(ca0,-ga);
	sb=w*z[J-1]+r*a[i]-cab;
	Ax[(J-1)*(5*I-2)-I+4*i]=0;
	Ax[(J-1)*(5*I-2)-I+4*i+1]=1/Delta+rho-centdiag[(J-1)*I+i];
	Ax[(J-1)*(5*I-2)-I+4*i+2]=0;
	Axt[(J-1)*(5*I-2)-I+4*i]=0;
	Axt[(J-1)*(5*I-2)-I+4*i+1]=centdiag[(J-1)*I+i];
	Axt[(J-1)*(5*I-2)-I+4*i+2]=0;
	if (sf>0) {
	  va_upwind+=vaf;
	  Ax[(J-1)*(5*I-2)-I+4*i+1]+=sf/da;
	  Ax[(J-1)*(5*I-2)-I+4*i+2]-=sf/da;
	  Axt[(J-1)*(5*I-2)-I+4*i+1]-=sf/da;
	  Axt[(J-1)*(5*I-2)-I+4*i+2]+=sf/da;	
	}
	if (sb<0) {
	  va_upwind+=vab;
	  Ax[(J-1)*(5*I-2)-I+4*i]=sb/da;
	  Ax[(J-1)*(5*I-2)-I+4*i+1]-=sb/da;
	  Axt[(J-1)*(5*I-2)-I+4*i]=-sb/da;
	  Axt[(J-1)*(5*I-2)-I+4*i+1]+=sb/da;
	}
	if ((sf<=0)&&(sb>=0)) va_upwind+=va0;
	c=pow(va_upwind,-1/ga);
	u=pow(c,(1-ga))/(1-ga);
	b[(J-1)*I+i] = u+v[(J-1)*I+i]/Delta;
      } // End of i Iteration //
      va_upwind=0;
      vab=(v[I*J-1]-v[I*J-2])/da;
      cab=pow(vab,-1/ga);
      sb=w*z[J-1]+r*a[I-1]-cab;
      ca0=w*z[J-1]+r*a[I-1];
      va0=pow(ca0,-ga);
      if (sb<0) {
	Ax[(J-1)*(5*I-2)-I+4*(I-1)]=sb/da;
	Ax[(J-1)*(5*I-2)-I+4*(I-1)+1]=1/Delta+rho-centdiag[J*I-1]-sb/da;
	Axt[(J-1)*(5*I-2)-I+4*(I-1)]=-sb/da;
	Axt[(J-1)*(5*I-2)-I+4*(I-1)+1]=centdiag[J*I-1]+sb/da;
	va_upwind+=vab;
      }
      else {
	va_upwind+=va0;
	Ax[(J-1)*(5*I-2)-I+4*(I-1)]=0;
	Ax[(J-1)*(5*I-2)-I+4*(I-1)+1]=1/Delta+rho-centdiag[J*I-1];
	Axt[(J-1)*(5*I-2)-I+4*(I-1)]=0;
	Axt[(J-1)*(5*I-2)-I+4*(I-1)+1]=centdiag[J*I-1];
      }
      c=pow(va_upwind,-1/ga);
      u=pow(c,(1-ga))/(1-ga);
      b[J*I-1]=u+v[J*I-1]/Delta;

      (void) umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,null,Info);
      (void) umfpack_di_solve(UMFPACK_At,Ap,Ai,Ax,v,b,Numeric,null,null);
      umfpack_di_free_numeric(&Numeric);

      tmpmax=0;
      for (int i=0;i<n;i++){
	if ((v[i]-vold[i])>tmpmax) tmpmax=v[i]-vold[i];
	if ((vold[i]-v[i])>tmpmax) tmpmax=vold[i]-v[i];
      }

      if (tmpmax<crit) {
	printf("Converged %e\n",tmpmax);
	break;
      }
    } // End of HBJ Iteration //



    for (int i=1;i<n;i++) {
      b[i]=0;
    }
    b[0]=0.1;
    Axt[0]=1; Axt[3]=0; Axt[4*I-2]=0;
    (void) umfpack_di_numeric(Ap,Ai,Axt,Symbolic,&Numeric,null,Info);
    (void) umfpack_di_solve(UMFPACK_A,Ap,Ai,Axt,gg,b,Numeric,null,Info);
    umfpack_di_free_numeric(&Numeric);

    g_sum=0;
    S=0;
    for (int j=0;j<J;j++) {
      for (int i=0;i<I;i++) {
	S+=a[i]*gg[j*I+i];
	g_sum+=gg[j*I+i];
      }
    }

    S=S/g_sum;
    if (S>K) tmpmax=S-K;
    else tmpmax=K-S;
    printf("K is %e and S is %e\n",K, S);

    if (tmpmax <critK){
      printf("Converged K %e\n",K);
      break;
    }

    K=relax*K+(1-relax)*S;
    r=alpha*pow(K,alpha-1)-delta;
    w=(1-alpha)*pow(K,alpha);

  } // End of Capital Iteration //

  //printf("%0.9f\n",Var);
  //(void) umfpack_di_symbolic (n,n,Ap,Ai,Ax, &Symbolic, null, null);
  //(void) umfpack_di_numeric (Ap,Ai,Ax, Symbolic, &Numeric,  null, null);
  //(void) umfpack_di_solve (UMFPACK_A,Ap,Ai,Ax,x,b,Numeric, null, null);
  //umfpack_di_free_symbolic (&Symbolic);
  //umfpack_di_free_numeric (&Numeric);
  /*for (int i=0;i<I*J;i++) {
    printf("x[%d]=%e\n",i,v[i]);
    }*/
  return (0);
}
