//#include "../CamposEquilibrio/CamposMagneticos.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>


#define PI 3.1415926535897932384626

#define n 1000
#define Rn sqrt(n*n+1)
#define Nx 1001
#define amplitud 0.0003

void B_Asdex(double rp,double zp,double *B,double *s_flux);

int main (void){

  //dimensiones del mapa pertubado
 double z0 = -0.9*2;
 double r0= 1.1*2;
 double z1 = 0.9*2;
 double r1= 2.2*2;
 //

double dr;
double dz;
int N_r=401;
int N_z=401;

 dr=(r1-r0)/(N_r-1.);
 dz=(z1-z0)/(N_z-1.);


 double B[3],s_flux;

	FILE *out;

	out=fopen("./Outputs/Perturbado/poincareBB33k12_a400.dat","w");//pertur amplitud=0.01, cambie jper,theta*,



	  FILE *in= fopen("CamposMagneticoPerturbado_A6_n401.txt", "r");
  	assert(in!=NULL);

//CAMPOS MAGNETICOS
	double Br=0,Bz=0,Bt=0;
	//double R0=1.72;
	double r=3.2,z=0,phii=0;
	double dphi=0.001;
	double hr[4]={0,0,0,0},hz[4]={0,0,0,0};
	double ri=0,zi=0;//,theta=0,rho=0;
	double an[4]={0.0,0.5,0.5,1.0};
	double A=0.2;
	double aux;
  //double omega = -0.0008;
//Para la interpolacion
double pi=0.0,qi=0.0;
//numero de pasos
double N=2.0*PI/dphi;
//contadores
int i=0,j=0,k=0,m=0;
 double fff;
 // ojo que cambie el orden de los B
 //double B1r1[N_z][N_r],B1z1[N_z][N_r],B1r2[N_z][N_r],B1z2[N_z][N_r];
 //double B1r1[N_r][N_z],B1z1[N_r][N_z],B1r2[N_r][N_z],B1z2[N_r][N_z];

int ii=0,jj=0;double b1ra[N_r*N_z],b1za[N_r*N_z],b1rb[N_r*N_z],b1zb[N_r*N_z];
double B1r=0.0,B1z=0.0; //dx=0.002,

//int l=0;
/* double **B1r1 = (double **)malloc(N_z*sizeof(double*)); */
/* double **B1z1 = (double **)malloc(N_z * sizeof(double*)); */
/* double **B1r2 = (double **)malloc(N_z * sizeof(double*)); */
/* double **B1z2 = (double **)malloc(N_z * sizeof(double*)); */


/* for(l = 0; l < N_z; l++){ 	 */
/* B1r1[l] = (double *)malloc(N_r * sizeof(double)); */
/* B1z1[l] = (double *)malloc(N_r * sizeof(double)); */
/* B1r2[l] = (double *)malloc(N_r * sizeof(double)); */
/* B1z2[l] = (double *)malloc(N_r * sizeof(double)); */
/* } */



/* for(j=0;j<N_z;j++){ */
/* 	for(i=0;i<(N_r);i++){ */

	  //			 fscanf(in,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&aux,&aux,&aux,&B1r1[i][j],&B1r2[i][j],&B1z1[i][j],&B1z2[i][j]);
	  //	  fscanf(in,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&aux,&aux,&B1r1[i][j],&B1r2[i][j],&B1z1[i][j],&B1z2[i][j]);
	for(i=0;i<N_z;i++){
	  for(j=0;j<N_r;j++){
	    fscanf(in,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&aux,&aux,&b1ra[i*N_r+j],&b1rb[i*N_r+j],&b1za[i*N_r+j],&b1zb[i*N_r+j]);
		//printf("%e \n",Br1[i][j]);

	}
}

printf("paso\n");

// exit(0);

//empiezo la integracion
for(m=0;m<21;m++){
  //for(m=0;m<18;m++){
  A=A+ 0.02;
  r= 3.46 + m*0.04;
  z= 0.246 ;

  //r=1.85+0.01*m;
  //z=0.16;

  	for(k=0;k<1000;k++){ //Original k<6000

  //	for(k=0;k<500;k++){
	  for(i=0; i<=N; i++ ){
			//empiezo runge kutta
			// phii=dphi*i;
      phii=dphi*i - 3.0*PI/2.0;
			for(j=0; j<4; j++){
				ri=r + an[j]*hr[j-1];
				zi=z + an[j]*hz[j-1];



				B_Asdex(ri,zi, &B[0], &s_flux);




				//	theta=Theta(ri,zi,R0);
				//	rho=Rho(ri,zi,R0);

				if(ri<r0 || ri>r1 || zi <z0 || zi>z1)
				  {
				    ii=0;
				    jj=0;
				    qi=0;
				    pi=0;
				    fff=0.0;
				  }
				else
				  {
				    jj=(ri-r0)/dr;
				    ii=(zi-z0)/dz;
				    pi=(ri-r0-jj*dr)/dr;
				    qi=(zi-z0-ii*dz)/dz;
				    fff=amplitud;
				  }

      double up = 1.0 - pi;
      double uq = 1.0 - qi;
      int nr = 401;
      int kk = jj+ii*nr;

      double cphi = cos(phii), sphi = sin(-phii);

				/* Br=B[0];				 */
				/* B1r=    (1.0-p)*(1.0-q)*( B1r1[ii][jj]*cos(phi)+B1r2[ii][jj]*sin(phi) )  + p*(1.0-q)*( B1r1[ii+1][jj]*cos(phi)+B1r2[ii+1][jj]*sin(phi) ) +  */
				/*  	 q*(1.0-p)*( B1r1[ii][jj+1]*cos(phi)+B1r2[ii][jj+1]*sin(phi) ) + p*q*( B1r1[ii+1][jj+1]*cos(phi)+B1r2[ii+1][jj+1]*sin(phi) ) ;  */



				/* Bz=B[2]; */
				/* 	B1z= (1.0-p)*(1.0-q)*( B1z1[ii][jj]*cos(phi)+B1z2[ii][jj]*sin(phi) )  + p*(1.0-q)*( B1z1[ii+1][jj]*cos(phi)+B1z2[ii+1][jj]*sin(phi) ) + */
				/* 	 q*(1.0-p)*( B1z1[ii][jj+1]*cos(phi)+B1z2[ii][jj+1]*sin(phi) ) + p*q*( B1z1[ii+1][jj+1]*cos(phi)+B1z2[ii+1][jj+1]*sin(phi) ) ; */


     B1r = up*uq*(b1ra[kk]*cphi+b1rb[kk]*sphi) + pi*uq*(b1ra[kk+1]*cphi+b1rb[kk+1]*sphi) + qi*up*(b1ra[kk+nr]*cphi+b1rb[kk+nr]*sphi) + pi*qi*(b1ra[kk+nr+1]*cphi+b1rb[kk+nr+1]*sphi);

     B1z = up*uq*(b1za[kk]*cphi+b1zb[kk]*sphi) + pi*uq*(b1za[kk+1]*cphi+b1zb[kk+1]*sphi) + qi*up*(b1za[kk+nr]*cphi+b1zb[kk+nr]*sphi) + pi*qi*(b1za[kk+nr+1]*cphi+b1zb[kk+nr+1]*sphi);
  // B1r=    (1.0-p)*(1.0-q)*( B1r1[ii*N_r + jj]*cos(phii)+B1r2[ii*N_r + jj]*sin(-phii) )  + q*(1.0-p)*( B1r1[(ii+1)*N_r+ jj]*cos(phii)+B1r2[(ii+1)*N_r+jj]*sin(-phii) ) +
  //   p*(1.0-q)*( B1r1[ii*N_r + jj+1]*cos(phii)+B1r2[ii*N_r + jj+1]*sin(-phii) ) + p*q*( B1r1[(ii+1)*N_r + jj+1]*cos(phii)+B1r2[(ii+1)*N_r + jj+1]*sin(-phii) ) ;
  //
  // B1z= (1.0-p)*(1.0-q)*( B1z1[ii*N_r + jj]*cos(phii)+B1z2[ii*N_r + jj]*sin(-phii) )  + q*(1.0-p)*( B1z1[(ii+1)*N_r+ jj]*cos(phii)+B1z2[(ii+1)*N_r+ jj]*sin(-phii) ) +
  //   p*(1.0-q)*( B1z1[ii*N_r + jj+1]*cos(phii)+B1z2[ii*N_r + jj+1]*sin(-phii) ) + p*q*( B1z1[(ii+1)*N_r + jj+1]*cos(phii)+B1z2[(ii+1)*N_r + jj+1]*sin(-phii) ) ;


				//if(fabs(ii)>=Nx || fabs(jj)>=Nx || ii==0 || jj==0 ){ B1r=0.0;B1z=0.0; }
				/*  Br=Br+ B1r*fff; */

				/* Bz=Bz + B1z*fff; */

  Br=B[0]+ B1r*fff;

  Bz=B[2] + B1z*fff;
  Bt=B[1];


				hr[j]=(	(	ri*Br	)/Bt	)*dphi;
				hz[j]=( (	ri*Bz	)/Bt	)*dphi;
			 }
			r= r + (hr[0]+2.0*hr[1]+2.0*hr[2]+hr[3])/6.0;
			z= z + (hz[0]+2.0*hz[1]+2.0*hz[2]+hz[3])/6.0;
	  }
	  fprintf(out,"%e %e %e %e\n",r,z,phii,s_flux);
	}
	fprintf(out,"\n");
 }


/* for( i=0; i<N_z; i++ ) { */
/* 	free( B1r1[i] ); */
/* 	free( B1z1[i] ); */
/* 	free( B1r2[i] ); */
/* 	free( B1z2[i] ); */


/* } */
/* free( B1r1 ); */
/* free( B1z1 ); */
/* free( B1r2 ); */
/* free( B1z2 ); */

	fclose(out);
	fclose(in);

	return 0;


	}

const double BT0 = 17.815757116271065;

void B_Asdex(double rp,double zp,double *B,double *s_flux) {
  double T_a=15.2329;
  double p_a=sqrt(T_a);
  double q_a=p_a/2.0;
  double nu_a=p_a*sqrt(3.0/4.0);

  double cc1=0.4733, cc2=-0.2164,cc3=0.0, cc4=0.0, cc5=0.0, cc6=0.0,cc7=-0.06830, cc8=0.01220, cc9=0.1687;
  double cc10=0.8635, cc11=-1.0682, cc12=0.02166,cc13=-0.002662, cc14=0.1178, cc15=1.4008, cc16=-0.2656,cc17=1.3770, cc18=0.2468;

	double r=rp*0.5;
	double z=zp*0.5;

	double csp=cos(p_a*z);
	double snp=sin(p_a*z);
	double csq=cos(q_a*z);
	double snq=sin(q_a*z);
	double csnu=cos(nu_a*z);
	double snnu=sin(nu_a*z);
	double jb1p=j1(p_a*r);
	double jb1q=j1(q_a*r);
	double jb1nu=j1(nu_a*r);
	double yb1q=y1(q_a*r);
	double yb1nu=y1(nu_a*r);
	double rho=sqrt(r*r+z*z);
	double Br=0.0;

  Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r))/(BT0);


	double jb0p=j0(p_a*r);
	double jb0q=j0(q_a*r);
	double jb0nu=j0(nu_a*r);
	double yb0q=y0(q_a*r);
	double yb0nu=y0(nu_a*r);
	double Bz=0.0;

  Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) )/(BT0);


	double Bt=0.0;
	double u_a=-(cc1*T_a);
	double F0_a= 30.4;
	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu);
	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r)/(BT0) ;


  B[0]=Br;
	B[1]=Bt;
	B[2]=Bz;
	s_flux[0]=Psi;

  if (rp >= -0.1 && rp < -0.099396) {
    printf("yb1q = %f c/ r = %f.\n", yb1q, rp);
  }

  return;
}

//
// #define T_a 15.2329
// #define p_a sqrt(T_a)
// #define q_a p_a/2.0
// #define nu_a p_a*sqrt(3.0/4.0)
// #define BT0asdex 17.815757116271065
// //Equilibrio Asdex
// void B_Asdex(double r,double z,double B_eq[],double *s_flux){
// double cc1=0.4733, cc2=-0.2164,cc3=0.0, cc4=0.0, cc5=0.0, cc6=0.0,cc7=-0.06830, cc8=0.01220, cc9=0.1687;
// double cc10=0.8635, cc11=-1.0682, cc12=0.02166,cc13=-0.002662, cc14=0.1178, cc15=1.4008, cc16=-0.2656,cc17=1.3770, cc18=0.2468;
//
//
//  r=r*0.5;
//  z=z*0.5;
//
// 	double csp=cos(p_a*z);
// 	double snp=sin(p_a*z);
// 	double csq=cos(q_a*z);
// 	double snq=sin(q_a*z);
// 	double csnu=cos(nu_a*z);
// 	double snnu=sin(nu_a*z);
// 	double jb1p=j1(p_a*r);
// 	double jb1q=j1(q_a*r);
// 	double jb1nu=j1(nu_a*r);
// 	double yb1q=y1(q_a*r);
// 	double yb1nu=y1(nu_a*r);
// 	double rho=sqrt(r*r+z*z);
// 	double Br=0.0;
// 	Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r))/(BT0asdex);
//
// 	double jb0p=j0(p_a*r);
// 	double jb0q=j0(q_a*r);
// 	double jb0nu=j0(nu_a*r);
// 	double yb0q=y0(q_a*r);
// 	double yb0nu=y0(nu_a*r);
// 	double Bz=0.0;
// 	Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) )/(BT0asdex) ;
//
// 	double Bt=0.0;
//
// 	double u_a=-(cc1*T_a);
// 	double F0_a= 30.4;
// 	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu)  ;
// 	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r)/(BT0asdex) ;
// 	B_eq[0]=Br;
// 	B_eq[1]=Bt;
// 	B_eq[2]=Bz;
// 	double Psi_max=1.0;
// 	//	*s_flux= 1.0-(Psi/Psi_max);
//
// 	*s_flux=Psi;
//
// }
