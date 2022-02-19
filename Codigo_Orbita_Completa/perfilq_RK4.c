//#include "../CamposEquilibrio/CamposMagneticos.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define PI 3.1415926535897932384626

#define n 1000

#define Nx 1001
#define amplitud 0.02
void B_Asdex(double r,double z,double B_eq[],double *s_flux);

int main (void){

  double Br=0,Bz=0,Bt=0;
  double R0=1.72;
  double   r,z,br,bq,bz,C,ep,fi,q0,dr,dz,r0,z0,hq,qs,q,rq,psi,A,B[3],sflux,RM,ZM,s_flux,flu_max;

  double dphi=0.001;
  double hr[4]={0,0,0,0},hz[4]={0,0,0,0};
  double ri=0,zi=0,theta=0,rho=0;
  double an[4]={0.0,0.5,0.5,1.0};
  double aux;

  //numero de pasos
  double N=2.0*PI/dphi;
  //contadores
  int i=0,j=0,k=0,m=0;


  double dx=0.002,B1r=0.0,B1z=0.0;

  int l=0,test;

	FILE *out;

	out=fopen("per_q2.dat","w");//pertur amplitud=0.006, cambie jper,theta*




 RM=1.717808;
 ZM=0.128492;
 flu_max = 1.410312;
 //	printf("busca el flujo maximo\n");
 //	find_Mag_axis(&RM,&ZM,&flu_max);

//empiezo la integracion
 for(i=0;i<210;i++){

   r0=RM+0.01+i*0.002;
   rq=r0;
   z0=ZM;
   hq=1.570796327e-4;
   dphi=1.570796327e-4;
   fi=0.0;
   q0=0.0;
   test=0;

   z=z0;
   q=q0;


   //empiezo runge kutta
   while(test<2)
     {

       //un paso RK4
       for(j=0; j<4; j++)
	 {
	   ri=r0 + an[j]*hr[j-1];
	   zi=z0 + an[j]*hz[j-1];

	   B_Asdex(ri,zi, &B[0], &s_flux);
	   Br=B[0];
	   Bz=B[2];
	   Bt=B[1];

	   hr[j]=((ri*Br)/Bt)*dphi;
	   hz[j]=((ri*Bz)/Bt)*dphi;
	 }
       r= r0 + (hr[0]+2.0*hr[1]+2.0*hr[2]+hr[3])/6.0;
       z= z0 + (hz[0]+2.0*hz[1]+2.0*hz[2]+hz[3])/6.0;
       q=q+dphi;

       if(((z0-ZM)/(z-ZM)) < 0.0)
	 test++;
       // printf("%d %f %f %f\n",i,r,q,z);
       z0=z;
       r0=r;
       q0=q;

     }
   //printf("%d %lf %lf %lf\n",i,r,z,q0);
   qs=q0/(2*PI);
   fprintf(out,"%e %e %e \n",rq,qs,s_flux);
 }

fclose(out);

return 0;
}


#define T_a 15.2329
#define p_a sqrt(T_a)
#define q_a p_a/2.0
#define nu_a p_a*sqrt(3.0/4.0)
#define BT0asdex 17.815757116271065
//Equilibrio Asdex
void B_Asdex(double r,double z,double B_eq[],double *s_flux){
double cc1=0.4733, cc2=-0.2164,cc3=0.0, cc4=0.0, cc5=0.0, cc6=0.0,cc7=-0.06830, cc8=0.01220, cc9=0.1687;
double cc10=0.8635, cc11=-1.0682, cc12=0.02166,cc13=-0.002662, cc14=0.1178, cc15=1.4008, cc16=-0.2656,cc17=1.3770, cc18=0.2468;


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
	Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r))/(BT0asdex);

	double jb0p=j0(p_a*r);
	double jb0q=j0(q_a*r);
	double jb0nu=j0(nu_a*r);
	double yb0q=y0(q_a*r);
	double yb0nu=y0(nu_a*r);
	double Bz=0.0;
	Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) )/(BT0asdex) ;

	double Bt=0.0;

	double u_a=-(cc1*T_a);
	double F0_a= 30.4;
	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu)  ;
	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r)/(BT0asdex) ;
	B_eq[0]=Br;
	B_eq[1]=Bt;
	B_eq[2]=Bz;
	double Psi_max=1.0;
	//	*s_flux= 1.0-(Psi/Psi_max);

	*s_flux=Psi;

}
