
#include <math.h>
#include<stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define PI 3.14159265358979323846
#define fase 0 
#define wresistivo 3.8e-4

#define T_a 15.2329 
#define p_a sqrt(T_a) 
#define q_a p_a/2.0
#define nu_a p_a*sqrt(3.0/4.0)
#define BT0asdex 17.815757116271065


void B_Asdex(double r,double z,double B_eq[],double *s_flux);
void  magnetic_field(double *B, double r, double qq, double zc);
void find_Mag_axis(double *RM,double *ZM,double *flu_max);
main()
{
  int i,j,test;
  double   r,z,br,bq,bz,C,ep,fi,q0,dr,dz,r0,z0,hq,qs,q,rq,psi,A,B[3],sflux,RM,ZM,s_flux,flu_max;
 FILE *qper,*bzper;

 qper=fopen("perfilqAl2.dat","w");
 // bzper=fopen("bzp.dat","w");


	printf("busca el flujo maximo\n");
	find_Mag_axis(&RM,&ZM,&flu_max);

	printf("eje magnetico RM=%f ZM=%f flu_max = %lf \n",RM,ZM,flu_max);



 A=3;

 for(i=0;i<210;i++)
   {
     //     r0=RM-0.48+i*0.0047;
     r0=RM+0.01+i*0.002;
     rq=r0;
     z0=ZM;
     hq=1.570796327e-4;
     fi=0.0;
     q0=0.0;

     test=0;

     z=z0;
     q=q0;
     while(test<2)
       {
	 //	 magnetic_field(&B[0],r0, q0,z0);

	  B_Asdex(r0,z0, &B[0], &s_flux);

	 //	 campo(r0,z0,br,bq,bz,psi);
	 //	 dr=r0*br*hq/bq;
	 // dz=r0*bz*hq/bq;
	 dr=r0*B[0]*hq/B[1];
	 dz=r0*B[2]*hq/B[1];
	 r=r0+dr;
	 z=z0+dz;
	 q=q0+hq;

	 if(((z0-ZM)/(z-ZM)) < 0.0)
	   test++;
	 //  fprintf(qper,"%f %f %f\n",r,q,z);
	 z0=z;
	 r0=r;
	 q0=q;
       }

     printf("%d %lf %lf %lf\n",i,r,z,q0);
     qs=q0/(2*PI);
     fprintf(qper,"%lf %lf %lf \n",rq,qs,s_flux);

   }

 fclose(qper);
}
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
	//	Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r))/(BT0asdex);
	Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r));



	double jb0p=j0(p_a*r);
	double jb0q=j0(q_a*r);	
	double jb0nu=j0(nu_a*r);
	double yb0q=y0(q_a*r);	
	double yb0nu=y0(nu_a*r);		
	double Bz=0.0;
	//Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) )/(BT0asdex) ;

	Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) ) ;
	
	double Bt=0.0;	
	
	double u_a=-(cc1*T_a);
	double F0_a= 30.4;	
	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu)  ;
	//	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r)/(BT0asdex) ;
	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r) ;

	B_eq[0]=Br;
	B_eq[1]=Bt;
	B_eq[2]=Bz;
	double Psi_max=1.410312;
	//*s_flux= 1.0-(Psi/Psi_max);
	*s_flux=Psi;
}


void  magnetic_field(double *B, double r, double qq, double zc){
// Coincide con el de Ricardo.
	// Necesito coord. toroidales:
	
	
  double fi;
  double A=3.;
	double ep = 1.0/A;
	double R0 = A;
	//	double r = sqrt(xc*xc +yc*yc);
	// Usamos x para la coord toroidal.
	double xx = sqrt((r-R0)*(r-R0) + zc*zc);
	double cfi = (r-R0)/xx;
	double sfi = zc/xx;
	
	//	int mm=1,nn=1;
	
	
	// constantes:
	double ak = 2.4048255577;		// 1er cero de la J_0.
	double p1 = 0.05;			// Ref: R. Farengo, Plasma Phys. Control. Fusion 54 (2012) 025007.
	double I0 = 1.0;
	double alpha = 4.0*p1/(ep*ep);
	double bw = 0.08;					// Campo poloidal en x = 1, theta = 0. Ref: Idem p1.
	double i1 = sqrt(0.25*ep*ep*ak*ak-p1);	// Corriente orden 1.
	double  j1k=j1(ak);
	// calculo la constante c:
	//	double aux = -0.5*(ep/(1.0+ep))*k*j1(k)*(1.0 + 0.5*ep*(1.0+2.0*af/(k*k)));
	//	double C = Bw/aux;		//corroborado con el codigo de Ricardo.
	double  C=-2.0*(1.0+ep)*bw/(ep*j1k*ak*(1.0
					 +ep*(1.0+2.0*alpha/(ak*ak))/2.0));
	 
	 

	
	double bfi;
	
	double bmod2,dbfidx;







	//	   xx=dsqrt((r-A)**2+z**2) !"radial" toroidal coordinate
	fi=asin(zc/xx);             // !poloidal angle
     if(r < A) 
	  fi=3.14159265359-fi;


	double xd=1.0-xx*xx; 
        double j0x=j0(ak*xx);
        double j1x=j1(ak*xx);
        double j0p=-j1x;
        double j1p=j0x-j1x/(ak*xx);

	double den=2.0*(1.0+ep*xx*cfi);
        double dp=1.0+ep*ep*xx*xx;

        double psi0=C*j0x;

       


	double aux=ak*xx+alpha*xx/ak+alpha/(ak*xx);
	

	bfi=-0.5*ep*C*ak*j1x;
	dbfidx=0.5*ep*C*ak*ak*(j1x/(ak*xx)-j0x);
       
        

	
	double psi1=cfi*C*0.5*(xx*j0x+alpha*j1x*xd/ak);

         double  br=-ep*C*(-ak*j1x*sfi+0.5*ep*sfi*cfi*(j1x*(-ak*xx
							    -2.0*alpha/(ak*xx))+j0x*alpha*xd))/den;
           
         double  bz=ep*C*(-cfi*ak*j1x+0.5*ep*(j0x*(1.0+alpha*xd*cfi*cfi)
		  +j1x*(alpha*xd*sfi*sfi/(ak*xx)
			-cfi*cfi*(ak*xx+alpha*(xx+1.0/xx)/ak))))/den;

	 double psi=psi0+ep*psi1;
	 

	  double  bq=2.0*sqrt(1.0+i1*i1*psi*psi)/den;
	  

	  //	 bqp is used to calcualte the perturbed fields. Set bqp=1.0d0 for "cylindrical case"

          B[0]=br;
          B[2]=bz;
          B[1]=bq;

	 
}



void find_Mag_axis(double *RM,double *ZM,double *flux_max)
{
  int i,j;
 
  double s_flux=0.0,s_flux_max=0,ri,zi,B_equilibrio[3];
  //  FILE *out;
 
  // out=fopen("flux.dat","w");

 double z0 = -1;
 double r0= 1;
 double z1 = 1;
 double r1= 2.2;

double N_r=220;
double N_z=180;
 double dr=(r1-r0)/(N_r-1.);
 double dz=(z1-z0)/(N_z-1);


 *RM=1;
 *ZM=1;
  

  for(i=0;i<N_z;i++)
    {
      for(j=0;j<N_r;j++)
	{
	  ri=r0+j*dr;
	  zi=z0+i*dz;
	  
	  B_Asdex(ri,zi, B_equilibrio, &s_flux);
	  //	  fprintf(out,"%f %f %f \n",ri,zi,s_flux);
	  if(s_flux>s_flux_max)
	    {
	      s_flux_max=s_flux;
	      *RM=ri;
	      *ZM=zi;
	    }
	}
    }

  *flux_max=s_flux_max;
  
  // fclose(out);
}
