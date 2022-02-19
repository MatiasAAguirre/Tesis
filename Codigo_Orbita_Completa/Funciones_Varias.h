#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cuentaD(double vpmod, double vmod2, double Bmod);
void B_Asdex(double r,double z,double *B,double *s_flux);
void magnetic_field(double *B, double *E, double r, double q, double z, double *s_flux);
double centro_giro2(double *r, double *rg, double *v, double *mu, double t0, double t);
void cond_i (double *t, double *r, double *v, double *rgc, double *vpar, double *mu);
void perturbacion(double *r, double t0, double t, double *b1, double *e1);
// double RHS_cil(double time, double *u, double *F);
// double RK4(double t, double *ri, double *vi, double *rs, double *vs);
double RHS_cil(int j, double t0, double time, double *u, double *F, double *Bmod2);
double RK46_NL(double t0, double t, double *ri, double *vi, double *rs, double *vs, double *Bmod2);
int integrador(double *t, double *r, double *v, double *rgc, double *vpar, double *mu, int *d);
void PROC(double t0, int nstepf, double *r, double *v, double *rgc, double *vpar, double *mu, int *d, double *rp, double *vp, double *rgcp, double *vparp, double *mup, double *tp, int *dp);

void cuentaD(double vpmod, double vmod2, double Bmod) {
  double m, q, ta;
  double Omega, a, B0, Rl;
  double tsim;

  B0 = 2.5; // Tesla
  // B0 = 2; // Tesla para P40 y 77.
  //R0 = 1.71, // radio mayor en metros
  //E = 70.e3*1.602e-19;    // 8keV Energy
  m = 2.0*1.6726e-27;           // D
  q = 1.0*1.6022e-19;          // D
  Omega = q*B0/m;       // Cyclotronic frequency (en rad/seg)
  ta = 1.0439E-8*2/B0;     // Período de Ciclotrón en segundos/rad, la inversa de Omega m/q*B0
  // ta = 1.0/((Omega*180)/3.14159); //Período de Ciclotrón en seg.
  //frec = 5000; //frecuencia del modo 5kHz
  a = 0.5;              // Minor radius
  Rl = a*gam*vpmod/Bmod; //Radio de Larmor en metros.
  tsim = nstep*dt*ta;

  printf("------Datos Varios------\n");
  printf("Gamma = %f\n", gam);
  printf("Omega = %.2e rad/s [%.2e Hz]\n", Omega, (Omega*180)/3.14159);
  printf("Radio de Larmor = %0.1f cm\n", Rl*100);
  printf("Tiempo de simulacion  %.2f ms \n",tsim*1000);
  printf("Cuantas vueltas da en un tiempo de simulacion = %.0f\n", tsim/ta);
  printf("Cuanto tarda en dar una vuelta: %.1e ms \n", ta*1000);
  printf("------------------------\n");

  return;
}

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

void magnetic_field(double *B, double *E, double r, double q, double z, double *s_flux) {
  //double *s_flux=0; //¡¡VER!! Siento que acá no importa meter psi como parametro ¡¡VER!!

   B_Asdex(r,z, B, s_flux);

   E[0]=0.0;
   E[1]=0.0;
   E[2]=0.0;

   return;
}

double centro_giro2(double *r, double *rg, double *v, double *mu, double t0, double t) {
  int i;
	double vpar, vpmod, rho, Bmod, vmod2;
	double vp[3];		// Velocidad perpendicular.
	double e[3];			// Vector unitario perpendicular a v y B.
	double B[3], E[3], dB[3], dE[3];
  double psi[1];



  for (i=0; i<3; i++) {
    dB[i] = 0.0;
    dE[i] = 0.0;
  }

	rg[0] = r[0];
	rg[1] = r[1];
	rg[2] = r[2];

  psi[0] = 0.0;
	magnetic_field(B, E, r[0], r[1], r[2], psi);

  perturbacion(r, t0, t, dB, dE);

  for (i=0; i<3; i++) {
    B[i] += dB[i];
  }

  Bmod = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);

  for(i=0;i<3;i++) {
		vp[i] = v[i]*B[i]/Bmod;
  }

  vpar = vp[0]+vp[1]+vp[2]; // Vel paralela.

  // for (i=0; i<3; i++) {
  //   vpar += v[i]*B[i];
  // }
  //
  // for (i=0; i<3; i++) {
  //   vp[i] = vpar*B[i]/Bmod;
  // }

	for(i=0;i<3;i++) {
		vp[i] = v[i] - vp[i];			// Vel perpendicular
  }

  vpmod = sqrt(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
  rho = gam*vpmod/Bmod;

	e[0] = v[1]*B[2] - v[2]*B[1];
	e[1] = v[2]*B[0] - v[0]*B[2];
	e[2] = v[0]*B[1] - v[1]*B[0];

  if (rho < 1.0E-8) {
    return 0; //Devuelve la posición de la partícula.
  }

	for(i=0;i<3;i++) {
    e[i] = e[i]/(vpmod*Bmod);
  }

	// Posición del centro de giro:
	for(i=0;i<3;i++) {
    rg[i] = r[i] + rho*e[i]; //Está con un + porque es v x B en vez de B x v.
  }

  if (t == t0) {
    vmod2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    cuentaD(vpmod, vmod2, Bmod);
    mu[0] = (gam*vpmod*vpmod)/(2.0*Bmod);
  }
  return vpar;
}

void cond_i(double *t, double *r, double *v, double *rgc, double *vpar, double *mu) {
  int i;
  double rr0, rq0, rz0, vr0, vq0, vz0, t0;
	double v0[3], r0[3], rg0[3];

	/* ------ Particula Pasante ------ */
  t0 = 0.0;

  rr0 = 4.022350;
	rq0 = 196.003;
	rz0 = 0.623175;

	vr0 = 0.075035;
	vq0 = 0.64219;
	vz0 = 0.762512;

	/* ------ Particula Atrapada ------ */
  // t0 = 0.0;
  //
  // rr0 = 3.41970;
	// rq0 = 21.2890;
	// rz0 = -0.781325;
  //
	// vr0 = 0.0181466;
	// vq0 = 0.0394052;
	// vz0 = 0.998925;

  /* ------ Particula 167 ------ */
  // t0 = 0.0;
  //
  // rr0 = 4.015180;
	// rq0 = 15.127500;
	// rz0 = 0.811112;
  //
	// vr0 = 0.816373;
	// vq0 = -0.053964;
	// vz0 = 0.574868;

  /* ------ Particula 64 ------ */
  // t0 = 0.0;
  //
  // rr0 = 3.642350;
	// rq0 = 20.581100;
	// rz0 = 1.178040;
  //
	// vr0 = 0.089354;
	// vq0 = 0.325686;
	// vz0 = 0.941076;

  /* ------ Particula 77 ------ */
  // t0 = 0.0;
  //
  // rr0 = 3.70682;
	// rq0 = 20.0959;
	// rz0 = 0.914920;
  //
	// vr0 = 0.918946;
	// vq0 = 0.0723818;
	// vz0 = 0.387400;

  /* ------ Particula 40 ------ */
  // t0 = 0.0;
  //
  // rr0 = 4.04292;
	// rq0 = 14.5101;
	// rz0 = 0.630314;
  //
	// vr0 = 0.631254;
	// vq0 = 0.254462;
	// vz0 = -0.732475;

  t[0] = t0;

	r0[0] = rr0;
	r0[1] = rq0;
	r0[2] = rz0;

	v0[0] = vr0;
	v0[1] = vq0;
	v0[2] = vz0;

	vpar[0] = centro_giro2(r0, rg0, v0, mu, t0, t0);


  for (i=0; i<3; i++) {
    r[i] = r0[i];
    v[i] = v0[i];
    rgc[i] = rg0[i];
  }

  return;
}

void perturbacion(double *r, double t0, double t, double *b1, double *e1) {
  int i, j, ii, jj, kk, ccb=0, cce=0, cib=0, cie=0;
  static double *b1ra, *b1rb, *b1za, *b1zb;
  static double *e1ra, *e1rb, *e1qa, *e1qb, *e1za, *e1zb;
  double rp = r[0]/2.0, zp = r[2]/2.0, pi, qi, up, uq, aux;

  if (t == t0) {
    FILE *Pert_Campo_Mag = fopen("./CamposMagneticoPerturbado_A6_n401.txt", "r");
    b1ra = (double *)malloc(nr*nr*sizeof(double));
    b1rb = (double *)malloc(nr*nr*sizeof(double));
    b1za = (double *)malloc(nr*nr*sizeof(double));
    b1zb = (double *)malloc(nr*nr*sizeof(double));

    FILE *Pert_Campo_Elec = fopen("./CamposElectricos4.txt", "r");
    e1ra = (double *)malloc(nr*nr*sizeof(double));
    e1rb = (double *)malloc(nr*nr*sizeof(double));
    e1za = (double *)malloc(nr*nr*sizeof(double));
    e1zb = (double *)malloc(nr*nr*sizeof(double));
    e1qa = (double *)malloc(nr*nr*sizeof(double));
    e1qb = (double *)malloc(nr*nr*sizeof(double));

    for (i=0; i<nr; i++) {
      for (j=0; j<nr; j++) {
        if(fscanf(Pert_Campo_Mag, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &aux, &aux, &b1ra[i*nr+j], &b1rb[i*nr+j], &b1za[i*nr+j], &b1zb[i*nr+j]) == 6) {
          ccb++;
        }
        else{
          cib++;
          printf("No se leyó bien el archivo de perturbación %d veces para la posición %d.\n", cib, i*nr+j);
          break;
        }
      }
    }

    for (i=0; i<nr; i++) {
      for (j=0; j<nr; j++) {
        if(fscanf(Pert_Campo_Elec, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &aux, &aux, &e1ra[i*nr+j], &e1rb[i*nr+j], &e1qa[i*nr+j], &e1qb[i*nr+j], &e1za[i*nr+j], &e1zb[i*nr+j], &aux, &aux) == 10) {
          cce++;
        }
        else{
          cie++;
          printf("No se leyó bien el archivo de perturbación %d veces para la posición %d.\n", cie, i*nr+j);
          break;
        }
      }
    }

    printf("\n");
    printf("Se leyeron %d lineas correctamente del archivo de perturbación Magnética y %d de perturbacion Eléctrica.\n", ccb, cce);
    printf("\n");

    for(i=0; i<160801; i++) {
      b1ra[i] = bscale*b1ra[i];
      b1rb[i] = bscale*b1rb[i];
      b1za[i] = bscale*b1za[i];
      b1zb[i] = bscale*b1zb[i];
      e1ra[i] = escale*e1ra[i];
      e1rb[i] = escale*e1rb[i];
      e1qa[i] = escale*e1qa[i];
      e1qb[i] = escale*e1qb[i];
      e1za[i] = escale*e1za[i];
      e1zb[i] = escale*e1zb[i];
    }

    fclose(Pert_Campo_Mag);
    fclose(Pert_Campo_Elec);

    if(rp<rp0 || rp>2*rp0 || zp<zp0 || zp>-zp0) {
        ii = 0;
        jj = 0;
        qi = 0.0;
        pi = 0.0;
      }

      else {
        ii = (int)((rp-rp0)/hr);
        jj = (int)((zp-zp0)/hz);
        pi = rp/hr - (float)ii - rp0/hr;
        qi = zp/hz - (float)jj - zp0/hz;
      }

    double phi = nmode * r[1] - omega * t;
    double cphi = cos(phi), sphi = sin(phi);

    kk = ii+jj*nr;

    up = 1.0 - pi;
    uq = 1.0 - qi;

    b1[0] = up*uq*(b1ra[kk]*cphi+b1rb[kk]*sphi) + pi*uq*(b1ra[kk+1]*cphi+b1rb[kk+1]*sphi) + qi*up*(b1ra[kk+nr]*cphi+b1rb[kk+nr]*sphi) + pi*qi*(b1ra[kk+nr+1]*cphi+b1rb[kk+nr+1]*sphi);
    b1[1] = 0.0;
    b1[2] = up*uq*(b1za[kk]*cphi+b1zb[kk]*sphi) + pi*uq*(b1za[kk+1]*cphi+b1zb[kk+1]*sphi) + qi*up*(b1za[kk+nr]*cphi+b1zb[kk+nr]*sphi) + pi*qi*(b1za[kk+nr+1]*cphi+b1zb[kk+nr+1]*sphi);

    e1[0] = up*uq*(e1ra[kk]*cphi+e1rb[kk]*sphi) + pi*uq*(e1ra[kk+1]*cphi+e1rb[kk+1]*sphi) + qi*up*(e1ra[kk+nr]*cphi+e1rb[kk+nr]*sphi) + pi*qi*(e1ra[kk+nr+1]*cphi+e1rb[kk+nr+1]*sphi);
    e1[1] = up*uq*(e1qa[kk]*cphi+e1qb[kk]*sphi) + pi*uq*(e1qa[kk+1]*cphi+e1qb[kk+1]*sphi) + qi*up*(e1qa[kk+nr]*cphi+e1qb[kk+nr]*sphi) + pi*qi*(e1qa[kk+nr+1]*cphi+e1qb[kk+nr+1]*sphi);
    e1[2] = up*uq*(e1za[kk]*cphi+e1zb[kk]*sphi) + pi*uq*(e1za[kk+1]*cphi+e1zb[kk+1]*sphi) + qi*up*(e1za[kk+nr]*cphi+e1zb[kk+nr]*sphi) + pi*qi*(e1za[kk+nr+1]*cphi+e1zb[kk+nr+1]*sphi);
  }

  else {

    if(rp<rp0 || rp>2*rp0 || zp<zp0 || zp>-zp0) {
        ii = 0;
        jj = 0;
        qi = 0.0;
        pi = 0.0;
      }

      else {
        ii = (int)((rp-rp0)/hr);
        jj = (int)((zp-zp0)/hz);
        pi = rp/hr - (float)ii - rp0/hr;
        qi = zp/hz - (float)jj - zp0/hz;
      }

    double phi = nmode * r[1] - omega * t;
    double cphi = cos(phi), sphi = sin(phi);

    kk = ii+jj*nr;

    up = 1.0 - pi;
    uq = 1.0 - qi;

    b1[0] = up*uq*(b1ra[kk]*cphi+b1rb[kk]*sphi) + pi*uq*(b1ra[kk+1]*cphi+b1rb[kk+1]*sphi) + qi*up*(b1ra[kk+nr]*cphi+b1rb[kk+nr]*sphi) + pi*qi*(b1ra[kk+nr+1]*cphi+b1rb[kk+nr+1]*sphi);
    b1[1] = 0.0;
    b1[2] = up*uq*(b1za[kk]*cphi+b1zb[kk]*sphi) + pi*uq*(b1za[kk+1]*cphi+b1zb[kk+1]*sphi) + qi*up*(b1za[kk+nr]*cphi+b1zb[kk+nr]*sphi) + pi*qi*(b1za[kk+nr+1]*cphi+b1zb[kk+nr+1]*sphi);

    e1[0] = up*uq*(e1ra[kk]*cphi+e1rb[kk]*sphi) + pi*uq*(e1ra[kk+1]*cphi+e1rb[kk+1]*sphi) + qi*up*(e1ra[kk+nr]*cphi+e1rb[kk+nr]*sphi) + pi*qi*(e1ra[kk+nr+1]*cphi+e1rb[kk+nr+1]*sphi);
    e1[1] = up*uq*(e1qa[kk]*cphi+e1qb[kk]*sphi) + pi*uq*(e1qa[kk+1]*cphi+e1qb[kk+1]*sphi) + qi*up*(e1qa[kk+nr]*cphi+e1qb[kk+nr]*sphi) + pi*qi*(e1qa[kk+nr+1]*cphi+e1qb[kk+nr+1]*sphi);
    e1[2] = up*uq*(e1za[kk]*cphi+e1zb[kk]*sphi) + pi*uq*(e1za[kk+1]*cphi+e1zb[kk+1]*sphi) + qi*up*(e1za[kk+nr]*cphi+e1zb[kk+nr]*sphi) + pi*qi*(e1za[kk+nr+1]*cphi+e1zb[kk+nr+1]*sphi);
  }

  return;
}

// double RHS_cil(double time, double *u, double *F) {
//   int i;
//   double Ba[3], dB[3], E[3], r[3];
//   double s_flux[1], c=1;
//
//   for (i=0; i<3; i++) {
//     Ba[i] = 0.0;
//     dB[i] = 0.0;
//     E[i] = 0.0;
//     r[i] = u[i];
//   }
//   s_flux[0] = 0.0;
//
//   magnetic_field(Ba, E, u[0], u[1], u[2], s_flux);
//
//   if (s_flux[0] < 0.01) {
//     printf("Partícula perdida.\n");
//     return c = 0;
//   }
//
//   perturbacion(r, time, dB);
//
//   for (i=0; i<3; i++) {
//     Ba[i] += dB[i];
//   }
//
//   F[0] = gam*u[3];
//   F[1] = gam*u[4]/u[0];
//   F[2] = gam*u[5];
//   F[3] = (gam*u[4]*u[4]/u[0]) + (u[4]*Ba[2]-u[5]*Ba[1]) + E[0]; //No pongo el q_Z, VER!! q_Z = 1.0 en lo de Hugo.
//   F[4] = -(gam*u[3]*u[4]/u[0]) + (u[5]*Ba[0]-u[3]*Ba[2]) + E[1]; //No pongo el q_Z, VER!!
//   F[5] = (u[3]*Ba[1]) - (u[4]*Ba[0]) + E[2]; //No pongo el q_Z, VER!!
//
// 	return c;
// }
//
// double RK4(double t, double *ri, double *vi, double *rs, double *vs) {
//   double F[6], k1[6], k2[6], k3[6], k4[6], u[6];  // F=dr;dv-k's-Pos y Vel intermedias.
//   double cc;
//   int i;
//
//   for(i=0;i<6;i++) { //Valores por defecto para las variables.
//     F[i] = 0.0;
//     k1[i] = 0.0;
//     k2[i] = 0.0;
//     k3[i] = 0.0;
//     k4[i] = 0.0;
//     u[i] = 0.0;
//   }
//
//
//   //Aplico las formulas para los k's, uno por cada elemento de drdt y dvdt.
//   //Para k1=h*f(t,u)
//   for(i=0;i<3;i++) {
//     u[i] = ri[i];
//     u[i+3] = vi[i];
//   }
//   cc = RHS_cil(t, u, F); //Es muy problable que la i ya no la necesite.
//   if (!cc) {
//     return cc;
//   }
//   for(i=0;i<6;i++) {
//     k1[i] = dt*F[i];
//   }
//
//   //Para k2=h*f(t+h/2,u+h*k1/2)
//   for(i=0;i<3;i++) {
//     u[i] = ri[i]+(k1[i]*0.5);
//     u[i+3] = vi[i]+(k1[i+3]*0.5);
//   }
//   cc = RHS_cil(t+(dt*0.5), u, F);
//   if (!cc) {
//     return cc;
//   }
//   for(i=0;i<6;i++) {
//     k2[i] = dt*F[i];
//   }
//
//   //Para k3=h*f(t+h/2,u+h*k2/2)
//   for(i=0;i<3;i++) {
//     u[i] = ri[i]+(k2[i]*0.5);
//     u[i+3] = vi[i]+(k2[i+3]*0.5);
//   }
//   cc = RHS_cil(t+(dt*0.5), u, F);
//   if (!cc) {
//     return cc;
//   }
//   for(i=0;i<6;i++) {
//     k3[i] = dt*F[i];
//   }
//
//   //Para k4=h*f(t+h,u+h*k3)
//   for(i=0;i<3;i++) {
//     u[i] = ri[i]+k3[i];
//     u[i+3] = vi[i]+k3[i+3];
//   }
//   cc = RHS_cil(t+dt, u, F);
//   if (!cc) {
//     return cc;
//   }
//   for(i=0;i<6;i++) {
//     k4[i] = dt*F[i];
//   }
//
//
//   //Evoluciono el sistema y=y0+1/6 * h * (k1+2*k2+2*k3+k4)
//   for(i=0;i<3;i++) {
//     rs[i] = ri[i] + (1.0/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
//     vs[i] = vi[i] + (1.0/6.0)*(k1[i+3]+2*k2[i+3]+2*k3[i+3]+k4[i+3]);
//   }
//
//   return cc;
// }

double RHS_cil(int j, double t0, double time, double *u, double *F, double *Bmod2) {
  int i;
  double Ba[3], dB[3], E[3], dE[3], r[3];
  double s_flux[1], c=1, EB=0.0;

  for (i=0; i<3; i++) {
    Ba[i] = 0.0;
    dB[i] = 0.0;
    E[i] = 0.0;
    dE[i] = 0.0;
    r[i] = u[6*j+i];
  }
  s_flux[0] = 0.0;

  magnetic_field(Ba, E, u[6*j], u[6*j+1], u[6*j+2], s_flux);

  if (s_flux[0] < 0.01) {
    printf("Partícula perdida.\n");
    return c = 0;
  }

  perturbacion(r, t0, time, dB, dE);

  for (i=0; i<3; i++) {
    Ba[i] += dB[i];
  }

  //----Resto la componente paralela de B a E----//
  Bmod2[0] = (Ba[0]*Ba[0]+Ba[1]*Ba[1]+Ba[2]*Ba[2]);

  EB = (dE[0]*Ba[0]+dE[1]*Ba[1]+dE[2]*Ba[2]);

  E[0] = dE[0] - (EB)*Ba[0]/Bmod2[0];
  E[1] = dE[1] - (EB)*Ba[1]/Bmod2[0];
  E[2] = dE[2] - (EB)*Ba[2]/Bmod2[0];

  E[0] /= gam;
  E[1] /= gam;
  E[2] /= gam;
  //---------------------------------------------//

  F[0] = gam*u[6*j+3];
  F[1] = gam*u[6*j+4]/u[6*j];
  F[2] = gam*u[6*j+5];
  F[3] = (gam*u[6*j+4]*u[6*j+4]/u[6*j]) + (u[6*j+4]*Ba[2]-u[6*j+5]*Ba[1]) + E[0]; //q_Z = 1.0 en lo de Hugo.
  F[4] = -(gam*u[6*j+3]*u[6*j+4]/u[6*j]) + (u[6*j+5]*Ba[0]-u[6*j+3]*Ba[2]) + E[1]; //q_Z = 1.0 en lo de Hugo.
  F[5] = (u[6*j+3]*Ba[1]) - (u[6*j+4]*Ba[0]) + E[2]; //q_Z = 1.0 en lo de Hugo.

	return c;
}

double RK46_NL(double t0, double t, double *ri, double *vi, double *rs, double *vs, double *Bmod2) {
  double u[42]; 	// velocidades y posiciones
  double w[42]; 	// vel./pos. intermedias
  double F[6];  // dr y dv
  double a[6],b[6],c[6]; // Constantes del método integrador
  double tw, cc;		// tw=tiempo intermedio;
  int i,j;

  for(i=0;i<6;i++) {
    F[i] = 0.0;
  }
  for(i=0;i<42;i++) {
    w[i] = 0.0;
  }

	// Valores de las constantes del método integrador
  a[0]=0.0;					      b[0]=0.032918605146;	c[0]=0.0;
  a[1]=-0.737101392796;		b[1]=0.823256998200;	c[1]=0.032918605146;
  a[2]=-1.634740794341;		b[2]=0.381530948900;	c[2]=0.249351723343;
  a[3]=-0.744739003780;		b[3]=0.200092213184;	c[3]=0.466911705055;
  a[4]=-1.469897351522;		b[4]=1.718581042715;	c[4]=0.582030414044;
  a[5]=-2.813971388035;		b[5]=0.27;				    c[5]=0.847252983783;

	// Condiciones iniciales:
  u[0] = ri[0]; // r_{\rho}
  u[1] = ri[1]; // r_{\theta}
  u[2] = ri[2]; // r_{z}
  u[3] = vi[0]; // v_{\rho}
  u[4] = vi[1]; // v_{\theta}
  u[5] = vi[2]; // v_{z}


  for(i=0;i<6;i++) {
    tw = t + c[i] * dt; //Etapas

    cc = RHS_cil(i, t0, tw, u, F, Bmod2);

    if (!cc) {
      return cc;
    }

    for(j=0;j<6;j++) {	// variables (pos./vel.)
        w[6*i+j+6] = a[i]*w[6*i+j] + dt * F[j];
        u[6*i+j+6] = u[6*i+j] + b[i] * w[6*i+j+6];
      }
  }

  // Actualizo Vel./pos. y tiempo:
  for(i=0;i<3;i++) {
    rs[i] = u[i+36];
    vs[i] = u[i+39];
  }

  //t[0] = t[0] + dt;

  return cc;
}

int integrador(double *t, double *r, double *v, double *rgc, double *vpar, double *mu, int *d) {
  int i, j, c, nstepf;
  double t0=t[0], vper2;
  double rs[3], ri[3], vs[3], vi[3], rg[3], Bmod2[1];

  printf("Integrador con nstep = %d.\n", nstep);

  for (i=0; i<3; i++) {
    ri[i] = r[i];
    vi[i] = v[i];
    rg[i] = 0.0;
  }

  Bmod2[0] = 0.0;

  for (i=1; i<nstep; i++) {
    t[0] = t0 + i * dt;

    // c = RK4(t, ri, vi, rs, vs);
    c = RK46_NL(t0, t[0], ri, vi, rs, vs, Bmod2);

    if (!c) {
      nstepf = i;
      return nstepf;
    }

    r[3*i] = rs[0];
    r[3*i+1] = rs[1];
    r[3*i+2] = rs[2];
    v[3*i] = vs[0];
    v[3*i+1] = vs[1];
    v[3*i+2] = vs[2];

    vpar[i] = centro_giro2(rs, rg, vs, mu, t0, t[0]);

    //----Calculo amu con Bmod usando vperi---//
    vper2 = (vs[0]*vs[0]+vs[1]*vs[1]+vs[2]*vs[2])-(vpar[i]*vpar[i]);
    mu[i] = (gam*vper2)/(2.0*sqrt(Bmod2[0]));

    if (vpar[i-1]*vpar[i] < 0.0) {
      d[i] = 1;
    }

    for (j=0; j<3; j++) {
      rgc[3*i+j] = rg[j];
    }

    for (j=0; j<3; j++) {
      ri[j] = rs[j];
      vi[j] = vs[j];
    }
  }


  return nstep;
}

void PROC(double t0, int nstepf, double *r, double *v,double *rgc, double *vpar,  double *mu, int *d, double *rp, double *vp, double *rgcp, double *vparp, double *mup, double *tp, int *dp) {
  int i,j,k,jstep=((float)nstep)/((float)nprint);
  j = 0;

  printf("PROC con jstep = %d.\n", jstep);

  for (i=0; i<nstepf; i++) {
    k = i%jstep;
    if (k == 0) {
      tp[j] = t0 + i*dt;
      rp[3*j] = r[3*i];
      rp[3*j+1] = r[3*i+1];
      rp[3*j+2] = r[3*i+2];
      vp[3*j] = v[3*i];
      vp[3*j+1] = v[3*i+1];
      vp[3*j+2] = v[3*i+2];
      rgcp[3*j] = rgc[3*i];
      rgcp[3*j+1] = *(rgc+3*i+1);
      rgcp[3*j+2] = rgc[3*i+2];
      vparp[j] = vpar[i];
      mup[j] = mu[i];
      j += 1;
    }
  }

  for (i=0; i<(int)(nstepf/jstep); i++) {
    if (vparp[i+1]*vparp[i]<0.0) {
      dp[i+1] = 1;
    }
  }

}
