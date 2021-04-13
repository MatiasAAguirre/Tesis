#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void B_Asdex(double r,double z,double *B,double *s_flux);
void magnetic_field(double *B, double *E, double r, double q, double z, double *s_flux);
void centro_giro(double *r, double *rg, double *vi, double *B);
double cond_i (double *rgc, double *vpar);
void perturbacion(double *rgc, double t, double *b1, double *e1, double *db1);
void campo(double *rgc, double t, double *Ba, double *Bmod, double *E, double *dBm, double *dbv);
double RHS_cil(int j, double time, double amu, double *rgce, double *drgc, double vpar);
double RK46_NL(double amu, double t, double *rgc, double *vpari, double *rgcs);
void integrador(double amu, double *rgc, double *vpar);
void PROC(double amu, double *rgc, double *vpar, double *rgcp, double *vparp, double *tp);

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

void centro_giro(double *r, double *rg, double *vi, double *B) {
  int i;
	double vpmod, rho, Bmod;
	double vp[3];		// Velocidad perpendicular.
	double e[3];			// Vector unitario perpendicular a v y B.
	double E[3];
  double psi[1];

	rg[0] = r[0];
	rg[1] = r[1];
	rg[2] = r[2];

  psi[0] = 0.0;
	magnetic_field(B, E, rg[0], rg[1], rg[2], psi);

	e[0] = vi[1]*B[2] - vi[2]*B[1];
	e[1] = vi[2]*B[0] - vi[0]*B[2];
	e[2] = vi[0]*B[1] - vi[1]*B[0];

	Bmod = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  for (i=0; i<3; i++) {
    vp[i] = e[i]/Bmod;
  }

	vpmod = sqrt(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);

	rho = gam*vpmod/Bmod;

	// Devuelve la posicion de la particula.
	if (rho < 1.0E-8) {
    return;
  }

	for(i=0;i<3;i++) {
    e[i] = e[i]/(vpmod*Bmod);
  }

	// Posición del centro de giro:
	for(i=0;i<3;i++) {
    rg[i] = r[i] + rho*e[i]; //Está con un + porque es v x B en vez de B x v.
  }

  return;
}

double cond_i(double *rgc, double *vpar) {
  int i;
  double rr0, rq0, rz0, vr0, vq0, vz0, vper2, Bmod, amu;
	double v0[3], r0[3], rg0[3], B[3];

	/* ------ Particula Pasante ------ */
	rr0 = 4.022350;
	rq0 = 196.003;
	rz0 = 0.623175;

	vr0 = 0.075035;
	vq0 = 0.64219;
	vz0 = 0.762512;

	/* ------ Particula Atrapada ------ */
	// rr0 = 3.41970;
	// rq0 = 21.2890;
	// rz0 = -0.781325;
  //
	// vr0 = 0.0181466;
	// vq0 = 0.0394052;
	// vz0 = 0.998925;

	r0[0] = rr0;
	r0[1] = rq0;
	r0[2] = rz0;

	v0[0] = vr0;
	v0[1] = vq0;
	v0[2] = vz0;

	centro_giro(r0, rg0, v0, B);

  for (i=0; i<3; i++) {
    rgc[i] = rg0[i];
  }

	Bmod = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	vpar[0] = (v0[0]*B[0]+v0[1]*B[1]+v0[2]*B[2])/Bmod;

	vper2 = 1.0 - vpar[0]*vpar[0];
	amu = (gam*vper2)/(2.0*Bmod);

  return amu;
}

void perturbacion(double *rgc, double t, double *b1, double *e1, double *db1) {
  int i, j, ii, jj, kk, cc=0, ci=0;
  static double *b1ra, *b1rb, *b1za, *b1zb;
  static double e1ra[160801], e1rb[160801], e1qa[160801], e1qb[160801], e1za[160801], e1zb[160801];
  double rp = rgc[0]/2.0, zp = rgc[2]/2.0, pi, qi, up, uq, aux;

  if (t == 0.0) {
    FILE *Pert_Campo_Mag = fopen("./CamposMagneticoPerturbado_A6_n401.txt", "r");
    b1ra = (double *)malloc(nr*nr*sizeof(double));
    b1rb = (double *)malloc(nr*nr*sizeof(double));
    b1za = (double *)malloc(nr*nr*sizeof(double));
    b1zb = (double *)malloc(nr*nr*sizeof(double));

    for (i=0; i<nr; i++) {
      for (j=0; j<nr; j++) {
        if(fscanf(Pert_Campo_Mag, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &aux, &aux, &b1ra[i*nr+j], &b1rb[i*nr+j], &b1za[i*nr+j], &b1zb[i*nr+j]) == 6) {
          cc++;
        }
        else{
          ci++;
          printf("No se leyó bien el archivo de perturbación %d veces para la posición %d.\n", ci, i*nr+j);
          break;
        }
      }
    }

    printf("Se leyeron %d lineas correctamente del archivo de perturbación.\n", cc);

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
  }

  else {
    double phi = nmode * rgc[1] + omega * t;
    double cphi = cos(phi), sphi = sin(phi);

    ii = (int)((rp-rp0)/hr);
    jj = (int)((zp-zp0)/hz);
    kk = ii+jj*nr;

    pi = rp/hr - (float)ii - rp0/hr;
    up = 1.0 - pi;
    qi = zp/hz - (float)jj - zp0/hz;
    uq = 1.0 - qi;

    b1[0] = up*uq*(b1ra[kk]*cphi+b1rb[kk]*sphi) + pi*uq*(b1ra[kk+1]*cphi+b1rb[kk+1]*sphi) + qi*up*(b1ra[kk+nr]*cphi+b1rb[kk+nr]*sphi) + pi*qi*(b1ra[kk+nr+1]*cphi+b1rb[kk+nr+1]*sphi);
    b1[1] = 0.0;
    b1[2] = up*uq*(b1za[kk]*cphi+b1zb[kk]*sphi) + pi*uq*(b1za[kk+1]*cphi+b1zb[kk+1]*sphi) + qi*up*(b1za[kk+nr]*cphi+b1zb[kk+nr]*sphi) + pi*qi*(b1za[kk+nr+1]*cphi+b1zb[kk+nr+1]*sphi);


    double db1rdrk = ((b1ra[kk+1]-b1ra[kk-1])*cphi + (b1rb[kk+1]-b1rb[kk-1])*sphi)/(2.0*hr);
    double db1rdrk1 = ((b1ra[kk+2]-b1ra[kk])*cphi + (b1rb[kk+2]-b1rb[kk])*sphi)/(2.0*hr);
    double db1rdrknr = ((b1ra[kk+nr+1]-b1ra[kk+nr-1])*cphi + (b1rb[kk+nr+1]-b1rb[kk+nr-1])*sphi)/(2.0*hr);
    double db1rdrknr1 = ((b1ra[kk+2+nr]-b1ra[kk+nr])*cphi + (b1rb[kk+2+nr]-b1rb[kk+nr])*sphi)/(2.0*hr);

    db1[0] = up*uq*db1rdrk + pi*uq*db1rdrk1 + qi*up*db1rdrknr + qi*pi*db1rdrknr1;


    double db1zdrk = ((b1za[kk+1]-b1za[kk-1])*cphi + (b1zb[kk+1]-b1zb[kk-1])*sphi)/(2.0*hr);
    double db1zdrk1=((b1za[kk+2]-b1za[kk])*cphi + (b1zb[kk+2]-b1zb[kk])*sphi)/(2.0*hr);
    double db1zdrknr=((b1za[kk+nr+1]-b1za[kk+nr-1])*cphi + (b1zb[kk+nr+1]-b1zb[kk+nr-1])*sphi)/(2.0*hr);
    double db1zdrknr1=((b1za[kk+2+nr]-b1za[kk+nr])*cphi + (b1zb[kk+2+nr]-b1zb[kk+nr])*sphi)/(2.0*hr);

    db1[1] = up*uq*db1zdrk + pi*uq*db1zdrk1 + qi*up*db1zdrknr + qi*pi*db1zdrknr1;


    db1[2] = nmode*(-b1ra[kk]*sphi + b1rb[kk]*cphi);


    db1[3] = nmode*(-b1za[kk]*sphi + b1zb[kk]*cphi);


    double db1rdzk=((b1ra[kk+nr]-b1ra[kk-nr])*cphi + (b1rb[kk+nr]-b1rb[kk-nr])*sphi)/(2.0*hz);
    double db1rdzk1=((b1ra[kk+1+nr]-b1ra[kk+1-nr])*cphi + (b1rb[kk+1+nr]-b1rb[kk+1-nr])*sphi)/(2.0*hz);
    double db1rdzknr=((b1ra[kk+nr+nr]-b1ra[kk])*cphi + (b1rb[kk+nr+nr]-b1rb[kk])*sphi)/(2.0*hz);
    double db1rdzknr1=((b1ra[kk+1+nr+nr]-b1ra[kk+1])*cphi + (b1rb[kk+1+nr+nr]-b1rb[kk+1])*sphi)/(2.0*hz);

    db1[4]=up*uq*db1rdzk+pi*uq*db1rdzk1+qi*up*db1rdzknr+qi*pi*db1rdzknr1;


    double db1zdzk = ((b1za[kk+nr]-b1za[kk-nr])*cphi + (b1zb[kk+nr]-b1zb[kk-nr])*sphi)/(2.0*hz);
    double db1zdzk1 = ((b1za[kk+1+nr]-b1za[kk+1-nr])*cphi + (b1zb[kk+1+nr]-b1zb[kk+1-nr])*sphi)/(2.0*hz);
    double db1zdzknr = ((b1za[kk+nr+nr]-b1za[kk])*cphi + (b1zb[kk+nr+nr]-b1zb[kk])*sphi)/(2.0*hz);
    double db1zdzknr1 = ((b1za[kk+1+nr+nr]-b1za[kk+1])*cphi + (b1zb[kk+1+nr+nr]-b1zb[kk+1])*sphi)/(2.0*hz);

    db1[5] = up*uq*db1zdzk + pi*uq*db1zdzk1 + qi*up*db1zdzknr + qi*pi*db1zdzknr1;


    db1[6] = omega*(-b1ra[kk]*sphi + b1rb[kk]*cphi);


    db1[7] = omega*(-b1za[kk]*sphi + b1zb[kk]*cphi);


    e1[0] = up*uq*(e1ra[kk]*cphi+e1rb[kk]*sphi) + pi*uq*(e1ra[kk+1]*cphi+e1rb[kk+1]*sphi) + qi*up*(e1ra[kk+nr]*cphi+e1rb[kk+nr]*sphi) + pi*qi*(e1ra[kk+nr+1]*cphi+e1rb[kk+nr+1]*sphi);

    e1[0] = e1[0] * omega;


    e1[1] = up*uq*(e1qa[kk]*cphi+e1qb[kk]*sphi) + pi*uq*(e1qa[kk+1]*cphi+e1qb[kk+1]*sphi) + qi*up*(e1qa[kk+nr]*cphi+e1qb[kk+nr]*sphi) + pi*qi*(e1qa[kk+nr+1]*cphi+e1qb[kk+nr+1]*sphi);

    e1[1] = e1[1] * omega;


    e1[2] = up*uq*(e1za[kk]*cphi+e1zb[kk]*sphi) + pi*uq*(e1za[kk+1]*cphi+e1zb[kk+1]*sphi) + qi*up*(e1za[kk+nr]*cphi+e1zb[kk+nr]*sphi) + pi*qi*(e1za[kk+nr+1]*cphi+e1zb[kk+nr+1]*sphi);

    e1[2] = e1[2] * omega;
  }

  return;
}

void campo(double *rgc, double t, double *Ba, double *Bmodv, double *E, double *dBm, double *dbv) {
  int i;
  double dBmdt, Bmod, ar, az;
  double B[3], s_flux[1], Bp[6], Bm[6], dB0[6], b1[3], e1[3], db1[8], b[3], dB[9];

  dBmdt = 0.0;
  s_flux[0] = 0.0;
  for (i=0; i<3; i++) {
    B[i] = 0.0;
    b1[i] = 0.0;
    e1[i] = 0.0;
    b[i] = 0.0;
  }
  for (i=0; i<6; i++) {
    Bp[i] = 0.0;
    Bm[i] = 0.0;
    dB0[i] = 0.0;
  }
  for (i=0; i<8; i++) {
    db1[i] = 0.0;
  }
  for (i=0; i<9; i++) {
    dB[i] = 0.0;
  }

  ar = rgc[0];
  az = rgc[2];
  B_Asdex(ar, az, B, s_flux);
  Ba[0] = B[0];
  Ba[1] = B[1];
  Ba[2] = B[2];

  if (s_flux[0]<0.01) {
    printf("Partícula perdida.\n");
    exit(0); //Ver como terminar el programa con esto, no me deja o no entendí el error.
  }

  ar = rgc[0] + dr;
  B_Asdex(ar, az, B, s_flux);
  Bp[0] = B[0];
  Bp[1] = B[1];
  Bp[2] = B[2];

  ar = rgc[0] - dr;
  B_Asdex(ar, az, B, s_flux);
  Bm[0] = B[0];
  Bm[1] = B[1];
  Bm[2] = B[2];

  ar = rgc[0];
  az = rgc[2] + dz;
  B_Asdex(ar, az, B, s_flux);
  Bp[3] = B[0];
  Bp[4] = B[1];
  Bp[5] = B[2];

  az = rgc[2] - dz;
  B_Asdex(ar, az, B, s_flux);
  Bm[3] = B[0];
  Bm[4] = B[1];
  Bm[5] = B[2];

  az = rgc[2];

  dB0[0]=(Bp[0]-Bm[0])/(2.0*dr); //dBrdr
  dB0[1]=(Bp[1]-Bm[1])/(2.0*dr); //dBqdr
  dB0[2]=(Bp[2]-Bm[2])/(2.0*dr); //dBzdr
  dB0[3]=(Bp[3]-Bm[3])/(2.0*dz); //dBrdz
  dB0[4]=(Bp[4]-Bm[4])/(2.0*dz); //dBqdz
  dB0[5]=(Bp[5]-Bm[5])/(2.0*dz); //dBzdz

  //Perturbaciones de los campos y deribadas.
  perturbacion(rgc, t, b1, e1, db1);

  b[0] = Ba[0]+b1[0];
  b[1] = Ba[1]+b1[1];
  b[2] = Ba[2]+b1[2];

  E[0] = e1[0];
  E[1] = e1[1];
  E[2] = e1[2];

  dB[0] = dB0[0] + db1[0]; //dBrdr c/ db1rdr
  dB[1] = dB0[1]; //dBqdr
  dB[2] = dB0[2] + db1[1]; //dBzdr c/ db1zdr
  dB[3] = db1[2]; //dBrdq c/ db1rdq
  dB[4] = 0.0; //dBqdq
  dB[5] = db1[3]; //dBzdq c/ db1zdq
  dB[6] = dB0[3] + db1[4]; //dBrdz c/ db1rdz
  dB[7] = dB0[4]; //dBqdz
  dB[8] = dB0[5] + db1[5]; //dBzdz c/ db1zdz

  Bmodv[0] = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
  Bmod = Bmodv[0];

  dBm[0] = (b[0]*dB[0]+b[2]*dB[2]+b[1]*dB[1])/Bmod; //dBmdr
  dBm[1] = (b[0]*dB[3]+b[2]*dB[5]+b[1]*dB[4])/Bmod; //dBmdq
  dBm[2] = (b[0]*dB[6]+b[2]*dB[8]+b[1]*dB[7])/Bmod; //dBmdz

  dBmdt = (b[0]*db1[6]+b[2]*db1[7])/Bmod;

  dbv[0] = (dB[1]-b[1]*dBm[0]/Bmod)/Bmod; //dbvqdr
  dbv[1] = (dB[2]-b[2]*dBm[0]/Bmod)/Bmod; //dbvzdr
  dbv[2] = (dB[3]-b[0]*dBm[1]/Bmod)/Bmod; //dbvrdq
  dbv[3] = (dB[5]-b[2]*dBm[1]/Bmod)/Bmod; //dbvzdq
  dbv[4] = (dB[6]-b[0]*dBm[2]/Bmod)/Bmod; //dbvrdz
  dbv[5] = (dB[7]-b[1]*dBm[2]/Bmod)/Bmod; //dbvqdz

  dbv[6] = (db1[6]-b[0]*dBmdt/Bmod)/Bmod; //dbvrdt
  dbv[7] = -b[1]*dBmdt/(Bmod*Bmod); //dbvqdt
  dbv[8] = (db1[7]-b[2]*dBmdt/Bmod)/Bmod; //dbvzdt

  Ba[0] = b[0];
  Ba[1] = b[1];
  Ba[2] = b[2];

  return;
}

double RHS_cil(int j, double time, double amu, double *rgce, double *drgc, double vpar) {
  int i;
  double Bmodv[1], bv[3], Ba[3], dBm[3], dbv[9], E[3], rgc[3];
  double Bmod, Bpars, dvpar;

  for (i=0; i<3; i++) {
    bv[i] = 0.0;
    Ba[i] = 0.0;
    dBm[i] = 0.0;
    E[i] = 0.0;
  }
  for (i=0; i<9; i++) {
    dbv[i] = 0.0;
  }
  Bmodv[0] = 0.0;

  for (i=0; i<3; i++) {
    rgc[i] = rgce[4*j+i];
  }

  campo(rgc, time, Ba, Bmodv, E, dBm, dbv);

  Bmod = Bmodv[0];
  bv[0] = Ba[0] / Bmod;
  bv[1] = Ba[1] / Bmod;
  bv[2] = Ba[2] / Bmod;

  Bpars = bv[0]*(Ba[0]+gam*vpar*(dbv[3]/rgc[0]-dbv[5]))+bv[1]*(Ba[1]+gam*vpar*(dbv[4]-dbv[1]))+bv[2]*(Ba[2]+gam*vpar*(bv[1]/rgc[0]+dbv[0]-dbv[2]/rgc[0]));
  drgc[0] = gam*(vpar*Ba[0]+gam*vpar*vpar*(dbv[3]/rgc[0]-dbv[5])+E[1]*bv[2]-E[2]*bv[1]-vpar*(dbv[7]*bv[2]-dbv[8]*bv[1])+amu*(bv[1]*dBm[2]-bv[2]*dBm[1]/rgc[0]))/Bpars;
  drgc[2] = gam*(vpar*Ba[2]+gam*vpar*vpar*(bv[1]/rgc[0]+dbv[0]-dbv[2]/rgc[0])+E[0]*bv[1]-E[1]*bv[1]-vpar*(dbv[6]*bv[1]-dbv[7]*bv[0])+amu*(dBm[1]*bv[0]/rgc[0]-dBm[0]*bv[1]))/Bpars;
  drgc[1] = gam*(vpar*Ba[1]+gam*vpar*vpar*(dbv[4]-dbv[1])+E[2]*bv[0]-E[0]*bv[2]-vpar*(dbv[8]*bv[0]-dbv[6]*bv[2])+amu*(dBm[0]*bv[2]-dBm[2]*bv[0]))/(Bpars*rgc[0]);
  dvpar = ((Ba[0]+gam*vpar*(dbv[3]/rgc[0]-dbv[5]))*(E[0]-amu*dBm[0])+(Ba[1]+gam*vpar*(dbv[4]-dbv[1]))*(E[1]-amu*dBm[1]/rgc[0])+(Ba[2]+gam*vpar*(bv[1]/rgc[0]+dbv[0]-dbv[2]/rgc[0]))*(E[2]-amu*dBm[2]))/Bpars;

	return dvpar;
}

double RK46_NL(double amu, double t, double *rgc, double *vpari, double *rgcs) {
  double u[28]; 	// velocidades y posiciones
  double w[28]; 	// vel./pos. intermedias
  double F[4];  // drgc
  double a[6],b[6],c[6]; // Constantes del método integrador
  double tw, dvpar;		// tiempo intermedio;
  int i,j;

  for(i=0;i<4;i++) {
    w[i] = 0.0;
    F[i] = 0.0;
  }

	// Valores de las constantes del método integrador
  a[0]=0.0;					      b[0]=0.032918605146;	c[0]=0.0;
  a[1]=-0.737101392796;		b[1]=0.823256998200;	c[1]=0.032918605146;
  a[2]=-1.634740794341;		b[2]=0.381530948900;	c[2]=0.249351723343;
  a[3]=-0.744739003780;		b[3]=0.200092213184;	c[3]=0.466911705055;
  a[4]=-1.469897351522;		b[4]=1.718581042715;	c[4]=0.582030414044;
  a[5]=-2.813971388035;		b[5]=0.27;				    c[5]=0.847252983783;

	// Condiciones iniciales:
  u[0] = rgc[0]; // r_{\rho}
  u[1] = rgc[1]; // r_{\theta}
  u[2] = rgc[2]; // r_{z}
  u[3] = vpari[0]; // vpar


  for(i=0;i<6;i++) {
    tw = t + c[i] * dt; //Etapas

    dvpar = RHS_cil(i, tw, amu, u, F, u[4*i+3]);

    for(j=0;j<3;j++) {	// variables (pos./vel.)
      w[4*i+j+4] = a[i]*w[4*i+j] + dt * F[j];
      u[4*i+j+4] = u[4*i+j] + b[i] * w[4*i+j+4];
    }
    w[4*i+7] = a[i]*w[4*i+3] + dt * dvpar;
    u[4*i+7] = u[4*i+3] + b[i] * w[4*i+7];
  }

  // Actualizo Vel./pos. y tiempo:
  for(i=0;i<3;i++) {
    rgcs[i] = u[i+24];
  }
  vpari[0] = u[27];
  t = t + dt;

  return t;
}

void integrador(double amu, double *rgc, double *vpar) {
  int i, j;
  double t;
  double rgcs[3], rgci[3], vpari[1];

  printf("Integrador con nstep = %d.\n", nstep);

  for (i=0; i<3; i++) {
    rgci[i] = rgc[i];
  }
  vpari[0] = vpar[0];

  for (i=1; i<nstep; i++) {
    t = i * dt;
    t = RK46_NL(amu, t, rgci, vpari, rgcs); //Creo que este t está de más, pero lo dejo porque así funciona.

    rgc[3*i] = rgcs[0];
    rgc[3*i+1] = rgcs[1];
    rgc[3*i+2] = rgcs[2];
    vpar[i] = vpari[0];

    for (j=0; j<3; j++) {
      rgci[j] = rgcs[j];
    }
  }


  return;
}

void PROC(double amu, double *rgc, double *vpar, double *rgcp, double *vparp, double *tp) {
  int i,j,k,jstep=((float)nstep)/((float)nprint);
  j = 0;

  printf("PROC con jstep = %d.\n", jstep);

  for (i=0; i<nstep; i++) {
    k = i%jstep;
    if (k == 0) {
      tp[j] = i*dt;
      rgcp[3*j] = rgc[3*i];
      rgcp[3*j+1] = rgc[3*i+1];
      rgcp[3*j+2] = rgc[3*i+2];
      vparp[j] = vpar[i];
      j = j + 1;
    }
  }
}
