#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cuentaD(double vpar0, double gam);
void B_Asdex(double r,double z,double *B,double *s_flux);
void magnetic_field(double *B, double *E, double r, double q, double z, double *s_flux);
double centro_giro2(double *r, double *rg, double *v, double t, double gam);
void cond_i (double *r, double *v, double *rgc, double *vpar, double gam);
void perturbacion(double *r, double t, double *b1);
void RHS_cil(double time, double *u, double *dr, double gam);
double RK46_NL(double t, double *ri, double *vi, double *rs, double *vs, double gam);
void integrador(double *r, double *v, double *rgc, double *vpar, double gam);
void PROC(double *r, double *v, double *rgc, double *vpar, double *rp, double *vp, double *rgcp, double *vparp, double *tp);

void cuentaD(double vpar0, double gam) {
  double E, m, q, ta;
  double Omega, a, B0, frec, R0, Rl;
  double Gamma, tsim;

  B0 = 2.5; // Tesla
  R0 = 1.71, // radio mayor en metros
  E = 70.e3*1.602e-19;    // 8keV Energy
  m = 2.0*1.6726e-27;           // D
  q = 1.0*1.6022e-19;          // D
  Omega = q*B0/m;       // Cyclotronic frequency
  ta = 1.0439E-8*2/B0;      // Período de Ciclotrón en segundos.
  frec = 5000; //frecuencia del modo 5kHz
  a = 0.5;              // Minor radius
  Rl = sqrt(1-vpar0*vpar0)*a; //Radio de Larmor en metros.
  tsim = nstep*dt*ta;

  printf("------Datos Varios------\n");
  printf("Gamma = %f\n", gam);
  printf("Omega = %.2e rad/s\n", Omega);
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

double centro_giro2(double *r, double *rg, double *v, double t, double gam) {
  int i;
	double vparmod, vpmod, rho, Bmod;
	double vp[3];		// Velocidad perpendicular.
	double e[3];			// Vector unitario perpendicular a v y B.
	double B[3], E[3], dB[3];
  double psi[1];

	rg[0] = r[0];
	rg[1] = r[1];
	rg[2] = r[2];

  psi[0] = 0.0;
	magnetic_field(B, E, r[0], r[1], r[2], psi);

  perturbacion(r, t, dB);

  for (i=0; i<3; i++) {
    B[i]=B[i]+dB[i];
  }

  Bmod = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);

  for(i=0;i<3;i++) {
		vp[i] = v[i]*B[i]/Bmod;			// Vel paralela.
  }

  vparmod = sqrt(vp[0]*vp[0]+vp[1]*vp[1]+vp[2]*vp[2]);

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

  return vparmod;
}

void cond_i(double *r, double *v, double *rgc, double *vpar, double gam) {
  int i;
  double rr0, rq0, rz0, vr0, vq0, vz0, vper2, Bmod, amu;
	double v0[3], r0[3], rg0[3];

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

	vpar[0] = centro_giro2(r0, rg0, v0, 0.0, gam); //Se supone que esto ya no va, igual preguntar.


  for (i=0; i<3; i++) {
    r[i] = r0[i]; //Según lo de Hugo, para r debería sumar la relación de aspecto 171.0/50.0, VER.
    v[i] = v0[i];
    rgc[i] = rg0[i];
  }

  return;
}

void perturbacion(double *r, double t, double *b1) {
  int i, j, ii, jj, kk, cc=0, ci=0;
  static double *b1ra, *b1rb, *b1za, *b1zb;
  static double e1ra[160801], e1rb[160801], e1qa[160801], e1qb[160801], e1za[160801], e1zb[160801];
  double rp = r[0]/2.0, zp = r[2]/2.0, pi, qi, up, uq, aux;

  if (t == 0.0) {
    FILE *Pert_Campo_Mag = fopen("./CamposMagneticoPerturbado_A6_n401.txt", "r");
    b1ra = (double *)malloc(nr*nr*sizeof(double));
    b1rb = (double *)malloc(nr*nr*sizeof(double));
    b1za = (double *)malloc(nr*nr*sizeof(double));
    b1zb = (double *)malloc(nr*nr*sizeof(double));

    for (i=0; i<nr; i++) {
      for (j=0; j<nr; j++) {
        if(fscanf(Pert_Campo_Mag, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &aux, &aux, &b1ra[i*nr+j], &b1rb[i*nr+j], &b1za[i*nr+j], &b1zb[i*nr+j]) == 6) {
          cc++; //Contador elementos correctamente leídos.
        }
        else{
          ci++; //Contador elementos incorrectamente leídos.
          printf("No se leyó bien el archivo de perturbación %d veces para la posición %d.\n", ci, i*nr+j);
          break;
        }
      }
    }

    printf("\n");
    printf("Se leyeron %d lineas correctamente del archivo de perturbación.\n", cc);
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
  }

  else {

    //Esta perte está en lo de Hugo, lo saco para ver si se me perturba.
    if(rp<rp0 || rp>2*rp0 || zp<zp0 || zp>-zp0) //Modificar los rp1 y zp1!!!!
      {
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

    double phi = nmode * r[1] + omega * t;
    double cphi = cos(phi), sphi = sin(phi);

    kk = ii+jj*nr;

    up = 1.0 - pi;
    uq = 1.0 - qi;

    b1[0] = up*uq*(b1ra[kk]*cphi+b1rb[kk]*sphi) + pi*uq*(b1ra[kk+1]*cphi+b1rb[kk+1]*sphi) + qi*up*(b1ra[kk+nr]*cphi+b1rb[kk+nr]*sphi) + pi*qi*(b1ra[kk+nr+1]*cphi+b1rb[kk+nr+1]*sphi);
    b1[1] = 0.0;
    b1[2] = up*uq*(b1za[kk]*cphi+b1zb[kk]*sphi) + pi*uq*(b1za[kk+1]*cphi+b1zb[kk+1]*sphi) + qi*up*(b1za[kk+nr]*cphi+b1zb[kk+nr]*sphi) + pi*qi*(b1za[kk+nr+1]*cphi+b1zb[kk+nr+1]*sphi);
  }

  return;
}

void RHS_cil(double time, double *u, double *F, double gam) {
  int i;
  double Ba[3], dB[3], E[3], r[3];
  double s_flux[1];

  for (i=0; i<3; i++) {
    Ba[i] = 0.0;
    E[i] = 0.0;
    r[i] = u[i];
  }
  s_flux[0] = 0.0;

  magnetic_field(Ba, E, u[0], u[1], u[2], s_flux);

  if (s_flux[0] < 0.01) {
    printf("Partícula perdida.\n");
    exit(0);
  }

  perturbacion(r, time, dB);

  for (i=0; i<3; i++) {
    Ba[i] = Ba[i] + dB[i];
  }

  F[0] = gam*u[3];
  F[1] = gam*u[4]/u[0];
  F[2] = gam*u[5];
  F[3] = gam*u[4]*u[4]/u[0] + (u[4]*Ba[2]-u[5]*Ba[1]) + E[0]; //No pongo el q_Z, VER!! q_Z = 1.0 en lo de Hugo.
  F[4] = -gam*u[3]*u[4]/u[0] + (u[5]*Ba[0]-u[3]*Ba[2]) + E[1]; //No pongo el q_Z, VER!!
  F[5] = u[3]*Ba[1] - u[4]*Ba[0] + E[2]; //No pongo el q_Z, VER!!

	return;
}

double RK46_NL(double t, double *ri, double *vi, double *rs, double *vs, double gam) {
  double u[6]; 	// velocidades y posiciones
  double w[6]; 	// vel./pos. intermedias
  double F[6];  // dr y dv
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
  u[0] = ri[0]; // r_{\rho}
  u[1] = ri[1]; // r_{\theta}
  u[2] = ri[2]; // r_{z}
  u[3] = vi[0]; // v_{\rho}
  u[4] = vi[1]; // v_{\theta}
  u[5] = vi[2]; // v_{z}


  for(i=0;i<6;i++) {
    tw = t + c[i] * dt; //Etapas

    RHS_cil(tw, u, F, gam);

  for(j=0;j<6;j++) {	// variables (pos./vel.)
      w[j] = a[i]*w[j] + dt * F[j];
      u[j] = u[j] + b[i] * w[j];
    }
  }

  // Actualizo Vel./pos. y tiempo:
  for(i=0;i<3;i++) {
    rs[i] = u[i];
    vs[i] = u[i+3];
  }

  t = t + dt;

  return t;
}

void integrador(double *r, double *v, double *rgc, double *vpar, double gam) {
  int i, j;
  double t;
  double rs[3], ri[3], vs[3], vi[3], rg[3];

  printf("Integrador con nstep = %d.\n", nstep);

  for (i=0; i<3; i++) {
    ri[i] = r[i];
    vi[i] = v[i];
    rg[i] = 0.0;
  }

  for (i=1; i<nstep; i++) {
    t = i * dt;

    RK46_NL(t, ri, vi, rs, vs, gam); //Saqué el t = RK46_NL() que estaba antes, veamos si sigue funcionando.

    r[3*i] = rs[0];
    r[3*i+1] = rs[1];
    r[3*i+2] = rs[2];
    v[3*i] = vs[0];
    v[3*i+1] = vs[1];
    v[3*i+2] = vs[2];

    vpar[i] = centro_giro2(rs, rg, vs, t, gam);
    for (j=0; j<3; j++) {
      rgc[3*i+j] = rg[j];
    }

    for (j=0; j<3; j++) {
      ri[j] = rs[j];
      vi[j] = vs[j];
    }
  }


  return;
}

void PROC(double *r, double *v, double *rgc, double *vpar, double *rp, double *vp, double *rgcp, double *vparp, double *tp) {
  int i,j,k,jstep=((float)nstep)/((float)nprint);
  j = 0;

  printf("PROC con jstep = %d.\n", jstep);

  for (i=0; i<nstep; i++) {
    k = i%jstep;
    if (k == 0) {
      tp[j] = i*dt;
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
      j = j + 1;
    }
  }
}
