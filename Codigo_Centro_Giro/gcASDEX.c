#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j, nstepf, jstep=((float)nstep)/((float)nprint), *dp;
  double B[3], s_flux[1], t[1], amu = 0.0, r, z, t0;
  //double rgc0[3];
  double *rgc, *vpar, *vper, *mu, *rgcp, *vparp, *vperp, *tp, *mup, r_B[1], r_E[1];
  //char part[3];
  FILE *frgc, *flux, *relas;

  clock_t begin = clock();

  //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
  //Caso sin perturbar:
  frgc = fopen("./Outputs/rgc_P40_dt=20_h=3_w=-8_es=0.01_90keV_1ms.out", "w"); //(estructura: rgc_particula_dt_h_w_cE_E0keV_t-computo.out)
  flux = fopen("./Outputs/sup_flujo_P.out", "w");
  relas = fopen("./Outputs/relaciones_ByE.out", "w");

  //Caso perturbado:
  // frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P40_cwsE.out", "w");
  // flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_PP.out", "w");

  rgc = (double *)malloc(3*nstep*sizeof(double));
  vpar = (double *)malloc(nstep*sizeof(double));
  vper = (double *)malloc(nstep*sizeof(double));
  mu = (double *)malloc(nstep*sizeof(double));

  rgcp = (double *)malloc(3*nprint*sizeof(double));
  vparp = (double *)malloc(nprint*sizeof(double));
  vperp = (double *)malloc(nprint*sizeof(double));
  mup = (double *)malloc(nprint*sizeof(double));

  tp = (double *)malloc(nprint*sizeof(double));
  dp = (int *)malloc(nprint*sizeof(int));

  t[0] = 0.0;
  r_B[0] = 0.0;
  r_E[0] = 0.0;

  for (i=0; i<3; i++) {
    B[i]=0.0;
    // rgc0[i]=0.0;
  }
  vpar[0]=0.0;
  s_flux[0]=0.0;

  for (i=0; i<nstep; i++) {
    rgc[3*i]=0.0;
    rgc[3*i+1]=0.0;
    rgc[3*i+2]=0.0;
    vpar[i]=0.0;
    vper[i]=0.0;
    mu[i]=0.0;
  }
  for (i=0; i<nprint; i++) {
    tp[i]=0.0;
    rgcp[3*i]=0.0;
    rgcp[3*i+1]=0.0;
    rgcp[3*i+2]=0.0;
    vparp[i]=0.0;
    vperp[i]=0.0;
    mup[i]=0.0;
    dp[i] = 0;
  }

  //-----------Cálculo de Trayectorias-------------//
  // printf("\n");
  // perturbacion(rgc0, t[0], t[0], b1, e1, db1);
  // printf("\n");

  amu = cond_i(t, rgc, vpar, vper, mu);

  t0 = t[0];

  printf("rgci = %.16f, qgci = %.16f, zgci = %.16f, vpari = %.16f y amu = %.16f.\n", rgc[0], rgc[1], rgc[2], vpar[0], amu);
  printf("\n");

  nstepf = integrador(amu, t, rgc, vpar, vper, mu, r_B, r_E);

  r_B[0] /= nstepf;
  r_E[0] /= nstepf;

  PROC(amu, t0, nstepf, rgc, vpar, vper, mu, rgcp, vparp, vperp, tp, dp, mup);
  printf("nstepf = %d, con r_B = %f y r_E = %f\n", nstepf, r_B[0], r_E[0]);

  for (i=1; i<(int)(nstepf/jstep); i++) {
    fprintf(frgc, "%f %f %f %f %f %f %f %d\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i], vperp[i], mup[i], dp[i]);
  }
  //----------------------------------------//

  //----------Superficie de Flujo----------//
  for(i=0; i<101; i++) {
    r = 2.0+3.0*(float)i/100;
    for(j=0; j<101; j++) {
      z = -2.0+4.0*(float)j/100;
      B_Asdex(r,z,B,s_flux);
      fprintf(flux, "%f %f %f\n", r, z, s_flux[0]);
    }
    fprintf(flux, "\n");
  }
  //---------------------------------------//

  //----------Longitudes de Campo----------//
  for(i=0; i<101; i++) {
    r = 2.0+3.0*(float)i/100;
    for(j=0; j<101; j++) {
      z = -2.0+4.0*(float)j/100;
      double rz[3];
      rz[0] = r;
      rz[1] = 0.0;
      rz[2] = z;
      double time = (2.0*3.14159/omega);
      double Ba[3], Bmodv[1], E[3], dBm[3], dEm[3], dbv[9];

      campo(rz, t0, time, Ba, Bmodv, E, dBm, dEm, dbv);

      double Bmod = Bmodv[0];
      double Emod = sqrt(E[0]*E[0]+E[1]*E[1]+E[2]*E[2]);

      double L_B = Bmod/sqrt(dBm[0]*dBm[0]+dBm[1]*dBm[1]+dBm[2]*dBm[2]);
      double L_E = Emod/sqrt(dEm[0]*dEm[0]+dEm[1]*dEm[1]+dEm[2]*dEm[2]);

      fprintf(relas, "%f %f %f %f\n", r, z, 0.032/L_B, 0.032/L_E);
    }
    fprintf(relas, "\n");
  }
  //---------------------------------------//

  free(rgc);
  free(vpar);
  free(vper);
  free(mu);
  free(rgcp);
  free(vparp);
  free(vperp);
  free(mup);
  free(tp);
  free(dp);

  fclose(frgc);
  fclose(flux);
  fclose(relas);

  clock_t end = clock();
  double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
  printf("Tiempo invertido en el cómputo: %f s.\n", time_spent);

  return 0;
}
