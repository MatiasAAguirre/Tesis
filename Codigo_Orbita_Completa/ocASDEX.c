#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j, nstepf, jstep=((float)nstep)/((float)nprint), *d, *dp;
  double B[3], s_flux[1], t[1], rr, zz, t0;
  double *r, *v, *vpar, *mu, *rgc, *rp, *vp, *rgcp, *vparp, *mup, *tp;
  FILE *fr, *frgc, *flux;

  clock_t begin = clock();

  //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
  //Caso sin perturbar:
  fr = fopen("./Outputs/r_PP_dt=0.2_1ms.out", "w"); //_h=0.0003_w=-3.75 _h=3_w=-8_cE_90keV _h=0.0003_w=-8
  frgc = fopen("./Outputs/rgc_PP_dt=0.2_1ms.out", "w"); //_h=0.0003_w=-3.75 _h=3_w=-8_cE_90keV _h=0.0003_w=-8
  flux = fopen("./Outputs/sup_flujo_A.out", "w");

  //Caso Perturbado:
  // fr = fopen("./Outputs/Perturbado/r_perturbado_P40_dt=0.2_cwsE.out", "w");
  // frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P40_dt=0.2_cwsE.out", "w");
  // flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_PP.out", "w");

  r = (double *)malloc(3*nstep*sizeof(double));
  v = (double *)malloc(3*nstep*sizeof(double));

  rgc = (double *)malloc(3*nstep*sizeof(double));
  vpar = (double *)malloc(nstep*sizeof(double));

  mu = (double *)malloc(nstep*sizeof(double));

  d = (int *)malloc(nstep*sizeof(int));

  rp = (double *)malloc(3*nprint*sizeof(double));
  vp = (double *)malloc(3*nprint*sizeof(double));

  rgcp = (double *)malloc(3*nprint*sizeof(double));
  vparp = (double *)malloc(nprint*sizeof(double));

  mup = (double *)malloc(nstep*sizeof(double));

  dp = (int *)malloc(nprint*sizeof(int));

  tp = (double *)malloc(nprint*sizeof(double));

  t[0] = 0.0;

  for (i=0; i<3; i++) {
    B[i]=0.0;
  }
  s_flux[0]=0.0;
  for (i=0; i<nstep; i++) {
    r[i]=0.0;
    r[i+1]=0.0;
    r[i+2]=0.0;
    v[i]=0.0;
    v[i+1]=0.0;
    v[i+2]=0.0;
    rgc[i]=0.0;
    rgc[i+1]=0.0;
    rgc[i+2]=0.0;
    vpar[i]=0.0;
  }
  for (i=0; i<nprint; i++) {
    rp[i]=0.0;
    rp[i+1]=0.0;
    rp[i+2]=0.0;
    vp[i]=0.0;
    vp[i+1]=0.0;
    vp[i+2]=0.0;
    rgcp[i]=0.0;
    rgcp[i+1]=0.0;
    rgcp[i+2]=0.0;
    vparp[i]=0.0;
    tp[i]=0.0;
    dp[i]=0;
  }

  cond_i(t, r, v, rgc, vpar, mu);

  printf("t0 = %.1f, ri = %.6f, qi = %.6f, zi = %.6f, vri = %.6f, vqi = %.6f y vzi = %.6f.\n", t[0], r[0], r[1], r[2], v[0], v[1], v[2]);
  printf("\n");

  t0 = t[0];

  nstepf = integrador(t, r, v, rgc, vpar, mu, d);
  PROC(t0, nstepf, r, v, rgc, vpar, mu, d, rp, vp, rgcp, vparp, mup, tp, dp);
  printf("\n");

  for (i=0; i<(int)(nstepf/jstep); i++) {
    fprintf(fr, "%f %f %f %f %f %f %f\n", tp[i], rp[3*i], rp[3*i+1], rp[3*i+2], vp[3*i], vp[3*i+1], vp[3*i+2]);
    fprintf(frgc, "%f %f %f %f %f %f %d\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i], mup[i], dp[i]);
  }


  for(i=0; i<101; i++) {
    rr = 2.0+3.0*(float)i/100;
    for(j=0; j<101; j++) {
      zz = -2.0+4.0*(float)j/100;
      B_Asdex(rr,zz,B,s_flux);
      fprintf(flux, "%f %f %f\n", rr, zz, s_flux[0]);
    }
    fprintf(flux, "\n");
  }

  free(r);
  free(v);
  free(rgc);
  free(vpar);
  free(rp);
  free(vp);
  free(rgcp);
  free(vparp);
  free(tp);

  fclose(fr);
  fclose(frgc);
  fclose(flux);

  clock_t end = clock();
  double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
  printf("Tiempo invertido en el cómputo: %f.\n", time_spent);

  return 0;
}
