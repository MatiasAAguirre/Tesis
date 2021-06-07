#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j, nstepf, jstep=((float)nstep)/((float)nprint);
  double B[3], s_flux[1], amu = 0.0, rr, zz;
  double *r, *v, *vpar, *rgc, *rp, *vp, *rgcp, *vparp, *tp;
  double t = 0, r0[3], b1[3], e1[3], db1[8]; //No se si ahora db1 tiene 8 componentes, VER.
  FILE *fr, *frgc, *flux;

  clock_t begin = clock();

  //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
  //Caso sin perturbar:
  // fr = fopen("./Outputs/r_P.out", "w");
  // frgc = fopen("./Outputs/rgc_P.out", "w");
  // flux = fopen("./Outputs/sup_flujo_P.out", "w");

  //Caso Perturbado:
  fr = fopen("./Outputs/Perturbado/r_perturbado_P_h=0.001.out", "w");
  frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P_h=0.001.out", "w");
  flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_P.out", "w");

  r = (double *)malloc(3*nstep*sizeof(double));
  v = (double *)malloc(3*nstep*sizeof(double));

  rgc = (double *)malloc(3*nstep*sizeof(double));
  vpar = (double *)malloc(nstep*sizeof(double));

  rp = (double *)malloc(3*nprint*sizeof(double));
  vp = (double *)malloc(3*nprint*sizeof(double));

  rgcp = (double *)malloc(3*nprint*sizeof(double));
  vparp = (double *)malloc(nprint*sizeof(double));

  tp = (double *)malloc(nprint*sizeof(double));


  for (i=0; i<3; i++) {
    B[i]=0.0;
    r0[i]=0.0;
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
  }

  cond_i(r, v, rgc, vpar); //CUIDADO, ahora no tengo el amu definido.

  printf("ri = %.6f, qi = %.6f, zi = %.6f, vri = %.6f, vqi = %.6f y vzi = %.6f.\n", r[0], r[1], r[2], v[0], v[1], v[2]);
  printf("\n");

  nstepf = integrador(r, v, rgc, vpar);
  PROC(nstepf, r, v, rgc, vpar, rp, vp, rgcp, vparp, tp);
  printf("\n");

  for (i=0; i<(int)(nstepf/jstep); i++) {
    fprintf(fr, "%f %f %f %f %f %f %f\n", tp[i], rp[3*i], rp[3*i+1], rp[3*i+2], vp[3*i], vp[3*i+1], vp[3*i+2]);
    fprintf(frgc, "%f %f %f %f %f\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i]);
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
