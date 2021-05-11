#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j;
  double B[3], s_flux[1], amu = 0.0, rr, zz;
  double *r, *v, *rgc, *vpar, *rp, *vp, *rgcp, *vparp, *tp;
  double t = 0, gam, r0[3], b1[3], e1[3], db1[8]; //No se si ahora db1 tiene 8 componentes, VER.
  char *frs, *frgcs;
  FILE *fr, *frgc, *flux;

  frs = malloc(11*sizeof("./Outputs/r_P_gam=0.083.out"));
  frgcs = malloc(11*sizeof("./Outputs/rgc_P_gam=0.083.out"));

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
  }
  for (i=0; i<nprint; i++) {
    rp[i]=0.0;
    rp[i+1]=0.0;
    rp[i+2]=0.0;
    rgc[i]=0.0;
    rgc[i+1]=0.0;
    rgc[i+2]=0.0;
    vp[i]=0.0;
    vp[i+1]=0.0;
    vp[i+2]=0.0;
    vpar[i]=0.0;
  }

  for (j=0; j<11; j++) {
    gam = 0.00327451 + 0.08*j/10; //Recorre los valores gam = [0.003;0.083]

    //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
    //Caso sin perturbar:
    sprintf(frs, "./Outputs/r_P_gam=%.3f.out", gam);
    sprintf(frgcs, "./Outputs/rgc_P_gam=%.3f.out", gam);

    fr = fopen(frs, "w");
    frgc = fopen(frgcs, "w");
    flux = fopen("./Outputs/sup_flujo_P.out", "w");

    //Caso Perturbado: (Falta modificar según valor de gam)
    //sprintf(frs, "./Outputs/r_P_gam=%.3d.out", gam);
    //sprintf(frgcs, "./Outputs/rgc_P_gam=%.3d.out", gam);

    // fr = fopen("./Outputs/Perturbado/r_perturbado_P.out", "w");
    // frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P.out", "w");
    // flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_P.out", "w");

    cond_i(r, v, rgc, vpar, gam); //CUIDADO, ahora no tengo el amu definido.

    cuentaD(vpar[0], gam);

    printf("ri = %.6f, qi = %.6f, zi = %.6f, vri = %.6f, vqi = %.6f, vzi = %.6f y gam = %.3f.\n", r[0], r[1], r[2], v[0], v[1], v[2], gam);
    printf("\n");

    integrador(r, v, rgc, vpar, gam);
    PROC(r, v, rgc, vpar, rp, vp, rgcp, vparp, tp);

    for (i=0; i<nprint; i++) {
      fprintf(fr, "%f %f %f %f %f %f %f\n", tp[i], rp[3*i], rp[3*i+1], rp[3*i+2], vp[3*i], vp[3*i+1], vp[3*i+2]);
      fprintf(frgc, "%f %f %f %f %f\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i]);
    }

    fclose(fr);
    fclose(frgc);
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

  free(frs);
  free(frgcs);
  free(r);
  free(v);
  free(rgc);
  free(vpar);
  free(rp);
  free(vp);
  free(rgcp);
  free(vparp);
  free(tp);

  fclose(flux);

  return 0;
}
