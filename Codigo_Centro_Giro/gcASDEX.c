#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j, nstepf, jstep=((float)nstep)/((float)nprint);
  double B[3], s_flux[1], amu = 0.0, r, z;
  double t = 0, rgc0[3], b1[3], e1[3], db1[8];
  double *rgc, *vpar, *rgcp, *vparp, *tp;
  char part[3];
  FILE *frgc, *flux;

  clock_t begin = clock();

  //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
  //Caso sin perturbar:
  // frgc = fopen("./Outputs/rgc_P.out", "w");
  // flux = fopen("./Outputs/sup_flujo_P.out", "w");

  //Caso perturbado:
  frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P_h=0.001_dt=0.1.out", "w");
  flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_P.out", "w");

  rgc = (double *)malloc(3*nstep*sizeof(double));
  vpar = (double *)malloc(nstep*sizeof(double));

  rgcp = (double *)malloc(3*nprint*sizeof(double));
  vparp = (double *)malloc(nprint*sizeof(double));

  tp = (double *)malloc(nprint*sizeof(double));

  for (i=0; i<3; i++) {
    B[i]=0.0;
    rgc0[i]=0.0;
  }
  vpar[0]=0.0;
  s_flux[0]=0.0;

  for (i=0; i<nstep; i++) {
    rgc[3*i]=0.0;
    rgc[3*i+1]=0.0;
    rgc[3*i+2]=0.0;
    vpar[i]=0.0;
  }
  for (i=0; i<nprint; i++) {
    tp[i]=0.0;
    rgcp[3*i]=0.0;
    rgcp[3*i+1]=0.0;
    rgcp[3*i+2]=0.0;
    vparp[i]=0.0;
  }

  printf("\n");
  perturbacion(rgc0, t, b1, e1, db1);
  printf("\n");

  amu = cond_i(rgc, vpar);

  printf("rgci = %.16f, qgci = %.16f, zgci = %.16f, vpari = %.16f y amu = %.16f.\n", rgc[0], rgc[1], rgc[2], vpar[0], amu);
  printf("\n");

  nstepf = integrador(amu, rgc, vpar);
  PROC(amu, nstepf, rgc, vpar, rgcp, vparp, tp);
  printf("nstepf = %d\n", nstepf);

  for (i=0; i<(int)(nstepf/jstep); i++) {
    fprintf(frgc, "%f %f %f %f %f\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i]);
  }

  for(i=0; i<101; i++) {
    r = 2.0+3.0*(float)i/100;
    for(j=0; j<101; j++) {
      z = -2.0+4.0*(float)j/100;
      B_Asdex(r,z,B,s_flux);
      fprintf(flux, "%f %f %f\n", r, z, s_flux[0]);
    }
    fprintf(flux, "\n");
  }

  free(rgc);
  free(vpar);
  free(rgcp);
  free(vparp);
  free(tp);

  fclose(frgc);
  fclose(flux);

  clock_t end = clock();
  double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
  printf("Tiempo invertido en el cómputo: %f s.\n", time_spent);

  return 0;
}
