#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "datos.h"
#include "Funciones_Varias.h"

int main(int argc, char const *argv[]) {
  int i, j;
  double B[3], s_flux[1], vpar[nstep], rgc[3*nstep], rgcp[3*nprint], vparp[nprint], tp[nprint], amu = 0.0, r, z;
  double t = 0, rgc0[3], b1[3], e1[3], db1[8];
  FILE *frgc, *flux;

  clock_t begin = clock();

  //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
  //Caso sin perturbar:
  frgc = fopen("./Outputs/rgc_P.out", "w");
  flux = fopen("./Outputs/sup_flujo_P.out", "w");

  //Caso perturbado:
  // frgc = fopen("./Outputs/Perturbado/rgc_perturbado_P.out", "w");
  // flux = fopen("./Outputs/Perturbado/sup_flujo_perturbado_P.out", "w");

  for (i=0; i<3; i++) {
    B[i]=0.0;
    rgc0[i]=0.0;
  }
  vpar[0]=0.0;
  s_flux[0]=0.0;

  for (i=0; i<nstep; i++) {
    rgc[i]=0.0;
    rgc[i+1]=0.0;
    rgc[i+2]=0.0;
    vpar[i]=0.0;
  }
  for (i=0; i<nprint; i++) {
    rgc[i]=0.0;
    rgc[i+1]=0.0;
    rgc[i+2]=0.0;
    vpar[i]=0.0;
  }

  printf("\n");
  perturbacion(rgc0, t, b1, e1, db1);
  printf("\n");

  amu = cond_i(rgc, vpar);

  cuentaD(vpar[0]);

  printf("rgci = %.16f, qgci = %.16f, zgci = %.16f, vpari = %.16f y amu = %.16f.\n", rgc[0], rgc[1], rgc[2], vpar[0], amu);
  printf("\n");

  integrador(amu, rgc, vpar);
  PROC(amu, rgc, vpar, rgcp, vparp, tp);
  printf("\n");

  for (i=0; i<nprint; i++) {
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

  fclose(frgc);
  fclose(flux);

  clock_t end = clock();
  double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
  printf("Tiempo invertido en el cómputo: %f.\n", time_spent);

  return 0;
}
