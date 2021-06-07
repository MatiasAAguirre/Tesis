#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "datos_CG.h"
#include "Funciones_Varias_CG.h"

int main(int argc, char const *argv[]) {
  int i, j;
  double B[3], s_flux[1], vpar[nstep], rgc[3*nstep], rgcp[3*nprint], vparp[nprint], tp[nprint], amu = 0.0, r, z;
  double t = 0, gam, rgc0[3], b1[3], e1[3], db1[8], pitch;
  char *frgcs;
  FILE *frgc, *flux;

  frgcs = malloc(11*sizeof("./Outputs/prog_rgc_P_gam=0.083.out"));

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

  for (j=0; j<11; j++) {
    gam = 0.00327451 + 0.08*j/10;

    pitch = 0;

    //RECORDAR CAMBIAR EL DIRECTORIO DE DESTINO DEPENDIENDO DE QUE ESTOY HACIENDO.
    sprintf(frgcs, "./Outputs/prog_rgc_P_gam=%.3f.out", gam);

    frgc = fopen(frgcs, "w");

    amu = cond_i(rgc, vpar, gam);

    printf("rgci = %.16f, qgci = %.16f, zgci = %.16f, vpari = %.16f y amu = %.16f.\n", rgc[0], rgc[1], rgc[2], vpar[0], amu);
    printf("\n");

    integrador(amu, rgc, vpar, gam);
    for (i=0; i<nstep; i++) {
      pitch = pitch + vpar[i];
    }
    printf("Valor medio pitch = %0.2f\n", pitch/nstep);
    PROC(amu, rgc, vpar, rgcp, vparp, tp);

    for (i=0; i<nprint; i++) {
      fprintf(frgc, "%f %f %f %f %f\n", tp[i], rgcp[3*i], rgcp[3*i+1], rgcp[3*i+2], vparp[i]);
    }

    fclose(frgc);
  }

  flux = fopen("./Outputs/sup_flujo_P.out", "w");

  for(i=0; i<101; i++) {
    r = 2.0+3.0*(float)i/100;
    for(j=0; j<101; j++) {
      z = -2.0+4.0*(float)j/100;
      B_Asdex(r,z,B,s_flux);
      fprintf(flux, "%f %f %f\n", r, z, s_flux[0]);
    }
    fprintf(flux, "\n");
  }

  free(frgcs);

  fclose(flux);

  return 0;
}
