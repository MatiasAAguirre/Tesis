
#include <math.h>
#include <stdio.h>

int main(void)
{
    int iter;
    double E, m, q;
    double Omega, v0, a,B0,frec,R0,vter,TkeV;
    double Gamma,tsim,frec_m;

    B0=2.5; // Tesla

    R0=1.71, // radio mayor en metros

           E=70.e3*1.602e-19;    // 8keV Energy
      // E=3.5e6*1.6022e-19;    // 3.5MeV Energy
    // m=58.7*1.66e-27;    // 58.7 u Mass Ni

    m=2.0*1.6726e-27;           // D
    //    q=26*1.6e-19;       // 26e Charge Ni particle
    
    //    q=10.*1.6e-19;       // 26e Charge Ni particle
    q=1.0*1.6022e-19;          // D
    Omega=q*B0/m;       // Cyclotronic frequency

     frec = 5000; //frecuencia del modo 5kHz
   
    //    v0= 2*3.14159*R0*frec; //velocity of an ion rotating at frec in the Rm pitch=1

     v0=sqrt(2*E/m);     // Initial speed
 
     // E = 0.5*m*v0*v0/(1e3*1.602e-19); //energy in keV

    a=0.5;              // Minor radius
  
    Gamma=v0/(Omega*a); // Nondimensional parameter

    //    TkeV=10.0;
    // vter= sqrt(2.*TkeV*1000.*1.602e-19/m);
    // printf("vtermica (m/s) %f  \n",vter);

    printf("frec %e \n",2*3.14159*frec/Omega);
    printf("gamma %f omega %e gamma/omega %e \n",Gamma, Omega, Gamma/Omega);
    printf("v0(m/s) %f E(keV) %f  %f \n",v0,E/1.602e-19/1000., v0*100.e-6);


    // tsim=0.05*Omega;
    //     tsim=2.27e+05;
    // tsim=2500000*0.16/Omega;  //3.3406 ms
    tsim=400000/Omega; 

    printf("tiempo de simulacion  %f ms \n",tsim*1000.);

    // tsim=10.975863*0.001*Omega;
     //frecuencia del modo
     // frec=3.8e-4*Omega;
    //    frec= 5.539218e+04; //(frec en rad/s)

       

    printf("radio de larmor  %e m \n",v0/Omega);
    printf("omega ciclotron %e 1/s frec modo (2,1) %e w %e \n",Omega,frec/Omega,2*3.14159*frec/Omega);
    printf("gamma= %lf \n",Gamma);
   
    
    printf("cuantas vueltas da en un tiempo de simulacion = %e  \n",v0*tsim/(2*3.14159*R0));
    
    printf("cuanto tarda en dar una vuelta: %e ms \n",((2*3.14159*(R0+0.5))/v0)*1000);
    
    //     printf("frecuancia modo (1,1) = %e \n",frec);    


    //printf("velocidad termica D  %e velocidad de rotacion del modo %e m \n",v0,frec*Omega*R0);
    //  printf("periodo de una part con  velocidad termica D  %e  \n",2*3.1415*R0/(v0));
    // printf("frec de una part con  velocidad termica D  %e  \n",v0/R0);


}
