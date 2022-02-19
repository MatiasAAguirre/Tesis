/* -------------------------------------
              Parametros Varios
--------------------------------------*/

const int nrnz = 160801;
const int npart = 1;
const int nr = 401;
const int nmode = -1;
const int nstep = 800000;//Valor Predeterminado: 1600000 800k=1ms-80k=0.1ms-8k=0.01ms (con dt=0.2)
const int nprint = 8000;//Valor Predeterminado: 1600

// const double gam = 0.02306869; //20keV con B0=2.5T
// const double gam = 0.02883586; //20keV con B0=2T
const double gam = 0.04327451; //70keV con B0=2.5T
// const double gam = 0.04893608; //90keV con B0=2.5T
// const double gam = 0.06234974; //90keV con B0=2T
const double dr = 0.001;
const double dq = 0.001;
const double dz = 0.001;
const double dt = 0.2; //Valor predeterminado 0.2
const double rp0 = 1.1;
const double zp0 = -0.9;
const double hr = 0.00275;
const double hz = 0.0045;
const double omega = 0.0;
// const double omega = -0.0008;
// const double omega = -0.000375;
const double BT0 = 17.815757116271065;
const double bscale = 0.0;
// const double bscale = 0.0003;
const double escale = 0.0;
// const double escale = 0.01*omega;
// const double escale = 0.001*omega;
