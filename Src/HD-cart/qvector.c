#include"main.h"
    
void funct_Q(double *q, double *u, double *x)
{
   double E, e;
   double rho, p, vx1=0, vx2=0, vx3=0;
   rho = u[0];
   p   = u[1];

#if DIM == 1
   vx1 = u[2];
#elif DIM == 2
   vx1 = u[2];
   vx2 = u[3];
#elif DIM == 3 || DIM == 4
   vx1 = u[2];
   vx2 = u[3];
   vx3 = u[4];
#endif

#if EOS == IDEAL
   EoS_Ideal(&e,u,x);
#endif

   E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho * e;

   q[0] = rho;
   q[1] = E;
   q[2] = rho*vx1;
   q[3] = rho*vx2;
   q[4] = rho*vx3;
}
