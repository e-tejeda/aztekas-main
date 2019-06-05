#include"main.h"
    
void Prim2FluxH(double *f, double *v, double *u, double *x)
{
   double E;
   eos_ eos;
   local_metric_ g;
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
   EoS_Ideal(&eos,u,x);
#endif

   Metric_Components(&g,x);

   E = 0.5 * rho * (vx1*vx1/g.gamma[0][0] \
                  + vx2*vx2/g.gamma[1][1] \
                  + vx3*vx3/g.gamma[2][2]) + rho * eos.e;

   f[0] = rho * vx3;
   f[1] = vx3 * (E + p);
   f[2] = rho * vx1 * vx3;
   f[3] = rho * vx2 * vx3;
   f[4] = rho * vx3 * vx3 + g.gamma[2][2]*p;

   v[0] = vx3 - eos.cs;
   v[1] = vx3 + eos.cs;
   v[2] = vx3;
}
