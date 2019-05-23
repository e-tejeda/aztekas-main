#include"main.h"
    
void Prim2Cons(double *q, double *u, double *x)
{
   eos_ eos;
   metric_ g;
   double W, VV = 0;
   double rho, p, vd[3];
   double D, tau, Sd[3];
   rho = u[0];
   p   = u[1];

#if DIM == 1
   vd[0] = u[2];
   vd[1] = 0;
   vd[2] = 0;
#elif DIM == 2
   vd[0] = u[2];
   vd[1] = u[3];
   vd[2] = 0;
#elif DIM == 3 || DIM == 4
   vd[0] = u[2];
   vd[1] = u[3];
   vd[2] = u[4];
#endif

#if EOS == IDEAL
   EoS_Ideal(&eos,u,x);
#endif

   /* Specific enthalpy */
   h = 1.0 + eos.e + p/rho;

   /* Lorentz factor */
   Metric_Components(&g,x)

   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      {
         VV += g.yu[i][j] * vu[i] * vu[j];
      }
   }

   W = 1.0/sqrt(1.0 - VV);

   D   = rho * W;
   tau = rho * h * W * W - p - D;
   for(i = 0; i < 3; i++)
   {
      Sd[i] = rho * h * W * W * vd[i];
   }

   q[0] = D;
   q[1] = tau;
   q[2] = Sd[0];
   q[3] = Sd[1];
   q[4] = Sd[2];
}
