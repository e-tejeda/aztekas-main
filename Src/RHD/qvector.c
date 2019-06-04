#include"main.h"
    
void Prim2Cons(double *q, double *u, double *x)
{
   int i, j;
   eos_ eos;
   metric_ m;
   double W = 1, VV = 0, h;
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
   Metric_Components(&m,x);

   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      {
         VV += m.yu[i][j] * vd[i] * vd[j];
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

   W = 1/sqrt((-pow(vd[2],2.0))-pow(vd[1],2.0)-pow(vd[0],2.0)+1);
   h = (K*p+(K-1)*rho)/((K-1)*rho);

   q[0] = W*rho;
   q[1] = (pow(W,2.0)*h-W)*rho-p;
   q[2] = pow(W,2.0)*h*rho*vd[0];
   q[3] = pow(W,2.0)*h*rho*vd[1];
   q[4] = pow(W,2.0)*h*rho*vd[2];
}
