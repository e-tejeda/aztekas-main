#include"main.h"
    
void Prim2FluxG(double *f, double *v, double *u, double *x)
{
   int i, j;
   eos_ eos;
   metric_ m;
   gauge_ g;
   double W = 1, VV = 0, V, h;
   double rho, p, vd[3], vu[3];
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
   Gauge_Components(&g,x);

   vu[0] = m.yu[0][0]*vd[0] + m.yu[0][1]*vd[1] + m.yu[0][2]*vd[2];
   vu[1] = m.yu[1][0]*vd[0] + m.yu[1][1]*vd[1] + m.yu[1][2]*vd[2];
   vu[2] = m.yu[2][0]*vd[0] + m.yu[2][1]*vd[1] + m.yu[2][2]*vd[2];
   
   for(i = 0; i < 3; i++)
   {
      VV += vu[i] * vd[i];
   }

   W = 1.0/sqrt(1.0 - VV);
   V = vu[1] - g.beta[1]/g.lapse;

   D   = rho * W;
   tau = rho * h * W * W - p - D;
   for(i = 0; i < 3; i++)
   {
      Sd[i] = rho * h * W * W * vd[i];
   }

   f[0] = D * V;
   f[1] = tau * V + p * vu[1];
   f[2] = Sd[0] * V;
   f[3] = Sd[1] * V + p;
   f[4] = Sd[2] * V;

   double term1, term2, term3, root;

   term1 = 1 - VV*pow(eos.cs,2.0);
   term2 = 1 - pow(eos.cs,2.0);
   term3 = g.lapse / term1;
   root  = (1 - VV)*(m.yu[1][1]*term1 - vu[1]*vu[1]*term2);

   v[0] = term3 * (vu[1]*term2 + eos.cs*sqrt(root)) - g.beta[1];
   v[1] = term3 * (vu[1]*term2 - eos.cs*sqrt(root)) - g.beta[1];
   v[2] = g.lapse * vu[1] - g.beta[1];
}
