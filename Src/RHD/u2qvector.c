#include"main.h"

void Prim2Cons_All(double *q, double *u)
{
   int i, j, k;
   int a, b;
   eos_ eos;
   metric_ m;
   double W = 1, VV = 0, h;
   double rho, p, vd[3];
   double D, tau, Sd[3];
   double x[4], P[eq+1];

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      rho   = u[c1(0,i)];
      p     = u[c1(1,i)];
      vd[0] = u[c1(2,i)];
      vd[1] = 0;
      vd[2] = 0;

      x[0] = time;
      x[1] = X1[i];
      x[2] = 0.0;
      x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      x[2] = M_PI_2;
      #endif

      P[0] = rho;
      P[1] = p;

      #if EOS == IDEAL
      EoS_Ideal(&eos,P,x);
      #endif

      /* Specific enthalpy */
      h = 1.0 + eos.e + p/rho;

      /* Lorentz factor */
      Metric_Components(&m,x);

      for(a = 0; a < 3; a++)
      {
         for(b = 0; b < 3; b++)
         {
            VV += m.yu[a][b] * vd[a] * vd[b];
         }
      }

      W  = 1.0/sqrt(1.0 - VV);
      VV = 0.0;

      D   = rho * W;
      tau = rho * h * W * W - p - D;
      for(a = 0; a < 3; a++)
      {
         Sd[a] = rho * h * W * W * vd[a];
      }

      q[c1(0,i)] = D;
      q[c1(1,i)] = tau;
      q[c1(2,i)] = Sd[0];

      W = 1/sqrt((-pow(vd[2],2.0))-pow(vd[1],2.0)-pow(vd[0],2.0)+1);
      h = (K*p+(K-1)*rho)/((K-1)*rho);

      q[c1(0,i)] = W*rho;
      q[c1(1,i)] = (pow(W,2.0)*h-W)*rho-p;
      q[c1(2,i)] = pow(W,2.0)*h*rho*vd[0];
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         rho   = u[c2(0,i,j)];
         p     = u[c2(1,i,j)];
         vd[0] = u[c2(2,i,j)];
         vd[1] = u[c2(3,i,j)];
         vd[2] = 0;

         x[0] = time;
         x[1] = X1[i];
         x[2] = X2[j];
         x[3] = 0.0;
         #if POLAR == TRUE
         x[2] = M_PI_2;
         #endif

         P[0] = rho;
         P[1] = p;

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif

         /* Specific enthalpy */
         h = 1.0 + eos.e + p/rho;

         /* Lorentz factor */
         Metric_Components(&m,x)

         for(a = 0; a < 3; a++)
         {
            for(b = 0; b < 3; b++)
            {
               VV += m.yu[a][b] * vd[a] * vd[b];
            }
         }

         W  = 1.0/sqrt(1.0 - VV);
         VV = 0.0;

         D   = rho * W;
         tau = rho * h * W * W - p - D;
         for(a = 0; a < 3; a++)
         {
            Sd[a] = rho * h * W * W * vd[a];
         }

         q[c2(0,i,j)] = D;
         q[c2(1,i,j)] = tau;
         q[c2(2,i,j)] = Sd[0];
         q[c2(3,i,j)] = Sd[1];
      }
   }
   
#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         rho   = u[c2(0,i,j)];
         p     = u[c2(1,i,j)];
         vd[0] = u[c2(2,i,j)];
         vd[1] = u[c2(3,i,j)];
         vd[2] = u[c2(4,i,j)];

         x[0] = time;
         x[1] = X1[i];
         x[2] = X2[j];
         x[3] = 0.0;

         P[0] = rho;
         P[1] = p;

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif

         /* Specific enthalpy */
         h = 1.0 + eos.e + p/rho;

         /* Lorentz factor */
         Metric_Components(&m,x)

         for(a = 0; a < 3; a++)
         {
            for(b = 0; b < 3; b++)
            {
               VV += m.yu[a][b] * vd[a] * vd[b];
            }
         }

         W  = 1.0/sqrt(1.0 - VV);
         VV = 0.0;

         D   = rho * W;
         tau = rho * h * W * W - p - D;
         for(a = 0; a < 3; a++)
         {
            Sd[a] = rho * h * W * W * vd[a];
         }

         q[c2(0,i,j)] = D;
         q[c2(1,i,j)] = tau;
         q[c2(2,i,j)] = Sd[0];
         q[c2(3,i,j)] = Sd[1];
         q[c2(4,i,j)] = Sd[2];
      }
   }
   
#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k+)
         {
            rho   = u[c2(0,i,j)];
            p     = u[c2(1,i,j)];
            vd[0] = u[c2(2,i,j)];
            vd[1] = u[c2(3,i,j)];
            vd[2] = u[c2(4,i,j)];
          
            x[0] = time;
            x[1] = X1[i];
            x[2] = X2[j];
            x[3] = X3[k];
          
            P[0] = rho;
            P[1] = p;
          
            #if EOS == IDEAL
            EoS_Ideal(&eos,P,x);
            #endif
          
            /* Specific enthalpy */
            h = 1.0 + eos.e + p/rho;
          
            /* Lorentz factor */
            Metric_Components(&m,x)
          
            for(a = 0; a < 3; a++)
            {
               for(b = 0; b < 3; b++)
               {
                  VV += m.yu[a][b] * vd[a] * vd[b];
               }
            }
          
            W = 1.0/sqrt(1.0 - VV);
          
            D   = rho * W;
            tau = rho * h * W * W - p - D;
            for(a = 0; a < 3; a++)
            {
               Sd[a] = rho * h * W * W * vd[a];
            }
          
            q[c3(0,i,j,k)] = D;
            q[c3(1,i,j,k)] = tau;
            q[c3(2,i,j,k)] = Sd[0];
            q[c3(3,i,j,k)] = Sd[1];
            q[c3(4,i,j,k)] = Sd[2];
         }
      }
   }

#endif
}
