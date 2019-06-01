#include"main.h"
    
void Prim2Cons_All(double *q, double *u)
{
   int i, j, k;
   double P[eq+1];
   eos_ eos;
   double x[4];
   
#if DIM == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      x[0] = time;
      x[1] = X1[i];
      x[2] = 0.0;
      x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      x[2] = M_PI_2;
      #endif
    
      rho = u[c1(0,i)];
      p   = u[c1(1,i)];
      vx1 = u[c1(2,i)];
      vx2 = 0.0;
      vx3 = 0.0;

      P[0] = rho;
      P[1] = p;

      #if EOS == IDEAL
      EoS_Ideal(&eos,P,x);
      #endif
 
      q[c1(0,i)] = rho;
      q[c1(1,i)] = 0.5*rho*(vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
      q[c1(2,i)] = rho*vx1;
      q[c1(3,i)] = rho*vx2;
      q[c1(4,i)] = rho*vx3;
   }

#elif DIM == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         x[0] = time;
         x[1] = X1[i];
         x[2] = X2[j];
         x[3] = 0.0;
         #if POLAR == TRUE
         x[2] = M_PI_2;
         #endif
       
         rho = u[c2(0,i,j)];
         p   = u[c2(1,i,j)];
         vx1 = u[c2(2,i,j)];
         vx2 = u[c2(3,i,j)];
         vx3 = 0.0;

         P[0] = rho;
         P[1] = p;

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif
 
         q[c2(0,i,j)] = rho;
         q[c2(1,i,j)] = 0.5*rho*(vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
         q[c2(2,i,j)] = rho*vx1;
         q[c2(3,i,j)] = rho*vx2;
         q[c2(4,i,j)] = rho*vx3;
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         x[0] = time;
         x[1] = X1[i];
         x[2] = X2[j];
         x[3] = 0.0;
       
         rho = u[c2(0,i,j)];
         p   = u[c2(1,i,j)];
         vx1 = u[c2(2,i,j)];
         vx2 = u[c2(3,i,j)];
         vx3 = u[c2(4,i,j)];

         P[0] = rho;
         P[1] = p;

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif
 
         q[c2(0,i,j)] = rho;
         q[c2(1,i,j)] = 0.5*rho*(vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
         q[c2(2,i,j)] = rho*vx1;
         q[c2(3,i,j)] = rho*vx2;
         q[c2(4,i,j)] = rho*vx3;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            x[0] = time;
            x[1] = X1[i];
            x[2] = X2[j];
            x[3] = X3[k];
        
            rho = u[c3(0,i,j,k)];
            p   = u[c3(1,i,j,k)];
            vx1 = u[c3(2,i,j,k)];
            vx2 = u[c3(3,i,j,k)];
            vx3 = u[c3(4,i,j,k)];

            P[0] = rho;
            P[1] = p;

            #if EOS == IDEAL
            EoS_Ideal(&eos,P,x);
            #endif
          
            q[c3(0,i,j,k)] = rho;
            q[c3(1,i,j,k)] = 0.5*rho*(vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
            q[c3(2,i,j,k)] = rho*vx1;
            q[c3(3,i,j,k)] = rho*vx2;
            q[c3(4,i,j,k)] = rho*vx3;
         }
      }
   }

#endif
}
