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
      #if SPHERICAL == TRUE
      x[2] = M_PI_2;
      #endif
    
      P[0] = u[c1(0,i)];
      P[1] = u[c1(1,i)];
      P[2] = u[c1(2,i)];
      P[3] = 0.0;
      P[4] = 0.0;

      #if EOS == IDEAL
      EoS_Ideal(&eos,P,x);
      #endif
 
      q[c1(0,i)] = P[0];
      q[c1(1,i)] = 0.5*P[0]*(P[2]*P[2] + P[3]*P[3] + P[4]*P[4]) + P[0]*eos.e;
      q[c1(2,i)] = P[0]*P[2];
      q[c1(3,i)] = P[0]*P[3];
      q[c1(4,i)] = P[0]*P[4];
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
       
         P[0] = u[c2(0,i,j)];
         P[1] = u[c2(1,i,j)];
         P[2] = u[c2(2,i,j)];
         P[3] = u[c2(3,i,j)];
         P[4] = 0.0;

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif
 
         q[c2(0,i,j)] = P[0];
         q[c2(1,i,j)] = 0.5*P[0]*(P[2]*P[2] + P[3]*P[3] + P[4]*P[4]) + P[0]*eos.e;
         q[c2(2,i,j)] = P[0]*P[2];
         q[c2(3,i,j)] = P[0]*P[3];
         q[c2(4,i,j)] = P[0]*P[4];
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
       
         P[0] = u[c2(0,i,j)];
         P[1] = u[c2(1,i,j)];
         P[2] = u[c2(2,i,j)];
         P[3] = u[c2(3,i,j)];
         P[4] = u[c2(4,i,j)];

         #if EOS == IDEAL
         EoS_Ideal(&eos,P,x);
         #endif
 
         q[c2(0,i,j)] = P[0];
         q[c2(1,i,j)] = 0.5*P[0]*(P[2]*P[2] + P[3]*P[3] + P[4]*P[4]) + P[0]*eos.e;
         q[c2(2,i,j)] = P[0]*P[2];
         q[c2(3,i,j)] = P[0]*P[3];
         q[c2(4,i,j)] = P[0]*P[4];
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
        
            P[0] = u[c2(0,i,j)];
            P[1] = u[c2(1,i,j)];
            P[2] = u[c2(2,i,j)];
            P[3] = u[c2(3,i,j)];
            P[4] = u[c2(4,i,j)];
          
            #if EOS == IDEAL
            EoS_Ideal(&eos,P,x);
            #endif
          
            q[c2(0,i,j)] = P[0];
            q[c2(1,i,j)] = 0.5*P[0]*(P[2]*P[2] + P[3]*P[3] + P[4]*P[4]) + P[0]*eos.e;
            q[c2(2,i,j)] = P[0]*P[2];
            q[c2(3,i,j)] = P[0]*P[3];
            q[c2(4,i,j)] = P[0]*P[4];
         }
      }
   }

#endif
}
