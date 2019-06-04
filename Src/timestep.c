/**
 * @file timestep.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Time-step calculation.
 */

//Do not erase any of these libraries//
#include"main.h"

double TimeStep()
{
   int i, j, k;
   double dtmin;
   double c, dt, cmax;
   double r;

   dtmin = 1000000;

#if DIM == 1
   for(i = gc; i <= Nx1-gc; i++)
   {
      c = sqrt(K*U(1,i) / (U(0,i)));
      dtmin = MIN(dx1/(fabs(U(2,i)) + fabs(c)),dtmin);
   }

#elif DIM == 2
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         c = sqrt(K*U(1,i,j) / (U(0,i,j)));
         dtmin = MIN(dx1/(fabs(U(2,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(3,i,j)) + fabs(c)),dtmin);
      }
   }
   
#elif DIM == 4
   
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         c = sqrt(K*U(1,i,j) / (U(0,i,j)));
         dtmin = MIN(dx1/(fabs(U(2,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx1/(fabs(U(4,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(3,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(4,i,j)) + fabs(c)),dtmin);
      }
   }
  
#elif DIM == 3
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            c = sqrt(K*U[c3(1,i,j,k)] / (U[c3(0,i,j,k)]));
            dtmin = MIN(dx1/(fabs(U[c3(2,i,j,k)]) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U[c3(3,i,j,k)]) + fabs(c)),dtmin);
            dtmin = MIN(dx3/(fabs(U[c3(4,i,j,k)]) + fabs(c)),dtmin);
         }
      }
   }

#endif

#if PHYSICS == HD //HD
   dt = cou*dtmin;
#elif PHYSICS == RHD //RHD
   #if DIM == 1
   dt = cou*MIN(dx1,1000);
   #elif DIM == 2 || DIM == 4
   dt = cou*MIN(dx1,dx2);
   #endif 
#endif

   return dt;
}
