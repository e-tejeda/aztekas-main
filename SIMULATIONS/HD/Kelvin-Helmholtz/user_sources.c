/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:13:33
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

void User_Source_Terms(double *a, double *u, int *I)
{
   double rho, p, vx1=0, vx2=0, vx3=0;
   double x[4];

   x[0] = grid.time;

#if DIM == 1

   x[1] = grid.X1[I[0]];
   x[2] = 0.0;
   x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   x[1] = grid.X1[I[0]];
   x[2] = grid.X2[I[1]];
   x[3] = 0.0;
   #if POLAR == TRUE
   x[2] = M_PI_2;
   #endif

#elif DIM == 3

   x[1] = grid.X1[I[0]];
   x[2] = grid.X2[I[1]];
   x[3] = grid.X3[I[2]];

#endif 

   rho = u[0];
   p = u[1];

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

   a[0] = 0.0;
   a[1] = 0.0;
   a[2] = 0.0;
   a[3] = 0.0;
   a[4] = 0.0;
}
