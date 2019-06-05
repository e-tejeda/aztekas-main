/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:50:45
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

void User_Source_Terms(double *a, double *uu, double *x)
{
   int i;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
#if DIM == 1
   u = uu[2];
#elif DIM == 2
   u = uu[2];
   v = uu[3];
#elif DIM == 3 || DIM == 4
   u = uu[2];
   v = uu[3];
   w = uu[4];
#endif

   a[0] = 0.0;
   a[1] = -n*v;
#if DIM == 1
   a[2] = 0.0;
#elif DIM == 2
   a[2] = 0.0;
   a[3] = -n;
#elif DIM == 3 || DIM == 4
   a[2] = 0.0;
   a[3] = -n;
   a[4] = 0.0;
#endif
}
