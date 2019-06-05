#include"main.h"

void Metric_Components(local_metric_ *g, double *x)
{ 
#if COORDINATES == CARTESIAN

   g->gamma[0][0] = 1.0;
   g->gamma[1][1] = 1.0;
   g->gamma[2][2] = 1.0;

#elif COORDINATES == CYLINDRICAL

   g->gamma[0][0] = 1.0;
   g->gamma[1][1] = 1.0;
   g->gamma[2][2] = 1.0 / (x[1]*x[1]);

#elif COORDINATES == SPHERICAL

   g->gamma[0][0] = 1.0;
   g->gamma[1][1] = 1.0 / (x[1]*x[1]);
   g->gamma[2][2] = 1.0 / (x[1]*x[1] * sin(x[2])*sin(x[2]));

#endif
}

double dety(double x1, double x2, double x3)
{
#if COORDINATES == CARTESIAN

   return 1;

#elif COORDINATES == CYLINDRICAL

   return x1;

#elif COORDINATES == SPHERICAL

   return x1*x1*sin(x2)*sin(x2);

#endif
}
