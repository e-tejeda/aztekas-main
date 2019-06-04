#include"main.h"

void Gauge_Components(gauge_ *g, double *x)
{
#if COORDINATES == CARTESIAN
   
   g->lapse = 1.0;
   g->dety  = 1.0;
   g->detg  = g->lapse * g->dety;

   g->beta[0] = 0.0;
   g->beta[0] = 0.0;
   g->beta[0] = 0.0;

#endif
}
