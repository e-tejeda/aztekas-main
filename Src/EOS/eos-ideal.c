#include"main.h"

void EoS_Ideal(double *e, double *u, double *x)
{
   double rho, p;
   rho = u[0];
   p   = u[1];

   *e = p / (rho * (K - 1.0));
}

void Sound_Speed(double *cs, double *u, double *x)
{
   double rho, p;
   rho = u[0];
   p   = u[1];

   *cs = sqrt(K * p / rho);
}
