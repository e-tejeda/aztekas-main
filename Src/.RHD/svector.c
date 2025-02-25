#include"main.h"
    
void Source_Terms(double *s, double *u, gauge_ local_grid)
{
   double rho, p, vx1=0, vx2=0, vx3=0;

   rho = u[0];
   p   = u[1];

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

#if COORDINATES == CARTESIAN

   double t = local_grid.x[0];
   double X = local_grid.x[1];
   double Y = local_grid.x[2];
   double Z = local_grid.x[3];

   s[0] = 0.0;
   s[1] = 0.0;
   s[2] = 0.0;
   s[3] = 0.0;
   s[4] = 0.0;

#elif COORDINATES == CYLINDRICAL

   double t   = local_grid.x[0];
   double R   = local_grid.x[1];
   double Z   = local_grid.x[2];
   double phi = local_grid.x[3];

   s[0] = 0;
   s[1] = 0;
   s[2] = (rho * vx3 * vx3 + p)/R;
   s[3] = 0;
   s[4] = - rho * vx1 * vx3/R;

#elif COORDINATES == SPHERICAL

   double E;
   eos_ eos;
   double t     = local_grid.x[0];
   double r     = local_grid.x[1];
   double theta = local_grid.x[2];
   double phi   = local_grid.x[3];

   EoS(&eos,u,local_grid);

   E = 0.5*rho*(vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

   s[0] = -2.0*vx1*rho/r - vx2*rho*cos(theta)/(r*sin(theta));
   s[1] = -2.0*vx1*(E+p)/r -vx2*(E+p)*cos(theta)/(r*sin(theta));
   s[2] = -2.0*vx1*(rho*vx1)/r - vx2*(rho*vx1)*cos(theta)/(r*sin(theta)) +\
          rho*(vx2*vx2 + vx3*vx3)/r;
   s[3] = -2.0*vx1*(rho*vx2)/r - vx2*(rho*vx2)*cos(theta)/(r*sin(theta)) -\
          rho*vx1*vx2/r + rho*vx3*vx3*cos(theta)/(r*sin(theta));
   s[4] = -2.0*vx1*(rho*vx3)/r - vx2*(rho*vx3)*cos(theta)/(r*sin(theta)) -\
          rho*vx1*vx3/r - rho*vx2*vx3*cos(theta)/(r*sin(theta));

#endif

}
