#include"main.h"
    
void Source_Terms(double *s, double *u, int *I)
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

   double t = x[0];
   double X = x[1];
   double Y = x[2];
   double Z = x[3];

   s[0] = 0;
   s[1] = 0;
   s[2] = 0;
   s[3] = 0;
   s[4] = 0;

#elif COORDINATES == CYLINDRICAL

   double t   = x[0];
   double R   = x[1];
   double Z   = x[2];
   double phi = x[3];

   double ChrRpp = -R;
   double ChrpRp = 1/R;

   double gpp = 1/(R*R);

   s[0] = 0;
   s[1] = 0;
   s[2] = - ChrRpp * (rho * vx3 * vx3 + gpp*p);
   s[3] = 0;
   s[4] = - 2 * ChrpRp * (rho * vx1 * vx3);

#elif COORDINATES == SPHERICAL

   double E;
   eos_ eos;
   double t     = x[0];
   double r     = x[1];
   double theta = x[2];
   double phi   = x[3];

   double Chrrtt = -r;
   double Chrrpp = -r*sin(theta)*sin(theta);
   double Chrtrt = 1/r;
   double Chrtpp = -sin(theta)*cos(theta);
   double Chrprp = 1/r;
   double Chrptp = cos(theta)/sin(theta);

   double grr = 1.0;
   double gtt = 1/(r*r);
   double gpp = 1/(r*r*sin(theta)*sin(theta));

   double dety   = r*r*sin(theta);
   double drdety = 2*r*sin(theta);
   double dtdety = r*r*cos(theta);
   double dpdety = 0.0;


   #if EOS == IDEAL
   EoS_Ideal(&eos,u,x);
   #endif

   E = 0.5*rho*(vx1*vx1/grr + vx2*vx2/gtt + vx3*vx3/gpp) + rho*eos.e;
/*
   s[2] = -2.0*(rho*vx1)*vx1/r - (rho*vx1)*vx2*cos(theta)/(r*sin(theta)) + rho*(vx2*vx2 + vx3*vx3)/r;
   s[3] = -2.0*(rho*vx2)*vx1/r - (rho*vx2)*vx2*cos(theta)/(r*sin(theta)) - rho*(vx2*vx1)/r + rho*vx3*vx3*cos(theta)/(r*sin(theta)) ;
   s[4] = -2.0*(rho*vx3)*vx1/r - (rho*vx3)*vx2*cos(theta)/(r*sin(theta)) - rho*(vx3*vx1)/r + rho*vx2*vx3*cos(theta)/(r*sin(theta)) ;
*/
/*
   s[0] = - (rho*vx1)*(drdety/dety) - (rho*vx2)*(dtdety/dety);
   s[1] = - (E+p)*vx1*(drdety/dety) - (E+p)*vx2*(dtdety/dety);
   s[2] = - Chrrtt * (rho * vx2 * vx2 + gtt*p)\
          - Chrrpp * (rho * vx3 * vx3 + gpp*p)\
          - (rho*vx1*vx1 + grr*p)*(drdety/dety)\
          - (rho*vx1*vx2)*(dtdety/dety);
   //s[3] = - 2 * Chrtrt * (rho * vx1 * vx2)\
          - Chrtpp * (rho * vx3 * vx3 + gpp*p)\
          - (rho*vx2*vx1)*(drdety/dety)\
          - (rho*vx2*vx2 + gtt*p)*(dtdety/dety);
   //s[4] = - 2 * Chrprp * (rho * vx1 * vx3)\
          - 2 * Chrptp * (rho * vx2 * vx3)\
          - (rho*vx3*vx1)*(drdety/dety)\
          - (rho*vx3*vx2)*(dtdety/dety);
*/          
   s[0] = 0.0;
   s[1] = 0.0;
   s[2] = - Chrrtt * (rho * vx2 * vx2 + gtt*p)\
          - Chrrpp * (rho * vx3 * vx3 + gpp*p);
   s[3] = - 2 * Chrtrt * (rho * vx1 * vx2)\
          - Chrtpp * (rho * vx3 * vx3 + gpp*p);
   s[4] = - 2 * Chrprp * (rho * vx1 * vx3)\
          - 2 * Chrptp * (rho * vx2 * vx3);

#endif

}
