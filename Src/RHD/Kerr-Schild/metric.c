#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   local_grid->lapse = 0.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 0.0; 
   local_grid->gamma_con[0][1] = 0.0; 
   local_grid->gamma_con[0][2] = 0.0; 
   local_grid->gamma_con[1][0] = 0.0; 
   local_grid->gamma_con[1][1] = 0.0; 
   local_grid->gamma_con[1][2] = 0.0; 
   local_grid->gamma_con[2][0] = 0.0; 
   local_grid->gamma_con[2][1] = 0.0; 
   local_grid->gamma_con[2][2] = 0.0; 

   local_grid->dety = 0.0;

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);

   local_grid->lapse = 0.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 0.0;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 0.0;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 0.0;

   local_grid->dety = 0.0;

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 
   double phi   = local_grid->x[3];

   double Sigma, Delta, rho2, a, M;
   double cos2, sin2;
   double r2, a2;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;

   cos2  = pow(cos(theta),2.0);
   sin2  = pow(sin(theta),2.0);
   r2    = pow(r,2.0);
   a2    = pow(a,2.0);

   Delta = r2 - 2.0*M*r + a2;
   rho2  = r2 + a2*cos2;
   Sigma = pow(r2 + a2,2.0) - Delta*a2*sin2;

   local_grid->lapse = 1.0/sqrt(1.0 + 2.0*M*r/rho2);

   local_grid->beta_con[0] = (2.0*M*r/rho2)*(1.0/(1.0 + 2.0*M*r/rho2));
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   #if POLAR == FALSE

   local_grid->gamma_con[0][0] = (a2*(rho2 + 2.0*M*r)*sin2 + pow(rho2,2.0))/(rho2*(rho2 + 2.0*M*r));
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = a/rho2;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/rho2;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = a/rho2;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(rho2*sin2);

   #elif POLAR == TRUE

   local_grid->gamma_con[0][0] = (a2*(1.0 + 2.0*M*r/rho2)*sin2 + rho2)/(rho2*(1.0 + 2.0*M*r/rho2));
   local_grid->gamma_con[0][1] = a/rho2;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = a/rho2;
   local_grid->gamma_con[1][1] = 1.0/(rho2*sin2);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/rho2;

   #endif

   //local_grid->dety = rho2*sin(theta);
   local_grid->dety = sqrt(rho2*sin2*(2.0*M*r + rho2));

#endif
}

void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   der->dlapse[0] = 0.0;
   der->dlapse[1] = 0.0;
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 0.0;
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 0.0;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 0.0;

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);

   der->dlapse[0] = 0.0;
   der->dlapse[1] = 0.0;
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;


   der->dgam[0][0][0] = 0.0;
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 0.0;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 0.0;

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 
   double phi   = local_grid->x[3];

   double Sigma, Delta, rho2, a, M;
   double cos2, sin2;
   double r2, a2;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;

   cos2  = pow(cos(theta),2.0);
   sin2  = pow(sin(theta),2.0);
   r2    = pow(r,2.0);
   a2    = pow(a,2.0);

   Delta = r2 - 2.0*M*r + a2;
   rho2  = r2 + a2*cos2;
   Sigma = pow(r2 + a2,2.0) - Delta*a2*sin2;

   #if POLAR == FALSE

   der->dlapse[0] = M*(r2 - a2*cos2)*sqrt(rho2)/(rho2*pow(2.0*M*r + rho2, 3.0/2.0));
   der->dlapse[1] = -2.0*M*r*a2*sin(theta)*cos(theta)/(pow(rho2, 2.0)*pow(2.0*M*r/rho2 + 1.0, 3.0/2.0));
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 2.0*M*(a2*cos2 - r2)/pow(2.0*M*r + rho2, 2.0);
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 2.0*M*r*a2*sin(2.0*theta)/pow(2.0*M*r + rho2, 2.0);
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 2.0*M*(a2*cos2 - r2)/pow(rho2, 2.0);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 2.0*M*a*(r2 - a2*cos2)*sin2/pow(rho2, 2.0);
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 2.0*M*a*(-a2*cos2 + r2)*sin2/pow(rho2, 2.0);
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*(-M*a2*(-a2*cos2 + r2)*sin2 + r*pow(rho2, 2.0))*sin2/pow(rho2, 2.0);

   der->dgam[1][0][0] = 4*M*a2*r*sin(theta)*cos(theta)/pow(rho2, 2.0);
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = -2.0*a*(2.0*M*r*a2*sin2 + rho2*(2.0*M*r + rho2))*sin(theta)*cos(theta)/pow(rho2, 2.0);
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = -a2*sin(2.0*theta);
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = -2.0*a*(2.0*M*r*a2*sin2 + rho2*(2.0*M*r + rho2))*sin(theta)*cos(theta)/pow(rho2, 2.0);
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 2.0*(2.0*M*r*a2*(a2 + r2)*sin2 + rho2*(a2*(2.0*M*r + rho2)*sin2 + pow(rho2, 2.0)))*sin(theta)*cos(theta)/pow(rho2, 2.0);

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

   #elif POLAR == TRUE

   #endif

#endif
}
