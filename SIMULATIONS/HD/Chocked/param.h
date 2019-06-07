#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}
#define outflow_x1max    TRUE
#define outflow_x1min    TRUE
#define reflective_x2max TRUE
#define reflective_x2min TRUE

//#define x1min_exc TRUE

#define RECONST         WENO5
#define FLUX            HLLC
#define GRID            UNIFORM //(UNIFORM, LOGMESH)

double gtheta(double th);

double density_0, pressure_0, velocity_0;
double rho_atm;
