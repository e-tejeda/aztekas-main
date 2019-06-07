#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define periodic_x1       TRUE
#define reflective_x2max  TRUE
#define reflective_x2min  TRUE

#define RECONST   WENO5
#define FLUX      HLL
#define GRID      UNIFORM

double nl, pl, vx1l, vx2l, vx3l;
double nr, pr, vx1r, vx2r, vx3r;
double x_0;
