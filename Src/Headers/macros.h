/* FUNCTIONS */
#define MIN(a,b) (((a)<(b))?(a):(b))                                            
#define MAX(a,b) (((a)>(b))?(a):(b)) 

/* LOGICAL */
#define TRUE  1
#define FALSE 0

/* PHYSICS */
#define HD    0
#define RHD   1

/* EOS */
#define IDEAL 0

/* COORDINATES */
#define CARTESIAN   0
#define CYLINDRICAL 1
#define SPHERICAL   2

/* GRID */
#define UNIFORM 0
#define LOGMESH 1

/* FLUX RECONSTRUCTOR */
#define HLL  0
#define HLLC 1

/* INTEGRATION METHOD */
#define STANDARD 0
#define PVRS     1

/* RECONSTRUCTION */
#define GODUNOV  0
#define MINMOD   1
#define MC       2
#define SUPERBEE 3
#define WENO5    4
