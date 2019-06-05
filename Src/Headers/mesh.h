typedef struct
{
   double time;
   double *X1;
   double *X1p, *X1m;
   double *X2;
   double *X2p, *X2m;
   double *X3;
   double *X3p, *X3m;
}grid_;

double dx1, dx2, dx3;
double dt;
double tmax, timefile, cou;

grid_ grid;

int Mesh(); 

double TimeStep();

