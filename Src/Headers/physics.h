typedef struct
{
   double e;
   double cs;
}eos_;

int funct_A(double *a, double *uu);
int funct_Q2U(double *q, double *u);

int GAUGE(double *a, double g1, double g2, double g3);

void Prim2Cons(double *q, double *u, double *x);
void Prim2Cons_All(double *u, double *q);
void Prim2FluxF(double *f, double *v, double *u, double *x);
void Prim2FluxG(double *f, double *v, double *u, double *x);
void Prim2FluxH(double *f, double *v, double *u, double *x);
void Source_Terms(double *s, double *uu, double *x);

void EoS_Ideal(eos_ *e, double *u, double *x);
