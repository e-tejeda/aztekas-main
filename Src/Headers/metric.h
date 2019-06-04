typedef struct
{
   double gd[4][4];
   double gu[4][4];

   double yd[3][3];
   double yu[3][3];
}metric_;

typedef struct
{
   double lapse;
   double dety;
   double detg;
   double beta[3];
}gauge_;

void Metric_Components(metric_ *g, double *x);
void Gauge_Components(gauge_ *g, double *x);
