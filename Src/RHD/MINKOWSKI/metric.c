void Metric_Components(metric_ *g, double *x)
{
#if COORDINATES == CARTESIAN
   
   g->gd[0][0] = -1.0;
   g->gd[1][1] =  1.0;
   g->gd[2][2] =  1.0;
   g->gd[3][3] =  1.0;

   g->gu[0][0] = -1.0;
   g->gu[1][1] =  1.0;
   g->gu[2][2] =  1.0;
   g->gu[3][3] =  1.0;

   g->yd[0][0] =  1.0;
   g->yd[1][1] =  1.0;
   g->yd[2][2] =  1.0;

   g->yu[0][0] =  1.0;
   g->yu[1][1] =  1.0;
   g->yu[2][2] =  1.0;

#elif COORDINATES == CYLINDRICAL

   double R   = x[1];

   g->gd[0][0] = -1.0;
   g->gd[1][1] =  1.0;
   g->gd[2][2] =  1.0;
   g->gd[3][3] =  R*R;

   g->gu[0][0] = -1.0;
   g->gu[1][1] =  1.0;
   g->gu[2][2] =  1.0;
   g->gu[3][3] =  1/(R*R);

   g->yd[0][0] =  1.0;
   g->yd[1][1] =  1.0;
   g->yd[2][2] =  R*R;

   g->yu[0][0] =  1.0;
   g->yu[1][1] =  1.0;
   g->yu[2][2] =  1/(R*R);
#elif COORDINATES == SPHERICAL

   double r     = x[1];
   double theta = x[2];

   g->gd[0][0] = -1.0;
   g->gd[1][1] =  1.0;
   g->gd[2][2] =  r*r;
   g->gd[3][3] =  r*r*sin(theta)*sin(theta);

   g->gd[0][0] = -1.0;
   g->gd[1][1] =  1.0;
   g->gd[2][2] =  1/(r*r);
   g->gd[3][3] =  1/(r*r*sin(theta)*sin(theta));

#endif
}
