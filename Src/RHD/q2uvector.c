#include"main.h"

int funct_Q2U(double *u, double *q)
{
   int i, j, k;
   int a, b;
   metric_ m;
   double D, t, sd[3];
   double x[4];
   double W = 1, SS = 0;
   double theta, theta_0;
   double f, derf, h, derh;

#if DIM == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      D     = q(0,i);
      t     = q(1,i);
      sd[0] = q(2,i);
      sd[1] = 0;
      sd[2] = 0;

      x[0] = grid.time;
      x[1] = grid.X1[i];
      x[2] = 0;
      x[3] = 0;
      #if COORDINATES == TRUE
      x[2] = M_PI_2;
      #endif

      Metric_Components(&m,x);

      for(a = 0; a < 3; a++)
      {
         for(b = 0; b < 3; b++)
         {
            SS += m.yu[a][b] * sd[a] * sd[b];
         }
      }

      theta_0 = U(1,i)/U(0,i);
      f = 1.0;      

      while (fabs(f) > 0.0000001)
      {
         #if EOS == IDEAL
         h    = 1.0 + (K / (K - 1.0)) * theta_0;
         derh = K / (K - 1.0);
         #endif
         W  = sqrt(1.0 + SS / pow(D*h,2.0));

         f    = h * W - (theta_0 / W) - (t / D) - 1.0;
         derf = (1.0 / W) * (derh - 1.0 - theta_0 * ((W * W - 1.0) / (W * W)) * (derh / h));

         theta = theta_0 - f / derf;
         theta_0 = theta;
      }
   
      #if EOS == IDEAL
      h = 1.0 + (K / (K - 1.0)) * theta_0;
      #endif
      W = sqrt(1.0 + SS / pow(D*h,2.0));
      SS = 0;

      u(0,i) = D / W;
      u(1,i) = D*h*W - t - D;
      u(2,i) = sd[0] / (D * h * W);
   }

#elif DIM == 2
      
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         D     = q(0,i,j);
         t     = q(1,i,j);
         sd[0] = q(2,i,j);
         sd[1] = q(3,i,j);
         sd[2] = 0;

         x[0] = grid.time;
         x[1] = grid.X1[i];
         x[2] = grid.X2[j];
         x[3] = 0;
         #if POLAR == TRUE
         x[2] = M_PI_2;
         #endif

         Metric_Components(&m,x);

         for(a = 0; a < 3; a++)
         {
            for(b = 0; b < 3; b++)
            {
               SS += m.yu[a][b] * sd[a] * sd[b];
            }
         }

         theta_0 = U(1,i,j)/U(0,i,j);
         f = 1.0;

         while (fabs(f) > 0.0000001)
         {
            #if EOS == IDEAL
            h    = 1.0 + (K / (K - 1.0)) * theta_0;
            derh = K / (K - 1.0);
            #endif
            W  = sqrt(1.0 + SS / pow(D*h,2.0));

            f    = h * W - (theta_0 / W) - (t / D) - 1.0;
            derf = (1.0 / W) * (derh - 1.0 - theta_0 * ((W * W - 1.0) / (W * W)) * (derh / h));

            theta = theta_0 - f / derf;
            theta_0 = theta;
         }
   
         #if EOS == IDEAL
         h = 1.0 + (K / (K - 1.0)) * theta_0;
         #endif
         W = sqrt(1.0 + SS / pow(D*h,2.0));
         SS = 0;

         u(0,i,j) = D / W;
         u(1,i,j) = D*h*W - t - D;
         u(2,i,j) = sd[0] / (D * h * W);
         u(3,i,j) = sd[1] / (D * h * W);
      }
   }

#elif DIM = 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         D     = q(0,i,j);
         t     = q(1,i,j);
         sd[0] = q(2,i,j);
         sd[1] = q(3,i,j);
         sd[2] = q(4,i,j);

         x[0] = grid.time;
         x[1] = grid.X1[i];
         x[2] = grid.X2[j];
         x[3] = 0;
         #if POLAR == TRUE
         x[2] = M_PI_2;
         #endif

         Metric_Components(&m,x);

         for(a = 0; a < 3; a++)
         {
            for(b = 0; b < 3; b++)
            {
               SS += m.yu[a][b] * sd[a] * sd[b];
            }
         }

         theta_0 = U(1,i,j)/U(0,i,j);
         f = 1.0;

         while (fabs(f) > 0.0000001)
         {
            #if EOS == IDEAL
            h    = 1.0 + (K / (K - 1.0)) * theta_0;
            derh = K / (K - 1.0);
            #endif
            W  = sqrt(1.0 + SS / pow(D*h,2.0));

            f    = h * W - (theta_0 / W) - (t / D) - 1.0;
            derf = (1.0 / W) * (derh - 1.0 - theta_0 * ((W * W - 1.0) / (W * W)) * (derh / h));

            theta = theta_0 - f / derf;
            theta_0 = theta;
         }
   
         #if EOS == IDEAL
         h = 1.0 + (K / (K - 1.0)) * theta_0;
         #endif
         W = sqrt(1.0 + SS / pow(D*h,2.0));
         SS = 0;

         u(0,i,j) = D / W;
         u(1,i,j) = D*h*W - t - D;
         u(2,i,j) = sd[0] / (D * h * W);
         u(3,i,j) = sd[1] / (D * h * W);
         u(4,i,j) = sd[2] / (D * h * W);
      }
   }

#elif DIM = 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            D     = q(0,i,j,k);
            t     = q(1,i,j,k);
            sd[0] = q(2,i,j,k);
            sd[1] = q(3,i,j,k);
            sd[2] = q(4,i,j,k);
          
            x[0] = grid.time;
            x[1] = grid.X1[i];
            x[2] = grid.X2[j];
            x[3] = grid.X3[k];
          
            Metric_Components(&m,x);
          
            for(a = 0; a < 3; a++)
            {
               for(b = 0; b < 3; b++)
               {
                  SS += m.yu[a][b] * sd[a] * sd[b];
               }
            }
          
            theta_0 = U(1,i,j,k)/U(0,i,j,k);
            f = 1.0;
          
            while (fabs(f) > 0.0000001)
            {
               #if EOS == IDEAL
               h    = 1.0 + (K / (K - 1.0)) * theta_0;
               derh = K / (K - 1.0);
               #endif
               W  = sqrt(1.0 + SS / pow(D*h,2.0));
          
               f    = h * W - (theta_0 / W) - (t / D) - 1.0;
               derf = (1.0 / W) * (derh - 1.0 - theta_0 * ((W * W - 1.0) / (W * W)) * (derh / h));
          
               theta = theta_0 - f / derf;
               theta_0 = theta;
            }
          
            #if EOS == IDEAL
            h = 1.0 + (K / (K - 1.0)) * theta_0;
            #endif
            W = sqrt(1.0 + SS / pow(D*h,2.0));
            SS = 0;
          
            u(0,i,j,k) = D / W;
            u(1,i,j,k) = D*h*W - t - D;
            u(2,i,j,k) = sd[0] / (D * h * W);
            u(3,i,j,k) = sd[1] / (D * h * W);
            u(4,i,j,k) = sd[2] / (D * h * W);
         }
      }
   }

#endif

   return 0;
}
