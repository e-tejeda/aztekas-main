#include"main.h"
    
int funct_Q2U(double *a, double *u)
{
   int i, j, k;
   double D, E, S1, S2, S3;
   
#if DIM == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      D  = u[c1(0,i)];
      E  = u[c1(1,i)];
      S1 = u[c1(2,i)];
      S2 = 0;
      S3 = 0;
 
      a[c1(0,i)] = D;
      #if EOS == IDEAL
      a[c1(1,i)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
      #endif
      a[c1(2,i)] = S1/D;
      a[c1(3,i)] = S2/D;
      a[c1(4,i)] = S3/D;
   }

#elif DIM == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = u[c2(0,i,j)];
         E  = u[c2(1,i,j)];
         S1 = u[c2(2,i,j)];
         S2 = u[c2(3,i,j)];
         S3 = 0;
 
         a[c2(0,i,j)] = D;
         #if EOS == IDEAL
         a[c2(1,i,j)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #endif
         a[c2(2,i,j)] = S1/D;
         a[c2(3,i,j)] = S2/D;
         a[c2(4,i,j)] = S3/D;
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = u[c2(0,i,j)];
         E  = u[c2(1,i,j)];
         S1 = u[c2(2,i,j)];
         S2 = u[c2(3,i,j)];
         S3 = u[c2(4,i,j)];
 
         a[c2(0,i,j)] = D;
         #if EOS == IDEAL
         a[c2(1,i,j)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #endif
         a[c2(2,i,j)] = S1/D;
         a[c2(3,i,j)] = S2/D;
         a[c2(4,i,j)] = S3/D;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            D  = u[c3(0,i,j,k)];
            E  = u[c3(1,i,j,k)];
            S1 = u[c3(2,i,j,k)];
            S2 = u[c3(3,i,j,k)];
            S3 = u[c3(4,i,j,k)];
 
            a[c3(0,i,j,k)] = D;
            #if EOS == IDEAL
            a[c3(1,i,j,k)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
            #endif
            a[c3(2,i,j,k)] = S1/D;
            a[c3(3,i,j,k)] = S2/D;
            a[c3(4,i,j,k)] = S3/D;
         }
      }
   }

#endif
   
   return 0;
}
