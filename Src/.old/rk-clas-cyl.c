/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"
#include"param.h"

int RK1D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i;
   int I[3];
   double Dx1 = dx1;
   double Dt  = dt;
   double L[eq+1], F[eq+1];
   double UU[eq+1];
   double g[6], g1p[6], g1m[6];
   vec_ v;
   lim_ l;

   for(i = gc; i <= Nx1-gc; i++)
   {
      I[0] = i;
      r = limiter;
 
      RECONST1D(u,r,&l,I);
      FLUX1D(&v,&l,I);
 
      for(n = 0; n < eq; n++)
      {
         F[n] = (X1p[i]*v.Fp[n] - X1m[i]*v.Fm[n])/(X1[i]*Dx1) - \
         v.S[n];
      }

      switch(order)
      {
         case 1:
            for(n = 0; n < eq; n++)
            {
               q1[c1(n,i)] = q[c1(n,i)] - (Dt)*(F[n]);
            }
         break;

         case 2:
            for(n = 0; n < eq; n++)
            {
               q2[c1(n,i)] = 0.5*(q1[c1(n,i)] + q[c1(n,i)] - (Dt)*F[n]);
            }
         break;
      }
   }

   return 0;
}

int RK2D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i, j;
   int I[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dt  = dt;
   double L[eq+1], F[eq+1];
   double g[6], g1p[6], g1m[6];
   double g2p[6], g2m[6];
   vec_ v;
   lim_ l;

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         I[0] = i;
         I[1] = j;
         r = limiter;

         RECONST2D(u,r,&l,I);
         FLUX2D(&v,&l,I);

         for(n = 0; n < eq; n++)
         {
            F[n] = (X1p[i]*v.Fp[n] - X1m[i]*v.Fm[n])/(X1[i]*Dx1) + \
            (v.Gp[n] - v.Gm[n])/(Dx2) - \
            v.S[n];
         }

         switch(order)
         {
            case 1:
               for(n = 0; n < eq; n++)
               {
                  q1[c2(n,i,j)] = q[c2(n,i,j)] - (Dt)*(F[n]);
               }
            break;

            case 2:
               for(n = 0; n < eq; n++)
               {
                  q2[c2(n,i,j)] = 0.5*(q1[c2(n,i,j)] + q[c2(n,i,j)] - (Dt)*F[n]);
               }
            break;
         }
      }
   }

   return 0;
}

int RK3D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i, j, k;
   int I[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dx3 = dx3;
   double Dt  = dt;
   double g[6], g1p[6], g1m[6];
   double g2p[6], g2m[6];
   double g3p[6], g3m[6];
   double L[eq+1], F[eq+1];
   vec_ v;
   lim_ l;

   #pragma omp parallel shared(u,Dt,Dx1,order,r) private(n,L,F,v,l,i,j,k) num_threads(4)
   {
      #pragma omp for
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            for(k = gc; k <= Nx3-gc; k++)
            {
               I[0] = i;
               I[1] = j;
               I[2] = k;
               r = limiter;
   
               RECONST3D(u,r,&l,I);
               FLUX3D(&v,&l,I);

               for(n = 0; n < eq; n++)
               {
                  F[n] = (X1p[i]*v.Fp[n] - X1m[i]*v.Fm[n])/(X1[i]*Dx1) + \
                  (v.Gp[n] - v.Gm[n])/(Dx2) + \
                  (v.Hp[n] - v.Hm[n])/(X1[i]*Dx3) - \
                  v.S[n];
               }

               switch(order)
               {
                  case 1:
                     for(n = 0; n < eq; n++)
                     {
                        q1[c3(n,i,j,k)] = q[c3(n,i,j,k)] - (Dt)*L[n];
                     }
                  break;

                  case 2:
                     for(n = 0; n < eq; n++)
                     {
                        q2[c3(n,i,j,k)] = 0.5*(q1[c3(n,i,j,k)] + \
                        q[c3(n,i,j,k)] - (Dt)*L[n]);
                     }
                  break;
               }
            }
         }
      }
   }

   return 0;
}
