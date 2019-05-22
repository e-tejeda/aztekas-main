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
#include"main.h"

int Sources(double *u, vec_ *v, int *I)
{
   int n;
   double x[4];
   double default_S[eq + 1], user_S[eq + 1];

   x[0] = time;

#if DIM == 1

   x[1] = X1[I[0]];
   x[2] = 0.0;
   x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   x[1] = X1[I[0]];
   x[2] = X2[I[1]];
   x[3] = 0.0;
   #if POLAR == TRUE
   x[2] = M_PI_2;
   #endif

#elif DIM == 3

   x[1] = X1[I[0]];
   x[2] = X2[I[1]];
   x[3] = X3[I[2]];

#endif 

   Source_Terms(default_S,u,x);
   User_Source_Terms(user_S,u,x);

#if integration == 1
   funct_A(v->A,u,x);
#endif

   for(n = 0; n < eq; n++)
   {
      v->S[n] = default_S[n] + user_S[n];
   }

   return 0;
}

///////////////////////////////////////////////////////////////////////////

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I)
{
   int n;
   double *u, lr, ll;
   double up[eq + 1];
   double um[eq + 1];
   double lp[3];
   double lm[3];
   double dup[eq + 1];
   double dum[eq + 1];
   double x[4];

   x[0] = time;

#if DIM == 1

   x[1] = X1[I[0]];
   x[2] = 0.0;
   x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   x[1] = X1[I[0]];
   x[2] = X2[I[1]];
   x[3] = 0.0;

#elif DIM == 3

   x[1] = X1[I[0]];
   x[2] = X2[I[1]];
   x[3] = X3[I[2]];

#endif 

   if(pm == 1)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1p;
            x[1] = X1p[I[0]];
         break;

         case 'g':
            u = l->ux2p;
            x[2] = X2p[I[1]];
         break;

         case 'h':
            u = l->ux3p;
            x[3] = X3p[I[2]];
         break;
      }
   }
   if(pm == 0)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1m;
            x[1] = X1m[I[0]];
         break;

         case 'g':
            u = l->ux2m;
            x[2] = X2m[I[1]];
         break;

         case 'h':
            u = l->ux3m;
            x[3] = X3m[I[2]];
         break;
      }
   }

   for(n = 0; n < eq; n++)
   {
      f->up[n] = u[1*eq + n];
      f->um[n] = u[0*eq + n];
   }

   Prim2Cons(f->qp,f->up,x);
   Prim2Cons(f->qm,f->um,x);

   switch(flux)
   {
      case 'f':
         Prim2FluxF(f->fp,lp,f->up,x);
         Prim2FluxF(f->fm,lm,f->um,x);
      break;

      case 'g':
         Prim2FluxG(f->fp,lp,f->up,x);
         Prim2FluxG(f->fm,lm,f->um,x);
      break;

      case 'h':
         Prim2FluxH(f->fp,lp,f->up,x);
         Prim2FluxH(f->fm,lm,f->um,x);
      break;
   }

   lr = MAX(lp[0],lp[1]);
   lr = MAX(lr,lp[2]);
   lr = MAX(0.0,lr);
   ll = MAX(lm[0],lm[1]);
   ll = MAX(ll,lm[2]);
   ll = MAX(0.0,ll);

   f->lp = MAX(lr,ll);

   lr = MIN(lp[0],lp[1]);
   lr = MIN(lr,lp[2]);
   lr = MIN(0.0,lr);
   ll = MIN(lm[0],lm[1]);
   ll = MIN(ll,lm[2]);
   ll = MIN(0.0,ll);

   f->lm = MIN(lr,ll);

   return 0;
}

///////////////////////////////////////////////////////////////////////////
