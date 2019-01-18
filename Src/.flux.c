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

int FLUX1D(vec_ *v, lim_ *l, int *I)
{
   int n;

   flx_ f;

   AMATRIX1D(l->ux,v,I);
#if riemann == 1 //HLL
   VECTOR(1,'f',l,&f,I);
   HLL(v->Fp,&f);
   VECTOR(0,'f',l,&f,I);
   HLL(v->Fm,&f);
#elif riemann == 2//HLLC
   VECTOR(1,'f',l,&f,I);
   HLLC(v->Fp,&f,1);
   VECTOR(0,'f',l,&f,I);
   HLLC(v->Fm,&f,1);
#endif

   return 0;
}

int FLUX2D(vec_ *v, lim_ *l, int *I)
{
   int n;

   flx_ f;

   AMATRIX2D(l->ux,v,I);
#if riemann == 1//HLL
   VECTOR(1,'f',l,&f,I);
   HLL(v->Fp,&f);
   VECTOR(1,'g',l,&f,I);
   HLL(v->Gp,&f);
   VECTOR(0,'f',l,&f,I);
   HLL(v->Fm,&f);
   VECTOR(0,'g',l,&f,I);
   HLL(v->Gm,&f);
#elif riemann == 2//HLLC
   VECTOR(1,'f',l,&f,I);
   HLLC(v->Fp,&f,1);
   VECTOR(1,'g',l,&f,I);
   HLLC(v->Gp,&f,2);
   VECTOR(0,'f',l,&f,I);
   HLLC(v->Fm,&f,1);
   VECTOR(0,'g',l,&f,I);
   HLLC(v->Gm,&f,2);
#endif

   return 0;
}

int FLUX3D(vec_ *v, lim_ *l, int *I)
{
   int n;

   flx_ f;

   AMATRIX3D(l->ux,v,I);
   VECTOR(1,'f',l,&f,I);
   HLL(v->Fp,&f);
   VECTOR(1,'g',l,&f,I);
   HLL(v->Gp,&f);
   VECTOR(1,'h',l,&f,I);
   HLL(v->Hp,&f);
   VECTOR(0,'f',l,&f,I);
   HLL(v->Fm,&f);
   VECTOR(0,'g',l,&f,I);
   HLL(v->Gm,&f);
   VECTOR(0,'h',l,&f,I);
   HLL(v->Hm,&f);

   return 0;
}

int HLL(double *F, flx_ *f)
{
   int n;
   double q;

   if(f->lm <= 0 && f->lp >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         q = (f->lp)*(f->lm)*(f->qp[n] - f->qm[n]);
         F[n] = ((f->lp)*(f->fm[n]) - (f->lm)*(f->fp[n]) + \
         q)/(f->lp - f->lm);
      }
   }
   else if(f->lm > 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fm[n];
      }
   }
   else if(f->lp < 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fp[n];
      }
   }

   return 0;
}

int HLLC(double *F, flx_ *f, int x)
{
   int n;
   double pstar, ustar, rhobar, abar;
   double rhoR, rhoL, pR, pL;
   double uR, uL, vR, vL, wR, wL;
   double aR, aL;
   double sR, s, sL;
   double qR, qL, eR, eL;
   double QR[10], QL[10];
   double FR[10], FL[10];
   double QsR[10], QsL[10];
   double D[10];
   double rsu;

   rhoR = f->up[0];
   rhoL = f->um[0];
   pR   = f->up[1];
   pL   = f->um[1];
   uR   = f->up[2];
   uL   = f->um[2];
   vR   = f->up[3];
   vL   = f->um[3];
   wR   = f->up[4];
   wL   = f->um[4];

   if(x == 1)
   {
      aR = sqrt(K*pR/rhoR);
      aL = sqrt(K*pL/rhoL);

      rhobar = 0.5*(rhoR + rhoL);
      abar   = 0.5*(aR + aL);

      pstar = 0.5*(pR + pL) + 0.5*(uL - uR)*(rhobar*abar);
      ustar = 0.5*(uR + uL) + 0.5*(pL - pR)/(rhobar*abar);

      if(pstar <= pL)
      {
         qL = 1.0;
      }
      else
      {
         qL = sqrt(1 + ((K + 1)/(2*K))*((pstar/pL) - 1.0));
      }

      if(pstar <= pR)
      {
         qR = 1.0;
      }
      else
      {
         qR = sqrt(1 + ((K + 1)/(2*K))*((pstar/pR) - 1.0));
      }

      sR = uR + aR*qR;
      sL = uL - aL*qL;
      ustar = (pR - pL + rhoL*uL*(sL - uL) - rhoR*uR*(sR - uR))/(rhoL*(sL - uL) - rhoR*(sR - uR));
      s = ustar;
      
      eR = 0.5*rhoR*(uR*uR + vR*vR + wR*wR) + pR/(K-1);
      eL = 0.5*rhoL*(uL*uL + vL*vL + wL*wL) + pL/(K-1);

      rsu = rhoR*((sR - uR)/(sR - s));

      QsR[0] = rsu;
      QsR[1] = rsu*(eR/rhoR + (s - uR)*(s + pR/(rhoR*(sR - uR))));
      QsR[2] = rsu*s;
      QsR[3] = rsu*vR;
      QsR[4] = rsu*wR;

      rsu = rhoL*((sL - uL)/(sL - s));

      QsL[0] = rsu;
      QsL[1] = rsu*(eL/rhoL + (s - uL)*(s + pL/(rhoL*(sL - uL))));
      QsL[2] = rsu*s;
      QsL[3] = rsu*vL;
      QsL[4] = rsu*wL;
   }

   if(x == 2)
   {
      aR = sqrt(K*pR/rhoR);
      aL = sqrt(K*pL/rhoL);

      rhobar = 0.5*(rhoR + rhoL);
      abar   = 0.5*(aR + aL);

      pstar = 0.5*(pR + pL) + 0.5*(vL - vR)*(rhobar*abar);
      ustar = 0.5*(vR + vL) + 0.5*(pL - pR)/(rhobar*abar);

      if(pstar <= pL)
      {
         qL = 1.0;
      }
      else
      {
         qL = sqrt(1 + ((K + 1)/(2*K))*((pstar/pL) - 1.0));
      }

      if(pstar <= pR)
      {
         qR = 1.0;
      }
      else
      {
         qR = sqrt(1 + ((K + 1)/(2*K))*((pstar/pR) - 1.0));
      }

      sR = vR + aR*qR;
      sL = vL - aL*qL;
      ustar = (pR - pL + rhoL*vL*(sL - vL) - rhoR*vR*(sR - vR))/(rhoL*(sL - vL) - rhoR*(sR - vR));
      s = ustar;
      
      eR = 0.5*rhoR*(uR*uR + vR*vR + wR*wR) + pR/(K-1);
      eL = 0.5*rhoL*(uL*uL + vL*vL + wL*wL) + pL/(K-1);

      rsu = rhoR*((sR - vR)/(sR - s));

      QsR[0] = rsu;
      QsR[1] = rsu*(eR/rhoR + (s - vR)*(s + pR/(rhoR*(sR - vR))));
      QsR[2] = rsu*uR;
      QsR[3] = rsu*s;
      QsR[4] = rsu*wR;

      rsu = rhoL*((sL - vL)/(sL - s));

      QsL[0] = rsu;
      QsL[1] = rsu*(eL/rhoL + (s - vL)*(s + pL/(rhoL*(sL - vL))));
      QsL[2] = rsu*uL;
      QsL[3] = rsu*s;
      QsL[4] = rsu*wL;
   }

   if(sL >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fm[n];
      }
   }
   else if(sL < 0 && s >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FL[n] = f->fm[n];
         QL[n] = f->qm[n];
         F[n] = FL[n] + sL*(QsL[n] - QL[n]);
      }
   }
   else if(s <= 0 && sR >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FR[n] = f->fp[n];
         QR[n] = f->qp[n];
         F[n] = FR[n] + sR*(QsR[n] - QR[n]);
      }
   }
   else if(sR <= 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fp[n];
      }
   }

   return 0;
}
