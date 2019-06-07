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

int Boundaries(double *B)
{
   int i, j, k, n, cell;

#if DIM == 1

   Outflow(B);

   for(cell = 0; cell <= gc; cell++)
   {
      U(0,Nx1-cell) = density_0;
      U(1,Nx1-cell) = pressure_0;
      U(2,Nx1-cell) = velocity_0;
   }

   for(i = 0; i <= Nx1; i++)
   {
      for(n = 0; n < eq; n++)
      {
         RoundGen(&B(n,i));
      }
   }

#elif DIM == 2

   Outflow(B);
   Reflection(B);

   for(j = 0; j <= Nx2; j++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         U(0,Nx1-cell,j) = density_0;
         U(1,Nx1-cell,j) = pressure_0;
         U(2,Nx1-cell,j) = velocity_0;
         U(3,Nx1-cell,j) = 0.0;
      }
   }

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(n = 0; n < eq; n++)
         {
            RoundGen(&B(n,i,j));
         }
      }
   }

#elif DIM == 4

   Outflow(B);
   Reflection(B);

   for(j = 0; j <= Nx2; j++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         U(0,Nx1-cell,j) = density_0;
         U(1,Nx1-cell,j) = pressure_0;
         U(2,Nx1-cell,j) = velocity_0;
         U(3,Nx1-cell,j) = 0.0;
         U(4,Nx1-cell,j) = 0.0;
      }
   }

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(n = 0; n < eq; n++)
         {
            RoundGen(&B(n,i,j));
         }
      }
   }

#endif

   return 0;
}
