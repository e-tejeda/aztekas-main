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

void Initial(double *dtprint)
{
   int n, i, j, k, cell;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   //////////////////////////////////
   // Kelvin-Helmholtz Instability //
   //////////////////////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(fabs(grid.X2[j]) >= x_0)
         {
            U(0,i,j) = nl;
            U(1,i,j) = pl;
            U(2,i,j) = vx1l*(1 + 0.01*cos(2*M_PI*grid.X1[i])*sin(2*M_PI*grid.X2[j])); 
            U(3,i,j) = vx2l*(1 + 0.01*cos(2*M_PI*grid.X1[i])*sin(2*M_PI*grid.X2[j]));
         }
         else if(fabs(grid.X2[j]) < x_0) 
         {
            U(0,i,j) = nr;
            U(1,i,j) = pr;
            U(2,i,j) = vx1r*(1 + 0.01*cos(2*M_PI*grid.X1[i])*sin(2*M_PI*grid.X2[j])); 
            U(3,i,j) = vx2r*(1 + 0.01*cos(2*M_PI*grid.X1[i])*sin(2*M_PI*grid.X2[j]));
         }
      }
   }
   //////////////////////////////////
}
