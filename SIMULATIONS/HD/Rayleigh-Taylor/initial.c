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

   /////////////////////////////////
   // Rayleigh-Taylor Instability //
   /////////////////////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X2[j] <= 0.0)
         {
            U(0,i,j) = 1.0;
            U(1,i,j) = 2.5 - U(0,i,j)*grid.X2[j];
            U(2,i,j) = 0.0;
            U(3,i,j) = -0.1*(1 + cos(2*M_PI*grid.X1[i]/0.5))*(1 + cos(2*M_PI*grid.X2[j]/1.5))/4;
         }
         if(grid.X2[j] > 0.0)
         {
            U(0,i,j) = 2.0;
            U(1,i,j) = 2.5 - U(0,i,j)*grid.X2[j];
            U(2,i,j) = 0.0;
            U(3,i,j) = -0.1*(1 + cos(2*M_PI*grid.X1[i]/0.5))*(1 + cos(2*M_PI*grid.X2[j]/1.5))/4;
         }
      }
   }
   //////////////////////////////////
}
