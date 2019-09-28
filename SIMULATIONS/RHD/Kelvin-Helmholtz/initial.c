/*
 * File Name : initial.c
 * Description : aztekas initial module for Relativistic Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:56:15
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
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
            U(RHO,i,j) = rhod;
            U(PRE,i,j) = pd;
            U(VX1,i,j) = vx1d*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j])); 
            U(VX2,i,j) = vx2d*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j]));
         }
         else if(fabs(grid.X2[j]) < x_0) 
         {
            U(RHO,i,j) = rhou;
            U(PRE,i,j) = pu;
            U(VX1,i,j) = vx1u*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j])); 
            U(VX2,i,j) = vx2u*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j]));
         }
      }
   }
   //////////////////////////////////
}
