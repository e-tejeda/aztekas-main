#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int main(int argc, char* argv[])
{
   FILE *fdata;
   char line[100];
   char file[100];
   int i, j;
   int Nx1, Nx2;
   int idum;
   double time, dum, acc, in, out;

   if(argc != 2)
   {  
      printf("%s\n", "Wrong number of arguments");
      printf("%s\n", "./acc file");
      exit(EXIT_FAILURE);
   }

   strcpy(file,argv[1]);

   fdata = fopen(file,"r");
   idum = fscanf(fdata,"%s\n",line);
   idum = fscanf(fdata,"%lf\n",&time); 
   idum = fscanf(fdata,"%d\n",&Nx1);
   idum = fscanf(fdata,"%d\n",&Nx2);
   idum = fscanf(fdata,"%s\n",line);

   double *acc_sum, *in_sum, *out_sum, *Radius, *Theta, *Density, *Velocity, *VelTheta;
   acc_sum = (double *)malloc(Nx1*Nx2*sizeof(double));
   in_sum = (double *)malloc(Nx1*Nx2*sizeof(double));
   out_sum = (double *)malloc(Nx1*Nx2*sizeof(double));
   Radius = (double *)malloc(Nx1*Nx2*sizeof(double));
   Theta = (double *)malloc(Nx1*Nx2*sizeof(double));
   Density = (double *)malloc(Nx1*Nx2*sizeof(double));
   Velocity = (double *)malloc(Nx1*Nx2*sizeof(double));
   VelTheta = (double *)malloc(Nx1*Nx2*sizeof(double));


   for(i = 0; i < Nx1; i++)
   {
      for(j = 0; j< Nx2; j++)
      {
         idum = fscanf(fdata,"%lf %lf %lf %lf %lf %lf %lf\n",&Radius[i*Nx2 + j],\
         &Theta[i*Nx2 + j],&Density[i*Nx2 + j],&dum,&Velocity[i*Nx2 + j],&VelTheta[i*Nx2 + j],&dum);
      }
   }

   dum = 0;
   in = 0.0;
   out = 0.0;
   double rho, W, VV, r, th, thm, vr, vt;
   double gamma = 1.0001;
   double Mb = 0.25*pow(2.0/(5.0 - 3.0*gamma),(5.0 - 3.0*gamma)/(2.0*(gamma - 1.0)));

   for(i = 10; i < Nx1-10; i++)
   {
      for(j = 0; j< Nx2; j++)
      {
         r = Radius[i*Nx2 + j];
         th = Theta[i*Nx2 + j];
         thm = Theta[i*Nx2 + j-1];
         vr = Velocity[i*Nx2 + j];
         vt = VelTheta[i*Nx2 + j];
         rho = Density[i*Nx2 + j];

         //acc = acc - rho*vr*r*r*sin(th)*(M_PI_2/(Nx2))/(Mb);
         acc = acc - rho*vr*r*r/(Mb*Nx2);

         if(i == Nx1-11 && vr <= 0)
         {
            in = in - rho*vr*r*r*sin(th)*(M_PI_2/(Nx2))/(Mb);
         }

         if(i == Nx1-11 && vr > 0)
         {
            out = out + rho*vr*r*r*sin(th)*(M_PI_2/(Nx2))/(Mb);
         }
      }
      acc_sum[i] = acc;
      acc = 0.0;
   }
   fclose(fdata);

   double acc_prom = 0.0;
   double in_prom  = 0.0;
   double out_prom = 0.0;

   for(i = 10; i < Nx1-10; i++)
   {
      acc_prom = acc_prom + acc_sum[i]/(Nx1-20);
   }

   if(out == 0){in = acc_prom;}

   printf("%e %e %e %e\n",time,acc_prom,in,in-acc_prom);


   return 0;
}
