#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"physics.h"

#include"metric.h"

#include"limiters.h"
#include"cond.h"

#include"const.h"
#include"macros.h"
#include"param.h"

/* Define pointers */
double *U, *U1, *U2, *U3;
double *Q, *Q1, *Q2, *Q3;
double *X1;
double *X1p, *X1m;
double *X2;
double *X2p, *X2m;
double *X3;
double *X3p, *X3m;


double x1, x2, x3;
double dx1, dx2, dx3;
double dt, time;
double tmax, timefile, cou;

double start, delta;
double K;

/* Define number of grids */
int Nx1, Nx2, Nx3;

/* Define domain */
double x1max, x2max, x3max;
double x1min, x2min, x3min;

//RIEMANN
double nl, pl, vx1l, vx2l, vx3l;
double nr, pr, vx1r, vx2r, vx3r;
double x_0;

//JET
double r_jet;
double z_jet;

double n_jet, p_jet, vx1_jet, vx2_jet, vx3_jet;
double n_atm, p_atm, vx1_atm, vx2_atm, vx3_atm;

double n_post, p_post, vx1_post, vx2_post, vx3_post;
double n_pre, p_pre, vx1_pre, vx2_pre, vx3_pre;

//Spherical Accretion
int rho_boundary;
double r_out, r_in;
double theta_0, delta_theta;
double density_0, pressure_0, velocity_0;

int binary;

//Paramfile
char paramfile_name[50], outputdirectory[50], outputfile[50];
char restartfile[50];
int read_parameters_file(char const *paramfile_name);
int restart_simulation, restart_filecount;

typedef struct
{
	double up[eq+1];
	double um[eq+1];
	double qp[eq+1];
	double qm[eq+1];
	double fp[eq+1];
	double fm[eq+1];
	double gp[eq+1];
	double gm[eq+1];
	double hp[eq+1];
	double hm[eq+1];
	double lp;
	double lm;
}flx_;

typedef struct
{
	double A[(eq+1)*(eq+1)];
   double Q[eq+1];
   double Q1[eq+1];
   double Q2[eq+1];
	double S[eq+1];
   double F[eq+1];
   double G[eq+1];
   double H[eq+1];
	double Fp[eq+1];
	double Fm[eq+1];
	double Gp[eq+1];
	double Gm[eq+1];
	double Hp[eq+1];
	double Hm[eq+1];
}vec_;

void Allocate_Array();

void New_Size();

int Mesh(); 

void Initial();

double TimeStep();

int Boundaries(double *B);

int PrintValues(double *tprint, double *dtprint, int *itprint);

int Output1(int *itprint);

int Output2(int *itprint);

int Output3(int *itprint);

int Output1_bin(int *itprint);

int Output2_bin(int *itprint);

int Output3_bin(int *itprint);

int Integration();

int RK1D(double *u, double *q, double *q1, double *q2, int order);

int RK2D(double *u, double *q, double *q1, double *q2, int order);

int RK3D(double *u, double *q, double *q1, double *q2, int order);

int Flux1D(vec_ *v, lim_ *l, int *I);
                                   
int Flux2D(vec_ *v, lim_ *l, int *I);
                                   
int Flux3D(vec_ *v, lim_ *l, int *I);

int Sources(double *u, vec_ *v, int *I);

int Hll(double *F, flx_ *f, int x);

int Hllc(double *F, flx_ *f, int x);

int AMATRIX1D(double *u, vec_ *v, int *I);
                                                            
int AMATRIX2D(double *u, vec_ *v, int *I);
                                              
int AMATRIX3D(double *u, vec_ *v, int *I);

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I);

int c1(int n, int i);

int c2(int n, int i, int j);

int c3(int n, int i, int j, int k);

int MxV(double *M, double *V, double *L);

void RoundGen(double *num);

#if DIM == 1

   #define  U(N,x)  U[(N)*(Nx1+1) + (x)]
   #define  Q(N,x)  Q[(N)*(Nx1+1) + (x)]
   #define  B(N,x)  B[(N)*(Nx1+1) + (x)]
   #define  u(N,x)  u[(N)*(Nx1+1) + (x)]
   #define  q(N,x)  q[(N)*(Nx1+1) + (x)]
   #define q1(N,x) q1[(N)*(Nx1+1) + (x)]
   #define q2(N,x) q2[(N)*(Nx1+1) + (x)]

#elif DIM == 2 || DIM == 4

   #define  U(N,x,y)  U[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  Q(N,x,y)  Q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  B(N,x,y)  B[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  u(N,x,y)  u[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  q(N,x,y)  q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q1(N,x,y) q1[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q2(N,x,y) q2[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

#elif DIM == 3

   #define  U(N,x,y,z)  U[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  Q(N,x,y,z)  Q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  B(N,x,y,z)  B[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  u(N,x,y,z)  u[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  q(N,x,y,z)  q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q1(N,x,y,z) q1[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q2(N,x,y,z) q2[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]

#endif

