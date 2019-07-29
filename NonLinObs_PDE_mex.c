/* NonLinObs_PDE_mex.c - Solves the minimal surface obstacle problem
 *
 *    -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
 *  
 *  with forcing f, subject to boundary condition u=g on partial U
 *  and obstacle constraints
 *
 *       ob1 <= u <= ob2
 *
 *  Usage from Matlab is 
 *
 *     u = NonLinObs_PDE_mex(f,ob1,ob2,ui,T,eps).
 *
 *  ui = initial condition (also encodes g, as ui=g on partial)
 *  T = max number of iterations
 *  eps = tolerance.
 *
 *  To compile: mex COPTIMFLAGS="-O3" NonLinObs_PDE_mex.c
 *
 *  The method used is PDE acceleration
 *
 *  J. Calder and A. Yezzi. PDE Acceleration: A convergence rate analysis and applications to obstacle problems. 2018
 *
 *  M. Benyamin, J. Calder, G. Sundaramoorthi, and A. Yezzi. Accelerated PDE’s for efficient solution of regularized inversion problems. 2018
 * 
 *  Author: Jeff Calder, 2018.
 *
 *
 */

#include "mex.h"
#include "time.h"
#include "math.h"

#define ABS(a) (((a)<0)?-(a):(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


/*Global variables*/
int n;
int m;

/*Allocate memory for a 2D array of ints*/
int** array_int(int val){

   int **ptr = mxMalloc(m*sizeof(int*));
   /*ptr[0] = mxMalloc((m*n+1)*sizeof(int));*/
   ptr[0] = mxMalloc(m*n*sizeof(int));
   int i,j;
   for(i=0;i<m;i++){
      ptr[i] = ptr[0] + n*i;
      for(j=0;j<n;j++){
         ptr[i][j] = val;
      }
   }
   return ptr;
}

/*Allocate memory for a 2D array of doubles*/
double** array_double(double val){

   double **ptr = mxMalloc(m*sizeof(double*));
   ptr[0] = mxMalloc(m*n*sizeof(double));
   int i,j;
   for(i=0;i<m;i++){
      ptr[i] = ptr[0] + n*i;
      for(j=0;j<n;j++){
         ptr[i][j] = val;
      }
   }
   return ptr;
}




void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){

   int i,j;
   double **uout, **f, **ui, **ob1, **ob2;
   double *u_ptr,*f_ptr,*ui_ptr, *ob1_ptr, *ob2_ptr, *Num_Iter;

   /*mexPrintf("A\n");mexEvalString("drawnow;");*/
   /* Check for proper number of arguments. */
   if(nrhs<6) {
      mexErrMsgIdAndTxt("MATLAB:affineflow_fm_mex:invalidNumInputs","Inputs: (f,ob1,ob2,ui,T,eps).");
   } 

   /* Retrieve input.*/
   n = mxGetM(prhs[0]);
   m = mxGetN(prhs[0]);
   int T = *mxGetPr(prhs[4]);
   double eps = *mxGetPr(prhs[5]);

   /* Create matrix for the return argument. */
   plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)m, mxREAL);
   plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

   /* Assign pointers to each input and output. */
   f_ptr = mxGetPr(prhs[0]);
   ob1_ptr = mxGetPr(prhs[1]);
   ob2_ptr = mxGetPr(prhs[2]);
   ui_ptr = mxGetPr(prhs[3]);
   u_ptr = mxGetPr(plhs[0]);
   Num_Iter = mxGetPr(plhs[1]);

   /*Initialize arrays. */
   uout = mxMalloc(n * sizeof(double*));
   f = mxMalloc(n * sizeof(double*));
   ob1 = mxMalloc(n * sizeof(double*));
   ob2 = mxMalloc(n * sizeof(double*));
   ui = mxMalloc(n * sizeof(double*));

   for(i=0;i<n;i++){
      uout[i] = u_ptr + n*i;
      f[i] = f_ptr + n*i;
      ob1[i] = ob1_ptr + n*i;
      ob2[i] = ob2_ptr + n*i;
      ui[i] = ui_ptr + n*i;
   }

   double **u = array_double(0);
   double **p1 = array_double(0);
   double **p2 = array_double(0);
   double **unext = array_double(0);
   double **uprev = array_double(0);
   
   double dx = 1.0/(n-1);
   double dy = 1.0/(m-1);
   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         u[i][j] = ui[i][j];
         unext[i][j] = ui[i][j];
         uprev[i][j] = ui[i][j];
      }
   }
   double dt = 0.8*sqrt(1.0/2.0)*dx;
   double a = 2*M_PI;

   double A = (2+a*dt)/(1+a*dt);
   double B = -1/(1+a*dt);
   double C = dt*dt/(1+a*dt);
   double idx = 1/dx;
   double idy = 1/dy;

   double uxf,uyf,N,F;
   double err = 1;
   int count = 0;
   while(err > eps && count < T){
      count++;
      err = 0;
      j = 0;
      for(i=0;i<n-1;i++){
         uxf = idx*(u[i+1][j] - u[i][j]);
         uyf = idy*(u[i][j+1] - u[i][j]);
         N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
         p1[i][j] = uxf*N;
         p2[i][j] = uyf*N;
      }
      i = 0;
      for(j=0;j<m-1;j++){
         uxf = idx*(u[i+1][j] - u[i][j]);
         uyf = idy*(u[i][j+1] - u[i][j]);
         N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
         p1[i][j] = uxf*N;
         p2[i][j] = uyf*N;
      }
      for(i=1;i<n-1;i++){
         for(j=1;j<m-1;j++){
            uxf = idx*(u[i+1][j] - u[i][j]);
            uyf = idy*(u[i][j+1] - u[i][j]);
            N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
            p1[i][j] = uxf*N;
            p2[i][j] = uyf*N;

            F = idx*(p1[i][j] - p1[i-1][j]) + idy*(p2[i][j] - p2[i][j-1]) + f[i][j]; 
            unext[i][j] = MIN(MAX(A*u[i][j] + B*uprev[i][j] + C*F,ob1[i][j]),ob2[i][j]);
            err = MAX(ABS(MIN(MAX(F,ob1[i][j] - u[i][j]),ob2[i][j] - u[i][j])),err); 
         }
      }
      double **t = uprev;
      uprev = u;
      u = unext;
      unext = t;
   }
   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         uout[i][j] = uprev[i][j];
      }
   }
   *Num_Iter = (double) count;
   printf("Number of Iterations = %d.\n",count);
}


