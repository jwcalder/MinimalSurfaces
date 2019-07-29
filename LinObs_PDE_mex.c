/* LinObs_PDE_mex.c - Solves the Dirichlet Problem
 *
 *    -div (ae*grad u) = f   in  U = (0,1)^2
 *  
 *  subject to boundary condition u=g on partial U
 *  and obstacle constraints
 *
 *       ob1 <= u <= ob2
 *
 *  Usage from Matlab is 
 *
 *     u = LinObs_PDE_mex(ae,f,ob1,ob2,ui,T,eps).
 *
 *  T = max number of iterations
 *  eps = tolerance.
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
#define TT 100


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
   double **uout, **f, **ae, **ob1, **ob2, **ui;
   double *uout_ptr, *f_ptr, *ae_ptr, *ob1_ptr, *ob2_ptr, *ui_ptr, *Num_Iter;

   /*mexPrintf("A\n");mexEvalString("drawnow;");*/
   /* Check for proper number of arguments. */
   if(nrhs<7) {
      mexErrMsgIdAndTxt("MATLAB:affineflow_fm_mex:invalidNumInputs","Inputs: (ae,f,ob1,ob2,ui,T,eps).");
   } 

   /* Retrieve input.*/
   n = mxGetM(prhs[0]);
   m = mxGetN(prhs[0]);
   int T = *mxGetPr(prhs[5]);
   double eps = *mxGetPr(prhs[6]);

   /* Create matrix for the return argument. */
   plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)m, mxREAL);
   plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

   /* Assign pointers to each input and output. */
   ae_ptr = mxGetPr(prhs[0]);
   f_ptr = mxGetPr(prhs[1]);
   ob1_ptr = mxGetPr(prhs[2]);
   ob2_ptr = mxGetPr(prhs[3]);
   ui_ptr = mxGetPr(prhs[4]);
   uout_ptr = mxGetPr(plhs[0]);
   Num_Iter = mxGetPr(plhs[1]);

   /*Initialize arrays. */
   uout = mxMalloc(n * sizeof(double*));
   ae = mxMalloc(n * sizeof(double*));
   f = mxMalloc(n * sizeof(double*));
   ob1 = mxMalloc(n * sizeof(double*));
   ob2 = mxMalloc(n * sizeof(double*));
   ui = mxMalloc(n * sizeof(double*));

   for(i=0;i<n;i++){
      uout[i] = uout_ptr + n*i;
      ae[i] = ae_ptr + n*i;
      f[i] = f_ptr + n*i;
      ob1[i] = ob1_ptr + n*i;
      ob2[i] = ob2_ptr + n*i;
      ui[i] = ui_ptr + n*i;
   }

   double **u = array_double(0);
   double **unext = array_double(0);
   double **uprev = array_double(0);
   double **axf = array_double(0);
   double **ayf = array_double(0);
   double **axb = array_double(0);
   double **ayb = array_double(0);
   
   double dx = 1.0/(n-1);
   double dy = 1.0/(m-1);
   double maxae = 0;
   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         u[i][j] = ui[i][j];
         unext[i][j] = ui[i][j];
         uprev[i][j] = ui[i][j];
         maxae = MAX(maxae,ae[i][j]);
      }
   }
   double dt = 0.8*sqrt(1.0/2.0/maxae)*dx;
   double a = 2*M_PI;
   if(nrhs == 8)
      a = *mxGetPr(prhs[7]);

   for(i=1;i<n-1;i++){
      for(j=1;j<m-1;j++){
         axf[i][j] = (ae[i][j] + ae[i+1][j])/(2*dx*dx);
         ayf[i][j] = (ae[i][j] + ae[i][j+1])/(2*dy*dy);
         axb[i][j] = (ae[i][j] + ae[i-1][j])/(2*dx*dx);
         ayb[i][j] = (ae[i][j] + ae[i][j-1])/(2*dy*dy);
      }
   }

   double A = (2+a*dt)/(1+a*dt);
   double B = -1/(1+a*dt);
   double C = dt*dt/(1+a*dt);

   double err = 1;
   int count = 0;
   double mu = 0;
   while(err > eps && count < T){
      count++;
      err = 0;
      mu = 0;
      for(i=1;i<n-1;i++){
         for(j=1;j<m-1;j++){
            mu = MAX(ABS(u[i][j]),mu);
            double uxf = u[i+1][j] - u[i][j];
            double uxb = u[i][j] - u[i-1][j];
            double uyf = u[i][j+1] - u[i][j];
            double uyb = u[i][j] - u[i][j-1];

            double F = axf[i][j]*uxf - axb[i][j]*uxb + ayf[i][j]*uyf - ayb[i][j]*uyb + f[i][j];
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


