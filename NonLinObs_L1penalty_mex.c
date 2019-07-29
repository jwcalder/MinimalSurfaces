//NonLinObs_primaldual_mex.c- Solves the minimal surface obstacle problem
//
//   -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
// 
// subject to boundary condition u=g on partial U
// and obstacle constraints
//
//      ob1 <= u
//
// Other argument:
//
//     ui = initial condition (also encodes g, as ui=g on partial U)
//     T = max number of iterations
//     tol = tolerance.
//
//  Uses the L1-penalty method from 
//
//  Tran, Giang, et al. "An L^1 Penalty Method for General Obstacle Problems." 
//  SIAM Journal on Applied Mathematics 75.4 (2015): 1424-1444.
//
// Author: Jeff Calder, 2019.
//
//
//

#include "mex.h"
#include "time.h"
#include "math.h"

#define ABS(a) (((a)<0)?-(a):(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ECHECK 10
//#define ERROR 


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
   double **u, **ui, **bi, **ob1, **ob2;
   double *u_ptr, *ui_ptr, *bi_ptr, *ob1_ptr, *ob2_ptr, *Num_Iter;

   /* Check for proper number of arguments. */
   if(nrhs<8) {
      mexErrMsgIdAndTxt("MATLAB:affineflow_fm_mex:invalidNumInputs","Inputs: (ob1,ob2,ui,bi,T,tol,lambda,mu).");
   } 

   /* Retrieve input.*/
   n = mxGetM(prhs[0]);
   m = mxGetN(prhs[0]);
   int T = *mxGetPr(prhs[4]);
   double tol = *mxGetPr(prhs[5]);
   double lambda = *mxGetPr(prhs[6]);
   double mu = *mxGetPr(prhs[7]);

   /* Create matrix for the return argument. */
   plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)m, mxREAL);
   plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

   /* Assign pointers to each input and output. */
   ob1_ptr = mxGetPr(prhs[0]);
   ob2_ptr = mxGetPr(prhs[1]);
   ui_ptr = mxGetPr(prhs[2]);
   bi_ptr = mxGetPr(prhs[3]);

   u_ptr = mxGetPr(plhs[0]);
   Num_Iter = mxGetPr(plhs[1]);

   /*Initialize arrays. */
   u = mxMalloc(n * sizeof(double*));
   ob1 = mxMalloc(n * sizeof(double*));
   ob2 = mxMalloc(n * sizeof(double*));
   ui = mxMalloc(n * sizeof(double*));
   bi = mxMalloc(n * sizeof(double*));

   for(i=0;i<n;i++){
      u[i] = u_ptr + n*i;
      ob1[i] = ob1_ptr + n*i;
      ob2[i] = ob2_ptr + n*i;
      ui[i] = ui_ptr + n*i;
      bi[i] = bi_ptr + n*i;
   }

   double **w = array_double(0);
   double **v = array_double(0);
   double **b = array_double(0);
   double **p1 = array_double(0);
   double **p2 = array_double(0);
   double **uprev = array_double(0);
   double **uerr = array_double(0);

   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         u[i][j] = ui[i][j];
         v[i][j] = ob1[i][j] - u[i][j];
         b[i][j] = bi[i][j];
         uprev[i][j] = u[i][j];
      }
   }

   double dx = 1.0/(n-1);
   double dy = 1.0/(m-1);
   double L = 2/(dx*dx);
   double dt = 1/L/4;
   double alpha = (sqrt(L) - sqrt(lambda))/(sqrt(L) + sqrt(lambda));
   double C = mu/lambda;
   double idx = 1/dx;
   double idy = 1/dy;

   double tol2 = dx*dx/100;

   double uxf,uyf,N,F;
   double err = 1;
   int count = 0;
   int outercount = 0;
   while(err > tol && count < T){
      outercount++;

      double err2 = 1;
      for(i=0;i<n;i++){
         for(j=0;j<m;j++){
            uprev[i][j] = u[i][j];
            uerr[i][j] = u[i][j];
         }
      }
      int count2 = 0;
      while(err2 > tol2 && count < T){
         count2++;
         count++;
         for(i=0;i<n;i++){
            for(j=0;j<m;j++){
               w[i][j] = u[i][j] + alpha*(u[i][j] - uprev[i][j]);
               uprev[i][j] = u[i][j];
            }
         }   
         j = 0;
         for(i=0;i<n-1;i++){
            uxf = idx*(w[i+1][j] - w[i][j]);
            uyf = idy*(w[i][j+1] - w[i][j]);
            N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
            p1[i][j] = uxf*N;
            p2[i][j] = uyf*N;
         }
         i = 0;
         for(j=0;j<m-1;j++){
            uxf = idx*(w[i+1][j] - w[i][j]);
            uyf = idy*(w[i][j+1] - w[i][j]);
            N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
            p1[i][j] = uxf*N;
            p2[i][j] = uyf*N;
         }
         err2 = 0;
         for(i=1;i<n-1;i++){
            for(j=1;j<m-1;j++){
               uxf = idx*(w[i+1][j] - w[i][j]);
               uyf = idy*(w[i][j+1] - w[i][j]);
               N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
               p1[i][j] = uxf*N;
               p2[i][j] = uyf*N;

               F = -idx*(p1[i][j] - p1[i-1][j]) - idy*(p2[i][j] - p2[i][j-1]); 
               u[i][j] = w[i][j] - dt*(F + lambda*(v[i][j] - ob1[i][j] + w[i][j] + b[i][j]));
               err2 = MAX(ABS(u[i][j] - uprev[i][j]),err2);
            }
         }
      }
#ifdef ERROR
      err = 0;
#else
      if(outercount % ECHECK == 0){
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
      }
#endif
      for(i=1;i<n-1;i++){
         for(j=1;j<m-1;j++){
            //Updates
            double z = ob1[i][j] - u[i][j] - b[i][j];
            v[i][j] = MAX(z-C,0) + MIN(z,0);
            b[i][j] = b[i][j] + u[i][j] + v[i][j] - ob1[i][j];

#ifdef ERROR
            err = MAX(ABS(uerr[i][j] - u[i][j]),err);
#else
            if(outercount % ECHECK == 0){
               //Compute error
               uxf = idx*(u[i+1][j] - u[i][j]);
               uyf = idy*(u[i][j+1] - u[i][j]);
               N = 1/sqrt(1 + uxf*uxf + uyf*uyf);
               p1[i][j] = uxf*N;
               p2[i][j] = uyf*N;

               F = idx*(p1[i][j] - p1[i-1][j]) + idy*(p2[i][j] - p2[i][j-1]); 
               err = MAX(ABS(MIN(MAX(F,ob1[i][j] - u[i][j]),ob2[i][j] - u[i][j])),err); 
            }
#endif
         }
      }
   }
   *Num_Iter = (double) count;
   printf("Number of Iterations = %d.\n",count);
   printf("Number of Outer Iterations = %d.\n",outercount);
}


