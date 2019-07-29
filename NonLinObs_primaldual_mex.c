/* NonLinObs_primaldual_mex.c - Solves the minimal surface obstacle problem
 *
 *    -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
 *  
 *  subject to boundary condition u=g on partial U
 *  and obstacle constraints
 *
 *       ob1 <= u <= ob2
 *
 *  Usage from Matlab is 
 *
 *     u = NonLinObs_primaldual_mex(f,ob1,ob2,ui,T,eps).
 *
 *  ui = initial condition (also encodes g, as ui=g on partial)
 *  T = max number of iterations
 *  eps = tolerance
 *
 *  To compile: mex COPTIMFLAGS="-Ofast" MS_primaldual_mex.c
 *
 *  The method used is the primal dual method from
 * 
 * Zosso, Dominique, et al. "An Efficient Primal-Dual Method for the Obstacle Problem." 
 * Journal of Scientific Computing 73.1 (2017): 416-437.
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
#define HEAVISIDE(a) (((a)>0)?1:0)
#define SIGN(a) (((a)>0)?1:(-1))
#define ECHECK 10


/*Global variables*/
int n;
int m;

static inline double sgn(double val) {
 return 1-2*signbit(val);
}

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
   double **uout, **ui, **ob1, **ob2, **f;
   double *u_ptr, *ui_ptr, *ob1_ptr, *ob2_ptr, *f_ptr, *Num_Iter;

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
   double **unext = array_double(0);
   double **ubar = array_double(0);
   double **ubarnext = array_double(0);
   double **p1 = array_double(0);
   double **p2 = array_double(0);
   double **q1 = array_double(0);
   double **q2 = array_double(0);
   double **alpha = array_double(0);
   double **alpha1 = array_double(0);
   double **alpha2 = array_double(0);
   double **N = array_double(0);

   double dx = 1.0/(n-1);
   double dy = 1.0/(m-1);
   int NB = (int) round(log(eps*dx*dx)/log(0.5)); /*Number of bisection iterations*/

   /* Check for proper number of arguments. */
   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         u[i][j] = ui[i][j];
         unext[i][j] = ui[i][j];
         ubar[i][j] = ui[i][j];
         ubarnext[i][j] = ui[i][j];
      }
   }

   double a = 4*M_PI*M_PI;
   double r2 = dx/sqrt(6*a);
   double r1 = a*r2;

   double C1 = r2/dx;
   double C2 = r2/dy;
   double idx = 1/dx;
   double idy = 1/dy;
   double ir1 = 1/r1;
   double r12 = r1*r1;

   double uxf,uyf,F,H,M,g,b,c,al;
   double err = 1;
   int count = 0;
   while(err > eps && count < T){
      count++;
      /*Initialization for bisection search*/
      for(i=0;i<n-1;i++){
         for(j=0;j<m-1;j++){
            uxf = idx*(ubar[i+1][j] - ubar[i][j]);
            uyf = idy*(ubar[i][j+1] - ubar[i][j]);
            b = p1[i][j] + r1*uxf;
            c = p2[i][j] + r1*uyf;
            N[i][j] = sqrt(b*b + c*c + 1E-10); 
            q1[i][j] = b/N[i][j]; 
            q2[i][j] = c/N[i][j];
            alpha1[i][j] = 0;
            alpha2[i][j] = MIN(N[i][j],1);
            alpha[i][j] = 0.5*(alpha1[i][j] + alpha2[i][j]);
         }
      }
      /*Bisection search*/
      for(int l=NB;l--;){
         for(i=0;i<n-1;i++){
            for(j=0;j<m-1;j++){
               al = alpha[i][j];
               c = al*al;	
               b = al-N[i][j];
               g = r12*c - (1-c)*b*b;
               switch(signbit(g)){
                  case 0 : alpha2[i][j] = alpha[i][j]; break;
                  case 1 : alpha1[i][j] = alpha[i][j]; break;
               }
               alpha[i][j] = 0.5*(alpha1[i][j] + alpha2[i][j]);
            }
         }
      }
      j = 0;
      for(i=0;i<n-1;i++){
         p1[i][j] = alpha[i][j]*q1[i][j];
         p2[i][j] = alpha[i][j]*q2[i][j];
      }
      i = 0;
      for(j=0;j<m-1;j++){
         p1[i][j] = alpha[i][j]*q1[i][j];
         p2[i][j] = alpha[i][j]*q2[i][j];
      }
      for(i=1;i<n-1;i++){
         for(j=1;j<m-1;j++){
            p1[i][j] = alpha[i][j]*q1[i][j];
            p2[i][j] = alpha[i][j]*q2[i][j];
            unext[i][j] = MIN(MAX(u[i][j] + C1*(p1[i][j] - p1[i-1][j]) + C2*(p2[i][j] - p2[i][j-1]) + r2*f[i][j],ob1[i][j]),ob2[i][j]);
            ubarnext[i][j] = 2*unext[i][j] - u[i][j];
         }
      }
      if(count % ECHECK == 0){
         err = 0;
         for(i=0;i<n-1;i++){
            for(j=0;j<m-1;j++){
               uxf = idx*(ubar[i+1][j] - ubar[i][j]);
               uyf = idy*(ubar[i][j+1] - ubar[i][j]);
               b = 1/sqrt(1 + uxf*uxf + uyf*uyf);
               q1[i][j] = uxf*b;
               q2[i][j] = uyf*b;

               if(i > 1 & j > 1)
                  err = MAX(ABS(MIN(MAX(idx*(q1[i][j] - q1[i-1][j]) + idy*(q2[i][j] - q2[i][j-1]) + f[i][j],ob1[i][j] - ubar[i][j]),ob2[i][j] - ubar[i][j])),err); 
            }
         }
      }
      double **t = u;
      u = unext;
      unext = t;
      t = ubar;
      ubar = ubarnext;
      ubarnext = t;
   }
   for(i=0;i<n;i++){
      for(j=0;j<m;j++){
         uout[i][j] = ubarnext[i][j];
      }
   }
   *Num_Iter = (double) count;
   printf("Number of Iterations = %d.\n",count);
}


