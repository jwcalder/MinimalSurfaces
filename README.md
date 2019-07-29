Minimal Surface Obstacle Problem Solvers
========================================

This code computes solutions of minimal surface obstacle problems using the PDE acceleration technique described in the paper

"PDE Acceleration: A convergence rate analysis and applications to obstacle problems", Jeff Calder, and Anthony Yezzi.  arXiv preprint, 2018. Available here: https://arxiv.org/abs/1810.01066

All simulations presented in the above paper can be reproduced with this code. We also provide code for the primal dual method from 

"An efficient primal-dual method for the obstacle problem", Dominique Zosso, Braxton Osting, Mandy Mengqi Xia, and Stanley J. Osher. Journal of Scientific Computing 73, no. 1 (2017): 416-437.

and the L1 penalty method from 

"An L^1 Penalty Method for General Obstacle Problems", Giang Tran, Hayden Schaeffer, William M. Feldman, and Stanley J. Osher. SIAM Journal on Applied Mathematics 75, no. 4 (2015): 1424-1444.

All code has a pure vectorized Matlab implementation, and a faster C code implementation via the mex interface. Run the mexmake.m file to compile the mex C code.

```
>> mexmake.m
```

There are three demo scripts: PoissonExample.m, HomogenizationExample.m, and MinimalSurfaceExample.m, that show how to run the code and reproduce experiments from the paper.
