
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
#include <stdio.h>
#include "simstruc.h"
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
double yPrev = 0;
double delta;
double kp = 1.612, ti = 5.5780, td = 1.3945;
double h = 0.01;
double alfa = 0.13945;
double a0 = 0, a1 = 0, a2 = 0, k = 0;
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void funkcjaTestowaInc_Start_wrapper(real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
/*
 * Custom Start code goes here.
 */
xD[0] = 0;
xD[1] = 0;
xD[2] = 0;
    xD[3] = 0;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void funkcjaTestowaInc_Outputs_wrapper(const real_T *eInc,
			real_T *yInc,
			const real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* This sample sets the output equal to the input
      y0[0] = u0[0]; 
 For complex signals use: y0[0].re = u0[0].re; 
      y0[0].im = u0[0].im;
      y1[0].re = u1[0].re;
      y1[0].im = u1[0].im;
 */

k = alfa*(alfa + h);
a0 = kp * ( 1 + h/ti + td/(h+alfa) );
a1 = kp * ( 1 + ((2-k)*td)/(h+alfa) );
a2 = kp * ( ( ((k-1)*(k-1)) * td )/(alfa+h) );

// printf("k = %f \n", k);
// printf("a0 = %f \n", a0);
// printf("a1 = %f \n", a1);
// printf("a2 = %f  \n", a2);
   printf("yPrev = %f  \n", xD[3]); 

yInc[0] = a0*xD[0] - a1*xD[1] + a2*xD[2] + xD[3];

/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Updates function
 *
 */
void funkcjaTestowaInc_Update_wrapper(const real_T *eInc,
			real_T *yInc,
			real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
xD[3] = yInc[0];
xD[2] = xD[1];        //ei-2
xD[1] = xD[0];        //ei-1
xD[0] = eInc[0];    //e1
// yPrev = yInc[0];
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}

