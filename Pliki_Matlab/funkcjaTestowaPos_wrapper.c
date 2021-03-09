
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
double alfa = 0.1394;
double p = 0, i = 0, d = 1;
double dPrev = 1;
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void funkcjaTestowaPos_Start_wrapper(real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
/*
 * Custom Start code goes here.
 */
xD[0] = 0;
xD[1] = 0;
xD[2] = 0;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void funkcjaTestowaPos_Outputs_wrapper(const real_T *ePos,
			real_T *yPos,
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
// delta = kp *  ((xD[0] - xD[1]) + (h/ti) * xD[0] + (td*alfa*(xD[0]-xD[1]))/(td+alfa)*h );
// delta = kp *  ((xD[0] - xD[1]) + (h/ti) * xD[0] + (td*(xD[0]-xD[1]))/h );
// delta = kp *  ((xD[0] - xD[1]) + (h/ti) * xD[0] + (td*(xD[0]-xD[1])/(1 + td/alfa)) );  
p = kp * xD[0];
i = kp/ti * h * xD[1];
d = kp * td * ((xD[0] - xD[2])/(h+alfa) + (alfa)/(h+alfa) * dPrev );
dPrev = d;
// d = kp*td * (1/h) * (xD[0] - xD[2]);
yPos[0]  = p + i ;
// printf("d = %f \n", d);
//     printf("dPrev = %f \n", dPrev);
yPrev = yPos[0];
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Updates function
 *
 */
void funkcjaTestowaPos_Update_wrapper(const real_T *ePos,
			real_T *yPos,
			real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
xD[2] = xD[0];
xD[1] += ePos[0];    // suma ei
xD[0] = ePos[0];    //ei
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}

