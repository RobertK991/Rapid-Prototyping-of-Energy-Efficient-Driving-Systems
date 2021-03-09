
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
double uPrev;
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void regulatorP_Start_wrapper(real_T *xD,
			const real_T *kp_Pid, const int_T p_width0,
			const real_T *ti_Pid, const int_T p_width1,
			const real_T *td_Pid, const int_T p_width2,
			const real_T *h, const int_T p_width3)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
//P
xD[0] = 0;
xD[1] = 0;
uPrev = 0;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void regulatorP_Outputs_wrapper(const real_T *eInc,
			const real_T *ePos,
			real_T *yInc,
			real_T *yPos,
			const real_T *xD,
			const real_T *kp_Pid, const int_T p_width0,
			const real_T *ti_Pid, const int_T p_width1,
			const real_T *td_Pid, const int_T p_width2,
			const real_T *h, const int_T p_width3)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* This sample sets the output equal to the input
      y0[0] = u0[0]; 
 For complex signals use: y0[0].re = u0[0].re; 
      y0[0].im = u0[0].im;
      y1[0].re = u1[0].re;
      y1[0].im = u1[0].im;
 */

printf("Nastawy: %f, %f, %f \n",kp_Pid[0], &td_Pid[0], &ti_Pid[0]);
double deltaP = kp_Pid[0] *(xD[0] - xD[1]);
double pPos = kp_Pid[0] * ePos[0];
yInc[0] = deltaP + uPrev;
yPos[0] = pPos;
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Updates function
 *
 */
void regulatorP_Update_wrapper(const real_T *eInc,
			const real_T *ePos,
			real_T *yInc,
			real_T *yPos,
			real_T *xD,
			const real_T *kp_Pid, const int_T p_width0,
			const real_T *ti_Pid, const int_T p_width1,
			const real_T *td_Pid, const int_T p_width2,
			const real_T *h, const int_T p_width3)
{
/* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
/*
 * Code example
 *   xD[0] = u0[0];
 */
xD[1] = xD[0];
xD[0] = eInc[0];
uPrev = yInc[0];
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}

