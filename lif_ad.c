/*=================================================================
 *
 * LIF_AD.C	C-version .MEX file simulating a LIF neuron with
 *	        Ca adaptation 
 *
 * The calling syntax is:
 *
 *		[t, x] = lif_ad(ie,rin,t0,vth,vreset,vk,tref,tca,ca_spk,gaphb)
 *
 *  The closest m-file code is in the function lif_adapt
 *
 * 
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	IE_IN	prhs[0] /*current (nA)*/
#define RIN_IN  prhs[1] /*input resistance MO*/
#define T0_IN   prhs[2] /*time constant (ms)*/
#define	VTH_IN	prhs[3] /*absolute spiking threshold (mV)*/
#define VRE_IN  prhs[4] /*reset potential (mV)*/
#define VK_IN   prhs[5] /*potassium reversal potential*/
#define TREF_IN prhs[6] /*refractory period (ms)*/
#define TCA_IN  prhs[7] /*calcium time constant (ms)*/
#define CASP_IN prhs[8] /*calcium entry per spike (muM)*/
#define GAHP_IN prhs[9] /*AHP conductance per unit [Ca] (nuS/muM)*/

/* Output Arguments */

#define	T_OUT	plhs[0]
#define X_OUT   plhs[1]

static	double	t_start = -50.0;
static	double	t_end = 700.0;
static  double  dt = 0.01; /* time step in msec */
static  double  vrest = -70.0;

double ie,rin,t0,vth,vreset,vk,tref,tca,caspk,gahpb,gahpn; 
double t0m1,c;

double f_dloc1(double x1, double x2)
{
  double xd1;
  
  xd1 = t0m1*(-(x1-vrest)-gahpn*x2*(x1-vk)+c);
  return xd1;
}


double f_dloc2(double x1, double x2)
{
  double xd2;
  
  xd2 = - x2/tca;

  return xd2;
}

static void lif_ad_comp(
		   double	tp[],
		   double	xp[]
		   )
{
    double curr_eff,c1,h;
    double k11, k12, k21, k22, k31, k32, k41, k42;
    int i, n;
    int c_start, c_end, tref_ind, tref_counter = 0;

    gahpn = gahpb*rin;

    /*fill in the time vector */
    n = (int) floor((t_end-t_start)/dt) + 1;
    tp[0] = t_start;
    for (i=1;i<n;i++){
      tp[i] = tp[i-1]+dt;
    }

    
    /* starting and ending index of current pulse */
    c_start = (int) (-t_start/dt);
    c_end = (int) ((500.0-t_start)/dt);

    /* initial conditions */
    xp[0] = vrest;
    xp[1] = 0;
    xp[2] = 0;

    curr_eff = rin*ie;
    c1 = 1/3.0;
    h = dt/2.0;
    t0m1 = 1/t0;
    
    tref_ind = (int) floor(tref/dt);
    if ( tref_ind <= 0 )
      tref_ind = 1;

    for (i=3;i<3*n-1;i=i+3){
      if ( (i >= 3*c_start) && (i <= 3*c_end) )
        c = curr_eff;
      else
        c = 0;
    
      
      k11 = h*f_dloc1(xp[i-3],xp[i-2]);
      k12 = h*f_dloc2(xp[i-3],xp[i-2]);
	    
      k21 = h*f_dloc1(xp[i-3]+k11,xp[i-2]+k12);
      k22 = h*f_dloc2(xp[i-3]+k11,xp[i-2]+k12);

      k31 = h*f_dloc1(xp[i-3]+k21,xp[i-2]+k22);
      k32 = h*f_dloc2(xp[i-3]+k21,xp[i-2]+k22);

      k41 = h*f_dloc1(xp[i-3]+2*k31,xp[i-2]+2*k32);
      k42 = h*f_dloc2(xp[i-3]+2*k31,xp[i-2]+2*k32);

      xp[i] = xp[i-3] + c1*(k11 + 2*k21 + 2*k31 + k41); 
      xp[i+1] = xp[i-2] + c1*(k12 + 2*k22 + 2*k32 + k42); 

      xp[i+2] = 0;

      if ( tref_counter > 0 ){
	tref_counter = tref_counter -1;
        
	/*override the vm value computed using Runge-Kutta*/
        xp[i] = vreset;
      }

      if ( xp[i] >= vth ){
        /*spike*/
        xp[i+2] = 1;
        xp[i+1] = xp[i-2] + caspk;
        tref_counter = tref_ind;
      }

	
    }

      return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *xp, *tp, *tmp; 
    unsigned int m,n; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 10) { 
	mexErrMsgTxt("Ten input arguments required."); 
    } else if (nlhs != 2) {
	mexErrMsgTxt("Two output arguments required."); 
    } 
    
    n = (int) floor((t_end-t_start)/dt) + 1;
    

    /* Create a matrix for the return arguments */ 
    T_OUT = mxCreateDoubleMatrix(1, n, mxREAL);
    X_OUT = mxCreateDoubleMatrix(3, n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 


    tp = mxGetPr(T_OUT);
    xp = mxGetPr(X_OUT);
   
    tmp =  mxGetPr(IE_IN); 
    ie = tmp[0];
    
    tmp = mxGetPr(RIN_IN);
    rin = tmp[0];

    tmp = mxGetPr(T0_IN); 
    t0 = tmp[0];

    tmp = mxGetPr(VTH_IN); 
    vth = tmp[0];

    tmp = mxGetPr(VRE_IN);
    vreset = tmp[0];

    tmp = mxGetPr(VK_IN);
    vk = tmp[0];

    tmp = mxGetPr(TREF_IN);
    tref = tmp[0];

    tmp = mxGetPr(TCA_IN);
    tca = tmp[0];

    tmp = mxGetPr(CASP_IN);
    caspk = tmp[0];

    tmp = mxGetPr(GAHP_IN);
    gahpb = tmp[0];


    /* Do the actual computations in a subroutine */
        
    lif_ad_comp(tp,xp); 
    
    return;
    
}




