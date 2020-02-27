#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define EPSILON 0.0000001

#ifndef enc_type
#define enc_type unsigned char
#endif
#define PI 3.14159265


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAIN /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    (void)nlhs;
    (void)nrhs;
    
    
    /* INPUTS */
    //RGB image
    int* mat_adj = (int*) mxGetPr(prhs[0]);
    float* SP_center = (float*) mxGetPr(prhs[1]);
    int SP_nbr = (int) mxGetScalar(prhs[2]);
    int SuperPatch_R = (int) mxGetScalar(prhs[3]);
    
    int k, l;
    float pos;
    float pos_k[2];
    float pos_l[2];
    
    int dims[2];
    dims[0] = SP_nbr;
    dims[1] = SP_nbr;
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int *mat_adj_out = (int*)mxGetPr(plhs[0]);
    for (k=0; k<SP_nbr*SP_nbr; k++) {
        mat_adj_out[k] = mat_adj[k];
    }
    
    plhs[1] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *mat_neigh_angle = (float*)mxGetPr(plhs[1]);
    
    
    for (k=0; k<SP_nbr; k++) {
        
        pos_k[0] = SP_center[k];
        pos_k[1] = SP_center[k+SP_nbr];
        
        for (l=0; l<SP_nbr; l++) {
            
            if (k != l) {
                
                pos_l[0] = SP_center[l];
                pos_l[1] = SP_center[l+SP_nbr];
                
                if (mat_adj_out[k + l*SP_nbr] == 1) {
                    float angle_l = atan2(-(pos_l[0]-pos_k[0]), pos_l[1]-pos_k[1])*180/PI; 
                    if (angle_l < 0)
                        angle_l += 360.0;
                    mat_neigh_angle[k + l*SP_nbr] = angle_l;
                }
                
                
                //Within the SuperPatch radius
                pos = (pos_k[0]-pos_l[0])*(pos_k[0]-pos_l[0])
                + (pos_k[1]-pos_l[1])*(pos_k[1]-pos_l[1]);
                
                if (pos < (SuperPatch_R*SuperPatch_R)) {
                    
                    //angle computation
                    if (mat_adj_out[k + l*SP_nbr] == 1)
                        mat_adj_out[k + l*SP_nbr] = 2;
                    
                    if (mat_adj_out[k + l*SP_nbr] == 0)
                        mat_adj_out[k + l*SP_nbr] = 3;
                    
                }
                
            }
        }
        
    }
    
    
}

