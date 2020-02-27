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



//////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAIN /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    (void)nlhs;
    (void)nrhs;
    
    
    /* INPUTS */
    //RGB image
    float* img = (float*) mxGetPr(prhs[0]);
    int SP_nbr = (int) mxGetScalar(prhs[1]);
    int* lab_map = (int*) mxGetPr(prhs[2]);
    int hist_bins = (int) mxGetScalar(prhs[3]);
    
    const int* img_dims = mxGetDimensions(prhs[0]);
    int h = img_dims[0];
    int w = img_dims[1];
    int hw = h*w;
    int hw2 = h*w*2;
    int pos;
    
    int h2 = hist_bins*hist_bins;
    int h3 = hist_bins*h2;
    
    int dims[2];
    dims[0] = SP_nbr;
    dims[1] = 2;
    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *SP_center = (float*)mxGetPr(plhs[0]);
    
    dims[1] = hist_bins*hist_bins*hist_bins;
    plhs[1] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *SP_hist_rgb3D = (float*)mxGetPr(plhs[1]);
    
    int * count_vect = (int *)calloc(SP_nbr, sizeof(int));
    
    
    for (int k=0; k< SP_nbr; k++) {
        //Init
        SP_center[k] = 0;
        SP_center[k+SP_nbr] = 0;
        for (int l=0; l<h3; l++)
            SP_hist_rgb3D[k + l*SP_nbr] = 0;
    }
    
    
    for (int i=0; i<h; i++) {
        for (int j=0; j<w; j++) {
            
            pos = i+j*h;
            int k = lab_map[pos] - 1;
            count_vect[k] += 1;
            
            //3D rgb hist
            int dr = (int) MIN(floor(img[pos]*hist_bins),hist_bins-1);
            int dg = (int) MIN(floor(img[pos+hw]*hist_bins),hist_bins-1);
            int db = (int) MIN(floor(img[pos+hw2]*hist_bins),hist_bins-1);
            SP_hist_rgb3D[k + (dr + dg*hist_bins + db*h2)*SP_nbr] += 1;
            
            
            //Centers
            SP_center[k] += i;
            SP_center[k+SP_nbr] += j;
            
            
        }
    }
    
    //Normalization
    int count;
    for (int k=0; k<SP_nbr; k++) {
        if (count_vect[k] > 0) {
            count = count_vect[k];
            
            SP_center[k] /= count;
            SP_center[k+SP_nbr] /= count;
            
            
            for (int ii=1; ii<hist_bins; ii++) {
                for (int jj=0; jj<hist_bins; jj++) {
                    for (int kk=0; kk<hist_bins; kk++) {
                        SP_hist_rgb3D[k + (ii + jj*hist_bins + kk*h2)*SP_nbr] += SP_hist_rgb3D[k + (ii-1 + jj*hist_bins + kk*h2)*SP_nbr];
                    }
                }
            }
            for (int jj=1; jj<hist_bins; jj++) {
                for (int ii=0; ii<hist_bins; ii++) {
                    for (int kk=0; kk<hist_bins; kk++) {
                        SP_hist_rgb3D[k + (ii + jj*hist_bins + kk*h2)*SP_nbr] += SP_hist_rgb3D[k + (ii + (jj-1)*hist_bins + kk*h2)*SP_nbr];
                    }
                }
            }
            for (int kk=1; kk<hist_bins; kk++) {
                for (int jj=0; jj<hist_bins; jj++) {
                    for (int ii=0; ii<hist_bins; ii++) {
                        SP_hist_rgb3D[k + (ii + jj*hist_bins + kk*h2)*SP_nbr] += SP_hist_rgb3D[k + (ii + jj*hist_bins + (kk-1)*h2)*SP_nbr];
                    }
                }
            }
            
            
            for (int i=0; i<h3; i++)
                SP_hist_rgb3D[k+i*SP_nbr] /= count;
            
            
        }
    }
    
    
    free(count_vect);
    
    
}

