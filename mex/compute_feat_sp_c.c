#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"
#include "pthread.h"

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



typedef struct{
    double *img ;
    int * lab_map;
    double sigma_color;
    double * cov_sp;
    double * feat_sp;
    int w;
    int h ;
    int ini ;
    int fin ;
    double * count;
}ct_struct;




void*
        optimal_core(void *arg)
{
    ct_struct inputs = *(ct_struct*) arg;
    double *img     = inputs.img;
    int * lab_map    = inputs.lab_map;
    double sigma_color = inputs.sigma_color;
    double * feat_sp = inputs.feat_sp;
    double * cov_sp = inputs.cov_sp;
    double * count = inputs.count;
    
    int wa             = inputs.w;
    int ha             = inputs.h;
    int ini             = inputs.ini;
    int fin             = inputs.fin;
    int ha_wa = ha*wa;
    
    double * tmp = (double*) calloc(5,sizeof(double));
    
    
    for (int j=ini; j<fin; j++) {
        for (int i=0; i<ha; i++) {
            
            int pos = i+ha*j;
            int lab = lab_map[pos]-1;
            
            tmp[0] = (float)i/(float)ha;
            tmp[1] = (float)j/(float)wa;
            tmp[2] = img[pos]*sigma_color;
            tmp[3] = img[pos+ha_wa]*sigma_color;
            tmp[4] = img[pos+ha_wa*2]*sigma_color;
            
            
            int feat_offset = lab*5;
            for (int k=0; k<5; k++)
                feat_sp[k+feat_offset] += tmp[k];
            
            
            int lab_offset = 25*lab;
            for (int k=0; k<5; k++) {
                for (int kk=k; kk<5; kk++) {
                    cov_sp[k + 5*kk + lab_offset] += tmp[k]*tmp[kk];
                }
            }
            count[lab] += 1;
            
        }
    }
    
    free(tmp);
    
    pthread_exit(0);
}





//////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAIN /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    
    //Inputs
    double* img = (double*) mxGetPr(prhs[0]);
    
    int idx = 1;
    int * lab_map = (int *) mxGetPr(prhs[idx++]);
    int sp_nbr = (int) mxGetScalar(prhs[idx++]);
    double sigma_color = (double) mxGetScalar(prhs[idx++]);
    
    
    const int *img_dims = mxGetDimensions(prhs[0]);
    int h = img_dims[0];
    int w = img_dims[1];
    
    int dims[3];
    dims[0] = 5;
    dims[1] = sp_nbr;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* feat_sp = (double*)mxGetPr(plhs[0]);
    for (int i=0; i<sp_nbr*5; i++)
        feat_sp[i] = 0;
    
    dims[0] = 5;
    dims[1] = 5;
    dims[2] = sp_nbr;
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* cov_sp = (double*)mxGetPr(plhs[1]);
    for (int i=0; i<5*5*sp_nbr; i++)
        cov_sp[i] = 0;
    
    
    double * count = (double *) calloc(sp_nbr,sizeof(double));
    
    
    //Thread stuff
    int thread_nbr = get_nprocs_conf();
    if (thread_nbr == 0)
        thread_nbr = 1;
    
    int step = floor(w/thread_nbr);
    int * slice_vect = (int *) calloc(thread_nbr+1, sizeof(int));
    
    slice_vect[0] = 0;
    for (int i=1; i<thread_nbr; i++)
        slice_vect[i] = i*step;
    slice_vect[thread_nbr] = w;
    
    
    //Thread argument structures
    pthread_t*   thread_list = (pthread_t*) calloc(thread_nbr, sizeof(pthread_t));
    ct_struct* thread_args = (ct_struct*)calloc(thread_nbr, sizeof(ct_struct));
    
    srand(time(NULL));
    //Launching of the THREADS
    for (int i=0; i < thread_nbr; i++) {
        
        //Thread arguments
        thread_args[i].img = img;
        thread_args[i].lab_map = lab_map;
        thread_args[i].h = h;
        thread_args[i].w = w;
        thread_args[i].sigma_color = sigma_color;
        thread_args[i].feat_sp = feat_sp;
        thread_args[i].cov_sp = cov_sp;
        thread_args[i].count = count;
        
        thread_args[i].ini = slice_vect[i];
        thread_args[i].fin = slice_vect[i+1];
        
        if (pthread_create(&thread_list[i], NULL, optimal_core, &thread_args[i]))
            printf("Error creating a thread!\n");
        
    }
    
    /*Wait for all threads to end*/
    for (int j=0; j<thread_nbr; j++) {
        pthread_join(thread_list[j],NULL);
    }
    
    
    //Normalization
    for (int i=0; i< sp_nbr; i++){
        int feat_offset = i*5;
        for (int k=0; k<5; k++)
            feat_sp[k+feat_offset] /= count[i];
        
        int lab_offset = 25*i;
        for (int k=0; k< 5; k++) {
            for (int kk=k; kk< 5; kk++) {
                cov_sp[k + 5*kk + lab_offset]  =  cov_sp[k + 5*kk + lab_offset]/count[i] - feat_sp[k+feat_offset]*feat_sp[kk+feat_offset];
            }
        }
        
        //Parallel copy covariance matrix
        for (int k=1; k< 5; k++) {
            for (int kk=k-1; kk>-1; kk--) {
                cov_sp[k + 5*kk + lab_offset] = cov_sp[kk + 5*k + lab_offset];
            }
        }
    }
    
    
    
    
    free(slice_vect);
    free(count);
    
    
    
    
}


