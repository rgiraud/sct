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



void prod_mat(float *a, float * b, int ni, int mi, int nj, int mj, float * out)
{
    
    for (int i=0; i<ni; i++) {
        for (int j=0; j<mj; j++)  {
            out[i+j*ni] = 0;
            for (int k=0; k<mi; k++)
                out[i+j*ni] += a[i+k*ni]*b[k+j*nj];
        }
    }
    
}

float prod_vect_mat_vect(float *vect, float * mat, int n)
{
    
    float * tmp_mat = (float *) malloc(n*sizeof(float));
    prod_mat(mat, vect, n, n, n, 1, tmp_mat);
    
    float ans = 0;
    for (int i=0; i<n; i++) {
        ans += vect[i]*tmp_mat[i];
    }
    
    free(tmp_mat);
    
    return ans;
    
}


typedef struct{
    float *img ;
    int * lab_mapa;
    float sigma_color;
    int sp_nbra;
    int * SP_nnf;
    float * Sj;
    int * mat_adj;
    float * cov_sp_a;
    float * Sj_b;
    float * a_c;
    float SP_radius;
    int sp_nbrb;
    int wa;
    int ha;
    int ini;
    int fin;
}ct_struct;




void*
        optimal_core(void *arg)
{
    ct_struct inputs = *(ct_struct*) arg;
    float *img     = inputs.img;
    int * lab_mapa    = inputs.lab_mapa;
    float sigma_color = inputs.sigma_color;
    int sp_nbra = inputs.sp_nbra;
    int * SP_nnf = inputs.SP_nnf;
    float * Sj = inputs.Sj;
    float * cov_sp_a = inputs.cov_sp_a;
    float * Sj_b = inputs.Sj_b;
    float SP_radius = inputs.SP_radius;
    float * a_c = inputs.a_c;
    int sp_nbrb = inputs.sp_nbrb;
    
    
    int wa             = inputs.wa;
    int ha             = inputs.ha;
    int ini             = inputs.ini;
    int fin             = inputs.fin;
    int ha_wa = ha*wa;
    
    float * p_i = (float *) malloc(5*sizeof(float));
    float * p_ij = (float *) malloc(5*sizeof(float));
    float * cov_k = (float *) malloc(5*5*sizeof(float));
    float * w_sp = (float *) malloc(sp_nbra*sizeof(float));
    
    for (int j=ini; j<fin; j++) {
        for (int i=0; i<ha; i++) {
            
            int pos = i+ha*j;
            int lab = lab_mapa[pos]-1;
            
            //Build vect
            p_i[0] = (float)i/(float)ha;
            p_i[1] = (float)j/(float)wa;
            p_i[2] = img[pos]*sigma_color;
            p_i[3] = img[pos+ha*wa]*sigma_color;
            p_i[4] = img[pos+ha*wa*2]*sigma_color;
            
            
            //Build cov matrix
            for (int m=0; m<5; m++) {
                for (int n=0; n<5; n++) {
                    cov_k[n+m*5] = cov_sp_a[n+m*5 + lab*25];
                }
            }
            
            
            float weight = 0;
            float min_wij = FLT_MAX;
            for (int k=0; k<sp_nbra; k++) {
                
                if ( sqrt((Sj[lab]-Sj[k])*(Sj[lab]-Sj[k]) + (Sj[lab+sp_nbra]-Sj[k+sp_nbra])*(Sj[lab+sp_nbra]-Sj[k+sp_nbra]) ) < SP_radius) {
                    
                    p_ij[0] = p_i[0] - Sj[k];
                    p_ij[1] = p_i[1] - Sj[k+sp_nbra];
                    p_ij[2] = p_i[2] - Sj[k+sp_nbra*2];
                    p_ij[3] = p_i[3] - Sj[k+sp_nbra*3];
                    p_ij[4] = p_i[4] - Sj[k+sp_nbra*4];
                    
                    float wij = prod_vect_mat_vect(p_ij, cov_k, 5);
                    
                    w_sp[k] = wij;
                    if (wij < min_wij)
                        min_wij = wij;
                    
                    
                }
            }
            
            for (int k=0; k<sp_nbra; k++) {
                
                if ( sqrt((Sj[lab]-Sj[k])*(Sj[lab]-Sj[k]) + (Sj[lab+sp_nbra]-Sj[k+sp_nbra])*(Sj[lab+sp_nbra]-Sj[k+sp_nbra]) ) <= SP_radius) {
                    
                    
                    float wij = exp(1-w_sp[k]/(min_wij+0.00001));
                    
                    int ind_b = SP_nnf[k] - 1;
                    int bt = SP_nnf[k+sp_nbra] - 1;
                    
                    int pos_b = ind_b + bt*sp_nbrb*5;
                    a_c[pos] += wij*Sj_b[pos_b + 2*sp_nbrb]/sigma_color;
                    a_c[pos+ha_wa] += wij*Sj_b[pos_b + 3*sp_nbrb]/sigma_color;
                    a_c[pos+2*ha_wa] += wij*Sj_b[pos_b + 4*sp_nbrb]/sigma_color;
                    
                    weight += wij;
                    
                }
            }
            
            //Normalization
            a_c[pos] /= weight;
            a_c[pos+ha_wa] /= weight;
            a_c[pos+2*ha_wa] /= weight;
            
            
            
        }
        
        
        
    }
    
    free(cov_k);
    free(p_i);
    free(p_ij);
    
    
    pthread_exit(0);
    
}






void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    
    //Inputs
    float* img = (float*) mxGetPr(prhs[0]);
    
    int idx = 1;
    int * lab_mapa = (int *) mxGetPr(prhs[idx++]);
    int * SP_nnf = (int *) mxGetPr(prhs[idx]);
    const int *SP_nnf_dims = mxGetDimensions(prhs[idx++]);
    int sp_nbra = SP_nnf_dims[0];
    float * Sj = (float *) mxGetPr(prhs[idx++]);
    float * cov_sp_a = (float *) mxGetPr(prhs[idx++]);
    float * Sj_b = (float *) mxGetPr(prhs[idx]);
    const int *SP_nnfb_dims = mxGetDimensions(prhs[idx++]);
    int sp_nbrb = SP_nnfb_dims[0];
    float sigma_color = (float) mxGetScalar(prhs[idx++]);
    float SP_radius = (float) mxGetScalar(prhs[idx++]);
    
    
    const int *img_dims = mxGetDimensions(prhs[0]);
    int ha = img_dims[0];
    int wa = img_dims[1];
    int za = img_dims[2];
    
    
    int dims[3];
    dims[0] = ha;
    dims[1] = wa;
    dims[2] = za;
    plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
    float* a_c = (float*)mxGetPr(plhs[0]);
    for (int i=0; i<ha*wa*za; i++)
        a_c[i] = 0;
    
    
    //Thread stuff
    int thread_nbr = get_nprocs_conf();
    if (thread_nbr == 0)
        thread_nbr = 1;
    
    int step = floor(wa/thread_nbr);
    int * slice_vect = (int *) calloc(thread_nbr+1, sizeof(int));
    
    slice_vect[0] = 0;
    for (int i=1; i<thread_nbr; i++)
        slice_vect[i] = i*step;
    slice_vect[thread_nbr] = wa;
    
    
    //Thread argument structures
    pthread_t*   thread_list = (pthread_t*) calloc(thread_nbr, sizeof(pthread_t));
    ct_struct* thread_args = (ct_struct*)calloc(thread_nbr, sizeof(ct_struct));
    
    srand(time(NULL));
    //Launching of the THREADS
    for (int i=0; i < thread_nbr; i++) {
        //Thread arguments
        
        thread_args[i].img = img;
        thread_args[i].lab_mapa = lab_mapa;
        thread_args[i].ha = ha;
        thread_args[i].wa = wa;
        thread_args[i].sigma_color = sigma_color;
        thread_args[i].sp_nbra = sp_nbra;
        thread_args[i].SP_nnf = SP_nnf;
        thread_args[i].Sj = Sj;
        thread_args[i].cov_sp_a = cov_sp_a;
        thread_args[i].Sj_b = Sj_b;
        thread_args[i].sp_nbrb = sp_nbrb;
        thread_args[i].SP_radius = SP_radius;
        thread_args[i].a_c = a_c;
        
        thread_args[i].ini = slice_vect[i];
        thread_args[i].fin = slice_vect[i+1];
        
        if (pthread_create(&thread_list[i], NULL, optimal_core, &thread_args[i]))
            printf("Error creating a thread!\n");
        
    }
    
    /*Wait for all threads to end*/
    for (int j=0; j<thread_nbr; j++) {
        pthread_join(thread_list[j],NULL);
    }
    
    free(slice_vect);
    
}


