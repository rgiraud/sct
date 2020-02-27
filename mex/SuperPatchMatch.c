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




float dist_angle(float a, float b) {
    return MIN(MIN(fabsf(a-b),fabsf(-(360-a)-b)),fabsf(-(360-b)-a));
}



void descriptor_extraction_a(float * desc_SP,  int i, int h_bins, int *mat_adj, float *SP_hist, float* SP_center, int SP_nbr) {
    
    int j, k, l;
    
    for (l=0; l<h_bins; l++)
        desc_SP[l] = SP_hist[i + l*SP_nbr];
    
    desc_SP[h_bins] = SP_center[i];
    desc_SP[h_bins+1] = SP_center[i + SP_nbr];
    
    j = 1;
    for (k=0; k<SP_nbr; k++) {
        if (mat_adj[i + k*SP_nbr] > 1) {
            
            for (l=0; l<h_bins; l++)
                desc_SP[l+j*(h_bins+2)] = SP_hist[k + l*SP_nbr];
            
            desc_SP[j*(h_bins+2) + h_bins] = SP_center[k];
            desc_SP[j*(h_bins+2) + h_bins+1] = SP_center[k + SP_nbr];
            
            j += 1;
        }
    }
    
    
}


void list_comp_rgb(float *cost, float *desc_SP_i, float * desc_SP_j, int j,  int t, float * SP_histb,
        int l_i, int* lab_mapb, int SP_nbrb, int h_bins, int h, int w, float sigma_spa, float *SP_centerb) {
    
    
    int ii,l,k; 
    float x, y;
    
    for (l=0; l<h_bins; l++)
        desc_SP_j[l] = SP_histb[j + SP_nbrb*(l + t*h_bins)];
    
    desc_SP_j[h_bins] = SP_centerb[j  + t*SP_nbrb*2];
    desc_SP_j[h_bins+1] = SP_centerb[j + SP_nbrb*(1 + t*2)];
    
    
    // Difference between SP
    int i_off; 
    float central_w = 1;
    
    float vij[2]; //ci - cj;
    vij[0] = desc_SP_i[h_bins] - desc_SP_j[h_bins];
    vij[1] = desc_SP_i[h_bins+1] - desc_SP_j[h_bins+1];
    
    float dist = 0;
    
    //Central SP
    for (l = 0; l<h_bins; l++) {
        dist += (desc_SP_i[l] - desc_SP_j[l])*(desc_SP_i[l] - desc_SP_j[l]);
    }
    dist *= central_w;
    
    //Approx distance to the closest match
    float dist_spa;
    float dist_tmp;
    for (ii=1; ii<l_i; ii++) {
        i_off = (h_bins+2)*ii;
        y = desc_SP_i[i_off+h_bins] - vij[0];
        x = desc_SP_i[i_off+h_bins+1] - vij[1];
        k = lab_mapb[(int) MIN(MAX(y,0),h-1) + (int) MIN(MAX(x,0),w-1)*h + t*h*w] -1;
        
        dist_tmp = 0;
        for (l = 0; l<h_bins; l++) {
            dist_tmp += (desc_SP_i[l+i_off] - SP_histb[k + SP_nbrb*(l + t*h_bins)])*(desc_SP_i[l+i_off] - SP_histb[k + SP_nbrb*(l + t*h_bins)]);
        }
        
        //Distance weight within superpatch
        dist_spa =  (desc_SP_i[i_off+h_bins] - desc_SP_i[h_bins])*(desc_SP_i[i_off+h_bins] - desc_SP_i[h_bins]) +
                (desc_SP_i[i_off+h_bins+1] - desc_SP_i[h_bins+1])*(desc_SP_i[i_off+h_bins+1] - desc_SP_i[h_bins+1]);
        dist_spa = exp(-sqrt(dist_spa)/(sigma_spa*sigma_spa));
        dist += dist_tmp*dist_spa;
        
        
    }
    if (l_i > 0)
        dist /= l_i;
    
    *cost = dist;
    
}





void association_decision(float *desc_SP_k, float * desc_SP_j, int j, int t_n, float * SP_hista, float * SP_histb,
        int *l_vecta, int* lab_mapb, int SP_nbra, int SP_nbrb_max, int h_bins, int hb, int wb, float sigma,
        int * mat_spa, float * SP_centera, float *SP_centerb,
        int i, int * nnf, float * nnfd, int offset_mp, int * M_use,
        int *bestj, int *bestt, float *bestd, float max_use_sp, float dist) {
    
    
    
    //Number of use of Bk
    float tmp = M_use[j + SP_nbrb_max*t_n];
    
    //If already too much use, switch possibility
    if (tmp >= max_use_sp) {
        
        int j_switch = 0;
        float max_dj = -10000000;
        float dist_dj = 0;
        for (int ll = 0; ll < SP_nbra; ll++) {
            if ( (ll != i) && (nnf[ll + 2*offset_mp] == j)) {  //SP de A qui pointe aussi vers SP j
                
                //cost of assigning SP in Aj to old match of Ai
                descriptor_extraction_a(desc_SP_k, ll, h_bins, mat_spa, SP_hista, SP_centera, SP_nbra);
                list_comp_rgb(&dist_dj, desc_SP_k, desc_SP_j, *bestj, *bestt, SP_histb, l_vecta[ll]+1, 
                        lab_mapb, SP_nbrb_max, h_bins, hb, wb, sigma, SP_centerb);
                
                if ( (nnfd[ll + offset_mp] - dist_dj) > max_dj) {
                    max_dj = nnfd[ll + offset_mp] - dist_dj;
                    j_switch = ll;
                }
                
                
            }
        }
        
        //If energy of switch is > 0
        if ( ((*bestd - dist) + (max_dj)) > 0 ) {
            
            //Aj -> old match of Ai
            nnf[j_switch  + offset_mp*2] = *bestj;
            nnf[j_switch  + SP_nbra + offset_mp*2] = *bestt;
            nnfd[j_switch + offset_mp] -= max_dj;
            
            *bestj = j;
            *bestt = t_n;
            *bestd = dist;
            
            nnf[i + offset_mp*2] = *bestj;
            nnf[i + SP_nbra + offset_mp*2] = *bestt;
            nnfd[i + offset_mp] = *bestd;
            
        }
        
        
        
    }
    else {  //new assignement
        
    
        //Fill M
        M_use[*bestj + SP_nbrb_max*(*bestt)] -= 1;
        M_use[j + SP_nbrb_max*t_n] += 1;
      
        *bestj = j;
        *bestt = t_n;
        *bestd = dist;
        nnf[i + offset_mp*2] = *bestj;
        nnf[i + offset_mp*2 +SP_nbra] = *bestt;
        nnfd[i + offset_mp] = *bestd;
        
        
    }
    
    
}




typedef struct{
    int * nnf;
    float * nnfd;
    float * nnfd_rand;
    int SP_nbra;
    int SP_nbrb_max;
    int * mat_spa;
    float * SP_hista;
    float * SP_histb;
    int h_bins;
    float * SP_centera;
    float * SP_centerb;
    float sigma;
    int * lab_mapb;
    int * l_vecta;
    int l_maxa;
    int l_maxb;
    int hb;
    int wb;
    int nt;
    int iter;
    int offset_mp;
    unsigned long next;
    int * lib_img_dims;
    int * tab_adja;
    int * tab_adjb;
    int * tab_adja_n;
    int * tab_adjb_n;
    float max_use_sp;
}pm_struct;



int
        pm_rand(unsigned long* next)
{
    (*next) = (*next) * 1103515245 + 12345;
    return (unsigned int)((*next)/65536) % 32768;
}


void*
        pm_hist_rgb(void *arg)
{
    pm_struct inputs    = *(pm_struct*) arg;
    int * nnf           = inputs.nnf;
    float * nnfd        = inputs.nnfd;
    int SP_nbra         = inputs.SP_nbra;
    int SP_nbrb_max     = inputs.SP_nbrb_max;
    int * mat_spa      = inputs.mat_spa;
    float * SP_hista    = inputs.SP_hista;
    float * SP_histb    = inputs.SP_histb;
    int h_bins          = inputs.h_bins;
    float * SP_centera  = inputs.SP_centera;
    float * SP_centerb  = inputs.SP_centerb;
    float sigma         = inputs.sigma;
    int * lab_mapb      = inputs.lab_mapb;
    int * l_vecta       = inputs.l_vecta;
    int l_maxa          = inputs.l_maxa;
    int l_maxb          = inputs.l_maxb;
    int hb              = inputs.hb;
    int wb              = inputs.wb;
    int nt              = inputs.nt;
    int iter            = inputs.iter;
    int offset_mp       = inputs.offset_mp;
    unsigned long next  = inputs.next;
    int * lib_img_dims  = inputs.lib_img_dims;
    int * tab_adja  = inputs.tab_adja;
    int * tab_adjb  = inputs.tab_adjb;
    int * tab_adja_n  = inputs.tab_adja_n;
    int * tab_adjb_n  = inputs.tab_adjb_n;
    float max_use_sp  = inputs.max_use_sp;
    float * nnfd_rand = inputs.nnfd_rand;
    
    float dist = 0;
    int i,j,ii,i_n,j_n,t_n,pos, pos_nnfd;
    float xy_i[2];
    int size_b = hb*wb;
    
    float * desc_SP_i = (float *)calloc((l_maxa+1)*(h_bins+2), sizeof(float));
    float * desc_SP_k = (float *)calloc((l_maxa+1)*(h_bins+2), sizeof(float));
    float * desc_SP_j = (float *)calloc((l_maxb+1)*(h_bins+2), sizeof(float));
    int * M_use = (int *)calloc(SP_nbrb_max*nt, sizeof(int));
    
    
    //int ay, ax;
    int xmin, ymin, xmax, ymax, bx, by, loop;
    
    //INIT
    for (i=0; i<SP_nbra; i++) {
        
        loop = 1;
        while (loop) {
            
            pos = i + offset_mp*2;
            pos_nnfd = i + offset_mp;
           
            int rand_val = pm_rand(&next);
            int t = rand_val%(nt);
            
            //Max dims of b
            int ymax_b = lib_img_dims[t];
            int xmax_b = lib_img_dims[t+nt];
            
            //Random match     
            xmin = 0;          
            xmax = xmax_b-1;   
            ymin = 0;          
            ymax = ymax_b-1;   
            
            rand_val = pm_rand(&next);
            bx = xmin+rand_val%(xmax-xmin);
            rand_val = pm_rand(&next);
            by = ymin+rand_val%(ymax-ymin);
            
            j = lab_mapb[(int) by + (int) bx*hb + t*size_b] - 1;
            
            if (M_use[j + SP_nbrb_max*t] < max_use_sp) {
                
                nnf[pos] = j;
                nnf[pos + SP_nbra] = t;
                
                //Fill M
                M_use[j + SP_nbrb_max*t] += 1;
                
                loop = 0;
            }
            
            
        }
    }
    
    // Distance computation
    for (i=0; i<SP_nbra; i++) {
        
        pos = i + offset_mp*2;
        pos_nnfd = i + offset_mp;
        
        int j = nnf[pos];
        int t = nnf[pos + SP_nbra];
        
        descriptor_extraction_a(desc_SP_i, i, h_bins, mat_spa, SP_hista, SP_centera, SP_nbra);
        list_comp_rgb(&dist, desc_SP_i, desc_SP_j, j, t, SP_histb, l_vecta[i]+1, lab_mapb,
                SP_nbrb_max, h_bins, hb, wb, sigma, SP_centerb);
        
        
        nnfd[pos_nnfd] = dist;
        nnfd_rand[pos_nnfd] = dist;
        
    }
    
    

    // MAIN LOOP
    int count = 0;
    while ((count < iter)) {
        
        int i_change = 1;
        int i_start = 0;
        int i_end = SP_nbra-1;
        
        if (count % 2 == 1) {
            i_start = SP_nbra-1;
            i_end = -1;
            i_change = -1;
            
        }
        
        
        for (int i=i_start; i!=i_end; i += i_change) {
            
            //Current position in nnf
            pos = i + offset_mp*2;
            pos_nnfd = i + offset_mp;
            
            descriptor_extraction_a(desc_SP_i, i, h_bins, mat_spa, SP_hista, SP_centera, SP_nbra);
            xy_i[0] = SP_centera[i]; //y
            xy_i[1] = SP_centera[i+SP_nbra]; //x
            
            int jbest = nnf[pos];
            int *bestj = &jbest;
            int tbest = nnf[pos+SP_nbra];
            int *bestt = &tbest;
            float dbest = nnfd[pos_nnfd];
            float *bestd = &dbest;
            
            
            // PROPAGATION: Improve current guess by trying
            //instead correspondences from left and above
            //(below and right on odd iterations).
            for (ii=0; ii<tab_adja_n[i]; ii++) {
                
                i_n = tab_adja[i+ii*SP_nbra];  //subscript of the neighbor
                
                if ( (i-i_n)*i_change > 0) {
                    
                    j_n = (int) nnf[i_n + offset_mp*2];  //subscript of the neighbor's match
                    t_n = (int) nnf[i_n + SP_nbra + offset_mp*2];  //Template of the neighbor's match
                    
                    //j is the index in t_n in B of the candidate for SP i
                    //Selection of a random candidate among the adjacent neighbors of j
                    j = tab_adjb[j_n + SP_nbrb_max*((int) pm_rand(&next)%tab_adjb_n[j_n+t_n*SP_nbrb_max])];
                    

                    //cost of association Ai -> Bk
                    list_comp_rgb(&dist, desc_SP_i, desc_SP_j, j, t_n, SP_histb, l_vecta[i]+1, lab_mapb, 
                            SP_nbrb_max, h_bins, hb, wb, sigma, SP_centerb);
                    
                    
                    if (dist < *bestd) {
                        association_decision(desc_SP_k, desc_SP_j, j, t_n, SP_hista, SP_histb, l_vecta, 
                                lab_mapb, SP_nbra, SP_nbrb_max, h_bins, hb, wb, sigma, mat_spa, SP_centera, SP_centerb,
                                i, nnf, nnfd, offset_mp, M_use, bestj, bestt, bestd, max_use_sp, dist);
                    }
                    
                    
                }
            }
            
            ////////////////////////////////////////////////////////
            //RANDOM SEARCH : fixed template
            //////////////////////////////////////////////////////
            
            int ybest = (int) SP_centerb[*bestj + *bestt*SP_nbrb_max*2]; //y
            int xbest = (int) SP_centerb[*bestj + SP_nbrb_max*(*bestt*2 + 1)];  //x
            tbest = *bestt;
            
            int ymax_b = lib_img_dims[tbest];
            int xmax_b = lib_img_dims[tbest+nt];
            
            // Sampling window
            for (int mag = MAX(ymax_b,xmax_b); mag >= 3; mag /= 2) {
                
                /* Box limits */
                int xmin = MAX(xbest-mag, 0);
                int xmax = MIN(xbest+mag+1, xmax_b);
                if(xmin == xmax) continue;
                
                int ymin = MAX(ybest-mag, 0);
                int ymax = MIN(ybest+mag+1, ymax_b);
                if(ymin == ymax) continue;
                
                int rand_val = pm_rand(&next);
                int bx = xmin+rand_val%(xmax-xmin);
                rand_val = pm_rand(&next);
                int by = ymin+rand_val%(ymax-ymin);
                
                j = lab_mapb[by + bx*hb + tbest*size_b] - 1;
                
                
                list_comp_rgb(&dist, desc_SP_i, desc_SP_j, j, tbest, SP_histb, l_vecta[i]+1, lab_mapb, 
                        SP_nbrb_max, h_bins, hb, wb, sigma, SP_centerb);
                
                
                if (dist < *bestd) {
                    
                    association_decision(desc_SP_k, desc_SP_j, j, tbest, SP_hista, SP_histb, l_vecta, 
                                lab_mapb, SP_nbra, SP_nbrb_max, h_bins, hb, wb, sigma, mat_spa, SP_centera, SP_centerb,
                                i, nnf, nnfd, offset_mp, M_use, bestj, bestt, bestd, max_use_sp, dist);
                    
                    
                }
                
            }     
            
        }
        
        
        count++;
    }
    
    
    free(desc_SP_i);
    free(desc_SP_j);
    free(desc_SP_k);
    free(M_use);
    
    
    pthread_exit(0);
}





//////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAIN /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    (void)nlhs;
    (void)nrhs;
    
    
    
    /* INPUTS */
    //RGB image
    int* mat_spa = (int*) mxGetPr(prhs[0]);
    int* mat_spb = (int*) mxGetPr(prhs[1]);
    
    int idx = 2;
    float* SP_hista = (float*) mxGetPr(prhs[idx++]);
    float* SP_histb = (float*) mxGetPr(prhs[idx++]);
    int h_bins = (int) mxGetScalar(prhs[idx++]);
    
    float* SP_centera = (float*) mxGetPr(prhs[idx++]);
    float* SP_centerb = (float*) mxGetPr(prhs[idx++]);
    
    int SP_nbra = (int) mxGetScalar(prhs[idx++]);
    int SP_nbrb_max = (int) mxGetScalar(prhs[idx++]);
    
    int* lab_mapb = (int*) mxGetPr(prhs[idx]);
    int b_dims_nbr = mxGetNumberOfDimensions(prhs[idx]);
    const int* b_dims = mxGetDimensions(prhs[idx++]);
    int hb = b_dims[0];
    int wb = b_dims[1];
    
    int nt;
    if (b_dims_nbr == 2)
        nt = 1;
    else
        nt = b_dims[2];
    
    float sigma = (float) mxGetScalar(prhs[idx++]);
    
    int iter = mxGetScalar(prhs[idx++]);
    int multiple_pm_nbr = mxGetScalar(prhs[idx++]);
    
    
    int* l_vecta = (int *) mxGetPr(prhs[idx++]);
    int l_maxa = (int) mxGetScalar(prhs[idx++]);
    int l_maxb = (int) mxGetScalar(prhs[idx++]);
    
    int* lib_img_dims = (int *) mxGetPr(prhs[idx++]);
    
    float max_use_sp = (float) mxGetScalar(prhs[idx++]);
    
    int i;
    
    
    //Thread argument structures
    pthread_t* thread_list = (pthread_t*) malloc(multiple_pm_nbr*sizeof(pthread_t));
    pm_struct* thread_args =(pm_struct*) malloc(multiple_pm_nbr*sizeof(pm_struct));
    
    
    int dims[3];
    dims[0] = SP_nbra;
    dims[1] = 2;
    dims[2] = multiple_pm_nbr;
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    int * nnf = (int*)mxGetPr(plhs[0]);
    
    
    dims[0] = SP_nbra;
    dims[1] = multiple_pm_nbr;
    plhs[1] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float * nnfd = (float*)mxGetPr(plhs[1]);
    
    for (int i=0; i<SP_nbra*multiple_pm_nbr; i++)
        nnfd[i] = FLT_MAX;
    
    plhs[2] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float * nnfd_rand = (float*)mxGetPr(plhs[2]);
    
    
    
    // Adjacency tabs
    int * tab_adja_n = (int *) calloc(SP_nbra, sizeof(int));
    
    int maxn = 0;
    for (int i=0; i < SP_nbra; i++) {
        int c = 0;
        for (int j=0; j < SP_nbra; j++) {
            if ( (mat_spa[i + j*SP_nbra] == 1) || (mat_spa[i + j*SP_nbra] == 2) )
                c += 1;
        }
        tab_adja_n[i] = c;
        maxn = MAX(maxn,c);
    }
    
    int * tab_adja = (int *) calloc(SP_nbra*maxn, sizeof(int));
    for (int i=0; i < SP_nbra; i++) {
        int c = 0;
        for (int j=0;j < SP_nbra; j++) {
            if ( (mat_spa[i + j*SP_nbra] == 1) || (mat_spa[i + j*SP_nbra] == 2) ) {
                tab_adja[i + SP_nbra*c] = j;
                c += 1;
            }
        }
    }
    
    int * tab_adjb_n = (int *) calloc(SP_nbrb_max*nt, sizeof(int));
    maxn = 0;
    for (int k=0; k < nt; k++) {
        for (int i=0; i < SP_nbrb_max; i++) {
            int c = 0;
            for (int j=0; j < SP_nbrb_max; j++) {
                if ( (mat_spb[i + SP_nbrb_max*(j + k*SP_nbrb_max)] == 1) || (mat_spb[i + SP_nbrb_max*(j + k*SP_nbrb_max)] == 2) )
                    c += 1;
            }
            tab_adjb_n[i+k*SP_nbrb_max] = c;
            maxn = MAX(maxn,c);
        }
    }
    
    int * tab_adjb = (int *) calloc(SP_nbrb_max*maxn*nt, sizeof(int));
    for (int k=0; k < nt; k++) {
        for (int i=0; i < SP_nbrb_max; i++) {
            int c = 0;
            for (int j=0;j < SP_nbrb_max; j++) {
                if ( (mat_spb[i + SP_nbrb_max*(j + k*SP_nbrb_max)] == 1) || (mat_spb[i + SP_nbrb_max*(j + k*SP_nbrb_max)] == 2) ) {
                    tab_adjb[i + SP_nbrb_max*(c + k*maxn)] = j;
                    c += 1;
                }
            }
        }
    }
    
    
    
    unsigned long next = 1789;  //seed for same (random) sequence
    
    //Launching of the THREADS
    for (int i=0; i < multiple_pm_nbr; i++) {
        
        //Thread arguments
        thread_args[i].SP_nbra = SP_nbra;
        thread_args[i].SP_nbrb_max = SP_nbrb_max;
        thread_args[i].mat_spa = mat_spa;
        thread_args[i].SP_centera = SP_centera;
        thread_args[i].SP_centerb = SP_centerb;
        thread_args[i].sigma = sigma;
        thread_args[i].lab_mapb = lab_mapb;
        
        thread_args[i].tab_adja = tab_adja;
        thread_args[i].tab_adjb = tab_adjb;
        thread_args[i].tab_adja_n = tab_adja_n;
        thread_args[i].tab_adjb_n = tab_adjb_n;
        
        thread_args[i].max_use_sp = max_use_sp;
        thread_args[i].lib_img_dims = lib_img_dims;
        
        thread_args[i].l_vecta = l_vecta;
        thread_args[i].l_maxa = l_maxa;
        thread_args[i].l_maxb = l_maxb;
        thread_args[i].hb = hb;
        thread_args[i].wb = wb;
        thread_args[i].nt = nt;
        
        thread_args[i].iter = iter;
        thread_args[i].offset_mp = (i%multiple_pm_nbr)*SP_nbra;
        
        next = pm_rand(&next);
        thread_args[i].next = next;
        
        thread_args[i].SP_hista = SP_hista;
        thread_args[i].SP_histb = SP_histb;
        thread_args[i].h_bins = h_bins;
        thread_args[i].nnf = nnf;
        thread_args[i].nnfd = nnfd;
        thread_args[i].nnfd_rand = nnfd_rand;
        
        
        if (pthread_create(&thread_list[i], NULL, pm_hist_rgb, &thread_args[i]))
            printf("Error creating a thread!\n");
        
    }
    
    for (i=0; i<multiple_pm_nbr; i++) {
        pthread_join(thread_list[i],NULL);
    }
    
    
    free(tab_adja);
    free(tab_adjb);
    free(tab_adja_n);
    free(tab_adjb_n);
    free(thread_args);
    free(thread_list);
    
    
    
}
