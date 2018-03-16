#include <vector>
#include <cmath>
#include <cassert>
#include <cstdio>

extern "C" {
#include "bicubic_interpolation.h"
}

#include "energy_structures.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////AUXILIAR FUNCTIONS NLTVL1///////////////////////////////
bool positive(int val) {
    return val >= 0;
}

float getsample_0(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return 0;
    return x[w*h*l + j*w + i]; // r g b
    //return x[(i+j*w)*pd + l]; // intrelazado
}

int validate_ap(int w, int h, int i, int j, int di, int dj) {
    const int r = j + dj; //Row
    const int c = i + di; //Colum
    if ( c < 0 || c >= w || r < 0 || r >= h){
        return -1;
    }
    return 0;
}

int validate_ap_patch(
        const int ii,     // initial column
        const int ij,     // initial row
        const int ei,     // end column
        const int ej,     // end row
        int i,
        int j
        ) {
    //ei <= w and ej <=h by set up.
    if ( i < ii || i >= ei || j < ij || j >= ej){
        return -1;
    }
    return 0;
}

float get_wspatial(
        int l,
        int k
        ) {
    float ws = NL_BETA;
    float w_tmp;
    float difS = 0.0;
    difS = hypot(l,k);
    w_tmp = exp(- difS/ws);
    //std::printf("Dif_S: %f W_S: %f\n",difS, w_tmp);

    return w_tmp;
}

float get_wcolor(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
        ) {
    float wi = NL_INTENSITY;
    float w_tmp;
    float difI = 0.0;
    for (int m = 0; m < pd; m++)
    {
        float aux = getsample_0(a, w, h, pd, i, j, m)
                - getsample_0(a, w, h, pd, i + l, j + k, m);
        difI += aux*aux;
        //std::printf("valI_%d:%f ",m, difI);
    }
    //std::printf("\n");
    difI = sqrt(difI);
    w_tmp = exp( -difI/wi);

    return w_tmp;
}

float get_weight(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
        ) {
    float wc = get_wcolor(a, w, h, i, j, l, k, pd);
    float ws = get_wspatial(l, k);
    assert(ws*wc >= 0);
    return wc*ws;
}


//////////////////////////////NLTVL1////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void nltv_ini_dual_variables(
        float *a,
        const int pd,
        const int w,
        const int h,
        const int n_d,
        const int radius,
        DualVariables *p,
        DualVariables *q
        ) {
    int size = w*h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++)
        {
            p[i].sc[j] = -2.0;
            p[i].wp[j] = -2.0;
            p[i].api[j] = -1; //Indicates that is out
            p[i].apj[j] = -1; //Indicate that is out.
            p[i].rp[j] = -1; //Indicates that it is out.

            q[i].sc[j] = -2.0;
            q[i].wp[j] = -2.0;
            q[i].api[j] = -1; //Indicates that is out
            q[i].apj[j] = -1; //Indicate that is out.
            q[i].rp[j] = -1; //Indicates that it is out.
            assert(std::isfinite(p[i].sc[j]));
            assert(std::isfinite(q[i].sc[j]));
        }

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int   it = 0;
            float wt = 0.0;
            const int pos = j*w + i;
            for (int k = -radius; k < (radius + 1); k++)
                for (int l = -radius; l < (radius + 1); l++) {
                    //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
                    if (!(k ==0 && k==l)) {

                        int ap = validate_ap(w, h, i, j, l, k);
                        if (ap == 0) {
                            p[pos].sc[it]  = q[pos].sc[it]  = 0.0;
                            p[pos].api[it] = q[pos].api[it] = i+l;
                            p[pos].apj[it] = q[pos].apj[it] = j+k;
                            p[pos].rp[it] = q[pos].rp[it] = n_d - (it + 1);
                            float wp = sqrt(get_weight(a, w, h, i, j, l, k, pd));
                            p[pos].wp[it] = q[pos].wp[it] = wp;
                            wt += wp;
                            assert(p[pos].api[it] == q[pos].api[it]);
                            assert(p[pos].apj[it] == q[pos].apj[it]);
                            assert(p[pos].rp[it] >= 0);
                            assert(q[pos].rp[it] >= 0);
                            assert(std::isfinite(p[pos].sc[it]));
                            assert(std::isfinite(q[pos].sc[it]));
                        }
                        it++;
                    }
                }
            //TODO: It is used to normalize
            p[pos].wt = wt;
            q[pos].wt = wt;
        }
    // std::printf(" Acaba\n");
}


void non_local_divergence(
        DualVariables *p,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int w,
        int n_d,
        float *div_p
        ) {
    //#pragma omp parallel for
    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            const int i = l*w + k;
            div_p[i] = 0.0;
            for (int j = 0; j < n_d; j++) {
                const int api  = p[i].api[j];
                const int apj  = p[i].apj[j];
                const int rp   = p[i].rp[j];
                const float wp = p[i].wp[j];
                const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);
                if (ap == 0) {
                    assert(ap>=0);
                    assert(rp>=0);
                    assert(wp>=0);
                    const float pxy = p[i].sc[j];
                    const float pyx = p[apj*w + api].sc[rp];
                    assert(std::isfinite(pxy));
                    assert(std::isfinite(pyx));
                    assert(std::isfinite(wp));
                    div_p[i] += wp*(pxy - pyx);
                }
            }
        }
}


///////////////////////CSAD AUXILIAR FUNCTIONS//////////////////////////////
void csad_ini_pos_nei(
        const int w,
        const int h,
        const int ndt,
        const int rdt,
        PosNei *pnei
        ) {
    int size = w*h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < ndt; j++) {
            pnei[i].api[j] = -1; //Indicates that is out
            pnei[i].apj[j] = -1;
            pnei[i].b[j] = 0.0;
        }


    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++){
            int   it = 0;
            float ne = 0.0;
            const int pos = j*w + i;
            for (int k = -rdt; k < (rdt+1); k++)
                for (int l = -rdt; l < (rdt+1); l++){
                    if (!(k ==0 && k==l)){

                        int ap = validate_ap(w, h, i, j, l, k);
                        if (ap == 0)
                        {
                            pnei[pos].api[it] = i + l;
                            pnei[pos].apj[it] = j + k;
                            assert(pnei[pos].api[it]>=0 && pnei[pos].apj[it]>=0);
                            ne++;
                        }
                        it++;
                    }
                }
            pnei[pos].n = ne;
            pnei[pos].ba.resize(2*ne + 1);
        }
}



