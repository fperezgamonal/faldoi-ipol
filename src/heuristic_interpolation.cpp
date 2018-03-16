#ifndef HEURISTIC_LOCAL_FALDOI
#define HEURISTIC_LOCAL_FALDOI

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <queue>

extern "C" {
#include "smapa.h"
#include "iio.h"
}


SMART_PARAMETER(MIN_DISTANCE,10000)
SMART_PARAMETER(THE_LAMBDA_FACTOR,40)

//Struct

struct SparseOF{
    int i; // column
    int j; // row
    int u; //x- optical flow component
    int v; //y- optical flow component
    float sim_node; //similarity measure for the actual pixel
    float sim_accu; //similarity measure for the accumulated path.
};

class CompareSparseOF {
public:
    bool operator()(SparseOF& e1, SparseOF& e2){
        return e1.sim_node > e2.sim_node;
    }
};

typedef std::priority_queue<SparseOF, std::vector<SparseOF>, CompareSparseOF> heap;


static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return INFINITY;
    return x[(i+j*w)*pd + l];
}

static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return NAN;
    return x[(i+j*w)*pd + l];
}


//Census tranform of a single image
static void pack_bits_into_bytes(unsigned char *out, int *bits, int nbits) {
    int nbytes = ceil(nbits/8.0);
    for (int i = 0; i < nbytes; i++)
    {
        out[i] = 0;
        for (int j = 0; j < 8; j++)
            out[i] = out[i] * 2 + bits[8*i+j];
    }
}

static void census_at(unsigned char *out, float *x, int w, int h, int pd,
                      int winradius, int p, int q) {
    int side = 2*winradius + 1;
    int nbits = pd * (side * side - 1);
    int *bits = new int [nbits];
    int cx = 0;
    for (int l = 0; l < pd; l++)
        for (int j = -winradius; j <= winradius; j++)
            for (int i = -winradius; i <= winradius; i++)
            {
                float a = getsample_nan(x, w, h, pd, p    , q    , l);
                float b = getsample_nan(x, w, h, pd, p + i, q + j, l);
                if (i || j)
                    bits[cx++] = a > b;
            }
    assert(cx == nbits);

    pack_bits_into_bytes(out, bits, nbits);

    delete [] bits;
}

static void color_census_transform(unsigned char *y, int opd,
                                   float *x, int w, int h, int pd, int winradius)
{

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
            census_at(y + opd * (w * j + i), x, w, h, pd, winradius, i, j);
}


float eval_displacement_by_census(float *a, float *b,
                                  int w, int h, int pd,
                                  int wrad, int i, int j, float d[2])
{
    assert(wrad == 0);
    int count_bits[256] = {
    #               define B2(n) n,     n+1,     n+1,     n+2
    #               define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
    #               define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
        B6(0), B6(1), B6(1), B6(2)
    }, r = 0;
    for (int l = 0; l < pd; l++)
    {
        unsigned char p= getsample_inf(a, w, h, pd, i       , j       , l);
        unsigned char q= getsample_inf(b, w, h, pd, i + d[0], j + d[1], l);
        if (std::isfinite(p) && std::isfinite(q)){
            r += count_bits[p ^ q];
        }
    }
    //Normalize the result between [0,1]
    return (r*1.0)/(8*pd);
}



//PATCH_METHOD (CENSUS)
static inline float matching_cost(float *a, float *b, int w, int h, int pd,
                                  int i, int j, int u, int v, int w_radio)
{
    float d[2] = {(float)u, (float )v};
    int side = 2 * w_radio + 1;
    int nbits = pd * (side * side - 1);
    int nbytes = std::ceil(nbits / 8.0);
    return eval_displacement_by_census(a, b, w, h, nbytes, 0, i, j, d);
}

static float energy(float *a, float *b, int w, int h, int pd,
                    int i, int j, int u, int v, int w_radio,
                    int *dx, int *dy)
{ 

    int neighborhood[9][2] = {
        {0,0}, {0,1}, {0,-1},{1,0},{-1,0},
        {1,1}, {1,-1},{-1,1}, {-1,-1}};
    int best = 0;
    float r = INFINITY;
    //Enforce the disparity gradient limit (k = 5)
    //Not Enforce (k = 9)
    for (int k = 0; k < 5; k++)
    {
        int dx = neighborhood[k][0];
        int dy = neighborhood[k][1];
        float rn = matching_cost(a, b, w, h, pd, i, j,
                                 u + dx, v + dy,w_radio);
        if (k != 0) rn *= THE_LAMBDA_FACTOR();
        if (rn < r) {
            r = rn;
            best = k;
        }
    }
    *dx = neighborhood[best][0];
    *dy = neighborhood[best][1];
    //if (r<0)
    // std::printf("%f\n",r);
    assert(r >= 0);
    assert(r <= 1);
    return r;
}

static  inline void add_neighboors(float *a, float *b, int w, int h, int pd,
                                   int i, int j, int u, int v, int w_radio,
                                   float *sim_val, float *sim_val_accu,
                                   int *fixed_points, heap *queue)
{
    int neighborhood[8][2] = {
        {0,1}, {0,-1},{1,0},{-1,0},
        {1,1}, {1,-1},{-1,1}, {-1,-1}};
    for (int k = 0; k < 8; k++){
        int px = i + neighborhood[k][0];
        int py = j + neighborhood[k][1];
        //fprintf(stderr, "px:%d py:%d\n",px,py);
        if (px >= 0 && px < w && py >=0 && py < h){
            if (!fixed_points[py*w +px]){
                int dx,dy;
                float ener_N = energy(a,b,w,h,pd,px,py,u,v,w_radio,&dx,&dy);
                assert(!std::isnan(ener_N));
                if (ener_N < sim_val[py*w + px] && ener_N < MIN_DISTANCE()){
                    sim_val[py*w + px] = ener_N;
                    sim_val_accu[py*w + px] = ener_N + sim_val_accu[j*w + i];
                    SparseOF element;
                    element.i = px; // column
                    element.j = py; // row
                    element.u = u + dx;
                    element.v = v + dy;
                    element.sim_node = ener_N;
                    element.sim_accu = ener_N + sim_val_accu[j*w + i];
                    queue->push(element);
                }
            }
        }
    }
}



static inline int not_fixed(
        heap *queue,
        float *a,
        float *b,
        float *sim_val,
        float *sim_val_accu,
        int *fixed_points,
        int u,
        int v,
        int i,
        int j,
        int w,
        int h,
        int pd,
        int w_radio)
{
    int neighborhood[8][2] = {
        {0,1}, {0,-1},{1,0},{-1,0},
        {1,1}, {1,-1},{-1,1}, {-1,-1}};
    for (int k = 0; k < 8; k++){
        int px = i + neighborhood[k][0];
        int py = j + neighborhood[k][1];
        //fprintf(stderr, "px:%d py:%d\n",px,py);
        if (px >= 0 && px < w && py >=0 && py < h){
            if (!fixed_points[py*w +px]){
                int dx,dy;
                float ener_N = energy(a,b,w,h,pd,px,py,u,v,w_radio,&dx,&dy);
                assert(!std::isnan(ener_N));
                if (ener_N < sim_val[py*w + px]){
                    sim_val[py*w + px] = ener_N;
                    sim_val_accu[py*w + px] = ener_N + sim_val_accu[j*w + i];
                    SparseOF element;
                    element.i = px; // column
                    element.j = py; // row
                    element.u = u + dx;
                    element.v = v + dy;
                    element.sim_node = ener_N;
                    element.sim_accu = ener_N + sim_val_accu[j*w + i];
                    queue->push(element);
                }
            }
        }
    }

}

void heuristic_interpolation(
        float *in,
        float *a,
        float *b,
        int w,
        int h,
        int pd,
        int w_radio,
        float *out
        )
{
    int size = w*h;
    int *fixed_points   = new int[size]();
    float *sim_val      = new float[size]();
    float *sim_val_accu  = new float[size]();

    std::printf("Census Block Matching\n");
    int side = 2 * w_radio + 1;
    int nbits = pd * (side * side - 1);
    int nbytes = std::ceil(nbits / 8.0);
    // allocate output image
    int size_C = w * h * nbytes;
    unsigned char *y = (unsigned char *) malloc(size_C);
    unsigned char *z = (unsigned char *) malloc(size_C);

    float *aC = new float[size_C];
    float *bC = new float[size_C];

    color_census_transform(y, nbytes, a, w, h, pd, w_radio);
    color_census_transform(z, nbytes, b, w, h, pd, w_radio);

    for (int i = 0; i< size_C; i ++)
    {
        aC[i] = y[i];
        bC[i] = z[i];
    }

    for (int i = 0; i < size; i++){
        sim_val_accu[i] = INFINITY;
        sim_val[i] = INFINITY;
        out[i] = NAN;
        if (std::isfinite(in[i]) && std::isfinite(in[i + size]))
        {
            //1-fixed 0-not
            fixed_points[i] = 1;
            sim_val_accu[i] = 0.0;
            sim_val[i] = 0.0;
            out[i] = in[i];
            out[size + i] = in[size +  i];
        }
    }


    heap queue;
    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            //check if some of this neighboor is not fixed,
            if (fixed_points[j*w + i] == 1)
            {
                int u = round(in[j*w + i]);
                int v = round(in[size + j*w + i]);
                not_fixed(&queue, aC, bC, sim_val, sim_val_accu, fixed_points,
                          u, v, i, j, w, h, pd, w_radio);
            }
        }
    }

    //Extract elements from the queue until is empty
    int nfixed = 0;
    while (! queue.empty()) {
        SparseOF element = queue.top();
        int i = element.i;
        int j = element.j;
        int u = element.u;
        int v = element.v;
        int pos = j*w + i;
        queue.pop();
        if (!fixed_points[pos]) {
            nfixed++;
            fixed_points[pos] = 1;

            // if (0 == nfixed%1)
            // {
            //   fprintf(stderr, "fixed %d/%d (%g%%)\n",
            //       nfixed, size, nfixed*100.0/size);
            //   fprintf(stderr, "queue size now = %d\n", (int)queue.size());
            //   char buf[1000];
            //   sprintf(buf, "%s_of_%08d.flo",GLOBAL_TMP_FILE, nfixed);
            //   iio_save_image_float_split(buf, out, w, h, 2);
            //   sprintf(buf, "%s_sim_node_%08d.tiff",GLOBAL_TMP_FILE, nfixed);
            //   iio_save_image_float(buf, sim_val, w, h);
            //   sprintf(buf, "%s_sim_accu_%08d.tiff",GLOBAL_TMP_FILE, nfixed);
            //   iio_save_image_float(buf, sim_val_accu, w, h);
            // }
            out[pos] = (float) u;
            out[w*h + pos] = (float) v;
            add_neighboors(aC, bC, w, h, pd, i, j, u, v, w_radio,
                           sim_val, sim_val_accu, fixed_points, &queue);
        }
    }

    free(y);
    free(z);
    delete [] fixed_points;
    delete [] sim_val;
    delete [] sim_val_accu;
    delete [] aC;
    delete [] bC;
}
#endif


// int main(int argc, char *argv[])
// {
//   if (argc != 5) {
//     fprintf(stderr, "usage:\n\t"
//     "%s  in.flo I0 I1 out.flo\n", *argv);
//     //0 1     2     3  4        5        6        
//     return 1;
//   }

//   char *tpd = mktemp(GLOBAL_TMP_FILE); 
//   if (!tpd)
//     return fprintf(stderr, "ERROR(%s) cannot create temp file \"%s\"\n",
//         *argv, GLOBAL_TMP_FILE);
//   std::printf("%s\n",GLOBAL_TMP_FILE );

//   char *filename_in = argv[1];
//   char *filename_i0 = argv[2];
//   char *filename_i1 = argv[3];
//   char *filename_out = argv[4];

//   int w[2], h[2], pd;
//   float *in = iio_read_image_float_split(filename_in, w, h, &pd);
//   float *a = iio_read_image_float_vec(filename_i0, w+1, h+1, &pd);
//   float *b = iio_read_image_float_vec(filename_i1, w+1, h+1, &pd);
//   if (w[0] != w[1] || h[0] != h[1])
//     return fprintf(stderr, "image and edge file size mismatch"
//         " %d %d != %d %d", w[0], h[0], w[1], h[1]);
//   int size = w[0]*h[0];
//   float *out = new float[size*2]();

//   int w_radio = 1;
//   heuristic_interpolation(in, a, b, w[0], h[0], pd, w_radio, out);
//   iio_save_image_float_split(filename_out, out, *w, *h, 2);

//   return 0;
// }
