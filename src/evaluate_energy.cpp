#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>


extern "C" {
#include "iio.h"
#include "mask.h"
#include "bicubic_interpolation.h"
}

#define M_TVL1       0
#define M_TVL1_W     1
#define M_NLTVL1     2 
#define M_NLTVL1_W   3 
#define M_TVCSAD     4
#define M_TVCSAD_W   5
#define M_NLTVCSAD   6
#define M_NLTVCSAD_W 7


void rgb_to_gray(float *in, int w, int h, float *out)
{
  int size = w*h;
  for (int i = 0; i < size; i++)
  {
    out[i] = .299*in[i] + .587*in[size + i] + .114*in[2*size + i];

  }

}

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
  if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
    return 0;
  return x[w*h*l + j*w + i]; // r g b
  //return x[(i+j*w)*pd + l]; // intrelazado
}

static float get_wspatial(
                float *a, 
                int w, 
                int h, 
                int i, 
                int j, 
                int l, 
                int k,
                int pd
 )
{
  float ws = 5;
  float w_tmp;
  float difS = 0.0;
  difS = hypot(l,k);
  w_tmp = exp(- difS/ws);
  //std::printf("Dif_S: %f W_S: %f\n",difS, w_tmp);

  return w_tmp; 
}

static float get_wcolor(
                float *a, 
                int w, 
                int h, 
                int i, 
                int j, 
                int l, 
                int k,
                int pd
 )
{
  float wi = 5;
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

static float get_weight(
                float *a, 
                int w, 
                int h, 
                int i, 
                int j, 
                int l, 
                int k,
                int pd
 )
{
  float wc = get_wcolor(a, w, h, i, j, l, k, pd);
  float ws = get_wspatial(a, w, h, i, j, l, k, pd);

  return wc*ws;
}

static int validate_ap(int w, int h, int i, int j, int di, int dj)
{
  const int r = j + dj; //Row
  const int c = i + di; //Colum
  if ( c < 0 || c >= w || r < 0 || r >= h){
    return -1;
  }
  return 0;
}

void eval_nltvcsad(
    float *a,         //Imagen de la que se sacan los pesos.
    float *I0,           // source image
    float *I1,           // target image
    float *u,
    float *v,
    int w, 
    int h,
    int pd,
    float lambda,
    int w_rd_nltv,
    int w_csad,
    float *ener_N
    )
{

  float *u1 = ofD->u1;
  float *u2 = ofD->u2;

  float *a = I0;
  //Columns and Rows
  const int w = ofD->w;
  const int h = ofD->h;

  float *I1w = new float[w*h];
  float ener = 0.0;
  
  int radius;

  bicubic_interpolation_warp(I1,  u, v, I1w, w, h, true);
  //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
  float ener = 0.0;
  for (int j = 0; j < h; j++){
  for (int i = 0; i < w; i++){
    const int center = j*w + i;

    //Regularization term
    float g = 0.0;
    radius = w_rd_nltv;
    float wt = 0.0;
    for (int k = -radius; k < (radius+1); k++)
    for (int l = -radius; l < (radius+1); l++)
    {           
        //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
        if (!(k ==0 && k==l))
        {
          int ap = validate_ap(w, h, i, j, l, k);
          if (ap == 0)
          {
            float wp = sqrt(get_weight(a, w, h, i, j, l, k, pd));
            wt += wp;
            int pos_n = (j + k)*w + (i+ l);
            g += fabs(u[center] - v[pos_n])*wp + fabs(u[center] - v[pos_n])*wp;
          }
          it++;
        }
    }
    g /=wt;//Normalize the regularization term.

    //DATA Term
    float dt = 0.0;
    radius = w_rd_csad;
    for (int k = -radius; k < (radius+1); k++)
    for (int l = -radius; l < (radius+1); l++)
    {           
        //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
        if (!(k ==0 && k==l))
        {
          int ap = validate_ap(w, h, i, j, l, k);
          if (ap == 0)
          {
            const int pos_n = (j + k)*w + (i+ l);
            dt +=  fabs(I0[center] - I0[pos_n] - I1w[center] + I1w[pos_n]);
          }
          it++;
        }
    }
    dt *=lambda;

    assert(dt>=0);

    ener +=dt + g;
    if (!std::isfinite(dt))
        std::printf("Corrupt data\n");
    if (!std::isfinite(g))
        std::printf("Corrupt regularization\n");
  }
  }
  assert(ener>=0);
  ener /=(w*h*1.0);
  (*ener_N) = ener;
  delete [] I1w;
}


void eval_tvl2coupled(
    float *I0,           // source image
    float *I1,           // target image
    float *u,
    float *v,
    const int w,
    const int h,
    const float lambda,
    float *ener_N
    )
{


  //Optical flow derivatives
  float *ux  = new float[w*h];
  float *uy  = new float[w*h];
  float *vx  = new float[w*h];
  float *vy  = new float[w*h];

  float *I1w = new float[w*h];

  float ener = 0.0;

  forward_gradient(u, ux, uy, w, h);
  forward_gradient(v, vx, vy, w, h); 
  bicubic_interpolation_warp(I1,  u, v, I1w, w, h, true);
  //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
  int m  = 0;
  for (int l = 0; l < w; l++){
  for (int k = 0; k < w; k++){
    const int i = l*wº + k;
    float dt = lambda*fabs(I1w[i]-I0[i]);
    float g1  = ux[i]*ux[i];
    float g12 = uy[i]*uy[i];
    float g21 = vx[i]*vx[i];
    float g2  = vy[i]*vy[i];
    float g  = sqrt(g1 + g12 + g21 + g2);
    if (!std::isfinite(dt)){
        std::printf("Corrupt data\n");
    }
    if (!std::isfinite(g)){
        std::printf("Corrupt regularization\n");
    }
    ener += dt + g;
  }
  }
  ener /=(w*h*1.0);
  (*ener_N) = ener;
  assert(ener >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////MAIN/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/**
 *
 *  Function to read images using the iio library
 *  It always returns an allocated the image.
 *
 */ 
static float *read_image(const char *filename, int *w, int *h)
{
  float *f = iio_read_image_float(filename, w, h);
  if (!f)
    fprintf(stderr, "ERROR: could not read image from file "
        "\"%s\"\n", filename);
  return f;
}



#include <cstdio>
#include <string.h>
// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d)
{
  int argc = *c;
  char **argv = *v;
  int id = d ? 1 : 0;
  for (int i = 0; i < argc - id; i++)
    if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
    {
      char *r = argv[i+id]+1-id;
      *c -= id+1;
      for (int j = i; j < argc - id; j++)
        (*v)[j] = (*v)[j+id+1];
      return r;
    }
  return d;
}
int main(int c, char *v[])
{
  // process input 
  char *windows_radio_reg = pick_option(&c,&v, (char *)"wr_dt",(char *) "2");
  char *windows_radio_data = pick_option(&c,&v, (char *)"wr_rt",(char *) "2");
  char *var_reg       = pick_option(&c,&v,(char *)"m",(char *)"0");
  if (c!=5) {
    fprintf(stderr, "usage %d :\n\t%s i0 i1 in0.flo method"
//                          0      1  2  3   4      
                    "\n",c, *v);  
    fprintf(stderr, "M_TVL1       0\n M_NLTVCSAD   6\n");
    return 1;
  }

  char *filename_i0  = v[1];
  char *filename_i1  = v[2];
  char *filename_flo  = v[3];
  int val_method = atoi(v[4]);

  // open input images
  int w[4], h[4], pd[4];
  float *a = iio_read_image_float_split(filename_i0, w + 0, h + 0, pd + 0);
  float *b = iio_read_image_float_split(filename_i1, w + 1, h + 1, pd + 1);
  float *go = iio_read_image_float_split(filename_flo, w + 2, h + 2, pd + 2);

  if (w[0] != w[1] || h[0] != h[1] || w[0] != w[2] || h[0] != h[2] || pd[2]!=2)
    return fprintf(stderr, "ERROR: input images size mismatch\n");

  float i0 = new float[w[0]*h[0]];
  float i1 = new float[w[0]*h[0]];
  float *u = go;
  float *v = go + w[0]*h[0];

  if (pd[0]==1){
    for (int i = 0; i < w[0]*h[0]; i++)
    {
      i0[i] = a[i];
      i1[i] = b[i];
    }
  }else{
    rgb_to_gray(a, w[0], h[0], i0);
    rgb_to_gray(b, w[0], h[0], i1);
  }




  float lambda;
  int w_rd_nltv = atoi(windows_radio_reg);
  int w_csad = atoi(windows_radio_data);
  float ener;


  switch(val_method)
  { 
    case M_NLTVCSAD: //NLTV-CSAD
      lambda = 3.47;
      fprintf(stderr,"NLTV-CSAD\n");
      eval_nltvcsad(a, i0, i1, u, v, w[0], h[0], pd[0],
                            lambda, w_rd_nltv, w_csad, &ener);
      break;
    default: //TV-l2 coupled
      fprintf(stderr,"TV-l2 coupled\n");
      lambda = 40;
      eval_tvl2coupled(i0, i1, u, v, w[0], h[0], lambda, &ener);
  }
  


  // cleanup and exit
  free(a);
  free(b);
  free(go);

  delete [] a;
  delete [] b;
  return 0;
}