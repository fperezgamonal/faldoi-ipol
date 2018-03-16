

////////////////////////////////////////////////////////////////////////////////
////////////////////////AUXILIAR FUNCTIONS NLTVL1///////////////////////////////
 bool positive(int val);
 float getsample_0(float *x, int w, int h, int pd, int i, int j, int l);

 int validate_ap(int w, int h, int i, int j, int di, int dj);

 int validate_ap_patch(
                const int ii,     // initial column
                const int ij,     // initial row
                const int ei,     // end column
                const int ej,     // end row
                int i, 
                int j
  );

 float get_wspatial(int l, int k);
 float get_wcolor(float *a, int w, int h, int i, int j, int l, int k,int pd);

 float get_weight(float *a, int w, int h, int i, int j, int l, int k,int pd);

 void nltv_ini_dual_variables(
                float *a,
                const int pd,
                const int w,
                const int h,
                const int n_d,
                const int radius,
                DualVariables *p,
                DualVariables *q
  );

 void non_local_divergence(
            DualVariables *p,
            const int ii, // initial column
            const int ij, // initial row
            const int ei, // end column
            const int ej, // end row
            const int w,
            int n_d,
            float *div_p
    );

 //////////////////////////CSAD///////////////////////////////////////////////
 void csad_ini_pos_nei(
                const int w,
                const int h,
                const int ndt,
                const int rdt,
                PosNei *pnei
  );

 float max(float a, float b);
 float min(float a, float b);
