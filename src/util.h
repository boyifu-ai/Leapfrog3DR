#ifndef gnx
#define gnx 128
#endif

#ifndef gny
#define gny 128
#endif

#ifndef gnz
#define gnz 128
#endif

#define iloop for (i = 1; i <= gnx; i++)
#define jloop for (j = 1; j <= gny; j++)
#define kloop for (k = 1; k <= gnz; k++)
#define ijkloop iloop jloop kloop
#define kjiloop kloop jloop iloop

#define iloopp for (i = 0; i <= gnx; i++)
#define jloopp for (j = 0; j <= gny; j++)
#define kloopp for (k = 0; k <= gnz; k++)
#define ijkloopp iloopp jloopp kloopp
#define kjiloopp kloopp jloopp iloopp

#define iloopt for (i = 1; i <= nxt; i++)
#define jloopt for (j = 1; j <= nyt; j++)
#define kloopt for (k = 1; k <= nzt; k++)
#define ijkloopt iloopt jloopt kloopt
#define kjiloopt kloopt jloopt iloopt

float ***cube(int xl, int xr, int yl, int yr, int zl, int zr);

void free_cube(float ***t, int xl, int xr, int yl, int yr, int zl, int zr);

float *dvector(long nl, long nh);

float **dmatrix(long nrl, long nrh, long ncl, long nch);

void free_dvector(float *v, long nl, long nh);

void free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch);

float sum_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr);

void zero_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr);

void zero_cube2(float ***a, float ***b,
                int xl, int xr, int yl, int yr, int zl, int zr);

void print_data(float ***phi);

void cube_add(float ***a, float ***b, float ***c,
              int xl, int xr, int yl, int yr, int zl, int zr);

void cube_add2(float ***a, float ***b, float ***c,
               float ***a2, float ***b2, float ***c2,
               int xl, int xr, int yl, int yr, int zl, int zr);

void cube_sub(float ***a, float ***b, float ***c,
              int xl, int xr, int yl, int yr, int zl, int zr);

void cube_sub2(float ***a, float ***b, float ***c,
               float ***a2, float ***b2, float ***c2,
               int xl, int xr, int yl, int yr, int zl, int zr);

void cube_copy(float ***a, float ***b,
               int xl, int xr, int yl, int yr, int zl, int zr);

void cube_copy2(float ***a, float ***b,
                float ***a2, float ***b2,
                int xl, int xr, int yl, int yr, int zl, int zr);

void augmenc(float ***c);

void augmen_phi(float ***c, int nxt, int nyt, int nzt);
