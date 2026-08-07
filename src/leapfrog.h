
#ifndef LEAPFROG_H
#define LEAPFROG_H

void initialization(float ***phi,float ***gf);

void Heat3d1(float ***c_old, float ***c_new, float ***gf);

void Heat3d2(float ***c_old, float ***cc_old, float ***c_new,float ***gf);

void vcycle1(float ***uf_new, float ***c_old, float ***su, float ***gf, int nxf, int nyf, int nzf, int ilevel);


void vcycle2(float ***uf_new, float ***c_old, float ***su, float ***gf, int nxf, int nyf, int nzf, int ilevel);

void relax1(float ***c, float ***oc, float ***sc,float ***gf, int ilevel, int nxt, int nyt, int nzt);


void relax2(float ***c, float ***oc, float ***sc,float ***gf, int ilevel, int nxt, int nyt, int nzt);

void restrict_grid(float ***uf, float ***uc, int nxt, int nyt, int nzt);

void prolong(float ***uc, float ***uf, int nxt, int nyt, int nzt);
float cube_max(float ***c, int nxt, int nyt, int nzt);
float error(float ***c_old, float ***c_new, int nxt, int nyt, int nzt);
void laplace_ch_gf(float ***a, float ***lap_a,float ***gf, int nxt, int nyt, int nzt);
float error3(float ***c_old, float ***c_new, int nxt, int nyt, int nzt);
void cahn1(float ***c_old, float ***c_new);
void print_data3( float ***c, int count);
void source(float ***c_old, float ***src_c);

#endif
