#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include "util.h"
#include "leapfrog.h"
#define NR_END 1
#define MAX_PATH_LEN 1024

static const char *nt_points_path = "examples/teapot/nt_points.m";
static const char *point_data_path = "examples/teapot/fun_data.m";
static const char *output_dir = "output";
static int write_phi_output = 1;

static void parse_args(int argc, char **argv);
static void print_usage(const char *prog);
static FILE *open_file_or_die(const char *path, const char *mode);
static void make_directory_or_die(const char *path);
static void join_path(char *dst, size_t dst_size, const char *dir, const char *name);
static int max_int(int a, int b);
static int min_int(int a, int b);

int nx, ny, nz, it, count, c_relax;
float gam, lam, dt, h, h2, xleft, xright, yleft, yright, zleft, zright, pi, Innmu,
    InnR, ***muu, oE, cE, mE, oIE, cahn, 
    ***ct, ***sc, ***intc, ***ooc, Sb, ***sorc, lamm, xxi, MM, C;
void print_vol(float ***aa, int yu, int yd, int xl, int xr, int zl, int zr, int num);

int main(int argc, char **argv)
{
   int i, j, k, max_it, ns;
   float ***oc, ***nc, ***gf;
   float elapsed, start, T, end;

   FILE *my, *myOE;
   char output_path[MAX_PATH_LEN];

   parse_args(argc, argv);
   make_directory_or_die(output_dir);

   start = clock();

   pi = 4.0 * atan(1.0);

   nx = gnx;
   ny = gny;
   nz = gnz;

   xleft = 0.0, xright = 1.0;
   yleft = 0.0, yright = xright * ny / (1.0 * nx);
   zleft = 0.0, zright = xright * nz / (1.0 * nx);
   c_relax = 5;
   count = 1;

   h = (xright - xleft) / (1.0 * nx);
   h2 = h * h;
   gam = 6.0 * h / (2.0 * sqrt(2) * atanh(0.9)); //6.0
   cahn = pow(gam, 2);
   dt = 2 * 1.25e-5;// 2
   T = 16 * dt;// 16
   max_it = (int)(T / dt + 0.5);
   ns = (int)(max_it / 16 + 0.5);// 16

   Sb = 2.0 / cahn;//2.0 / cahn
   lamm = 0.0;

   C = 0.0;

   printf("nx = %d , ny = %d , nz = %d\n", nx, ny, nz);
   printf("dt      = %f\n", dt);
   printf("max_it  = %d\n", max_it);
   printf("ns      = %d\n", ns);
   printf("gam     = %f\n", gam);
   printf("h       = %f\n", h);
   join_path(output_path, sizeof(output_path), output_dir, "remarks.m");
   my = open_file_or_die(output_path, "w");
   fprintf(my, "\n\n");
   fprintf(my, "nx      = %d\n", nx);
   fprintf(my, "ny      = %d\n", ny);
   fprintf(my, "nz      = %d\n", nz);
   fprintf(my, "dt      = %f\n", dt);
   fprintf(my, "max_it  = %d\n", max_it);
   fprintf(my, "ns      = %d\n", ns);
   fprintf(my, "gam     = %f\n", gam);
   fprintf(my, "h       = %f\n", h);
   fclose(my);

   muu = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
   sorc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
   ooc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);

   gf = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
   oc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
   nc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
   ct = cube(1, nx, 1, ny, 1, nz);
   sc = cube(1, nx, 1, ny, 1, nz);

   initialization(oc, gf);

   cube_copy(nc, oc, 1, nx, 1, ny, 1, nz);
   cube_copy(ooc, oc, 1, nx, 1, ny, 1, nz);

   augmen_phi(oc, nx, ny, nz);

   MM = 0.0;
   oE = 0.0;
   ijkloop
   {
      MM = MM + 0.25 * pow(pow(oc[i][j][k], 2) - 1.0, 2) / (gam * gam) - 0.25 * lamm * (0.2 * pow(oc[i][j][k], 2) + oc[i][j][k] - (2.0 / 3) * pow(oc[i][j][k], 3));
      oE = oE + 0.5 * pow(oc[i + 1][j][k] - oc[i][j][k], 2) / (h * h) + 0.5 * pow(oc[i][j + 1][k] - oc[i][j][k], 2) / (h * h) + 0.5 * pow(oc[i][j][k + 1] - oc[i][j][k], 2) / (h * h) + 0.25 * Sb * pow(oc[i][j][k] - oc[i][j][k], 2);
   }

   MM = MM * h * h * h + C;

   oE = oE * h * h * h + MM - C;
   cE = oE;
   mE = oE;

   print_data3(oc, count);
   count++;
   augmen_phi(gf, nx, ny, nz);

   join_path(output_path, sizeof(output_path), output_dir, "Oene.m");
   myOE = open_file_or_die(output_path, "w");


   fprintf(myOE, "%16.12f \n", oE);
 

   for (it = 1; it <= max_it; it++)
   {
      
/*
      if (it == 1)
      {
         Heat3d1(oc, nc, gf);
        

         augmen_phi(nc, nx, ny, nz);
         augmen_phi(oc, nx, ny, nz);


         oE = 0.0;

         ijkloop{

            oE = oE + 0.25*( pow(nc[i+1][j][k]-nc[i][j][k],2)/h2 + pow(nc[i][j+1][k]-nc[i][j][k],2)/h2 + pow(nc[i][j][k+1]-nc[i][j][k],2)/h2
            + pow(oc[i+1][j][k]-oc[i][j][k],2)/h2 + pow(oc[i][j+1][k]-oc[i][j][k],2)/h2 + pow(oc[i][j][k+1]-oc[i][j][k],2)/h2 )
            + (0.25/cahn)*pow(nc[i][j][k]*oc[i][j][k]-1.0,2) + 0.5*Sb*pow(nc[i][j][k]-oc[i][j][k],2);
         }

         oE = oE*h*h*h;

        
           cube_copy(oc, nc, 1, nx, 1, ny, 1, nz);
      }
      else
      {*/
         Heat3d2(oc, ooc, nc, gf);
         

         augmen_phi(nc, nx, ny, nz);
         augmen_phi(oc, nx, ny, nz);


         oE = 0.0;

         ijkloop{

            oE = oE + 0.25*( pow(nc[i+1][j][k]-nc[i][j][k],2)/h2 + pow(nc[i][j+1][k]-nc[i][j][k],2)/h2 + pow(nc[i][j][k+1]-nc[i][j][k],2)/h2
            + pow(oc[i+1][j][k]-oc[i][j][k],2)/h2 + pow(oc[i][j+1][k]-oc[i][j][k],2)/h2 + pow(oc[i][j][k+1]-oc[i][j][k],2)/h2 )
            + (0.25/cahn)*pow(nc[i][j][k]*oc[i][j][k]-1.0,2) + 0.5*Sb*pow(nc[i][j][k]-oc[i][j][k],2);
         }

         oE = oE*h*h*h;

        
         cube_copy(ooc, oc, 1, nx, 1, ny, 1, nz);
         cube_copy(oc, nc, 1, nx, 1, ny, 1, nz);
      //}
     

      printf("iteration=   %d \n", it);

      if ((it % ns == 0))
      {
         print_data3(oc, count);
         count++;

         fprintf(myOE, "%16.12f \n", oE);
        

         printf("\n counts=   %d \n", count);
         end = clock();
         elapsed = ((float)(end - start)) / CLOCKS_PER_SEC;
         printf("Time elapsed %f\n", elapsed);
      }
   }

   end = clock();
   elapsed = ((float)(end - start)) / CLOCKS_PER_SEC;
   fclose(myOE);

   return 0;
}

static void parse_args(int argc, char **argv)
{
   int argi;

   for (argi = 1; argi < argc; argi++)
   {
      if (strcmp(argv[argi], "--nt-points") == 0 && argi + 1 < argc)
      {
         nt_points_path = argv[++argi];
      }
      else if (strcmp(argv[argi], "--points") == 0 && argi + 1 < argc)
      {
         point_data_path = argv[++argi];
      }
      else if (strcmp(argv[argi], "--output-dir") == 0 && argi + 1 < argc)
      {
         output_dir = argv[++argi];
      }
      else if (strcmp(argv[argi], "--no-phi") == 0)
      {
         write_phi_output = 0;
      }
      else if (strcmp(argv[argi], "--help") == 0 || strcmp(argv[argi], "-h") == 0)
      {
         print_usage(argv[0]);
         exit(0);
      }
      else
      {
         fprintf(stderr, "Unknown or incomplete option: %s\n", argv[argi]);
         print_usage(argv[0]);
         exit(2);
      }
   }
}

static void print_usage(const char *prog)
{
   fprintf(stderr,
           "Usage: %s --nt-points PATH --points PATH --output-dir DIR [--no-phi]\n"
           "\n"
           "Inputs:\n"
           "  --nt-points PATH   text file containing the number of point-cloud samples\n"
           "  --points PATH      text file with x y z point-cloud coordinates, one value per line\n"
           "  --output-dir DIR   directory for remarks.m, Oene.m, and phi*.m outputs\n"
           "  --no-phi           write energy only; skip large phi*.m output files\n",
           prog);
}

static FILE *open_file_or_die(const char *path, const char *mode)
{
   FILE *file = fopen(path, mode);
   if (file == NULL)
   {
      fprintf(stderr, "Error: cannot open %s with mode %s\n", path, mode);
      exit(1);
   }
   return file;
}

static void make_directory_or_die(const char *path)
{
#ifdef _WIN32
   int rc = _mkdir(path);
#else
   int rc = mkdir(path, 0777);
#endif
   if (rc != 0 && errno != EEXIST)
   {
      fprintf(stderr, "Error: cannot create output directory %s\n", path);
      exit(1);
   }
}

static void join_path(char *dst, size_t dst_size, const char *dir, const char *name)
{
   size_t len = strlen(dir);
   const char *sep = (len > 0 && (dir[len - 1] == '/' || dir[len - 1] == '\\')) ? "" : "/";
   int written = snprintf(dst, dst_size, "%s%s%s", dir, sep, name);
   if (written < 0 || (size_t)written >= dst_size)
   {
      fprintf(stderr, "Error: path is too long: %s/%s\n", dir, name);
      exit(1);
   }
}

static int max_int(int a, int b)
{
   return a > b ? a : b;
}

static int min_int(int a, int b)
{
   return a < b ? a : b;
}

void initialization(float ***phi, float ***gf)
{
   extern float h;
   float gam2;
   double valuexyz;
   FILE *fp;
   int i, j, k, ik, npo, ia, ja, ka, msize;
   int i_start, i_end, j_start, j_end, k_start, k_end;
   float x, y, z, d;
   gam2 = 3.0 * h / (2.0 * sqrt(2) * atanh(0.9));

   fp = open_file_or_die(nt_points_path, "r");
   if (fscanf(fp, "%lf", &valuexyz) != 1)
   {
      fprintf(stderr, "Error: cannot read point count from %s\n", nt_points_path);
      exit(1);
   }
   npo = (int)(valuexyz + 0.01);
   fclose(fp);
   fp = open_file_or_die(point_data_path, "r");

   ijkloop
   {
      phi[i][j][k] = 1.0;
   }

   msize = (int)(gam2 / h * 6 + 0.1);
   for (ik = 1; ik <= npo; ik++)
   {

      if (fscanf(fp, "%lf", &valuexyz) != 1)
      {
         fprintf(stderr, "Error: cannot read x coordinate for point %d from %s\n", ik, point_data_path);
         exit(1);
      }
      x = valuexyz;
      if (fscanf(fp, "%lf", &valuexyz) != 1)
      {
         fprintf(stderr, "Error: cannot read y coordinate for point %d from %s\n", ik, point_data_path);
         exit(1);
      }
      y = valuexyz;
      if (fscanf(fp, "%lf", &valuexyz) != 1)
      {
         fprintf(stderr, "Error: cannot read z coordinate for point %d from %s\n", ik, point_data_path);
         exit(1);
      }
      z = valuexyz;
      ia = (int)(x / h + 0.499999) - 1;
      ja = (int)(y / h + 0.499999) - 1;
      ka = (int)(z / h + 0.499999) - 1;
      i_start = max_int(1, ia - msize + 1);
      i_end = min_int(nx, ia + msize);
      j_start = max_int(1, ja - msize + 1);
      j_end = min_int(ny, ja + msize);
      k_start = max_int(1, ka - msize + 1);
      k_end = min_int(nz, ka + msize);

      for (i = i_start; i <= i_end; i++)
      {

         for (j = j_start; j <= j_end; j++)
         {

            for (k = k_start; k <= k_end; k++)
            {

               d = sqrt(pow((i - 0.5) * h - x, 2) + pow((j - 0.5) * h - y, 2) + pow((k - 0.5) * h - z, 2));
               if (d < phi[i][j][k])
               {
                  phi[i][j][k] = d;
               }
            }
         }
      }
   }

   fclose(fp);
   ijkloop
   {
      phi[i][j][k] = -tanh((phi[i][j][k] - sqrt(2) * gam2 - 2.0 * h) / (sqrt(2) * gam2));
      gf[i][j][k] = 1.0 - phi[i][j][k] * phi[i][j][k];
   }
}

void Heat3d1(float ***c_old, float ***c_new, float ***gf)
{
   extern int nx, ny, nz;
   extern float gam, lamm, h, Sb, dt, cahn, ***sorc;

   int i, j, k, max_it_CH = 100, it_mg = 1;
   float tol = 1.0e-6, resid = 1.0;

   float ***sor1;
   sor1 = cube(1, nx, 1, ny, 1, nz);

   cube_copy(sor1, c_new, 1, nx, 1, ny, 1, nz);

   augmen_phi(c_old, nx, ny, nz);

   ijkloop
   {
      sorc[i][j][k] = c_old[i][j][k] / dt + gf[i][j][k] * c_old[i][j][k] / cahn + gf[i][j][k] * Sb * c_old[i][j][k];
   }

   while (it_mg <= max_it_CH && resid > tol)
   {
      vcycle1(c_new, c_old, sorc, gf, nx, ny, nz, 1);
      resid = error3(sor1, c_new, nx, ny, nz);
      cube_copy(sor1, c_new, 1, nx, 1, ny, 1, nz);
      it_mg++;
   }

   printf("Vel: Heat equation iteration = %d  residual = %10.7f \n", it_mg, resid);

   free_cube(sor1, 1, nx, 1, ny, 1, nz);

   return;
}

void Heat3d2(float ***c_old, float ***cc_old, float ***c_new, float ***gf)
{
   extern int nx, ny, nz;
   extern float gam, lamm, h, Sb, dt, cahn, ***sorc;

   int i, j, k, max_it_CH = 100, it_mg = 1;
   float tol = 1.0e-6, resid = 1.0;

   float ***sor1;
   sor1 = cube(1, nx, 1, ny, 1, nz);

   cube_copy(sor1, c_new, 1, nx, 1, ny, 1, nz);

   augmen_phi(cc_old, nx, ny, nz);

   ijkloop
   {
      sorc[i][j][k] = cc_old[i][j][k] / (2.0 * dt) + gf[i][j][k] * c_old[i][j][k] / cahn - gf[i][j][k] * cc_old[i][j][k]*pow(c_old[i][j][k], 2) / (2.0 * cahn) - gf[i][j][k] * Sb * (-2.0 * c_old[i][j][k] + cc_old[i][j][k]) + 0.5 * gf[i][j][k] * (cc_old[i + 1][j][k] + cc_old[i - 1][j][k] + cc_old[i][j + 1][k] + cc_old[i][j - 1][k] + cc_old[i][j][k + 1] + cc_old[i][j][k - 1] - 6.0 * cc_old[i][j][k]) / h2;
   }

   while (it_mg <= max_it_CH && resid > tol)
   {
      vcycle2(c_new, c_old, sorc, gf, nx, ny, nz, 1);
      resid = error3(sor1, c_new, nx, ny, nz);
      cube_copy(sor1, c_new, 1, nx, 1, ny, 1, nz);
      it_mg++;
   }

   printf("Vel: Heat equation iteration = %d  residual = %10.7f \n", it_mg, resid);

   free_cube(sor1, 1, nx, 1, ny, 1, nz);

   return;
}

void vcycle1(float ***uf_new, float ***c_old, float ***su, float ***gf, int nxf, int nyf, int nzf, int ilevel)
{
   relax1(uf_new, c_old, su, gf, ilevel, nxf, nyf, nzf);
}

void vcycle2(float ***uf_new, float ***c_old, float ***su, float ***gf, int nxf, int nyf, int nzf, int ilevel)
{
   relax2(uf_new, c_old, su, gf, ilevel, nxf, nyf, nzf);
}

void relax1(float ***c, float ***oc, float ***sc, float ***gf, int ilevel, int nxt, int nyt, int nzt)
{

   extern float dt, Sb, xleft, xright, cahn;
   extern int c_relax;
   int i, j, k, ik;
   float ht2, x_fac, y_fac, z_fac, a[1], f[1], xs, xl, ys, yl, zs, zl;

   (void)ilevel;
   ht2 = pow((xright - xleft) / (float)nxt, 2);

   xs = 1.0;
   xl = 1.0;
   ys = 1.0;
   yl = 1.0;
   zs = 1.0;
   zl = 1.0;
   for (ik = 0; ik < c_relax; ik++)
      ijkloopt
      {

         f[0] = 0.0;

         if (i == 1)
         {
            x_fac = xl;
            f[0] += c[i + 1][j][k] * xl;
         }
         else if (i == nxt)
         {
            x_fac = xs;
            f[0] += c[i - 1][j][k] * xs;
         }
         else
         {
            x_fac = xs + xl;
            f[0] += c[i - 1][j][k] * xs + c[i + 1][j][k] * xl;
         }

         if (j == 1)
         {
            y_fac = yl;
            f[0] += c[i][j + 1][k] * yl;
         }
         else if (j == nyt)
         {
            y_fac = ys;
            f[0] += c[i][j - 1][k] * ys;
         }
         else
         {
            y_fac = ys + yl;
            f[0] += c[i][j - 1][k] * ys + c[i][j + 1][k] * yl;
         }

         if (k == 1)
         {
            z_fac = zl;
            f[0] += c[i][j][k + 1] * zl;
         }
         else if (k == nzt)
         {
            z_fac = zs;
            f[0] += c[i][j][k - 1] * zs;
         }
         else
         {
            z_fac = zs + zl;
            f[0] += c[i][j][k - 1] * zs + c[i][j][k + 1] * zl;
         }

         a[0] = 1.0 / dt + pow(oc[i][j][k], 2) / cahn + Sb * gf[i][j][k] + (x_fac + y_fac + z_fac) / ht2 * gf[i][j][k];
         f[0] = f[0] / ht2 * gf[i][j][k] + sc[i][j][k];

         c[i][j][k] = f[0] / a[0];
      }
}

void relax2(float ***c, float ***oc, float ***sc, float ***gf, int ilevel, int nxt, int nyt, int nzt)
{

   extern float dt, Sb, cahn, xleft, xright;
   extern int c_relax;
   int i, j, k, ik;
   float ht2, x_fac, y_fac, z_fac, a[1], f[1], xs, xl, ys, yl, zs, zl;

   (void)ilevel;
   ht2 = pow((xright - xleft) / (float)nxt, 2);

   xs = 1.0;
   xl = 1.0;
   ys = 1.0;
   yl = 1.0;
   zs = 1.0;
   zl = 1.0;
   for (ik = 0; ik < c_relax; ik++)
      ijkloopt
      {

         f[0] = 0.0;

         if (i == 1)
         {
            x_fac = xl;
            f[0] += c[i + 1][j][k] * xl;
         }
         else if (i == nxt)
         {
            x_fac = xs;
            f[0] += c[i - 1][j][k] * xs;
         }
         else
         {
            x_fac = xs + xl;
            f[0] += c[i - 1][j][k] * xs + c[i + 1][j][k] * xl;
         }

         if (j == 1)
         {
            y_fac = yl;
            f[0] += c[i][j + 1][k] * yl;
         }
         else if (j == nyt)
         {
            y_fac = ys;
            f[0] += c[i][j - 1][k] * ys;
         }
         else
         {
            y_fac = ys + yl;
            f[0] += c[i][j - 1][k] * ys + c[i][j + 1][k] * yl;
         }

         if (k == 1)
         {
            z_fac = zl;
            f[0] += c[i][j][k + 1] * zl;
         }
         else if (k == nzt)
         {
            z_fac = zs;
            f[0] += c[i][j][k - 1] * zs;
         }
         else
         {
            z_fac = zs + zl;
            f[0] += c[i][j][k - 1] * zs + c[i][j][k + 1] * zl;
         }

         a[0] = 1.0 / (2.0 * dt) + gf[i][j][k] * pow(oc[i][j][k], 2) / (2.0 * cahn) + Sb * gf[i][j][k] + 0.5 * (x_fac + y_fac + z_fac) / ht2 * gf[i][j][k];
         f[0] = 0.5 * f[0] / ht2 * gf[i][j][k] + sc[i][j][k];

         c[i][j][k] = f[0] / a[0];
      }
}

void restrict_grid(float ***uf, float ***uc, int nxt, int nyt, int nzt)
{
   int i, j, k;

   ijkloopt
   {

      uc[i][j][k] = (uf[2 * i][2 * j][2 * k] + uf[2 * i - 1][2 * j][2 * k] + uf[2 * i][2 * j - 1][2 * k] + uf[2 * i][2 * j][2 * k - 1] + uf[2 * i - 1][2 * j - 1][2 * k] + uf[2 * i - 1][2 * j][2 * k - 1] + uf[2 * i][2 * j - 1][2 * k - 1] + uf[2 * i - 1][2 * j - 1][2 * k - 1]) / 8.0;
   }
}

void prolong(float ***uc, float ***uf, int nxt, int nyt, int nzt)
{
   int i, j, k;

   ijkloopt
   {

      uf[2 * i][2 * j][2 * k] = uc[i][j][k];
      uf[2 * i - 1][2 * j][2 * k] = uc[i][j][k];
      uf[2 * i][2 * j - 1][2 * k] = uc[i][j][k];
      uf[2 * i][2 * j][2 * k - 1] = uc[i][j][k];

      uf[2 * i - 1][2 * j - 1][2 * k] = uc[i][j][k];
      uf[2 * i - 1][2 * j][2 * k - 1] = uc[i][j][k];
      uf[2 * i][2 * j - 1][2 * k - 1] = uc[i][j][k];
      uf[2 * i - 1][2 * j - 1][2 * k - 1] = uc[i][j][k];
   }
}

float cube_max(float ***c, int nxt, int nyt, int nzt)
{
   int i, j, k;
   float value = 0.0;

   ijkloopt
   {

      if (fabs(c[i][j][k]) > value)
         value = fabs(c[i][j][k]);
   }

   return value;
}

void laplace_ch_gf(float ***a, float ***lap_a, float ***gf, int nxt, int nyt, int nzt)
{
   extern float xright;

   int i, j, k;
   float dxt2, dadx_L, dadx_R, dady_B, dady_F, dadz_D, dadz_U;

   dxt2 = pow(xright / (float)nxt, 2);

   ijkloopt
   {

      if (i > 1)
         dadx_L = a[i][j][k] - a[i - 1][j][k];
      else
         dadx_L = 0.0;

      if (i < nxt)
         dadx_R = a[i + 1][j][k] - a[i][j][k];
      else
         dadx_R = 0.0;

      if (j > 1)
         dady_B = a[i][j][k] - a[i][j - 1][k];
      else
         dady_B = 0.0;

      if (j < nyt)
         dady_F = a[i][j + 1][k] - a[i][j][k];
      else
         dady_F = 0.0;

      if (k > 1)
         dadz_D = a[i][j][k] - a[i][j][k - 1];
      else
         dadz_D = 0.0;

      if (k < nzt)
         dadz_U = a[i][j][k + 1] - a[i][j][k];
      else
         dadz_U = 0.0;
      /*

            lap_a[i][j][k] = 0.5*((gf[i+1][j][k]+gf[i][j][k])*dadx_R
                      -(gf[i-1][j][k]+gf[i][j][k])*dadx_L + (gf[i][j+1][k]+gf[i][j][k])*dady_F
                      - (gf[i][j-1][k]+gf[i][j][k])*dady_B + (gf[i][j][k+1]+gf[i][j][k])*dadz_U
                      - (gf[i][j][k-1]+gf[i][j][k])*dadz_D) / dxt2;
      */
      lap_a[i][j][k] = gf[i][j][k] * (dadx_R - dadx_L + dady_F - dady_B + dadz_U - dadz_D) / dxt2;
   }
}

float error3(float ***c_old, float ***c_new, int nxt, int nyt, int nzt)
{
   int i, j, k;
   float value;

   value = 0.0;

   ijkloopt
   {
      value += pow(c_old[i][j][k] - c_new[i][j][k], 2);
   }

   value = sqrt(value / (1.0 * nxt * nyt * nzt));

   return value;
}

void print_data3(float ***c, int count)
{
   char filename[64];
   char buffer[MAX_PATH_LEN];
   FILE *fc;
   int i, j, k;

   if (!write_phi_output)
   {
      return;
   }

   snprintf(filename, sizeof(filename), "phi%d.m", count);
   join_path(buffer, sizeof(buffer), output_dir, filename);
   fc = open_file_or_die(buffer, "w");
   ijkloop
       fprintf(fc, " %10.7f\n", c[i][j][k]);
   fclose(fc);

   return;
}

float error(float ***c_old, float ***c_new, int nxt, int nyt, int nzt)
{
   int i, j, k;
   float value;

   value = 0.0;

   ijkloopt
   {
      value += pow(c_old[i][j][k] - c_new[i][j][k], 2);
   }

   value = sqrt(value / (1.0 * nxt * nyt * nzt));

   return value;
}

/************* util *****************/

float ***cube(int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, nrow = xr - xl + 1, ncol = yr - yl + 1, ndep = zr - zl + 1;
   float ***t;

   t = (float ***)malloc(((nrow + 1) * sizeof(float **)));
   t += 1;
   t -= xl;

   t[xl] = (float **)malloc(((nrow * ncol + 1) * sizeof(float *)));
   t[xl] += 1;
   t[xl] -= yl;

   t[xl][yl] = (float *)malloc(((nrow * ncol * ndep + 1) * sizeof(float)));
   t[xl][yl] += 1;
   t[xl][yl] -= zl;

   for (j = yl + 1; j <= yr; j++)
      t[xl][j] = t[xl][j - 1] + ndep;
   for (i = xl + 1; i <= xr; i++)
   {
      t[i] = t[i - 1] + ncol;
      t[i][yl] = t[i - 1][yl] + ncol * ndep;
      for (j = yl + 1; j <= yr; j++)
         t[i][j] = t[i][j - 1] + ndep;
   }

   return t;
}

void free_cube(float ***t, int xl, int xr, int yl, int yr, int zl, int zr)
{
   (void)xr;
   (void)yr;
   (void)zr;
   free((char *)(t[xl][yl] + zl - 1));
   free((char *)(t[xl] + yl - 1));
   free((char *)(t + xl - 1));
}

float *dvector(long nl, long nh)
/* allocate a float  vector with subscript range v[nl..nh] */
{
   float *v;

   v = (float *)malloc((nh - nl + 1 + NR_END) * sizeof(float));
   return v - nl + NR_END;
}

float **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   float **m;
   long i, nrow = nrh - nrl + 1 + NR_END, ncol = nch - ncl + 1 + NR_END;

   m = (float **)malloc((nrow) * sizeof(float *));
   m += NR_END;
   m -= nrl;

   m[nrl] = (float *)malloc((nrow * ncol) * sizeof(float));
   m[nrl] += NR_END;

   m[nrl] -= ncl;

   for (i = nrl + 1; i <= nrh; i++)
      m[i] = m[i - 1] + ncol;

   return m;
}

void free_dvector(float *v, long nl, long nh)
{
   (void)nh;
   free(v + nl - NR_END);

   return;
}

void free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch)
{
   (void)nrh;
   (void)nch;
   free(m[nrl] + ncl - NR_END);
   free(m + nrl - NR_END);

   return;
}

float sum_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;
   float sum = 0.0;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
            sum += a[i][j][k];

   return sum;
}

void zero_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = 0.0;
         }

   return;
}

void zero_cube2(float ***a, float ***b,
                int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k] = 0.0;
         }

   return;
}

void cube_add(float ***a, float ***b, float ***c,
              int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k] + c[i][j][k];
         }

   return;
}

void cube_add2(float ***a, float ***b, float ***c,
               float ***a2, float ***b2, float ***c2,
               int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k] + c[i][j][k];
            a2[i][j][k] = b2[i][j][k] + c2[i][j][k];
         }

   return;
}

void cube_sub(float ***a, float ***b, float ***c,
              int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k] - c[i][j][k];
         }

   return;
}

void cube_sub2(float ***a, float ***b, float ***c,
               float ***a2, float ***b2, float ***c2,
               int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k] - c[i][j][k];
            a2[i][j][k] = b2[i][j][k] - c2[i][j][k];
         }

   return;
}

void cube_copy(float ***a, float ***b,
               int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k];
         }

   return;
}

void cube_copy2(float ***a, float ***b,
                float ***a2, float ***b2,
                int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i = xl; i <= xr; i++)
      for (j = yl; j <= yr; j++)
         for (k = zl; k <= zr; k++)
         {

            a[i][j][k] = b[i][j][k];
            a2[i][j][k] = b2[i][j][k];
         }

   return;
}

void augmen_phi(float ***c, int nxt, int nyt, int nzt)
{

   int i, j, k;

   for (j = 1; j <= nyt; j++)
      for (k = 1; k <= nzt; k++)
      {

         c[0][j][k] = c[1][j][k];
         c[nxt + 1][j][k] = c[nxt][j][k];
      }

   for (i = 0; i <= nxt + 1; i++)
      for (k = 1; k <= nzt; k++)
      {

         c[i][0][k] = c[i][1][k];
         c[i][nyt + 1][k] = c[i][nyt][k];
      }

   for (i = 0; i <= nxt + 1; i++)
      for (j = 0; j <= nyt + 1; j++)
      {

         c[i][j][0] = c[i][j][1];
         c[i][j][nzt + 1] = c[i][j][nzt];
      }
}
