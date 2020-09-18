#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>


#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000
#define PI 3.14159265
#define TEND 50.

int N;
FILE *gp;

float dt;

void initialize(float x[N], float y[N], float z[N], float t[N]){
  t[0] = 0;
  x[0] = 1;
  y[0] = 0;
  z[0] = 0;
}

void RKint(float x[N], float y[N], float z[N], float t[N], float s, float b, float r){
  float kx1, ky1, kz1, mx, my, mz, kx2, ky2, kz2;
  for (int i=0; i<N;i++){
    t[i+1] = t[i] + dt;
    //Runge Kutta method
    //estimate the change in x,y,z using standard Euler integration
    kx1 = (s*(y[i]-x[i]))*dt;
    ky1 = (-x[i]*z[i] + r*x[i] - y[i])*dt;
    kz1 = (x[i]*y[i] - b * z[i])*dt;

    mx = x[i] + kx1/2.; //estimate x,y,z at the center of the interval
    my = y[i] + ky1/2.;
    mz = z[i] + kz1/2.;

    //now estimate the change in x,y,z at the center of the interval
    kx2 = (s*(my-mx))*dt;
    ky2 = (-mx*mz + r*mx-my)*dt;
    kz2 = (mx*my - b*mz)*dt;

    //use the slope at the center of the interval
    x[i+1] = x[i] + kx2;
    y[i+1] = y[i] + ky2;
    z[i+1] = z[i] + kz2;
  }
}

void sv(float x[N], float y[N], float z[N], float t[N], char* fname){
    FILE *ptr;
    ptr = fopen(fname,"w");
    if (ptr ==NULL){
      printf("Error opening file");
      exit(0);
    }
    char gpcmd[GPCMDLEN];
    memset(gpcmd,0,sizeof(gpcmd));

    strcat(gpcmd, "x,y,z,t\n");
    char *row;
    for (int i=0; i<N-1; i++){
      if (0 > asprintf(&row, "%lf,%lf,%lf,%lf\n", x[i], y[i], z[i], t[i])) exit(0);
      if (strlen(row)+strlen(gpcmd)>=GPCMDLEN){
        //write the current buffer and create new buffer
        fprintf(ptr, "%s", gpcmd);
        memset(gpcmd,0,sizeof(gpcmd));
      }
      strcat(gpcmd, row);
      free(row);
    }
    fprintf(ptr, "%s", gpcmd);

    fclose(ptr);
}

void sv_poincare(float x[N], float y[N], float z[N], float t[N], char* fname, float r){
  FILE *ptr;
  ptr = fopen(fname,"w");
  if (ptr ==NULL){
    printf("Error opening file");
    exit(0);
  }
  char gpcmd[GPCMDLEN];
  memset(gpcmd,0,sizeof(gpcmd));

  // strcat(gpcmd, "x,y,z,t\n");

  int sl = 0;
  char *row;
  for (int i=0; i<N-2; i++){
    if (t[i]>48){ //transients have gone away
      //for the bifurcation diagram, we want to plot the extrema
      // this is based on the assumption that we are not in the chaotic regime
      //  ie. the extrema repeat
      if (z[i-1]<z[i] && z[i+1]<z[i]){ //we are at a maxima
        //Usually theres only 5 or so maxima in the last 2 seconds
        if (0 > asprintf(&row, "%lf, %lf,%lf,%lf,%lf\n", r, x[i],y[i],z[i],t[i])) exit(0);
        sl = strlen(row)+strlen(gpcmd);
        if (sl>=GPCMDLEN){
          //write the current buffer and create new buffer
          fprintf(ptr, "%s", gpcmd);
          memset(gpcmd,0,sizeof(gpcmd));
        }
        strcat(gpcmd, row);
        free(row);
      }
    }
  }
  fprintf(ptr, "%s", gpcmd);

  fclose(ptr);
  free(x);
  free(y);
  free(z);
  free(t);
}

char* simulate(float r){
  float *x, *y, *z, *t, b, s;

  s = 10;
  b = 8./3.;

  x = (float *) malloc(N*sizeof(float));
  y = (float *) malloc(N*sizeof(float));
  z = (float *) malloc(N*sizeof(float));
  t = (float *) malloc(N*sizeof(float));

  char *f, *ff;

  initialize(x,y,z,t);
  RKint(x,y,z,t,s,b,r); //Runge Kutta integrate
  int np, id;
  np = omp_get_num_threads();
  id = omp_get_thread_num();
  // printf("Simulated for r = %lf from thread %d of %d\n", r, id, np);


  asprintf(&f, "lor_dat/r-%1.2lf.dat", r);
  asprintf(&ff, "\"%s\" ",f);

  if (r==160. || r==150. || r== 137.){
    sv(x, y, z, t, f);
  }
  sv_poincare(x,y,z,t,f,r);

  return ff;
}

int main(int argc, char **argv){
  dt = 0.0001;
  N = (int)TEND / dt;
  char files[1000000];
  char* ff;

  float rs[] = {5, 10, 25};

  #pragma omp parallel for num_threads ( 6 )
  for (int i=0; i<16000; i++){
  // for (int i=0; i<3; i++){
    ff = simulate(40+0.01*i);
    // ff = simulate(rs[i]);
    #pragma omp critical
    strcat(files,ff);
  }

  char* cmd;
  // asprintf(&cmd, "python -i lorenz_bifurc.py --files %s", files);

  asprintf(&cmd, "echo \"r,x,y,z,t\" > \"header\"; cat \"header\" %s > lorenz_bifurc.dat; python -i lorenz_bifurc.py --files lorenz_bifurc.dat", files);


  // asprintf(&cmd, "python plot_lorenz.py --files %s", files);
  printf("%s",cmd);
  system(cmd);
}
