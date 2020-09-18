#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
// #include "mercury_aux.h"

#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000
#define PI 3.1 4159265
#define TEND 1. //years

int N;
FILE *gp;

float dt;
void initialize(float x[N], float y[N], float t[N], float x0, float y0){
  t[0] = 0;
  x[0] = x0;
  y[0] = y0;
}

void sv(float x[N], float y[N], float t[N], char* fname){
  FILE *ptr;
  ptr = fopen(fname,"w");
  if (ptr ==NULL){
    printf("Error opening file");
    exit(0);
  }
  char gpcmd[GPCMDLEN];
  memset(gpcmd,0,sizeof(gpcmd));

  strcat(gpcmd, "x,y,t\n");

  int sl = 0;
  char *row;
  for (int i=0; i<N-2; i++){
    if (0 > asprintf(&row, "%lf,%lf,%lf\n", x[i],y[i],t[i])) exit(0);
    sl = strlen(row)+strlen(gpcmd);
    if (sl>=GPCMDLEN){
      fprintf(ptr, "%s", gpcmd);
      memset(gpcmd,0,sizeof(gpcmd));
    }
    strcat(gpcmd, row);
    free(row);
  }
  fprintf(ptr, "%s", gpcmd);

  fclose(ptr);
  free(x);
  free(y);
  free(t);
}

void Eint(float x[N], float y[N], float t[N], float vx0, float vy0){
  float r;
  float alpha = 0.01;
  float vx = vx0;
  float vy = vy0;
  float x1,y1, x2, y2;
  for (int i=0; i<N;i++){
    t[i+1] = t[i] + dt;
    r = sqrt(pow(x[i],2)+pow(y[i],2));
    vx = vx - (4* pow(PI,2) *x[i]*dt) *(1/pow(r,3) + alpha/pow(r,4));
    vy = vy - (4* pow(PI,2) *y[i]*dt) *(1/pow(r,3) + alpha/pow(r,4)) ;
    //using the new velocity value, rather than the old one
    // this is the Euler-Kromer method
    // keeps the energy of the system stable
    x[i+1] = x[i] + vx*dt;
    y[i+1] = y[i] + vy*dt;

    x2 = x[i]; y2 = y[i];
    x1 = x[i+1]; y1 = y[i+1];

    printf("%lf\n", pow((x1*x2+y1*y2)*x2-x1/sqrt(pow(x2,2)+pow(y2,2)),2) + pow((x1*x2+y1*y2)*y2-y1/sqrt(pow(x2,2)+pow(y2,2)),2));
  }
}

char* simulate(float vx0, float vy0, float x0, float y0){
  float *x, *y, *t;

  x = (float *) malloc(N*sizeof(float));
  y = (float *) malloc(N*sizeof(float));
  t = (float *) malloc(N*sizeof(float));

  char *f, *ff;

  initialize(x,y, t, x0, y0);
  Eint(x, y, t, vx0, vy0); //Euler Kromer integrate

  asprintf(&f, "mercury_traj.dat");
  asprintf(&ff, "%s",f);
  sv(x,y,t,ff);

  return ff;
}

int main(int argc, char **argv){
  //dt should be < 1% of the characteristic timescale of system (the period, 1 year)
  dt = 0.0001;
  N = (int)TEND / dt;
  simulate(0, 8.2, 0.47, 0);
  system("python plot_planetary_traj.py --f \"mercury_traj.dat\"");
}
