//baseball, including spin effect on drag force
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000

int Tend, N, NT;
float dt;
FILE *gp;

void calculate(float x[N], float y[N], float z[N], float dt, float v_init, float theta, float B_m, float w, float s0m){
  for (int i =0 ; i < N; i++){
    x[i] = 0;
    y[i] = 0;
    z[i] = 0;
  }
  float vx, vy, vz,v, f, a, delta;
  vx = v_init*cos(theta);
  vy = v_init*sin(theta);
  vz = 0;
  int i;
  float vd = 35.;
  delta = 5.;
  for(i=1; i<N; i++){
    x[i] = x[i-1] + vx*dt;
    y[i] = y[i-1] + vy*dt;
    z[i] = z[i-1] + vz*dt;
    v = sqrt(vx*vx + vy*vy + vz*vz);
    B_m = 0.0039 + 0.0058/(1+exp(((v-vd)/delta)));
    f = B_m * v;
    vy = vy - 9.8*dt;
    vx = vx - f*vx *dt;
    vz = vz - s0m * vx * w * dt;
  }
}

void load_str(char gpcmd[GPCMDLEN], float x[N], float y[N], float z[N]){
  char *row;
  for (int i=0; i<N; i++){
    if (x[i] > 60) break;
    if (0 > asprintf(&row, "%lf %lf\n", x[i], y[i])) exit(0);
    strcat(gpcmd, row);
    free(row);
  }
  strcat(gpcmd, "e\n");
  for (int i=0; i<N; i++){
    // if
    if (0 > asprintf(&row, "%lf %lf\n", x[i], z[i])) exit(0);
    strcat(gpcmd, row);
    free(row);
  }
  strcat(gpcmd, "e\n");
}

int main(int argc, char **argv){
  Tend = 5;
  N = 1000;
  dt = (float)Tend/N;
  float v0, B_m, s0m, w;
  v0 = 70; // 70mph
  v0*=0.4470; //meters per second
  B_m = 4e-5;
  // S0/m, ratio of ball drag force to mass
  s0m = 4.1e-4; // m=149g
  w = 30.; // 30 revolutions per second

  float x[N];
  float y[N];
  float z[N];

  gp = popen(GNUPLOT, "w");
  if (gp ==NULL){
          printf("Error opening Gnuplot");
          exit(0);
  }
  char gpcmd[GPCMDLEN];
  strcat(gpcmd, "set title 'Trajectory of cannon shell'\n");
  strcat(gpcmd, "set xlabel 'x(feet)'\n");
  strcat(gpcmd, "set ylabel 'y,z(feet)'\n");
  strcat(gpcmd, "set xrange [0:60]\n");

  int thetas[] = {5};

  char g[1000];
  strcat(gpcmd,  "plot ");
  sprintf(g,  "'-' with line ls %d title 'y',", 1);
  strcat(gpcmd, g);
  sprintf(g,  "'-' with line ls %d title 'z'", 2);
  strcat(gpcmd, g);

  // strcat(gpcmd,",");

  strcat(gpcmd,"\n");
  calculate(x, y, z, dt, v0, thetas[0]*3.1415/180.0, B_m, w, s0m);
  load_str(gpcmd, x, y, z);
  fprintf(gp, "%s", gpcmd);
}
