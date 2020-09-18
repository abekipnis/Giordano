#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000

int Tend, N, NT;
float dt;
FILE *gp;

void calculate(float x[N], float y[N], float dt, float v_init, float theta, float B_m){
  for (int i =0 ; i < N; i++){
    x[i] = 0;
    y[i] = 0;
  }
  float vx, vy, f, a;
  vx = v_init*cos(theta);
  vy = v_init*sin(theta);
  int i;
  for(i=1; i<N; i++){
    x[i] = x[i-1] + vx*dt;
    y[i] = y[i-1] + vy*dt;
    f = B_m * sqrt(vx*vx + vy*vy) * exp(-y[i]/1e4);
    vy = vy - 9.8*dt -f*vy*dt;
    vx = vx - f*vx *dt;
    if (y[i] < 0) break;
  }
  a = -y[i] / y[i-1];
  x[i] = (x[i] +a*x[i-1]) /(a+1);
  y[i] = 0;
}

void load_str(char gpcmd[GPCMDLEN], float x[N], float y[N]){
  char *row;
  for (int i=0; i<N; i++){
    if (0 > asprintf(&row, "%lf %lf\n", x[i]/1000., y[i]/1000.)) exit(0);
    strcat(gpcmd, row);
    free(row);
  }
  strcat(gpcmd, "e\n");
}

int main(int argc, char **argv){
  Tend = 200;
  N = 100;
  dt = (float)Tend/N;
  float v0; v0 = 700;
  float B_m = 4e-5;

  float x[N];
  float y[N];

  gp = popen(GNUPLOT, "w");
  if (gp ==NULL){
          printf("Error opening Gnuplot");
          exit(0);
  }
  char gpcmd[GPCMDLEN];
  strcat(gpcmd, "set title 'Trajectory of cannon shell, corrected for air density altitude dependent drag force'\n");
  strcat(gpcmd, "set xlabel 'x(km)'\n");
  strcat(gpcmd, "set ylabel 'y(km)'\n");

  int thetas[] = {30, 35, 40, 45, 50, 55};

  char g[1000];
  strcat(gpcmd,  "plot ");
  for (int i=0; i<6; i++){
    sprintf(g,  "'-' with line ls %d title 'Numerical solution, angle %d B2/m=4e-5/m'", i, thetas[i]);
    strcat(gpcmd, g);
    if (i<5){
      strcat(gpcmd,",");
    }
  }
  strcat(gpcmd,"\n");
  for (int i=0; i<6; i++){
    calculate(x, y, dt, v0,  thetas[i]*3.1415/180.0, B_m);
    load_str(gpcmd, x, y);
  }
  fprintf(gp, "%s", gpcmd);
}
