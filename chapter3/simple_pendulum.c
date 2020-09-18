#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000
#define PI 3.14159265
#define TEND 10000.
int N;
FILE *gp;

float length, dt, drive_frequency;
float *theta, *omega, *t, *E;

void initialize(float th0){//, float* theta, float omega, float t, float E){

  dt = 0.04;
  N = TEND / dt; //N 16000
  theta = (float *) malloc(N*sizeof(float));
  t = (float *) malloc(N*sizeof(float));
  omega = (float *) malloc(N*sizeof(float));
  E = (float *) malloc(N*sizeof(float));

  t[0] = 0;

  theta[0] = th0; // radians
  omega[0] = 0; //radians/s
  length = 9.8; // meters
  E[0] = .5*pow(omega[0],2) + pow(theta[0]*length,2);
}

void calculate_euler(float theta[N], float omega[N], float t[N], float E[N]){
  float g = 9.8;
  float period = 2*3.14159 * sqrt(g/length);
  for (int i=0; i<N-1;i++){
    t[i+1] = t[i] + dt;
    omega[i+1] = omega[i] - (g/length)*theta[i]*dt;
    theta[i+1] = theta[i] + omega[i]*dt;
    E[i+1] = .5*pow(omega[i+1],2) + pow(theta[i+1]*length,2);
  }
}

void calculate_euler_cromer(float theta[N], float omega[N], float t[N], float E[N]){
  float g = 9.8;
  float drive_force = 1.2;
  drive_frequency = 2/3.;
  float q = 1/2.; //damping
  float period = 2*PI / sqrt(g/length); // natural period of pendulum (s)
  for (int i=0; i<N-1;i++){
    t[i+1] = t[i] + dt;
    omega[i+1] = omega[i]
                  - (g/length)*sin(theta[i])*dt
                  - q*omega[i]*dt
                  + drive_force * sin(drive_frequency * t[i]) * dt;
    theta[i+1] = theta[i] + omega[i+1]*dt;

    if (theta[i+1] < PI) {
      theta[i+1]+=2*PI;
    }
    if (theta[i+1] >PI){
      theta[i+1]-=2*PI;
    }
    E[i+1] = .5*pow(omega[i+1],2) + pow(theta[i+1]*length,2);
  }
}

void save(float theta[N], float omega[N], float t[N], char* filename){
    FILE *ptr;
    ptr = fopen(filename,"w");
    if (ptr ==NULL){
            printf("Error opening file");
            exit(0);
    }
    char gpcmd[GPCMDLEN];
    memset(gpcmd,0,sizeof(gpcmd));
    strcat(gpcmd, "time,theta,omega\n");
    char *row;
    for (int i=0; i<N-1; i++){
      if (0 > asprintf(&row, "%lf, %lf, %lf\n", t[i], theta[i], omega[i])) exit(0);
      strcat(gpcmd, row);
      free(row);
    }
    fprintf(ptr, "%s", gpcmd);
    fclose(ptr);
}

//Poincare section - ie only plot points in phase with driving force
// where t ~= n*pi/driving_freq |(t-n*pi/driving_freq)| < dt/2
// if we want to plot points out of phase w/ driving force
//    we plot points where t  ~=
void sv_strb(float th[N], float om[N], float t[N], char* fname, float str_ph){
    printf("saving poincaré section");
    FILE *ptr;
    ptr = fopen(fname,"w");
    if (ptr ==NULL){
      printf("Error opening file");
      exit(0);
    }
    char gpcmd[GPCMDLEN];
    memset(gpcmd,0,sizeof(gpcmd));
    strcat(gpcmd, "time,theta,omega\n");
    char *row;
    int nmax = TEND * drive_frequency / PI; // for poincare section
    for (int i=0; i<N-1; i++){
      for (int n=(int)drive_frequency*t[i]/(PI); n<nmax; n++){
        if (fabs(t[i] - n*(PI)/drive_frequency + str_ph) < dt/2){
          if (0 > asprintf(&row, "%lf, %lf, %lf\n", t[i], th[i], om[i])) exit(0);
          strcat(gpcmd, row);
          free(row);
          break;
        }
      }
    }
    fprintf(ptr, "%s", gpcmd);
    fclose(ptr);
}

int main(int argc, char **argv){
  float th0 = 0.02;
  char *f;
  char files[1000];
  char* ff;
  float strobe_phases[] = {PI/4, PI/3, PI/2, 0}; //PI over
  for (int i=0; i<4; i++){
    initialize(th0);
    calculate_euler_cromer(theta, omega, t, E);
    asprintf(&f, "Poincare_phi-f_d-%lf.dat", strobe_phases[i]);
    asprintf(&ff, "\"%s\" ",f);
    strcat(files,ff);
    sv_strb(theta, omega, t, f, strobe_phases[i]);
  }
  char* cmd;
  asprintf(&cmd, "python plot_pendulum.py --files %s --driving_freq=%lf",
            files, drive_frequency);
  printf("%s",cmd);
  system(cmd);
}
