#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>


#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000
#define PI 3.14159265
#define TEND 1000.

int N;
FILE *gp;

float length, dt, drive_frequency;

void initialize(float th0, float* theta, float* omega, float* t, float* E){
  t[0] = 0;
  theta[0] = th0; // radians
  omega[0] = 0; //radians/s
  length = 9.8; // meters
  E[0] = .5*pow(omega[0],2) + pow(theta[0]*length,2);
}


void calculate_euler_cromer(float theta[N], float omega[N], float t[N], float E[N], float drive_force){
  float g = 9.8;

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

//Poincare section - ie only plot points in phase with driving force
// where t ~= n*pi/driving_freq |(t-n*pi/driving_freq)| < dt/2
// if we want to plot points out of phase w/ driving force
//    we plot points where t  ~=
void sv_strb(float th[N], float om[N], float t[N], char* fname){
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
      if (i*dt > 300*drive_frequency){
        for (int n=(int)drive_frequency*t[i]/(PI); n<nmax; n++){
          if (fabs(t[i] - n*(PI)/drive_frequency + 0) < dt/2){
            if (0 > asprintf(&row, "%lf, %lf, %lf\n", t[i], th[i], om[i])) exit(0);
            strcat(gpcmd, row);
            free(row);
            break;
          }
        }
      }
    }
    fprintf(ptr, "%s", gpcmd);
    fclose(ptr);
}


char* simulate(int i, int j){
  float th0s[] = {0.02, 0.05, 0.10, 0.2, 0.3, 0.4,0.5,0.6,0.785, 0.9, 1, 1.1, 1.2, 1.3, 1.4};

  // float th0 = 0.02;
  float th0 = th0s[j];
  float *theta, *omega, *t, *E;
  float drive_force = 1.35+ i*0.0005;

  theta = (float *) malloc(N*sizeof(float));
  t = (float *) malloc(N*sizeof(float));
  omega = (float *) malloc(N*sizeof(float));
  E = (float *) malloc(N*sizeof(float));
  char *f, *ff;

  initialize(th0, theta, omega, t, E);
  calculate_euler_cromer(theta, omega, t, E, drive_force);
  asprintf(&f, "bif_dat/th0-%1.2lf_phi-f_d-%1.3lf.dat", th0, drive_force);
  asprintf(&ff, "\"%s\" ",f);
  sv_strb(theta, omega, t, f);

  free(theta);
  free(t);
  free(omega);
  free(E);

  return ff;
}

int main(int argc, char **argv){
  dt = 0.01;
  N = (int)TEND / dt;
  char files[1000000];
  char* ff;

  for (int j=0;j<16; j++){
    #pragma omp parallel for num_threads ( 2 )
      for (int i=0; i<300; i++){
        ff = simulate(i,j);
        #pragma omp critical
        strcat(files,ff);
    }
  }

  char* cmd;
  asprintf(&cmd, "python plot_bifurc.py --files %s --driving_freq=%lf",
            files, drive_frequency);
  printf("%s",cmd);
  system(cmd);
}
