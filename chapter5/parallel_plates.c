#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
// #include "mercury_aux.h"

#define GNUPLOT "gnuplot -persist"
int Nx = 50;
int Ny = 50;
float **v1, **v2;

void update(float **v1, float **v2, int Nx, int Ny, float diff){
  float tmp;
  for (int i=1;i<Nx;i++){
    for (int j=1;j<Ny;j++){
      tmp = v2[i][j];
      v2[i][j] = (v1[i-1][j]+v1[i+1][j] + v1[i][j+1] + v1[i][j-1])/4;
      diff+=abs(tmp - v2[i][j]);
    }
  }
  diff = diff/(Nx*Ny); //average diff btw current & next iteration
}

void initialize(){
  v1 = (float **) malloc(Nx*sizeof(float *));
  v2 = (float **) malloc(Nx*sizeof(float *));
  for (int i =0 ;i< Nx; i++) {
    v1[i] = (float *) malloc(Ny*sizeof(float));
    v2[i] = (float *) malloc(Ny*sizeof(float));
  }

  //initialize two parallel plates
  // plate on the left at -1V
  // plate on the right at +1V
  // top & bottom lines between plates aea linear fcn btw +-1V
  for (int j=0;j<Ny; j++){
    v1[0][j] = -1;
    v1[Nx][j] = 1;
  }
  float val = 0;
  for (int j=0;j<Ny; j++){
    val = 2.0/Ny*j;
    v1[j][0] =val;
    v1[j][Ny] =val;
  }
}

void calculate(){
  initialize();
  int i = 0;
  float diff = 10;
  while(diff > 1e-5){
    update(v1, v2, Nx, Ny, diff);
    update(v2,v1, Nx, Ny, diff);
  }
}

int main(int argc, char **argv){
  calculate();
}
