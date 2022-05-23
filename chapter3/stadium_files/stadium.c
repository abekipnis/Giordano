#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#define GNUPLOT "gnuplot -persist"
#define GPCMDLEN 1000000
#define PI 3.14159265
#define TEND 300.

int N;
FILE *gp;

float dt;

void initialize(float x[N], float y[N], float t[N], float x0, float y0){
  t[0] = 0;
  x[0] = x0;
  y[0] = y0;
}

void Eint(float x[N], float y[N], float t[N], float a, float vx0, float vy0, float r){
  float kx1, ky1, kz1, vx, vy, dx, dy, vperpx, vperpy, vparx, vpary, vidn;
  vx = vx0;
  vy = vy0;

  for (int i=0; i<N;i++){
    t[i+1] = t[i] + dt;

    //estimate the change in x,y,z using standard Euler integration
    kx1 = (vx)*dt;
    ky1 = (vy)*dt;

    float xn, yn; // the next point in xy space
    xn = x[i] + kx1;
    yn = y[i] + ky1;

    /*
    check for collision
    */
    float delt = 0;
    // if (yn>-a*r && yn<a*r){ // in the rectangle of the stadium
    //   //only have to check for collision in x
    if(xn<-r){// collision
      //backtrack
      //start@ x[i], y[i], integrate 100x smaller step until another collision
      xn = x[i];
      yn = y[i];

      while(xn>-r){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01; // keep track of how much time has passed
      }
      xn += (vx)*dt*0.01; //take col. pt as pt after col. detected
      yn += (vy)*dt*0.01;
      delt += dt*0.01;

      vx*=-1; //reflection off vertical - only x component reflects

      //integrate until timestep completed
      while(delt<dt){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01;
      }
      x[i+1] = xn;
      y[i+1] = yn;
    }
    else if(xn>r){// collision
      //backtrack
      //start@ x[i], y[i], integrate 100x smaller step until another collision
      xn = x[i];
      yn = y[i];

      while(xn<r){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01; // keep track of how much time has passed
      }
      xn += (vx)*dt*0.01; //take col. pt as pt after col. detected
      yn += (vy)*dt*0.01;
      delt += dt*0.01;

      vx*=-1; //reflection off vertical - only x component reflects

      //integrate until timestep completed
      while(delt<dt){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01;
      }
      x[i+1] = xn;
      y[i+1] = yn;
    }
    //equation for a circle centered at x=0, y=-a*r
    else if(pow(xn,2)+pow((yn+a*r),2)>pow(r,2) && yn<-a){ //handle collision
      //backtrack
      //start@ x[i], y[i], integrate 100x smaller step until another collision
      xn = x[i];
      yn = y[i];
      while (pow(xn,2)+pow((yn+a*r),2)<pow(r,2) && yn<-a){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01; // keep track of time
      }
      xn += (vx)*dt*0.01; //take col. pt as pt after col. detected
      yn += (vy)*dt*0.01;
      delt += dt*0.01;

      // calculate the normal to the wall at the collision point
      /*
        can get the x and y components of the normal vector
        by getting the slope of the line between the edge and the center
        ie getting dx and dy between the center and the edge
        and dividing by dx**2 + dy**2 to normalize
      */
      dx = xn; //origin of the circle is always at x=0
      dy = yn + a;

      float mag = sqrt(pow(dx,2) + pow(dy,2));
      dx/=mag;
      dy/=mag;
      //get velocity components perpendicular to wall so we can reflect them
      // v_(i, perp) = (v_i*\hat{n})\hat{n})
      // dot product of initial velocity with unit normal to wall
      vidn = vx*dx + vy*dy;
      vperpx = vidn * dx;
      vperpy = vidn * dy;

      //velocity components parallel to the wall
      vparx = vx - vperpx;
      vpary = vy - vperpy;

      //the perpendicular components reflect! parallel components stay the same
      vperpx*=-1;
      vperpy*=-1;

      //the new velocity components
      vx = vperpx + vparx;
      vy = vperpy + vpary;

      //integrate until the timestep is completed
      while(delt<dt){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01;
      }
      x[i+1] = xn;
      y[i+1] = yn;
    }
    //upper hemisphere
    else if(pow(xn,2)+pow((yn-a*r),2)>pow(r,2) && yn>a){ //handle collision
      //backtrack
      //start@ x[i], y[i], integrate 100x smaller step until another collision
      xn = x[i];
      yn = y[i];
      while (pow(xn,2)+pow((yn-a*r),2)<pow(r,2) && yn>a){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01; // keep track of time
      }
      xn += (vx)*dt*0.01; //take col. pt as pt after col. detected
      yn += (vy)*dt*0.01;
      delt += dt*0.01;

      // calculate the normal to the wall at the collision point
      /*
        can get the x and y components of the normal vector
        by getting the slope of the line between the edge and the center
        ie getting dx and dy between the center and the edge
        and dividing by dx**2 + dy**2 to normalize
      */
      dx = xn; //origin of the circle is at x=0
      dy = yn - a;

      float mag = sqrt(pow(dx,2) + pow(dy,2));
      dx/=mag;
      dy/=mag;

      //get velocity components perpendicular to wall so we can reflect them
      // v_(i, perp) = (v_i*\hat{n})\hat{n})
      // dot product of initial velocity with unit normal to wall
      vidn = vx*dx + vy*dy;
      vperpx = vidn * dx;
      vperpy = vidn * dy;

      //velocity components parallel to the wall
      vparx = vx - vperpx;
      vpary = vy - vperpy;

      //the perpendicular components reflect! parallel components stay the same
      vperpx*=-1;
      vperpy*=-1;

      //the new velocity components
      vx = vperpx + vparx;
      vy = vperpy + vpary;

      //integrate until the timestep is completed
      while(delt<dt){
        xn += (vx)*dt*0.01;
        yn += (vy)*dt*0.01;
        delt += dt*0.01;
      }
      x[i+1] = xn;
      y[i+1] = yn;
    }
    else{
      //no collision, business as usual
      kx1 = (vx)*dt;
      ky1 = (vy)*dt;

      xn = x[i] + kx1;
      yn = y[i] + ky1;

      x[i+1] = xn;
      y[i+1] = yn;
    }
    // printf("Energy: %lf\n", pow(vx,2)+pow(vy,2));

  }
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
    char *row;
    for (int i=0; i<N-1; i++){
      if (0 > asprintf(&row, "%lf,%lf,%lf\n", x[i], y[i], t[i])) exit(0);
      if (strlen(row)+strlen(gpcmd)>=GPCMDLEN){
        fprintf(ptr, "%s", gpcmd);
        memset(gpcmd,0,sizeof(gpcmd));
      }
      strcat(gpcmd, row);
      free(row);
    }
    fprintf(ptr, "%s", gpcmd);

    fclose(ptr);
}

void sv_poincare(float x[N], float y[N], float t[N], char* fname, float r){
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

char* simulate(float a, float r, float vx0, float vy0, float x0, float y0, int n, int nit){
  float *x, *y, *t;

  x = (float *) malloc(N*sizeof(float));
  y = (float *) malloc(N*sizeof(float));
  t = (float *) malloc(N*sizeof(float));

  char *f, *ff;

  initialize(x,y, t, x0, y0);
  Eint(x, y, t, a, vx0, vy0, r); //Runge Kutta integrate
  int np, id;
  np = omp_get_num_threads();
  id = omp_get_thread_num();

  asprintf(&f, "sta_dat/traj/x0-%1.2lf_y0-%1.2lf_r-%1.2lf_a-%1.2lf_%d_nit-%d.dat",  x0, y0, r, a, n, nit);
  asprintf(&ff, "\"%s\" ",f);
  sv_poincare(x,y,t,f,r);

  return ff;
}

int main(int argc, char **argv){
  dt = 0.01;
  N = (int)TEND / dt;
  printf("N %d\n",N );
  char files[1000000];
  char* ff, *ff2;

  char* cmd;

  //grid of initial conditions
  int nx, ny, nit; //nit = number of tries with random changes
  nx = ny = 4; //we have n*n grid points
  nit = 4;

  float x0, y0, vx0, vy0;
  float alpha = 0.01;
  float radius = 1;

  float alphas[10];
  float lambdas[10];
  for (int i=0; i<10; i++){
    alphas[i] = pow(10, i-8);
  }
  clock_t t;
  t = clock();
  for(int a=0; a<10; a++){
    alpha = alphas[a];
    printf("testing alpha=%.10lf\n", alpha);
    for(int i=1; i<nx; i++){
      // #pragma omp parallel for num_threads ( 2 ) private(vx0, vy0, cmd)
      for (int j=1; j<ny; j++){
        vx0 = (float)rand()/RAND_MAX;
        vy0 = (float)rand()/RAND_MAX;

        #pragma omp parallel for num_threads ( 2 ) private(files, ff)
        for (int n=0; n<nit; n++){
          char* incmd, *ff, *ff2;

          //start randomly on circle of radius r
          float r =  (float)rand()/RAND_MAX;
          float theta = (float)rand()/RAND_MAX*2*PI;
          x0 = r*cos(theta);
          y0 = r*sin(theta);

          ff = simulate(alpha, radius, vx0, vy0, x0, y0, 1, n);
          #pragma omp critical
          strcat(files,ff);

          float dx = (float)rand()/RAND_MAX*1e-3;
          float dy = (float)rand()/RAND_MAX*1e-3;

          printf("%lf, %lf : %lf, %lf\n", x0, y0, x0 + dx, y0 + dy);

          ff2 = simulate(alpha, radius, vx0, vy0, x0 + dx, y0 + dy, 2, n);
          #pragma omp critical
          strcat(files,ff2);

          // get distance between two trajectories
          asprintf(&incmd, "paste -d ',' %s | "
                          "awk -F',' 'FNR>2  {"
                          " printf(\"%%s, %%lf\\n \", $3, sqrt(($1-$4)**2+($2-$5)**2)); }'"
              ">  sta_dat/diff/x-%lf_y-%lf_n-%d_i-%d_j-%d.dat", files, x0, y0, n, i, j);
          system(incmd);

          // asprintf(&incmd, "python plot_stadium.py --files %s", files);
          // system(incmd);


          #pragma omp critical
          asprintf(&incmd, "rm %s", ff);
          system(incmd);

          #pragma omp critical
          asprintf(&incmd, "rm %s", ff2);
          system(incmd);

          #pragma omp critical
          memset(&files[0], 0, sizeof(files));
        }

        //average files for each x0, y0 test of n random samples
        asprintf(&cmd,  "paste -d ',' sta_dat/diff/x-%lf_y-%lf_n-?_i-%d_j-%d.dat | "
                        "awk  -F',' '{ for(i=2;i<=NF;i+=2) array[$1]+=$i; if (i = NF) printf(\"%%s, %%lf\\n \", $1, array[$1]/NF*2); }' > "
                        " sta_dat/avg_diff/x-%lf_y-%lf.dat", x0, y0, i, j, x0, y0);
        system(cmd);

        asprintf(&cmd, "rm sta_dat/diff/x-%lf_y-%lf_n-*_i-%d_j-%d.dat",  x0, y0, i, j);
        system(cmd);
        memset(&cmd[0], 0, sizeof(cmd));
      }
    }

    // paste -d ',' sta_dat/avg_diff_*.dat | awk  -F ',' '{ for(i=2;i<=NF;i+=2) array[$1]+=$(i+1); if (i = NF) printf("%lf, %lf\n", $2, array[$1]/NF*2); }'
    asprintf(&cmd,  "paste -d ',' sta_dat/avg_diff/*.dat | "
                    "awk  -F ', ' '{ for(i=2;i<=NF;i+=2) array[$1]+=$(i); if (i = NF) printf(\"%%lf, %%lf\\n\", $1, array[$1]/NF*2); }' > "
                    " sta_dat/avg_diff_over_init_conds_alpha-%.10lf.dat", alpha);
    printf("%s\n", cmd);
    system(cmd);


    //run the Lyanpunov exponent code in Python and get the return value lambda
    asprintf(&cmd, "python plot_stad_lyapunov.py --files \"sta_dat/avg_diff_over_init_conds_alpha-%.10lf.dat\"", alpha);
    FILE *fp;
    char path[1035];
    fp = popen(cmd, "r");
    if (fp == NULL) {
      printf("Failed to run command\n" );
      exit(1);
    }
    float lambda;
    while (fgets(path, sizeof(path), fp) != NULL) {
      lambda = atof(path);
    }
    /* close */
    pclose(fp);
    // printf("%lf", lambda);
    lambdas[a] = lambda;
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("took %lf seconds to calculate l for alpha=%lf\n", time_taken, alpha);
  }

  FILE* lva;
  lva = fopen("lambda_v_alpha.dat", "w");
  fprintf(lva, "alpha, lambda\n");
  char *row;
  for (int i =0;i<10; i++) {
    asprintf(&row, "%lf, %lf\n", alphas[i], lambdas[i]);
    fprintf(lva, "%s", row);
  }
}
