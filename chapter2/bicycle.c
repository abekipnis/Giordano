//bicycle (page 19)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define GNUPLOT "gnuplot -persist"
#define MEMORY_BLOCK 100

int plot_data(int N, int npwr, int ns, float power, float vd[N], float sls[ns], float powers[npwr], float mass, float v[N], float t[N],  float vdp[npwr][N], float vg[ns][N]){

    FILE *gp;
    gp = popen(GNUPLOT, "w");
    if (gp ==NULL){
            printf("Error opening Gnuplot");
            exit(0);
    }
    char gpcmd[1000000];
    // printf("here");
    strcat(gpcmd, "set title 'Biker velocity'\n");
    strcat(gpcmd, "set yrange [0:35]\n");
    strcat(gpcmd, "plot '-' w l ls 1 title 'Numerical solution, no air resistance', ");
    strcat(gpcmd, "'-' with line title 'Numerical solution, with air resistance' lc rgb 'red', ");
    // strcat(gpcmd, "'-' w p ls 3 title 'Numerical solution, with 30pct less air resistance', ");

    char g[2000000];
    for (int i=0; i<npwr; i++){
      sprintf(g, "'-' w l ls %d title 'Numerical solution, power=%d 30pct less SA',",i+2, (int)powers[i]);
      strcat(gpcmd, g);
    }
    for (int i=0; i<ns; i++){
      sprintf(g, "'-' w l ls %d title 'Numerical solution, slope=%1.2lf 30pct SA',",i+2, sls[i]);
      strcat(gpcmd, g);
    }
    sprintf(g, "sqrt(%lf**2 + 2*%lf*x/%lf) title 'Analytic solution, no air resistance'\n", v[0], power, mass);
    strcat(gpcmd, g);



    /*
      no air resistance
    */
    char *row;
    for (int i=0; i<N; i++){
      if (0 > asprintf(&row, "%lf %lf\n", t[i], v[i])) exit(0);
      strcat(gpcmd, row);
      free(row);
    }
    strcat(gpcmd, "e\n");



    /*
      w/resistance
    */
    for (int i=0; i<N; i++){
      if (0 > asprintf(&row, "%lf %lf\n", t[i], vd[i])) exit(0);
      strcat(gpcmd, row);
      free(row);
    }
    strcat(gpcmd, "e\n");


    /*
     with reduced power, but 30% less surface area
    */
    for (int i=0; i<npwr; i++){
      for (int j=0; j<N; j++){
        if (0 > asprintf(&row, "%lf %lf\n", t[j], vdp[i][j])) exit(0);
        strcat(gpcmd, row);
        free(row);
      }
      strcat(gpcmd,"e\n");
    }

    /*
     30% less surface area, with additional slope
    */
    for (int i=0; i<ns; i++){
      for (int j=0; j<N; j++){
        if (0 > asprintf(&row, "%lf %lf\n", t[j], vg[i][j])) exit(0);
        strcat(gpcmd, row);
        free(row);
      }
      strcat(gpcmd,"e\n");
    }

    //print the whole gnuplot command to gnuplot
    fprintf(gp, "%s", gpcmd);
  return 0;
}

int main(int argc, char **argv){
  int Tend = 200;
  int N = 100;
  float v0; v0 = 4; // initial velocity, m/s
  float dt = (float)Tend/N;

  float g; g = 9.8; //gravitational constant

  float power; power = 404; //watts (kg m^2 s^(-3))
  float mass; mass = 70; //kg

  float v[N]; //velocity array
  float vd[N]; //velocity with drag
  float vf[N]; //velocity with reduced frontal area (ie middle of peloton)
  float t[N]; //time

  //investigating the effects of drag on input power in peloton
  int npwr = 5; // testing npwr different number of powers
  float vdp[npwr][N];
  float powers[npwr];
  for (int i=0; i < npwr; i++){
    powers[i] = power-i*11;
    vdp[i][0] = v0;
  }

  //investigating the effects of mountainous terrain
  int ns = 5; // number of slopes
  float vg[ns][N]; // testing what slope is required for vterm of 70 mph = 30 m/s
  float sls[ns]; // slopes
  for (int i=0; i<ns; i++){
    sls[i] = 0.01+0.015*i;
    vg[i][0] = v0;
  }

  //initial time and velocity
  t[0] = 0;
  v[0] = v0; // m/s
  vd[0] = v0; // + drag
  float a; a = 0.33; //frontal area of rider, m^2
  float rho; rho = 1.225; //kg/m^3
  float c; c = 0.5; //drag coefficient

  for (int i=1; i<N; i++){
    t[i] = t[i-1] + dt;
    v[i] = v[i-1] + power*dt / (mass * v[i-1]); //Euler method
    vd[i] = vd[i-1] + power*dt / (mass * vd[i-1]) - c*rho*a*pow(vd[i-1],2)/mass*dt;
    for (int j=0; j<npwr; j++){
      vdp[j][i] = vdp[j][i-1] + powers[j]*dt / (mass * vdp[j][i-1]) - c*rho*a*0.7*pow(vdp[j][i-1],2)/mass*dt;
    }
    for (int j=0; j<ns; j++){
      vg[j][i] = vg[j][i-1] + power*dt / (mass * vg[j][i-1]) - c*rho*a*0.7*pow(vg[j][i-1],2)/mass*dt + g*sin(sls[j])*dt;
    }
  }

  // printf("Final velocity of person @ front of peloton: %lf\n");

  //finding the power a rider at the back of the peloton has to exert to keep up w/ front
  float split = 10;
  float eps = 0.005;
  float alpha = 0.5;

  float vt = vd[N-1]; //velocity target (peloton)
  float p0 = powers[npwr-1];
  float vfp = vdp[npwr-1][N-1]; // initial guess

  int steps = 0;


  while (fabsf(vt-vfp) > eps){
    if (vfp > vt){
      while (vfp > vt){
        p0-=split;
        for (int i=1; i<N; i++){
          vdp[npwr-1][i] = 0; //resetting the DE solver
        }
        for (int i=1; i<N; i++){
          vdp[npwr-1][i] = vdp[npwr-1][i-1] + p0*dt / (mass * vdp[npwr-1][i-1]) - c*rho*a*0.7*pow(vdp[npwr-1][i-1],2)/mass*dt;
        }
        vfp = vdp[npwr-1][N-1];
      }
    }
    else if (vfp < vt){
      while (vfp < vt) {
        p0 += split;
        for (int i=1; i<N; i++){
          vdp[npwr-1][i] = 0; //resetting the DE solver
        }
        for (int i=1; i<N; i++){
          vdp[npwr-1][i] = vdp[npwr-1][i-1] + p0*dt / (mass * vdp[npwr-1][i-1]) - c*rho*a*0.7*pow(vdp[npwr-1][i-1],2)/mass*dt;
        }
        vfp = vdp[npwr-1][N-1];
      }
    }
    split*=alpha;
    steps+=1;
    printf("split: %lf\n", split);
  }
  printf("converged on solution in %d steps\n", steps);
  printf("power for rider in middle of peloton, 70 pct air resistance as front: %lf\n", p0);
  printf("power for rider in front of peloton, 100pct air resistance: %lf\n", power);
  powers[npwr-1] = p0;

  //finding the slope a rider needs to reach 70 mph
  vt = 31.2929; //velocity target (peloton)
  split = 0.1;
  float s0 = sls[0]; // initial guess for slope
  p0 = 400; //Watts
  vfp = vg[ns-1][N-1]; // initial guess

  steps = 0;
  while (fabsf(vt-vfp) > eps){
    if (vfp > vt){
      while (vfp > vt){
        s0-=split;
        for (int i=1; i<N; i++){
          vdp[ns-1][i] = 0; //resetting the DE solver
        }
        for (int i=1; i<N; i++){
          vg[ns-1][i] = vg[ns-1][i-1] + p0*dt / (mass * vg[ns-1][i-1]) - c*rho*a*0.7*pow(vg[ns-1][i-1],2)/mass*dt + g*sin(s0)*dt;
        }
        vfp = vg[ns-1][N-1];
        steps+=1;
      }
    }
    else if (vfp < vt){
      while (vfp < vt) {
        s0 += split;
        for (int i=1; i<N; i++){
          vdp[ns-1][i] = 0; //resetting the DE solver
        }
        for (int i=1; i<N; i++){
          vg[ns-1][i] = vg[ns-1][i-1] + p0*dt / (mass * vg[ns-1][i-1]) - c*rho*a*0.7*pow(vg[ns-1][i-1],2)/mass*dt + g*sin(s0)*dt;
        }
        vfp = vg[ns-1][N-1];
        steps+=1;
      }
    }
    split*=alpha;
    steps+=1;
    printf("split: %lf\n", split);
  }
  printf("converged on solution in %d steps\n", steps);
  printf("slope for rider w 70 pct air resistance as front to get to 70mph: %lf\n", s0);
  sls[ns-1] = s0;




  //creating gnuplot command in a string to pass all at once using fprintf
  plot_data(N, npwr, ns, power, vd, sls, powers, mass, v, t, vdp, vg);
}



  /*
  PSEUDOCODE FOR ROOT FINDING
  TO GET THE POWER REQUIRED FOR A RIDER IN THE MIDDLE OF THE PELOTON TO KEEP UP

  find v_f(target)
  find the power (from the ones already calculated) that give the closest v_f to v_f(target)
      or, take a guess. call this power p

  while abs(vtarg-vfp) > eps:
    get v_f(p)
    if v_f(p) > vtarg:
      while v_f(p) > vtarg
       subtract split from p
       get v_f(p)
    else if vfp < vtarg:
      while v_f(p) < vtarg
        add split to p
        get vfp
    split*=alpha

  */
