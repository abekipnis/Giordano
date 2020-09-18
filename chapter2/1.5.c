#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define GNUPLOT "gnuplot -persist"
/*void plot(float *x, float *y, int N, char *label){

for (int i=0; i<N; i++){
		fprintf(dat, "%lf %lf\n", x[N], y[N]);
	}
	fclose(dat);

       fprintf(gp, "set style data filledcurves; rep '-' lt rgb 'red'\n");
      fclose(gp);
}*/

int main(int argc, char **argv)
{
	float tau = 50;
	int Tend = 100;
	int N = 100;
	float dt = (float)Tend/N;
	float Na[N], Nb[N], t[N];
	Na[0] = 10;
	Nb[0] = 5;
	t[0] = 0;

	for (int i=1; i<N; i++) {
		//first order Eulers method
		Na[i] = Na[i-1] + (Nb[i-1]/tau - Na[i-1]/tau) * dt;
		Nb[i] = Nb[i-1] + (Na[i-1]/tau - Nb[i-1]/tau) * dt;
		t[i] = t[i-1] + dt;
	}
 	FILE *gp;
  gp = popen(GNUPLOT, "w");
  if (gp ==NULL){
          printf("Error opening Gnuplot");
          exit(0);
  }
	fprintf(gp, "set title 'dual nuclei decay pathway'\n");

	//have to plot all the points/lines in the same command
	fprintf(gp, "plot '-' w p ls 1 title 'Nb',\
			  '-' w p ls 2 title 'Na', \
		  	%lf*exp(-%lf*x)+%lf,\
			  %lf*exp(-%lf*x)+%lf\n",\
			  Na[0]/2.-Nb[0]/2, 2./tau, 0.5*(Nb[0]+Na[0]), -Na[0]/2.+Nb[0]/2., 2./tau, .5*(Na[0]+Nb[0]));

	for (int i=0; i<N; i++){
                fprintf(gp, "%g %g\n", t[i], Nb[i]);
        }
	fprintf(gp, "e\n"); //to mark end of this data
 	for (int i=0; i<N; i++){
                fprintf(gp, "%g %g\n", t[i], Na[i]);
      }
	fprintf(gp, "e\n");
	fprintf(gp, "set style line 1 lc rgb 'black'\n");
	fprintf(gp, "set style line 2 lc rgb 'red'\n");

	fclose(gp);
}
