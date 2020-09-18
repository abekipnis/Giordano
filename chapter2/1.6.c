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
	int Tend = 10;
	int N = 10000;
	float dt = (float)Tend/N;
	float Na[N],Nb[N], t[N];

	Na[0] = 50;
	Nb[0] = 50;
	t[0] = 0;

	float a, b, a1, b1;
	b = 0;
	a = 10;

	//testing another parameter family
	b1 = 3;
	a1 = 10;
	for (int i=1; i<N; i++) {
		//first order Eulers method
		Na[i] = Na[i-1] + (a*Na[i-1] - b*Na[i-1]*Na[i-1]) * dt;
		Nb[i] = Nb[i-1] + (a1*Nb[i-1] - b1*Nb[i-1]*Nb[i-1]) * dt;

		t[i] = t[i-1] + dt;
	}
 	FILE *gp;
  gp = popen(GNUPLOT, "w");
  if (gp ==NULL){
          printf("Error opening Gnuplot");
          exit(0);
  }
	fprintf(gp, "set title 'Population dynamics dn/dt=an-bn^2'\n");
	fprintf(gp, "set yrange [0:10]\n");

	//have to plot all the points/lines in the same command
	fprintf(gp, "plot '-' w p ls 1 title 'No death',\
							'-' w p ls 2 title 'Including N^2 death',\
		  				%lf*exp(%lf*x)\n",\
			  			Na[0], a);

	for (int i=0; i<N; i++){
  	fprintf(gp, "%g %g\n", t[i], Na[i]);
  }
	fprintf(gp, "e\n");
	for (int i=0; i<N; i++){
		fprintf(gp, "%g %g\n", t[i], Nb[i]);
	}
	fprintf(gp, "e\n");
	fprintf(gp, "set style line 1 lc rgb 'black'\n");
	fclose(gp);
}
