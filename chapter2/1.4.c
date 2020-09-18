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
	float ta = 1;
	float tb = 1.1;
	int Tend = 10;
	int N = 100;
	float dt = (float)Tend/N;
	float Na[N], Nb[N], t[N];
	Na[0] = 10;
	Nb[0] = 5;
	t[0] = 0;

	for (int i=1; i<N; i++) {
		//first order Eulers method
		Na[i] = Na[i-1] - Na[i-1]/ta*dt;
		Nb[i] = Nb[i-1] + (Na[i-1]/ta - Nb[i-1]/tb) * dt;
		t[i] = t[i-1] + dt;
	}
 	FILE *gp;
        gp = popen(GNUPLOT, "w");
        if (gp ==NULL){
                printf("Error opening Gnuplot");
                exit(0);
        }
	fprintf(gp, "set title 'two nuclei decay pathway'\n");

	fprintf(gp, "plot %1.2lf*exp(-%1.2lf*x) with lines linestyle 1\n", Na[0], 1.0/ta);
//	fprintf(gp, "%lf*exp(-%lf*x)+%lf*exp(-%lf*x)\n", Nb[0]-Na[0]*(ta/tb-1.0),tb, Na[0]*(ta/tb-1), 1/ta);


	float A0; A0 = Na[0];
	float x; x = 1.0/(ta/tb-1);
	fprintf(gp, "plot '-' w p ls 1 title 'Nb',\
			 '-' w p ls 2 title 'Na', \
			%1.2lf*exp(-%1.2lf*x) with lines linestyle 1,\
			(%lf*exp(-%lf*x)+%lf*exp(-%lf*x))\n", A0, 1.0/ta, Nb[0]-x*Na[0], 1.0/tb, x*Na[0], 1.0/ta);

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
//	plot(t, Na, N, "Na");
//	plot(t, Nb, N, "Nb");
}
