#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define GNUPLOT "gnuplot -persist"
void plot(float *x, float *y){
	
}


int main(int argc, char **argv)
{
	float a = 10; //applied force 
	float b = 1; //damping, friction
	int N = 1000; //number of points
	float v[N], t[N];
	float dt = 0.01;
	v[0] = 0; 
	t[0] = 0; 
	for (int i=1; i<N; i++){
		//Eulers:  dv/dt = a - bv
		v[i] = v[i-1] + (a - b * v[i-1])*dt;   
		t[i] = t[i-1] + dt;
	}

	
	//plotting	
	FILE *gp; 
	gp = popen(GNUPLOT, "w");
	if (gp ==NULL){
		printf("Error opening Gnuplot");
		exit(0);
	}
	fprintf(gp, "set samples %d\n", N);
	fprintf(gp, "set xlabel 'time'\n");
	fprintf(gp, "set ylabel 'velocity'\n");
	
	fprintf(gp, "plot '-'\n");
//	fprintf(gp, "set style data filledcurves; rep '-' lt rgb 'red'\n");
	for (int i=0; i<N; i++){
		fprintf(gp, "%g %g\n", t[i], v[i]); 
	}
	fclose(gp);
}
