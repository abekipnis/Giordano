void sv_poincare(float x[N], float y[N], float t[N], char* fname, float r);
char* simulate(float a, float r, float vx0, float vy0, float x0, float y0, int n, int nit);
void initialize(float x[N], float y[N], float t[N], float x0, float y0);
void Eint(float x[N], float y[N], float t[N], float a, float vx0, float vy0, float r);

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
