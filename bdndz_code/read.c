/* Reads the 2MASS redshift distribution */
void get_2mass(double **dist) {
  FILE *fp;
  int iz, sample;
  double j1,j2;
  double total[]={0,0,0,0};

  /* Read data */
  fp = fopen("data/2mass_b.dat", "r");
  for(iz=0; iz<40; iz++) {
    fscanf(fp, "%lg %lg %lg %lg %lg %lg", &j1, &j2, dist[0]+iz, dist[1]+iz, dist[2]+iz, dist[3]+iz);
    for(sample=0; sample<4; sample++)
      total[sample] += dist[sample][iz];
  }
  fclose(fp);

  /* Normalize */
  for(sample=0; sample<4; sample++)
    for(iz=0; iz<40; iz++)
      dist[sample][iz] /= ZSTEP*total[sample];
}

/* Reads the LRG redshift distribution */
void get_lrg(double **dist) {
  FILE *fp;
  int iz, sample;
  double j1, j2;
  double total[]={0,0};

  /* Read data */
  fp = fopen("data/lrg_lo.dat", "r");
  for(iz=0; iz<75; iz++) {
    fscanf(fp, "%lg %lg %lg", &j1, &j2, dist[0]+iz);
    total[0] += dist[0][iz];
  }
  fclose(fp);
  fp = fopen("data/lrg_hi.dat", "r");
  for(iz=0; iz<75; iz++) {
    fscanf(fp, "%lg %lg %lg", &j1, &j2, dist[1]+iz);
    total[1] += dist[1][iz];
  }
  fclose(fp);

  /* Normalize */
  for(sample=0; sample<2; sample++)
    for(iz=0; iz<75; iz++)
      dist[sample][iz] /= ZSTEP*total[sample];
}

/* Reads quasar redshift distribution */
void get_qso(double **dist) {
  int iz;
  FILE *fp;

  fp = fopen("data/dndz-qso0.dat", "r");
  for(iz=0; iz<300; iz++) {
    fscanf(fp, "%lg", dist[0]+iz);
  }
  fclose(fp);

  fp = fopen("data/dndz-qso1.dat", "r");
  for(iz=0; iz<300; iz++) {
    fscanf(fp, "%lg", dist[1]+iz);
  }
  fclose(fp);
}
