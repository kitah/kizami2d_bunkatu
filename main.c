#include <stdio.h>
#include <stdlib.h>
#include <imageio.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include "lpf.h"
#include "bpf.h"
#include "func.h"

/* go[NZ][NX][NY]  
 * | NHP |  NHP  | NHP |
 * | NHP | NX*NY | NHP | 
 * | NHP |  NHP  | NHP |
 *
 * g_sff[NZ][NX][NY]
 * | NX * NY |
 * 
 * y[2*NZ][BSN+4*NHP][BSM+4*NHP]
 * | NHP | NHP |   NHP   | NHP | NHP |
 * | NHP | NHP |   NHP   | NHP | NHP |
 * | NHP | NHP | BSN*BSM | NHP | NHP |
 * | NHP | NHP |   NHP   | NHP | NHP |
 * | NHP | NHP |   NHP   | NHP | NHP |
 * 
 * dy01, dy12, dy20[BSN][BSM]
 * | BSN*BSM |
 * 
 * p[zn][2*NHP+1][2*NHP+1]
 * | NHP | NHP | NHP |
 * | NHP |  2  | NHP |
 * | NHP | NHP | NHP |
 * 
 *  d
 * | NHP*K | NHLPF | NHP*K |   NHP*K   | NHP*K | NHLPF | NHP*K |
 * | NHP*K | NHLPF | NHP*K |   NHP*K   | NHP*K | NHLPF | NHP*K |
 * | NHP*K | NHLPF | NHP*K |  NX*NY*K  | NHP*K | NHLPF | NHP*K |
 * | NHP*K | NHLPF | NHP*K |   NHP*K   | NHP*K | NHLPF | NHP*K |
 * | NHP*K | NHLPF | NHP*K |   NHP*K   | NHP*K | NHLPF | NHP*K |
 * 
 * d_model, d_est, k_est
 * | NB * MB | 
 * 
 * E[NB][MB][DZ/DD]
 */


//---------------------------------------------------------------------------
int main(void){
  FILE *fp;
  IMAGE img_go;
  char flnm[256];
  int zn, n, m, k, l, bn, bm, b0[2], b1[2], bc[2];
  int nr, nw, hbsn = (BSNN-1) / 2, hbsm = (BSMM-1) / 2;
  double ***go, ***g_sff, ***y, **dy01, **dy12, **dy20, ***p, ***go_t;
  double **d, **d_model, **d_est, **k_est, **d_sff, **k_sff, ***E;
  double **dm_est, **dp_est, ***ydm, ***yd, ***ydp, **dp_model, **dm_model;
  //float data[(NHP+NX+NHP) * (NHP+NY+NHP)];
  double data[(NHP+NX+NHP) * (NHP+NY+NHP)];
  double dmin, kmin, Emin, dtest, ktest, z, sgm, d_tmp;
  double EE[3], sum_EE = 0.0, tmp, RMSE;
  time_t t1, t2;
  t1 = time(NULL);

  double alpha = 0.0, dmax[NNB], d_min[NNB];
  
#if SRAND_FLG == 1
  srandom(time(NULL));
#endif

  //---------------------------------------------------------------------------
  // 領域確保 
  // go[zn][NX][NY]
  if((go = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((go[zn] = (double **)malloc(sizeof(double *) * (NX+2*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    go[zn] += NHP;
    for(n = -NHP ; n < NX + NHP ; n++){
      if((go[zn][n] = (double *)malloc(sizeof(double ) * (NY+2*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      go[zn][n] += NHP;
    }
  }

  // go_t[zn][15*20][15*13]
  if((go_t = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((go_t[zn] = (double **)malloc(sizeof(double *) * (NNX+2*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    go_t[zn] += NHP;
    for(n = -NHP ; n < NNX + NHP ; n++){
      if((go_t[zn][n] = (double *)malloc(sizeof(double ) * (NNY+2*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      go_t[zn][n] += NHP;
    }
  }

  // g_sff[zn][NX][NY]
  if((g_sff = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((g_sff[zn] = (double **)malloc(sizeof(double *) * NX)) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    for(n = 0 ; n < NX ; n++){
      if((g_sff[zn][n] = (double *)malloc(sizeof(double ) * NY)) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
    }
  }

  // y[2*NZ][BSNN+4*NHP][BSMM+4*NHP]
  if((y = (double ***)malloc(sizeof(double **) * NZ*2)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < 2*NZ ; zn++){
    if((y[zn] = (double **)malloc(sizeof(double *) * (BSNN+4*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    y[zn] += 2*NHP;
    for(n = -2*NHP ; n < BSNN + 2*NHP ; n++){
      if((y[zn][n] = (double *)malloc(sizeof(double ) * (BSMM+4*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      y[zn][n] += 2*NHP;
    }
  }

  // ydm, yd, ydp[][][]
  if((ydm = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL ||
     (yd = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL ||
      (ydp = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < 2*NZ ; zn++){
    if((ydm[zn] = (double **)malloc(sizeof(double *) * (INBS + 4*NHP))) == NULL ||
       (yd[zn] = (double **)malloc(sizeof(double *) * (INBS + EXBS*2 + 4*NHP))) == NULL ||
       (ydp[zn] = (double **)malloc(sizeof(double *) * (INBS + 4*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);
    }
    ydm[zn] += 2*NHP;
    yd[zn] += 2*NHP;
    ydp[zn] += 2*NHP;
    for(n = -2*NHP ; n < INBS + EXBS*2 + 2*NHP ; n++){
      if((ydm[zn][n] = (double *)malloc(sizeof(double ) * (BSMM + 4*NHP))) == NULL ||
	 (yd[zn][n] = (double *)malloc(sizeof(double ) * (INBS + EXBS*2 + 4*NHP))) == NULL ||
	 (ydp[zn][n] = (double *)malloc(sizeof(double ) * (BSMM + 4*NHP))) == NULL){
	fprintf(stderr, "\nError : malloc\n\n");exit(1);
      }
      ydm[zn][n] += 2*NHP;
      yd[zn][n] += 2*NHP;
      ydp[zn][n] += 2*NHP;
    }
  }
  

  // dy01, dy12, dy20[BSNN][BSMM]
  if ((dy01  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL ||
      (dy12  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL ||
      (dy20  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL) {
        fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  dy01 += NHP;
  dy12 += NHP;
  dy20 += NHP;
  for(n = -NHP ; n < BSNN  + NHP; n++) {
    if ((dy01[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL ||
        (dy12[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL ||
        (dy20[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    dy01[n] += NHP;
    dy12[n] += NHP;
    dy20[n] += NHP;
  }

  // p[zn][2*NHP+1][2*NHP+1]
  if((p = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((p[zn] = (double **)malloc(sizeof(double *) * (2*NHP+1))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);
    }
    p[zn] += NHP;
    for(n = -NHP ; n < 1+NHP ; n++){
      if((p[zn][n] = (double *)malloc(sizeof(double ) * (2*NHP+1))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      p[zn][n] += NHP;
    }
  }

  // d
  if ( (d = (double **)malloc(sizeof(double *) * (NNX*K + 2*NHLPF + 4*NHP*K)) ) == NULL) {
    fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  d += NHLPF + 2*NHP*K;
  for(n = -(NHLPF + 2*NHP*K) ; n < NNX*K + NHLPF + 2*NHP*K ; n++) {
    if ( (d[n] = (double *)malloc(sizeof(double) * (NNY*K + 2*NHLPF + 4*NHP*K)) ) == NULL) {
      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    d[n] += NHLPF + 2*NHP*K;
  }

  // d_model, d_est, k_est[NB][MB], dm,dp_est
  if ((d_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dm_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dp_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (d_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dm_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dp_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (k_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL) {
        fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  for(n = 0 ; n < NNB ; n++) {
    if ((d_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dm_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dp_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
        (d_est[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dm_est[n]  = (double *)malloc(sizeof(double ) * MMB) ) == NULL ||
	(dp_est[n]  = (double *)malloc(sizeof(double ) * MMB) ) == NULL ||
        (k_est[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
  }

  // d_sff, k_sff[NB][MB]
  if ((d_sff  = (double **)malloc(sizeof(double *) * NB) ) == NULL ||
      (k_sff  = (double **)malloc(sizeof(double *) * NB) ) == NULL) {
       fprintf(stderr, "\nError : malloc\n\n"); exit(1); 
  }
  for(n = 0 ; n < NB ; n++) {
    if ((d_sff[n]  = (double *)malloc(sizeof(double) * MB) ) == NULL ||
        (k_sff[n]  = (double *)malloc(sizeof(double) * MB) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
  }

  // E
  if ( (E = (double ***)malloc(sizeof(double **) * NNB) ) == NULL) {
    fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  for(n = 0 ; n < NNB ; n++) {
    if ( (E[n] = (double **)malloc(sizeof(double *) * MMB) ) == NULL) {
      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    for(m = 0 ; m < MMB ; m++){    
      if ( (E[n][m] = (double *)malloc(sizeof(double) * (int)(DZ / DD)) ) == NULL) {
	      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
      }
    }
  }

  //---------------------------------------------------------------------------
  // d[]に対象物体のdepth形状を設定
#if HEIMEN == 0
  for(n = -(2*NHP*K + NHLPF) ; n < (NNX*K + 2*NHP*K + NHLPF) ; n++){
    for(m = -(2*NHP*K + NHLPF) ; m < (NNY*K + 2*NHP*K + NHLPF) ; m++){    
      d[n][m] = D0 + (D1-D0) / (double)(NX*K - 1) * n;
    }
  }
#endif

#if HEIMEN == 1
  for(n = -(2*NHP*K + NHLPF) ; n < (NNX*K + 2*NHP*K + NHLPF) ; n++){
    for(m = -(2*NHP*K + NHLPF) ; m < (NNY*K + 2*NHP*K + NHLPF) ; m++){    
      d[n][m] = 0.9;
    }
  }
#endif
  
  sprintf(flnm, "./d.dat");
  fp = fopen(flnm, "w");
  //for(n = -(2*NHP*K + NHLPF) ; n < (NX*K + 2*NHP*K + NHLPF) ; n++){
  //for(m = -(2*NHP*K + NHLPF) ; m < (NY*K + 2*NHP*K + NHLPF) ; m++){}
  for(n = 0 ; n < NNX ; n++){
    for(m = 0 ; m < NNY ; m++){
      fprintf(fp, "%d %d %f\n", n, m, d[n*K][m*K]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  for(n = 0 ; n < NNB ; n++){
    dmax[n] = d[n*BSNN*K + BSNN*K-K][0];
  }
  sprintf(flnm, "./dmax.dat");
  fp = fopen(flnm, "w");
  for(n = 0 ; n < NNB ; n++){
      fprintf(fp, "%d %f\n", n, dmax[n]);
  }
  fclose(fp);

  for(n = 0 ; n < NNB ; n++){
    d_min[n] = d[n*BSNN*K][0];
  }
  sprintf(flnm, "./d_min.dat");
  fp = fopen(flnm, "w");
  for(n = 0 ; n < NNB ; n++){
      fprintf(fp, "%d %f\n", n, d_min[n]);
  }
  fclose(fp);
  
  
  //---------------------------------------------------------------------------  
  // freedする data -> goへ
#if HEIMEN == 0
  for(zn = 0 ; zn < NZ ; zn++){
    //sprintf(flnm, "./zn3_0515/d0515_bi_g%d", zn);
    sprintf(flnm, "./../obs_img/k16/d0515/d0515_ra_0%d.bin", zn);
    fp = fopen(flnm, "r");
    //nr = fread(data, sizeof(float), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    nr = fread(data, sizeof(double), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    printf("nr = %f\n", nr);
    nw = (NHP+NX+NHP) * (NHP+NY+NHP);
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }

    for(m = -NHP ; m < NY + NHP ; m++){
      for(n = -NHP ; n < NX + NHP ; n++){
        go[zn][n][m] = (double)data[(m+NHP) * (NHP+NX+NHP) + (n+NHP)];
      }
    }
  }
#endif

  // 平面
#if HEIMEN == 1
  for(zn = 0 ; zn < NZ ; zn++){
    sprintf(flnm, "./d09/d09_bi_g%d", zn);
    fp = fopen(flnm, "r");
    nr = fread(data, sizeof(float), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    nw = (NHP+NX+NHP) * (NHP+NY+NHP);
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }
    for(m = -NHP ; m < NY + NHP ; m++){
      for(n = -NHP ; n < NX + NHP ; n++){
	go[zn][n][m] = (double)data[(m+NHP) * (NHP+NX+NHP) + (n+NHP)];
      }
    }
  }
#endif

  //---------------------------------------------------------------------------
  // go[zn][]に観測雑音を付与、量子化してgo[zn][]に
  for(zn = 0 ; zn < NZ ; zn++){
    for(n = -NHP ; n < NX + NHP ; n++){
      for(m = -NHP ; m < NY + NHP ; m++){
        
#if ONOISE_FLG == 1
        go[zn][n][m] = go[zn][n][m] + ON_SGM * gauss();
#endif

#if QUANT_FLG == 1
        go[zn][n][m] = floor(go[zn][n][m] + 0.5);
#endif

#if LIMIT_FLG == 1
        if (go[zn][n][m] > 255.0) {
          go[zn][n][m] = 255.0;
        } else if (go[zn][n][m] < 0.0) {
          go[zn][n][m] = 0.0;
        } else {
          go[zn][n][m] = floor(go[zn][n][m] + 0.5);
        }
#endif
      }
    }

    //---------------------------------------------------------------------------
    // go[zn][]の保存
    sprintf(flnm, "./model/go%d.dat", zn);
    cIMAGE(NX + 2*NHP, NY + 2*NHP, &img_go, MONO);
    fp = fopen(flnm, "w");
    for(n = -NHP ; n < NX + NHP ; n++){
      for(m = -NHP ; m < NY + NHP ; m++){
        fprintf(fp, "%4d %4d  %f\n", n, m, go[zn][n][m]);
        D(img_go, n+NHP, m+NHP) = (unsigned char)(go[zn][n][m] + 0.5);
      }
    }
    fclose(fp);
    sprintf(flnm, "./model/go%d.bmp", zn);
    wIMAGE(img_go, flnm);
  }//zn
  
  for(zn = 0 ; zn < NZ ; zn++){
    sprintf(flnm, "./model/bs15_go%d.dat", zn);
    //cIMAGE(NNX + 2*NHP, NNY + 2*NHP, &img_go, MONO);
    fp = fopen(flnm, "w");
    for(n = -NHP ; n < NNX + NHP; n++){
      for(m = -NHP ; m < NNY + NHP; m++){
	go_t[zn][n][m] = go[zn][n][m];
      }
    }
    /*
    for(n = 0 ; n < NHP ; n++){
      for(m = 0 ; m < NHP ; m++){
	go_t[zn][n + NNX][m + NNY] = go[zn][n + NX][m + NY];
      }
    }
    
    for(n = 0 ; n < NHP ; n++){
      for(m = -NHP ; m < NNY ; m++){
	go_t[zn][n + NNX][m] = go[zn][n + NX][m];
      }
    }

    for(n = -NHP ; n < NNX ; n++){
      for(m = 0 ; m < NHP ; m++){
	go_t[zn][n][m + NNY] = go[zn][n][m + NY];
      }
      }*/
    
    for(n = -NHP ; n < NNX + NHP ; n++){
      for(m = -NHP ; m < NNY + NHP ; m++){
	fprintf(fp, "%4d %4d  %f\n", n, m, go_t[zn][n][m]);
      }
    }
    fclose(fp);
  }
  
  //---------------------------------------------------------------------------
  // go[zn][]に逆のボケを付与してy[]へ、その後誤差測定
  // 観測画像のブロック1つずつy[]を生成してる
  printf("---start---\n");
  for(bn = 0 ; bn < NNB ; bn++){ // 20 * 13
    printf("bn = %d\n", bn);
    b0[0] = bn * BSNN;       //  0, 13, 26...
    b1[0] = b0[0] + BSNN;    // 13, 26, 39...
    bc[0] = b0[0] + hbsn;
    for(bm = 0 ; bm < MMB ; bm++){
      printf("  bm = %d\n", bm);
      b0[1] = bm * BSMM;     //  0, 13, 26...
      b1[1] = b0[1] + BSMM;  // 13, 26, 39...
      bc[1] = b0[1] + hbsm;

      dmin = 0.0;
      kmin = 0.0;
      Emin = LRG;
      sprintf(flnm, "./Edk/Edk%02d_%02d.dat", bn, bm);
      fp = fopen(flnm, "w");

      // d,kに対するEの最小値探索
      // DD ni sita
      for(dtest = 1.0 - DD_F*DDN ; dtest <= 1.0 + DD_F*DDN ; dtest+=DD_F){
          // p0, p1, p2を生成
          for(zn = 0 ; zn < NZ ; zn++){
            z = ZA + DZ * (double)zn;
            sgm = LENS_K * fabs(z - dtest);
            hsamp_gauss(2*NHP + 1, sgm, p[zn]);
          }
	  
	  // y[0]~y[5]をクリア
          /*for(n = -2*NHP ; n < BSNN + 2*NHP ; n++){
            for(m = -2*NHP ; m < BSMM + 2*NHP ; m++){
	      y[0][n][m] = 0.0;
              y[1][n][m] = 0.0;
              y[2][n][m] = 0.0;
              y[3][n][m] = 0.0;
	      y[4][n][m] = 0.0;
              y[5][n][m] = 0.0;
	    }
          }*/
	  
	  for(n = -2*NHP ; n < INBS + EXBS*2 + 2*NHP ; n++){
	    for(m = -2*NHP ; m < INBS + EXBS*2 + 2*NHP ; m++){
	      yd[0][n][m] = 0.0;
	      yd[1][n][m] = 0.0;
	      yd[2][n][m] = 0.0;
	      yd[3][n][m] = 0.0;
	      yd[4][n][m] = 0.0;
	      yd[5][n][m] = 0.0;
	    }
	  }
	  
          // y[0]~y[5]を生成
          /*for(n = -NHP ; n < BSNN + NHP ; n++){
            for(m = -NHP ; m < BSMM + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
                  // 0-1
                  y[0][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[1][k][l];
                  y[1][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[0][k][l];
                  // 1-2
                  y[2][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[3][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[1][k][l];
                  // 2-0
                  y[4][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[5][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[0][k][l];
                }
              }
            }
	    }*/
	  for(n = -NHP ; n < INBS + EXBS*2 + NHP ; n++){
            for(m = -NHP ; m < INBS + EXBS + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
                  // 0-1
                  yd[0][n + k][m + l] += go_t[0][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[1][k][l];
                  yd[1][n + k][m + l] += go_t[1][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[0][k][l];
                  // 1-2
                  yd[2][n + k][m + l] += go_t[1][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[2][k][l];
                  yd[3][n + k][m + l] += go_t[2][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[1][k][l];
                  // 2-0
                  yd[4][n + k][m + l] += go_t[0][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[2][k][l];
                  yd[5][n + k][m + l] += go_t[2][n + b0[0] + INBS - EXBS][m + b0[1] + INBS] * p[0][k][l];
                }
              }
            }
	  }

	  // EXBSを用いて差分を取る際のブロックを可変に
          // 差分をとる
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	  for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              dy01[n][m] = y[0][n][m] - y[1][n][m];
              dy12[n][m] = y[2][n][m] - y[3][n][m];
              dy20[n][m] = y[4][n][m] - y[5][n][m];
            }
          }
  
          // 誤差EEの評価
          EE[0] = EE[1] = EE[2] = 0.0;
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	  for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              EE[0] += dy01[n][m] * dy01[n][m];
              EE[1] += dy12[n][m] * dy12[n][m];
              EE[2] += dy20[n][m] * dy20[n][m];              
            }
	  }
	  
	  /*
	  // 差分をとる
          for(n = 0 ; n < INBS ; n++){
            for(m = 0 ; m < INBS ; m++){
              dy01[n][m] = yd[0][n][m] - yd[1][n][m];
              dy12[n][m] = yd[2][n][m] - yd[3][n][m];
              dy20[n][m] = yd[4][n][m] - yd[5][n][m];
            }
          }
  
          // 誤差EEの評価
          EE[0] = EE[1] = EE[2] = 0.0;
          for(n = 0 ; n < INBS ; n++){
            for(m = 0 ; m < INBS ; m++){
              EE[0] += dy01[n][m] * dy01[n][m];
              EE[1] += dy12[n][m] * dy12[n][m];
              EE[2] += dy20[n][m] * dy20[n][m];              
            }
	    }*/
	    
          // 誤差EEを最小にするd, k
	  sum_EE = 0.0;
          sum_EE = EE[0] + EE[1] + EE[2];
	  //printf("   sumEE = %f\n", sum_EE);
          if(Emin > sum_EE){
            Emin = sum_EE;
            dmin = dtest;
	    kmin = LENS_K;
	  }
           fprintf(fp, "%f  %f  %f\n", LENS_K, dtest, sum_EE);
	   fprintf(fp, "\n\n");
      }//d

      
      // 2回目
      for(dtest = dmin - DD_F, Emin = LRG ; dtest <= dmin + DD_F ; dtest+=DD_S){
          //p0, p1, p2を生成
          for(zn = 0 ; zn < NZ ; zn++){
            z = ZA + DZ * (double)zn;
            sgm = LENS_K * fabs(z - dtest);
            hsamp_gauss(2 * NHP + 1, sgm, p[zn]);
          }

          // y[0]~y[5]の初期化
          for(n = -2*NHP ; n < BSNN + 2*NHP ; n++){
            for(m = -2*NHP ; m < BSMM + 2*NHP ; m++){
              y[0][n][m] = 0.0;
              y[1][n][m] = 0.0;
              y[2][n][m] = 0.0;
              y[3][n][m] = 0.0;
              y[4][n][m] = 0.0;
              y[5][n][m] = 0.0;
            }
          }

          // y[0]~y[5]の生成
          for(n = -NHP ; n < BSNN + NHP ; n++){
            for(m = -NHP ; m < BSMM + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
		   // 0-1
                  y[0][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[1][k][l];
                  y[1][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[0][k][l];
                  // 1-2
                  y[2][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[3][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[1][k][l];
                  // 2-0
                  y[4][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[5][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[0][k][l];
		}
              }
            }
          }

    	  // EXBSを用いて差分を取る際のブロックを可変に
          // 差分をとる
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	    for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              dy01[n][m] = y[0][n][m] - y[1][n][m];
              dy12[n][m] = y[2][n][m] - y[3][n][m];
              dy20[n][m] = y[4][n][m] - y[5][n][m];
            }
          }
  
          // 誤差EEの評価
          EE[0] = EE[1] = EE[2] = 0.0;
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	    for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              EE[0] += dy01[n][m] * dy01[n][m];
              EE[1] += dy12[n][m] * dy12[n][m];
              EE[2] += dy20[n][m] * dy20[n][m];              
            }
	  }
	  
          sum_EE = 0.0;
          sum_EE = EE[0] + EE[1] + EE[2];

          // 誤差EEを最小にするd, kの最小値を探索
          if(Emin > sum_EE){
	    Emin = sum_EE;
            dmin = dtest;
            kmin = LENS_K;
          }
          fprintf(fp, "%f  %f  %f\n", LENS_K, dtest, sum_EE);          
	  fprintf(fp, "\n\n");
      }//d

      // 3回目
      for(dtest = dmin - DD_S, Emin = LRG ; dtest <= dmin + DD_S ; dtest+=DD){
          // p0, p1, p2
          for(zn = 0 ; zn < NZ ; zn++){
            z = ZA + DZ * (double)zn;
            sgm = LENS_K * fabs(z - dtest);
            hsamp_gauss(2*NHP + 1, sgm, p[zn]);
          }

          // y[0]~y[5]の初期化
          for(n = -2*NHP ; n < BSNN + 2*NHP ; n++){
            for(m = -2*NHP ; m < BSMM + 2*NHP ; m++){
              y[0][n][m] = 0.0;
              y[1][n][m] = 0.0;
              y[2][n][m] = 0.0;
              y[3][n][m] = 0.0;
              y[4][n][m] = 0.0;
              y[5][n][m] = 0.0;
            }
          }

          // y[0]~y[5]の生成
          for(n = -NHP ; n < BSNN + NHP ; n++){
            for(m = -NHP ; m < BSMM + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
		  // 0-1
                  y[0][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[1][k][l];
                  y[1][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[0][k][l];
                  // 1-2
                  y[2][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[3][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[1][k][l];
                  // 2-0
                  y[4][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[2][k][l];
                  y[5][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[0][k][l];
                }
              }
            }
          }

	  // EXBSを用いて差分を取る際のブロックを可変に
          // 差分をとる
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	    for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              dy01[n][m] = y[0][n][m] - y[1][n][m];
              dy12[n][m] = y[2][n][m] - y[3][n][m];
              dy20[n][m] = y[4][n][m] - y[5][n][m];
            }
          }
  
          // 誤差EEの評価
          EE[0] = EE[1] = EE[2] = 0.0;
          //for(n = 0 ; n < BSNN ; n++){
	  //for(m = 0 ; m < BSMM ; m++){
	  for(n = INBS - EXBS ; n < INBS*2 + EXBS ; n++){
	    for(m = INBS - EXBS ; m < INBS*2 + EXBS ; m++){
              EE[0] += dy01[n][m] * dy01[n][m];
              EE[1] += dy12[n][m] * dy12[n][m];
              EE[2] += dy20[n][m] * dy20[n][m];              
            }
	  }
	  
          sum_EE = 0.0;
          sum_EE = EE[0] + EE[1] + EE[2];

          // 誤差EEを最小にするd, kの最小値を探索
          if(Emin > sum_EE){
            Emin = sum_EE;
            dmin = dtest;
            kmin = LENS_K;
	  }
          fprintf(fp, "%f  %f  %f\n", LENS_K, dtest, sum_EE); 
	  fprintf(fp, "\n\n");
      }//d
      
      fclose(fp);
      //printf("dmin = %f\n", dmin);
      d_est[bn][bm] = dmin;
      printf("    d_est = %f\n", d_est[bn][bm]);
      k_est[bn][bm] = kmin;
    }
  }
  printf("\n");

  //-----------------------------------------------------------------------------
  // d0を求めたので、dpを求める
  printf("---dp start---\n");
  for(bn = 0 ; bn < NNB ; bn++){ // 20 * 13
    printf("bn = %d\n", bn);
    b0[0] = bn * BSNN;       //  0, 13, 26...
    b1[0] = b0[0] + BSNN;    // 13, 26, 39...
    bc[0] = b0[0] + hbsn;
    for(bm = 0 ; bm < MMB ; bm++){
      printf("  bm = %d\n", bm);
      b0[1] = bm * BSMM;     //  0, 13, 26...
      b1[1] = b0[1] + BSMM;  // 13, 26, 39...
      bc[1] = b0[1] + hbsm;
      
      alpha = (dmax[bn] - d_est[bn][bm])/hbsn;
      //printf(" alpha = %f\n", alpha);
      dmin = 0.0;
      kmin = 0.0;
      Emin = LRG;
      d_tmp = 0.0;
      sprintf(flnm, "./Edk/Edpk%02d_%02d.dat", bn, bm);
      fp = fopen(flnm, "w");
      d_tmp = d_est[bn][bm] + INBS*alpha;
 
      // d,kに対するEの最小値探索
      // DD ni sita
      for(dtest =  d_tmp - 0.01 ; dtest <= d_tmp + 0.01 ; dtest+=DD){
          // p0, p1, p2を生成
	for(zn = 0 ; zn < NZ ; zn++){
	  z = ZA + DZ * (double)zn;
	  sgm = LENS_K * fabs(z - dtest);
	  hsamp_gauss(2*NHP + 1, sgm, p[zn]);
	}

	// y[0]~y[5]をクリア
	for(n = -2*NHP ; n < INBS + 2*NHP ; n++){
	  //for(m = -2*NHP ; m < INBS + 2*NHP ; m++){
	    for(m = -2*NHP ; m < BSMM + 2*NHP ; m++){
	    ydp[0][n][m] = 0.0;
	    ydp[1][n][m] = 0.0;
	    ydp[2][n][m] = 0.0;
	    ydp[3][n][m] = 0.0;
	    ydp[4][n][m] = 0.0;
	    ydp[5][n][m] = 0.0;
	  }
	}

	// y[0]~y[5]を生成
	/*for(n = -NHP ; n < INBS + NHP ; n++){
	  for(m = -NHP ; m < INBS + NHP ; m++){
	    for(k = -NHP ; k <= NHP ; k++){
	      for(l = -NHP ; l <= NHP ; l++){
		// 0-1
		ydp[0][n + k][m + l] += go_t[0][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[1][k][l];
		ydp[1][n + k][m + l] += go_t[1][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[0][k][l];
		// 1-2
		ydp[2][n + k][m + l] += go_t[1][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[2][k][l];
		ydp[3][n + k][m + l] += go_t[2][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[1][k][l];
		// 2-0
		ydp[4][n + k][m + l] += go_t[0][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[2][k][l];
		ydp[5][n + k][m + l] += go_t[2][n + b0[0] + INBS*2][m + b0[1] + INBS] * p[0][k][l];
	      }
	    }
	  }
	  }*/
	
	for(n = -NHP ; n < INBS + NHP ; n++){
	  for(m = -NHP ; m < BSMM + NHP ; m++){
	    for(k = -NHP ; k <= NHP ; k++){
	      for(l = -NHP ; l <= NHP ; l++){
		// 0-1
		ydp[0][n + k][m + l] += go_t[0][n + b0[0] + INBS*2][m + b0[1]] * p[1][k][l];
		ydp[1][n + k][m + l] += go_t[1][n + b0[0] + INBS*2][m + b0[1]] * p[0][k][l];
		// 1-2
		ydp[2][n + k][m + l] += go_t[1][n + b0[0] + INBS*2][m + b0[1]] * p[2][k][l];
		ydp[3][n + k][m + l] += go_t[2][n + b0[0] + INBS*2][m + b0[1]] * p[1][k][l];
		// 2-0
		ydp[4][n + k][m + l] += go_t[0][n + b0[0] + INBS*2][m + b0[1]] * p[2][k][l];
		ydp[5][n + k][m + l] += go_t[2][n + b0[0] + INBS*2][m + b0[1]] * p[0][k][l];
	      }
	    }
	  }
	}
	
	// 差分をとる
	for(n = 0 ; n < INBS ; n++){
	  for(m = 0 ; m < BSMM ; m++){
	    //for(m = 0 ; m < INBS ; m++){
	    dy01[n][m] = ydp[0][n][m] - ydp[1][n][m];
	    dy12[n][m] = ydp[2][n][m] - ydp[3][n][m];
	    dy20[n][m] = ydp[4][n][m] - ydp[5][n][m];
	  }
	}
	
	// 誤差EEの評価
	EE[0] = EE[1] = EE[2] = 0.0;
	for(n = 0 ; n < INBS ; n++){
	  for(m = 0 ; m < BSMM ; m++){
	    //for(m = 0 ; m < INBS ; m++){
	    EE[0] += dy01[n][m] * dy01[n][m];
	    EE[1] += dy12[n][m] * dy12[n][m];
	    EE[2] += dy20[n][m] * dy20[n][m];              
	  }
	}
	  
	// 誤差EEを最小にするd, k
	sum_EE = EE[0] + EE[1] + EE[2];
	//printf("   sumEE = %f\n", sum_EE);
	if(Emin > sum_EE){
	  Emin = sum_EE;
	  dmin = dtest;
	  kmin = LENS_K;
	}
	fprintf(fp, "%f  %f  %f\n", LENS_K, dtest, sum_EE);
	fprintf(fp, "\n\n");
      }//d
      //dp_est[bn][bm] = d_est[bn][bm] + 5*alpha;
      dp_est[bn][bm]  = dmin;
      printf("    dp_est = %f\n", dp_est[bn][bm]);
    }
  }
  printf("\n");

  //-----------------------------------------------------------------------------
  // d0を求めたので、dmを求める
  printf("---dm start---\n");
  for(bn = 0 ; bn < NNB ; bn++){ // 20 * 13
    printf("bn = %d\n", bn);
    b0[0] = bn * BSNN;       //  0, 13, 26...
    b1[0] = b0[0] + BSNN;    // 13, 26, 39...
    bc[0] = b0[0] + hbsn;
    for(bm = 0 ; bm < MMB ; bm++){
      printf("  bm = %d\n", bm);
      b0[1] = bm * BSMM;     //  0, 13, 26...
      b1[1] = b0[1] + BSMM;  // 13, 26, 39...
      bc[1] = b0[1] + hbsm;

      alpha = (dmax[bn] - d_est[bn][bm])/hbsn;
      //alpha = (d_est[bn][bm] - d_min[bn])/7;
      //printf(" alpha = %f\n", alpha);
      dmin = 0.0;
      kmin = 0.0;
      d_tmp = 0.0;
      Emin = LRG;
      sprintf(flnm, "./Edk/Edmk%02d_%02d.dat", bn, bm);
      fp = fopen(flnm, "w");

      d_tmp = d_est[bn][bm] - INBS*alpha;
      printf(" d_tmp = %f\n", d_tmp);
      
      // d,kに対するEの最小値探索
      // DD ni sita
      for(dtest = d_tmp - 0.01 ; dtest <= d_tmp + 0.01 ; dtest+=DD){
          // p0, p1, p2を生成
      for(zn = 0 ; zn < NZ ; zn++){
	z = ZA + DZ * (double)zn;
	sgm = LENS_K * fabs(z - dtest);
	hsamp_gauss(2*NHP + 1, sgm, p[zn]);
      }

	  // y[0]~y[5]をクリア
          for(n = -2*NHP ; n < INBS + 2*NHP ; n++){
	    //for(m = -2*NHP ; m < INBS + 2*NHP ; m++){
	      for(m = -2*NHP ; m < BSMM + 2*NHP ; m++){
	      ydm[0][n][m] = 0.0;
	      ydm[1][n][m] = 0.0;
	      ydm[2][n][m] = 0.0;
	      ydm[3][n][m] = 0.0;
	      ydm[4][n][m] = 0.0;
	      ydm[5][n][m] = 0.0;
	    }
	  }
	  
          // y[0]~y[5]を生成
          /*for(n = -NHP ; n < INBS + NHP ; n++){
            for(m = -NHP ; m < INBS + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
		  // 0-1
		  ydm[0][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1] + INBS] * p[1][k][l];
		  ydm[1][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1] + INBS] * p[0][k][l];
		  // 1-2
		  ydm[2][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1] + INBS] * p[2][k][l];
		  ydm[3][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1] + INBS] * p[1][k][l];
		  // 2-0
		  ydm[4][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1] + INBS] * p[2][k][l];
		  ydm[5][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1] + INBS] * p[0][k][l];
                }
              }
            }
	    }*/
	  for(n = -NHP ; n < INBS + NHP ; n++){
	    for(m = -NHP ; m < BSMM + NHP ; m++){
              for(k = -NHP ; k <= NHP ; k++){
                for(l = -NHP ; l <= NHP ; l++){
                  // 0-1
                  ydm[0][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[1][k][l];
                  ydm[1][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[0][k][l];
                  // 1-2
                  ydm[2][n + k][m + l] += go_t[1][n + b0[0]][m + b0[1]] * p[2][k][l];
                  ydm[3][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[1][k][l];
                  // 2-0
                  ydm[4][n + k][m + l] += go_t[0][n + b0[0]][m + b0[1]] * p[2][k][l];
                  ydm[5][n + k][m + l] += go_t[2][n + b0[0]][m + b0[1]] * p[0][k][l];
                }
              }
            }
	  }
  
	  // 差分をとる
          for(n = 0 ; n < INBS ; n++){
            //for(m = 0 ; m < INBS ; m++){
	      for(m = 0 ; m < BSMM ; m++){
              dy01[n][m] = ydm[0][n][m] - ydm[1][n][m];
              dy12[n][m] = ydm[2][n][m] - ydm[3][n][m];
              dy20[n][m] = ydm[4][n][m] - ydm[5][n][m];
            }
          }
  
          // 誤差EEの評価
          EE[0] = EE[1] = EE[2] = 0.0;
          for(n = 0 ; n < INBS ; n++){
	    for(m = 0 ; m < BSMM ; m++){
	      //for(m = 0 ; m < INBS ; m++){
              EE[0] += dy01[n][m] * dy01[n][m];
              EE[1] += dy12[n][m] * dy12[n][m];
              EE[2] += dy20[n][m] * dy20[n][m];              
            }
	  }
	  
          // 誤差EEを最小にするd, k
          sum_EE = EE[0] + EE[1] + EE[2];
	  //printf("   sumEE = %f\n", sum_EE);
          if(Emin > sum_EE){
            Emin = sum_EE;
            dmin = dtest;
	    kmin = LENS_K;
	  }
           fprintf(fp, "%f  %f  %f\n", LENS_K, dtest, sum_EE);
	   fprintf(fp, "\n\n");
      }//d
      // dm_est[bn][bm] = d_est[bn][bm] - 5*alpha;
      dm_est[bn][bm] = dmin;
       printf("    dm_est = %f\n", dm_est[bn][bm]);
    }
  }
  printf("\n");

  //---------------------------------------------------------------------------  
  // ブロック単位のモデル形状
  for(bn = 0 ; bn < NNB ; bn++){
    b0[0] = bn * BSNN;    // 0, 15, 30, ...
    b1[0] = b0[0] + BSNN; // 15, 30, ...
    for(bm = 0 ; bm < MMB ; bm++){
      b0[1] = bm * BSMM;
      b1[1] = b0[1] + BSMM;

      tmp = 0.0;
      for(n = b0[0] ; n < b1[0] ; n++){
        for(m = b0[1] ; m < b1[1] ; m++){
          tmp += d[n * K][m * K];
        }
      }

      d_model[bn][bm] = tmp / (double)(BSNN * BSMM);
    }
  }
  sprintf(flnm, "./d_model.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*BSNN + hbsn, bm*BSMM + hbsm, d_model[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  // dp_model
  for(n = 0 ; n < NNB ; n++){
    for(m = 0 ; m < MMB ; m++){
      dp_model[n][m] = d[n*BSNN*K + ((BSNN-1)/2 + INBS)*K][0];
    }
  }
  sprintf(flnm, "./dp_model.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*BSNN + hbsn + INBS, bm*BSMM + hbsm, dp_model[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
   // dm_model
  for(n = 0 ; n < NNB ; n++){
    for(m = 0 ; m < MMB ; m++){
      dm_model[n][m] = d[n*BSNN*K + ((BSNN-1)/2 - INBS)*K][0];
    }
  }
  sprintf(flnm, "./dm_model.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*BSNN + hbsn - INBS, bm*BSMM + hbsm, dm_model[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  //---------------------------------------------------------------------------  
  // 結果保存
  sprintf(flnm, "./dk_est.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f  %f  %f\n",
              bn*BSNN + hbsn, bm*BSMM + hbsm, d_est[bn][bm], d_model[bn][bm], k_est[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sprintf(flnm, "./d_est.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n",
              bn*BSNN + hbsn, bm*BSMM + hbsm, d_est[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  // dp_est
  sprintf(flnm, "./dp_est.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n",
              bn*BSNN + hbsn + INBS, bm*BSMM + hbsm, dp_est[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  // dm_est
  sprintf(flnm, "./dm_est.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n",
              bn*BSNN + hbsn - INBS, bm*BSMM + hbsm, dm_est[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  //---------------------------------------------------------------------------  
  // 誤差保存
  sprintf(flnm, "./e.dat");
  fp = fopen(flnm, "w");
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      fprintf(fp, "%4d %4d  %f\n",
              bn*BSNN + hbsn, bm*BSMM + hbsm, d_est[bn][bm] - d_model[bn][bm]);
    }
    fprintf(fp, "\n");
  }

  //---------------------------------------------------------------------------  
  // RMSE
  printf("------------\n");
  printf("RMSE output\n");
  printf("------------\n");
  // dm
  RMSE = 0.0;
  //tmp = 0.0;
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      //printf("dm_model[%d][%d] = %f\n", bn, bm, dm_model[bn][bm]);
      tmp = dm_model[bn][bm] - dm_est[bn][bm];
      RMSE += tmp * tmp;
    }
    //printf("\n");
  }
  RMSE = sqrt(RMSE / (double)(NNB*MMB/3));
  printf("\n");
  printf("RMSE(dm) = %f\n\n", RMSE);
  
  // d0
  RMSE = 0.0;
  //tmp = 0.0;
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      //printf("d_model[%d][%d] = %f\n", bn, bm, d_model[bn][bm]);
      tmp = d_model[bn][bm] - d_est[bn][bm];
      RMSE += tmp * tmp;
    }
    //printf("\n");
  }
  RMSE = sqrt(RMSE / (double)(NNB*MMB));
  printf("\n");
  printf("RMSE(d0) = %f\n\n", RMSE);

   // dp
  RMSE = 0.0;
  //tmp = 0.0;
  for(bn = 0 ; bn < NNB ; bn++){
    for(bm = 0 ; bm < MMB ; bm++){
      //printf("dp_model[%d][%d] = %f\n", bn, bm, dp_model[bn][bm]);
      tmp = dp_model[bn][bm] - dp_est[bn][bm];
      RMSE += tmp * tmp;
    }
    //printf("\n");
  }
  //  RMSE = sqrt(RMSE / (double)(NNB*MMB));
  RMSE = sqrt(RMSE / (double)(NNB*MMB/3));
  printf("\n");
  printf("RMSE(dp) = %f\n\n", RMSE);
  
  t2 = time(NULL);
  printf("\n");
  printf("[time] : %d[minutes]\n\n", (int)((t2 - t1)/60));
}//main
