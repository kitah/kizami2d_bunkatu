#ifndef __FUNC_H__
#define __FUNC_H__

// フラグ
#define QUANT_FLG   0
#define ONOISE_FLG  0
#define LIMIT_FLG   0
#define SRAND_FLG   0
#define HEIMEN      0

// LPFなどの範囲
#define EPS    1.0E-8
#define LRG    1.0E+8
#define NHP        21//11 //9
#define NHLPF      30//25
#define NHBPF       3

// レンズ係数
#define LENS_K    1.6

// ブロックサイズBSN, BSM ブロック数NB, MB 
#define BSN        13 //31
#define BSM        13 //31
#define BSNN       15//15
#define BSMM       15//15
#define NB         24 //8
#define MB         16 //6
#define NX   (BSN*NB)
#define NY   (BSM*MB)
#define NNB        21//20
#define MMB        16//13
#define NNX  (BSNN*NNB) // 15*20
#define NNY  (BSMM*MMB) // 15*13
#define INBS        5
#define EXBS        0

#define SGMC      0.87 //0.87

// 拡大率K D0-D1の範囲に物体 観測画像NZ枚 合焦面ZAからDZずつ 
#define K           8
#define D0        0.5
#define D1        1.5 //D0
#define NZ          3
#define ZA        0.0
#define DZ        1.0

// 観測雑音のやつかな
#define ON_SGM    0.5

// d, kの推定には3回刻むのでそのためのやつ
#define DD_F    0.025
#define DD_S    0.005
#define DD      0.001
#define DDN        20

#define KDD_F    0.25
#define KDD_S    0.05
#define KD       0.01
#define KDN         5

#define QQ        100 // quant


int hsamp_gauss(int nf, double sgm, double **h);
double gauss(void);

#endif
