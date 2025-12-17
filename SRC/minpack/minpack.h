/* SPDX-License-Identifier: BSD-3-Clause */
/*
 * Copyright (c) 1980-1999 The University of Chicago, as Operator of
 * Argonne National Laboratory. All rights reserved.
 * Copyright (c) 2024 SciPy developers.
 *
 * This file is a C translation of the Fortran MINPACK-1 library written by
 * Jorge J. Moré, Burton S. Garbow, and Kenneth E. Hillstrom.
 *
 * ---------- OpenSees integration note ----------
 * Adapted in 2025 for OpenSees 
 * No external dependencies or Fortran runtime are introduced.
 */
// minpack_core.h
// 纯净的MINPACK核心算法C语言接口，供C++项目使用
// 此文件基于提供的 __minpack.c 源码生成，签名准确
#ifndef MINPACK_CORE_H
#define MINPACK_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

/* 1. CHKDER - 检查导数的正确性 */
void CHKDER(const int m, const int n, double* x, const double* fvec, const double* fjac,
            const int ldfjac, double* xp, const double* fvecp, const int mode, double* err);

/* 2. HYBRD - 使用混合方法求解非线性方程组（无需提供雅可比矩阵）*/
void HYBRD(int(*fcn)(int* n, double* x, double* fvec, int* iflag),
           const int n,
           double* x,
           double* fvec,
           const double xtol,
           const int maxfev,
           const int ml,
           const int mu,
           const double epsfcn,
           double* diag,
           const int mode,
           const double factor,
           const int nprint,
           int* info,
           int* nfev,
           double* fjac,
           const int ldfjac,
           double* r,
           const int lr,
           double* qtf,
           double* wa1,
           double* wa2,
           double* wa3,
           double* wa4);

/* 3. HYBRJ - 使用混合方法求解非线性方程组（需提供雅可比矩阵）*/
void HYBRJ(int(*fcn)(int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag),
           const int n,
           double* x,
           double* fvec,
           double* fjac,
           const int ldfjac,
           const double xtol,
           const int maxfev,
           double* diag,
           const int mode,
           const double factor,
           const int nprint,
           int* info,
           int* nfev,
           int* njev,
           double* r,
           const int ldr,
           double* qtf,
           double* wa1,
           double* wa2,
           double* wa3,
           double* wa4);

/* 4. LMDIF - 使用Levenberg-Marquardt算法求解最小二乘问题（无需提供雅可比矩阵）*/
void LMDIF(int(*fcn)(int* m, int* n, double* x, double* fvec, int* iflag),
           const int m,
           const int n,
           double* x,
           double* fvec,
           const double ftol,
           const double xtol,
           const double gtol,
           const int maxfev,
           const double epsfcn,
           double* diag,
           const int mode,
           const double factor,
           const int nprint,
           int* info,
           int* nfev,
           double* fjac,
           const int ldfjac,
           int* ipvt,
           double* qtf,
           double* wa1,
           double* wa2,
           double* wa3,
           double* wa4);

/* 5. LMDER - 使用Levenberg-Marquardt算法求解最小二乘问题（需提供雅可比矩阵）*/
void LMDER(int(*fcn)(int* m, int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag),
           const int m,
           const int n,
           double* x,
           double* fvec,
           double* fjac,
           const int ldfjac,
           const double ftol,
           const double xtol,
           const double gtol,
           const int maxfev,
           double* diag,
           const int mode,
           const double factor,
           const int nprint,
           int* info,
           int* nfev,
           int* njev,
           int* ipvt,
           double* qtf,
           double* wa1,
           double* wa2,
           double* wa3,
           double* wa4);

/* 注：LMSTR 等其他函数可按需后续添加 */

#ifdef __cplusplus
}
#endif

#endif // MINPACK_CORE_H