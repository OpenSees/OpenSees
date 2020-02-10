#if !defined(_petsc_blockmatmult_h)
#define _petsc_blockmatmult_h

#include <petscsys.h>

#define PetscKernel_v_gets_A_times_w_1_exp(v,A,w,exp) \
do {                                                  \
  v[0] exp A[0]*w[0];                                 \
} while (0)

#define PetscKernel_v_gets_A_times_w_2_exp(v,A,w,exp) \
do {                                                  \
  v[0] exp A[0]*w[0] + A[2]*w[1];                      \
  v[1] exp A[1]*w[0] + A[3]*w[1];                      \
} while (0)

#define PetscKernel_v_gets_A_times_w_3_exp(v,A,w,exp) \
do {                                                  \
  v[0] exp A[0]*w[0] + A[3]*w[1] + A[6]*w[2];         \
  v[1] exp A[1]*w[0] + A[4]*w[1] + A[7]*w[2];         \
  v[2] exp A[2]*w[0] + A[5]*w[1] + A[8]*w[2];         \
} while (0)

#define PetscKernel_v_gets_A_times_w_4_exp(v,A,w,exp)       \
do {                                                        \
  v[0] exp A[0]*w[0] + A[4]*w[1] + A[8] *w[2] + A[12]*w[3]; \
  v[1] exp A[1]*w[0] + A[5]*w[1] + A[9] *w[2] + A[13]*w[3]; \
  v[2] exp A[2]*w[0] + A[6]*w[1] + A[10]*w[2] + A[14]*w[3]; \
  v[3] exp A[3]*w[0] + A[7]*w[1] + A[11]*w[2] + A[15]*w[3]; \
} while (0)

#define PetscKernel_v_gets_A_times_w_4_exp(v,A,w,exp)       \
do {                                                        \
  v[0] exp A[0]*w[0] + A[4]*w[1] + A[8] *w[2] + A[12]*w[3]; \
  v[1] exp A[1]*w[0] + A[5]*w[1] + A[9] *w[2] + A[13]*w[3]; \
  v[2] exp A[2]*w[0] + A[6]*w[1] + A[10]*w[2] + A[14]*w[3]; \
  v[3] exp A[3]*w[0] + A[7]*w[1] + A[11]*w[2] + A[15]*w[3]; \
} while (0)

#define PetscKernel_v_gets_A_times_w_5_exp(v,A,w,exp)                    \
do {                                                                     \
  v[0] exp A[0]*w[0] + A[5]*w[1] + A[10]*w[2] + A[15]*w[3] + A[20]*w[4]; \
  v[1] exp A[1]*w[0] + A[6]*w[1] + A[11]*w[2] + A[16]*w[3] + A[21]*w[4]; \
  v[2] exp A[2]*w[0] + A[7]*w[1] + A[12]*w[2] + A[17]*w[3] + A[22]*w[4]; \
  v[3] exp A[3]*w[0] + A[8]*w[1] + A[13]*w[2] + A[18]*w[3] + A[23]*w[4]; \
  v[4] exp A[4]*w[0] + A[9]*w[1] + A[14]*w[2] + A[19]*w[3] + A[24]*w[4]; \
} while (0)

#define PetscKernel_v_gets_A_times_w_6_exp(v,A,w,exp)                                  \
do {                                                                                   \
  v[0] exp A[0]*w[0] + A[6] *w[1] + A[12]*w[2] + A[18]*w[3] + A[24]*w[4] + A[30]*w[5]; \
  v[1] exp A[1]*w[0] + A[7] *w[1] + A[13]*w[2] + A[19]*w[3] + A[25]*w[4] + A[31]*w[5]; \
  v[2] exp A[2]*w[0] + A[8] *w[1] + A[14]*w[2] + A[20]*w[3] + A[26]*w[4] + A[32]*w[5]; \
  v[3] exp A[3]*w[0] + A[9] *w[1] + A[15]*w[2] + A[21]*w[3] + A[27]*w[4] + A[33]*w[5]; \
  v[4] exp A[4]*w[0] + A[10]*w[1] + A[16]*w[2] + A[22]*w[3] + A[28]*w[4] + A[34]*w[5]; \
  v[5] exp A[5]*w[0] + A[11]*w[1] + A[17]*w[2] + A[23]*w[3] + A[29]*w[4] + A[35]*w[5]; \
} while (0)

#define PetscKernel_v_gets_A_times_w_7_exp(v,A,w,exp)                                               \
do {                                                                                                \
  v[0] exp A[0]*w[0] + A[7] *w[1] + A[14]*w[2] + A[21]*w[3] + A[28]*w[4] + A[35]*w[5] + A[42]*w[6]; \
  v[1] exp A[1]*w[0] + A[8] *w[1] + A[15]*w[2] + A[22]*w[3] + A[29]*w[4] + A[36]*w[5] + A[43]*w[6]; \
  v[2] exp A[2]*w[0] + A[9] *w[1] + A[16]*w[2] + A[23]*w[3] + A[30]*w[4] + A[37]*w[5] + A[44]*w[6]; \
  v[3] exp A[3]*w[0] + A[10]*w[1] + A[17]*w[2] + A[24]*w[3] + A[31]*w[4] + A[38]*w[5] + A[45]*w[6]; \
  v[4] exp A[4]*w[0] + A[11]*w[1] + A[18]*w[2] + A[25]*w[3] + A[32]*w[4] + A[39]*w[5] + A[46]*w[6]; \
  v[5] exp A[5]*w[0] + A[12]*w[1] + A[19]*w[2] + A[26]*w[3] + A[33]*w[4] + A[40]*w[5] + A[47]*w[6]; \
  v[6] exp A[6]*w[0] + A[13]*w[1] + A[20]*w[2] + A[27]*w[3] + A[34]*w[4] + A[41]*w[5] + A[48]*w[6]; \
} while (0)

#define PetscKernel_v_gets_A_times_w_1(v,A,w) PetscKernel_v_gets_A_times_w_1_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_2(v,A,w) PetscKernel_v_gets_A_times_w_2_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_3(v,A,w) PetscKernel_v_gets_A_times_w_3_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_4(v,A,w) PetscKernel_v_gets_A_times_w_4_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_5(v,A,w) PetscKernel_v_gets_A_times_w_5_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_6(v,A,w) PetscKernel_v_gets_A_times_w_6_exp(v,A,w,=)
#define PetscKernel_v_gets_A_times_w_7(v,A,w) PetscKernel_v_gets_A_times_w_7_exp(v,A,w,=)
#define PetscKernel_v_gets_v_plus_A_times_w_1(v,A,w) PetscKernel_v_gets_A_times_w_1_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_2(v,A,w) PetscKernel_v_gets_A_times_w_2_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_3(v,A,w) PetscKernel_v_gets_A_times_w_3_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_4(v,A,w) PetscKernel_v_gets_A_times_w_4_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_5(v,A,w) PetscKernel_v_gets_A_times_w_5_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_6(v,A,w) PetscKernel_v_gets_A_times_w_6_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_plus_A_times_w_7(v,A,w) PetscKernel_v_gets_A_times_w_7_exp(v,A,w,+=)
#define PetscKernel_v_gets_v_minus_A_times_w_1(v,A,w) PetscKernel_v_gets_A_times_w_1_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_2(v,A,w) PetscKernel_v_gets_A_times_w_2_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_3(v,A,w) PetscKernel_v_gets_A_times_w_3_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_4(v,A,w) PetscKernel_v_gets_A_times_w_4_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_5(v,A,w) PetscKernel_v_gets_A_times_w_5_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_6(v,A,w) PetscKernel_v_gets_A_times_w_6_exp(v,A,w,-=)
#define PetscKernel_v_gets_v_minus_A_times_w_7(v,A,w) PetscKernel_v_gets_A_times_w_7_exp(v,A,w,-=)

#endif
