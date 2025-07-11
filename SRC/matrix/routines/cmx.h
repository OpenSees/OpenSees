#ifndef cmx_h
#define cmx_h

#ifdef __cplusplus
extern "C" {
#endif

int cmx_inv2   (const double *a, double *ainv, int*ok_flag);
int cmx_inv3   (const double *a, double *ainv, int*ok_flag);
int cmx_inv4   (const double *a, double *ainv, int*ok_flag);
int cmx_inv5   (const double *a, double *ainv, int*ok_flag);
int cmx_inv6   (const double *a, double *ainv, int*ok_flag);
int cmx_inv6_v2(const double *a, double *ainv, int*ok_flag);
int cmx_inv6_v3(const double *a, double *ainv, int*ok_flag);

#define cmx_inv(a, inv, flag) _Generic(&(a), \
    double(*)[ 4]: cmx_inv2, \
    double(*)[ 9]: cmx_inv3, \
    double(*)[16]: cmx_inv4, \
    double(*)[25]: cmx_inv5, \
    double(*)[36]: cmx_inv6  \
    )(a, inv, flag)

#ifdef __cplusplus
}
#endif

#endif // cmx_h
