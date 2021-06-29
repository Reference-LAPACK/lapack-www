  /* Enumerated types */

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_trans_type {
            blas_no_trans   = 111,
            blas_trans      = 112,
            blas_conj_trans = 113 };

enum blas_uplo_type  {
            blas_upper = 121,
            blas_lower = 122 };

enum blas_diag_type {
            blas_non_unit_diag = 131,
            blas_unit_diag     = 132 };

enum blas_side_type {
            blas_left_side  = 141,
            blas_right_side = 142 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

enum blas_sort_type {
            blas_increasing_order = 181,
            blas_decreasing_order = 182 };

enum blas_conj_type {
            blas_conj    = 191,
            blas_no_conj = 192 };

enum blas_jrot_type {
            blas_jrot_inner  = 201,
            blas_jrot_outer  = 202,
            blas_jrot_sorted = 203 };

enum blas_prec_type {
            blas_single     = 211,
            blas_double     = 212,
            blas_indigenous = 213,
            blas_extra      = 214 };

enum blas_base_type {
            blas_zero_base = 221,
            blas_one_base  = 222 };

enum blas_symmetry_type {
            blas_general          = 231,
            blas_symmetric        = 232,
            blas_hermitian        = 233,
            blas_triangular       = 234,
            blas_lower_triangular = 235,
            blas_upper_triangular = 236 };

enum blas_field_type {
            blas_complex          = 241,
            blas_real             = 242,
            blas_double_precision = 243,
            blas_single_precision = 244  };

enum blas_size_type {
            blas_num_rows      = 251,
            blas_num_cols      = 252,
            blas_num_nonzeros  = 253  };

enum blas_valid_handle_type {
            blas_valid_handle = 261 };

enum blas_sparsity_optimization_type {
            blas_regular       = 271,
            blas_irregular     = 272,
            blas_block         = 273,
            blas_unassembled   = 274 };

  /* Prototypes for Chapter 4 */

  /* Reduction Operations */

#ifndef __CBLAS_H_
#define __CBLAS_H_

#include <malloc.h>

enum type_prec {
    real_S = 1,
    real_D = 2,
    real_I = 3,
    real_E = 4,
    complex_S = 5,
    complex_D = 6,
    complex_I = 7,
    complex_E = 8
};

/* constants */

#define BITS_S  24
#define BITS_D  53
#define BITS_E  106

/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
#define split   (134217729.0)

/* macros */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* prototypes */

extern void *blas_malloc(size_t);
extern void blas_free(void *);
extern void *blas_realloc(void *, size_t);

extern void BLAS_ddot_s_s(enum blas_conj_type conj,int n, double alpha,
                        const float* x, int incx, double beta,
                        const float* y, int incy, double* r);

extern void BLAS_ddot_s_d(enum blas_conj_type conj,int n, double alpha,
                        const float* x, int incx, double beta,
                        const double* y, int incy, double* r);

extern void BLAS_ddot_d_s(enum blas_conj_type conj,int n, double alpha,
                        const double* x, int incx, double beta,
                        const float* y, int incy, double* r);

extern void BLAS_zdot_c_c(enum blas_conj_type conj,int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r);

extern void BLAS_zdot_c_z(enum blas_conj_type conj,int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r);

extern void BLAS_zdot_z_c(enum blas_conj_type conj,int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r);

extern void BLAS_cdot_s_s(enum blas_conj_type conj,int n, void* alpha,
                        const float* x, int incx, void* beta,
                        const float* y, int incy, void* r);

extern void BLAS_cdot_s_c(enum blas_conj_type conj,int n, void* alpha,
                        const float* x, int incx, void* beta,
                        const void* y, int incy, void* r);

extern void BLAS_cdot_c_s(enum blas_conj_type conj,int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const float* y, int incy, void* r);

extern void BLAS_zdot_d_d(enum blas_conj_type conj,int n, void* alpha,
                        const double* x, int incx, void* beta,
                        const double* y, int incy, void* r);

extern void BLAS_zdot_z_d(enum blas_conj_type conj,int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const double* y, int incy, void* r);

extern void BLAS_zdot_d_z(enum blas_conj_type conj,int n, void* alpha,
                        const double* x, int incx, void* beta,
                        const void* y, int incy, void* r);

extern void BLAS_sdot_x(enum blas_conj_type conj, int n, float alpha,
                        const float* x, int incx, float beta,
                        const float* y, int incy, float* r,
                        enum blas_prec_type prec);

extern void BLAS_ddot_x(enum blas_conj_type conj, int n, double alpha,
                        const double* x, int incx, double beta,
                        const double* y, int incy, double* r,
                        enum blas_prec_type prec);

extern void BLAS_cdot_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_ddot_s_s_x(enum blas_conj_type conj, int n, double alpha,
                        const float* x, int incx, double beta,
                        const float* y, int incy, double* r,
                        enum blas_prec_type prec);

extern void BLAS_ddot_s_d_x(enum blas_conj_type conj, int n, double alpha,
                        const float* x, int incx, double beta,
                        const double* y, int incy, double* r,
                        enum blas_prec_type prec);

extern void BLAS_ddot_d_s_x(enum blas_conj_type conj, int n, double alpha,
                        const double* x, int incx, double beta,
                        const float* y, int incy, double* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_c_c_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_c_z_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_z_c_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_cdot_s_s_x(enum blas_conj_type conj, int n, void* alpha,
                        const float* x, int incx, void* beta,
                        const float* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_cdot_s_c_x(enum blas_conj_type conj, int n, void* alpha,
                        const float* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_cdot_c_s_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const float* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_d_d_x(enum blas_conj_type conj, int n, void* alpha,
                        const double* x, int incx, void* beta,
                        const double* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_d_z_x(enum blas_conj_type conj, int n, void* alpha,
                        const double* x, int incx, void* beta,
                        const void* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_zdot_z_d_x(enum blas_conj_type conj, int n, void* alpha,
                        const void* x, int incx, void* beta,
                        const double* y, int incy, void* r,
                        enum blas_prec_type prec);

extern void BLAS_csum_x( int n, const void* x, int incx, void* sum,
                      enum blas_prec_type prec );
extern void BLAS_dsum_x( int n, const double* x, int incx, double* sum,
                      enum blas_prec_type prec );
extern void BLAS_ssum_x( int n, const float* x, int incx, float* sum,
                      enum blas_prec_type prec );
extern void BLAS_zsum_x( int n, const void* x, int incx, void* sum,
                      enum blas_prec_type prec );

extern void BLAS_saxpby_x(int n, float alpha, const float *x, int incx,
                       float beta, float *y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_daxpby_x(int n, double alpha, const double *x, int incx,
                       double beta, double *y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_caxpby_x(int n, void *alpha, const void *x, int incx,
                       void *beta, void *y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zaxpby_x(int n, void *alpha, const void *x, int incx,
                       void *beta, void *y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_daxpby_s_x(int n, double alpha, const float *x, int incx,
                         double beta, double *y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_zaxpby_c_x(int n, void *alpha, const void *x, int incx,
                         void *beta, void *y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_caxpby_s_x(int n, void *alpha, const float *x, int incx,
                         void *beta, void *y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_zaxpby_d_x(int n, void *alpha, const double *x, int incx,
                         void *beta, void *y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_daxpby_s(int n, double alpha, const float *x, int incx,
                       double beta, double *y, int incy);

extern void BLAS_caxpby_s(int n, void *alpha, const float *x, int incx,
                       void *beta, void *y, int incy);

extern void BLAS_zaxpby_c(int n, void *alpha, const void *x, int incx,
                       void *beta, void *y, int incy);

extern void BLAS_zaxpby_d(int n, void *alpha, const double *x, int incx,
                       void *beta, void *y, int incy);

extern void BLAS_swaxpby_x(int n, float alpha, const float *x, int incx,
                        float beta, const float *y, int incy, float *w,
                        int incw, enum blas_prec_type prec);

extern void BLAS_dwaxpby_x(int n, double alpha, const double *x, int incx,
                        double beta, const double *y, int incy, double *w,
                        int incw, enum blas_prec_type prec);

extern void BLAS_cwaxpby_x(int n, void *alpha, const void *x, int incx,
                        void *beta, const void *y, int incy, void *w,
                        int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_x(int n, void *alpha, const void *x, int incx,
                        void *beta, const void *y, int incy, void *w,
                        int incw, enum blas_prec_type prec);

extern void BLAS_dwaxpby_s_s_x(int n, double alpha, const float *x, int incx,
                            double beta, const float *y, int incy, double *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_dwaxpby_s_d_x(int n, double alpha, const float *x, int incx,
                            double beta, const double *y, int incy, double *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_dwaxpby_d_s_x(int n, double alpha, const double *x, int incx,
                            double beta, const float *y, int incy, double *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_cwaxpby_s_s_x(int n, void *alpha, const float *x, int incx,
                            void *beta, const float *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_cwaxpby_c_s_x(int n, void *alpha, const void *x, int incx,
                            void *beta, const float *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_cwaxpby_s_c_x(int n, void *alpha, const float *x, int incx,
                            void *beta, const void *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_c_c_x(int n, void *alpha, const void *x, int incx,
                            void *beta, const void *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_z_c_x(int n, void *alpha, const void *x, int incx,
                            void *beta, const void *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_c_z_x(int n, void *alpha, const void *x, int incx,
                            void *beta, const void *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_d_d_x(int n, void *alpha, const double *x, int incx,
                            void *beta, const double *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_z_d_x(int n, void *alpha, const void *x, int incx,
                            void *beta, const double *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_zwaxpby_d_z_x(int n, void *alpha, const double *x, int incx,
                            void *beta, const void *y, int incy, void *w,
                            int incw, enum blas_prec_type prec);

extern void BLAS_dwaxpby_s_s(int n, double alpha, const float *x, int incx,
                          double beta, const float *y, int incy, double *w,
                          int incw);

extern void BLAS_dwaxpby_s_d(int n, double alpha, const float *x, int incx,
                          double beta, const double *y, int incy, double *w,
                          int incw);

extern void BLAS_dwaxpby_d_s(int n, double alpha, const double *x, int incx,
                          double beta, const float *y, int incy, double *w,
                          int incw);

extern void BLAS_cwaxpby_s_s(int n, void *alpha, const float *x, int incx,
                          void *beta, const float *y, int incy, void *w,
                          int incw);

extern void BLAS_cwaxpby_c_s(int n, void *alpha, const void *x, int incx,
                          void *beta, const float *y, int incy, void *w,
                          int incw);

extern void BLAS_cwaxpby_s_c(int n, void *alpha, const float *x, int incx,
                          void *beta, const void *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_c_c(int n, void *alpha, const void *x, int incx,
                          void *beta, const void *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_z_c(int n, void *alpha, const void *x, int incx,
                          void *beta, const void *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_c_z(int n, void *alpha, const void *x, int incx,
                          void *beta, const void *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_d_d(int n, void *alpha, const double *x, int incx,
                          void *beta, const double *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_z_d(int n, void *alpha, const void *x, int incx,
                          void *beta, const double *y, int incy, void *w,
                          int incw);

extern void BLAS_zwaxpby_d_z(int n, void *alpha, const double *x, int incx,
                          void *beta, const void *y, int incy, void *w,
                          int incw);

extern void BLAS_sspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
                      int n, float alpha, float* AP, float* x, int incx,
                      float beta, float* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_dspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
                      int n, double alpha, double* AP, double* x, int incx,
                      double beta, double* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_cspmv(enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, void* alpha, void* AP, void* x, int incx,
                    void* beta, void* y, int incy);

extern void BLAS_cspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
                      int n, void* alpha, void* AP, void* x, int incx,
                      void* beta, void* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_zspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
                      int n, void* alpha, void* AP, void* x, int incx,
                      void* beta, void* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_dspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, double alpha, float* AP, float* x,
                     int incx, double beta, double* y, int incy);

extern void BLAS_dspmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, double alpha, float* AP, float* x,
                       int incx, double beta, double* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_dspmv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, double alpha, float* AP, double* x,
                     int incx, double beta, double* y, int incy);

extern void BLAS_dspmv_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, double alpha, float* AP, double* x,
                       int incx, double beta, double* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_dspmv_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, double alpha, double* AP, float* x,
                     int incx, double beta, double* y, int incy);

extern void BLAS_dspmv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, double alpha, double* AP, float* x,
                       int incx, double beta, double* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, void* AP, void* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, void* AP, void* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, void* AP, void* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, void* AP, void* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, void* AP, void* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, void* AP, void* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_cspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, float* AP, float* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_cspmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, float* AP, float* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_cspmv_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, float* AP, void* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_cspmv_s_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, float* AP, void* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_cspmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, void* AP, float* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_cspmv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, void* AP, float* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, double* AP, double* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, double* AP, double* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, void* AP, double* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, void* AP, double* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_zspmv_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
                     int n, void* alpha, double* AP, void* x,
                     int incx, void* beta, void* y, int incy);

extern void BLAS_zspmv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, void* alpha, double* AP, void* x,
                       int incx, void* beta, void* y, int incy,
                       enum blas_prec_type prec);

extern void BLAS_dgbmv_s_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        double alpha, float* a, int lda,
                        float* x, int incx, double beta,
                        double* y, int incy);
extern void BLAS_dgbmv_s_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        double alpha, float* a, int lda,
                        double* x, int incx, double beta,
                        double* y, int incy);
extern void BLAS_dgbmv_d_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        double alpha, double* a, int lda,
                        float* x, int incx, double beta,
                        double* y, int incy);

extern void BLAS_zgbmv_c_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, void* a, int lda,
                        void* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_zgbmv_c_z(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, void* a, int lda,
                        void* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_zgbmv_z_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, void* a, int lda,
                        void* x, int incx, void* beta,
                        void* y, int incy);

extern void BLAS_cgbmv_s_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, float* a, int lda,
                        float* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_cgbmv_s_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, float* a, int lda,
                        void* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_cgbmv_c_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, void* a, int lda,
                        float* x, int incx, void* beta,
                        void* y, int incy);

extern void BLAS_zgbmv_d_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, double* a, int lda,
                        double* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_zgbmv_z_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, void* a, int lda,
                        double* x, int incx, void* beta,
                        void* y, int incy);
extern void BLAS_zgbmv_d_z(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, int kl, int ku,
                        void* alpha, double* a, int lda,
                        void* x, int incx, void* beta,
                        void* y, int incy);

extern void BLAS_sgbmv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, int kl, int ku,
                      float alpha, float* a, int lda,
                      float* x, int incx, float beta,
                      float* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_dgbmv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, int kl, int ku,
                      double alpha, double* a, int lda,
                      double* x, int incx, double beta,
                      double* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_cgbmv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, int kl, int ku,
                      void* alpha, void* a, int lda,
                      void* x, int incx, void* beta,
                      void* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_zgbmv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, int kl, int ku,
                      void* alpha, void* a, int lda,
                      void* x, int incx, void* beta,
                      void* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_dgbmv_s_s_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha, float* a, int lda,
                         float* x, int incx, double beta,
                         double* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_dgbmv_s_d_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha, float* a, int lda,
                         double* x, int incx, double beta,
                         double* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_dgbmv_d_s_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha, double* a, int lda,
                         float* x, int incx, double beta,
                         double* y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_zgbmv_c_c_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, void* a, int lda,
                         void* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_zgbmv_c_z_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, void* a, int lda,
                         void* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_zgbmv_z_c_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, void* a, int lda,
                         void* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_cgbmv_s_s_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, float* a, int lda,
                         float* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_cgbmv_s_c_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, float* a, int lda,
                         void* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_cgbmv_c_s_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, void* a, int lda,
                         float* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_zgbmv_d_d_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, double* a, int lda,
                         double* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_zgbmv_d_z_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, double* a, int lda,
                         void* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);
extern void BLAS_zgbmv_z_d_x(enum blas_order_type order,
                         enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         void* alpha, void* a, int lda,
                         double* x, int incx, void* beta,
                         void* y, int incy,
                         enum blas_prec_type prec);

extern void BLAS_strsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
             float alpha, const float* T, int ldt, float* x, int incx,
             enum blas_prec_type prec);

extern void BLAS_dtrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
             double alpha, const double* T, int ldt, double* x, int incx,
             enum blas_prec_type prec);

extern void BLAS_dtrsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
             double alpha, const float* T, int ldt, double* x, int incx,
             enum blas_prec_type prec);

extern void BLAS_dtrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
             double alpha, const float* T, int ldt, double* x, int incx);

extern void BLAS_sgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, float alpha, float* a, int lda,
                      float* x, int incx, float beta, float* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_dgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, double alpha, double* a, int lda,
                      double* x, int incx, double beta, double* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_zgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, void* alpha, void* a, int lda,
                      void* x, int incx, void* beta, void* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_cgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, void* alpha, void* a, int lda,
                      void* x, int incx, void* beta, void* y, int incy,
                      enum blas_prec_type prec);

void BLAS_dgemv_s_s_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, double alpha, float* a, int lda,
                   float* x, int incx, double beta, double* y, int incy,
                   enum blas_prec_type prec);
void BLAS_dgemv_s_d_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, double alpha, float* a, int lda,
                   double* x, int incx, double beta, double* y, int incy,
                   enum blas_prec_type prec);
void BLAS_dgemv_d_s_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, double alpha, double* a, int lda,
                   float* x, int incx, double beta, double* y, int incy,
                   enum blas_prec_type prec);

void BLAS_zgemv_c_c_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, void* a, int lda,
                   void* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_zgemv_c_z_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, void* a, int lda,
                   void* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_zgemv_z_c_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, void* a, int lda,
                   void* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);

void BLAS_cgemv_s_s_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, float* a, int lda,
                   float* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_cgemv_s_c_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, float* a, int lda,
                   void* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_cgemv_c_s_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, void* a, int lda,
                   float* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);

void BLAS_zgemv_d_d_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, double* a, int lda,
                   double* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_zgemv_d_z_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, double* a, int lda,
                   void* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);
void BLAS_zgemv_z_d_x(enum blas_order_type order,
                   enum blas_trans_type trans,
                   int m, int n, void* alpha, void* a, int lda,
                   double* x, int incx, void* beta, void* y, int incy,
                   enum blas_prec_type prec);

extern void BLAS_sgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, float alpha, float* a, int lda,
                      float* x, int incx, float beta, float* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_dgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, double alpha, double* a, int lda,
                      double* x, int incx, double beta, double* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_zgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, void* alpha, void* a, int lda,
                      void* x, int incx, void* beta, void* y, int incy,
                      enum blas_prec_type prec);
extern void BLAS_cgemv_x(enum blas_order_type order,
                      enum blas_trans_type trans,
                      int m, int n, void* alpha, void* a, int lda,
                      void* x, int incx, void* beta, void* y, int incy,
                      enum blas_prec_type prec);

extern void BLAS_sgemv(enum blas_order_type order,
                    enum blas_trans_type trans,
                    int m, int n, float alpha, float* a, int lda,
                    float* x, int incx, float beta, float* y, int incy);
extern void BLAS_dgemv(enum blas_order_type order,
                    enum blas_trans_type trans,
                    int m, int n, double alpha, double* a, int lda,
                    double* x, int incx, double beta, double* y, int incy);
extern void BLAS_zgemv(enum blas_order_type order,
                    enum blas_trans_type trans,
                    int m, int n, void* alpha, void* a, int lda,
                    void* x, int incx, void* beta, void* y, int incy);
extern void BLAS_cgemv(enum blas_order_type order,
                    enum blas_trans_type trans,
                    int m, int n, void* alpha, void* a, int lda,
                    void* x, int incx, void* beta, void* y, int incy);

extern void BLAS_dgemv_s_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, double alpha, float* a, int lda,
                        float* x, int incx, double beta, double* y, int incy);
extern void BLAS_dgemv_s_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, double alpha, float* a, int lda,
                        double* x, int incx, double beta, double* y, int incy);
extern void BLAS_dgemv_d_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, double alpha, double* a, int lda,
                        float* x, int incx, double beta, double* y, int incy);

extern void BLAS_zgemv_c_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, void* a, int lda,
                        void* x, int incx, void* beta, void* y, int incy);
extern void BLAS_zgemv_c_z(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, void* a, int lda,
                        void* x, int incx, void* beta, void* y, int incy);
extern void BLAS_zgemv_z_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, void* a, int lda,
                        void* x, int incx, void* beta, void* y, int incy);

extern void BLAS_cgemv_s_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, float* a, int lda,
                        float* x, int incx, void* beta, void* y, int incy);
extern void BLAS_cgemv_s_c(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, float* a, int lda,
                        void* x, int incx, void* beta, void* y, int incy);
extern void BLAS_cgemv_c_s(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, void* a, int lda,
                        float* x, int incx, void* beta, void* y, int incy);

extern void BLAS_zgemv_d_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, double* a, int lda,
                        double* x, int incx, void* beta, void* y, int incy);
extern void BLAS_zgemv_d_z(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, double* a, int lda,
                        void* x, int incx, void* beta, void* y, int incy);
extern void BLAS_zgemv_z_d(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n, void* alpha, void* a, int lda,
                        double* x, int incx, void* beta, void* y, int incy);

extern void BLAS_dgemm_d_s(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        double alpha, double* a, int lda, float* b, int ldb,
                        double beta, double* c, int ldc);
extern void BLAS_dgemm_s_d(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        double alpha, float* a, int lda, double* b, int ldb,
                        double beta, double* c, int ldc);
extern void BLAS_dgemm_s_s(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        double alpha, float* a, int lda, float* b, int ldb,
                        double beta, double* c, int ldc);
extern void BLAS_zgemm_z_c(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, void* a, int lda, void* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_zgemm_c_z(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, void* a, int lda, void* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_zgemm_c_c(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, void* a, int lda, void* b, int ldb,
                        void* beta, void* c, int ldc );
extern void BLAS_cgemm_c_s(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, void* a, int lda, float* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_cgemm_s_c(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, float* a, int lda, void* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_cgemm_s_s(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, float* a, int lda, float* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_zgemm_z_d(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, void* a, int lda, double* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_zgemm_d_z(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, double* a, int lda, void* b, int ldb,
                        void* beta, void* c, int ldc);
extern void BLAS_zgemm_d_d(enum blas_order_type order,
                        enum blas_trans_type transa,
                        enum blas_trans_type transb, int m, int n, int k,
                        void* alpha, double* a, int lda, double* b, int ldb,
                        void* beta, void* c, int ldc);

extern void BLAS_sgemm_x(enum blas_order_type order, enum blas_trans_type transa,
                      enum blas_trans_type transb, int m, int n, int k,
                      float alpha, float* a, int lda, float* b, int ldb,
                      float beta, float* c, int ldc, enum blas_prec_type prec);
extern void BLAS_dgemm_x(enum blas_order_type order, enum blas_trans_type transa,
                      enum blas_trans_type transb, int m, int n, int k,
                      double alpha, double* a, int lda, double* b, int ldb,
                      double beta, double* c, int ldc,
                      enum blas_prec_type prec);
extern void BLAS_cgemm_x(enum blas_order_type order, enum blas_trans_type transa,
                      enum blas_trans_type transb, int m, int n, int k,
                      void* alpha, void* a, int lda, void* b, int ldb,
                      void* beta, void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zgemm_x(enum blas_order_type order, enum blas_trans_type transa,
                      enum blas_trans_type transb, int m, int n, int k,
                      void* alpha, void* a, int lda, void* b, int ldb,
                      void* beta, void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_dgemm_d_s_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          double alpha, double* a, int lda, float* b, int ldb,
                          double beta, double* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_dgemm_s_d_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          double alpha, float* a, int lda, double* b, int ldb,
                          double beta, double* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_dgemm_s_s_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          double alpha, float* a, int lda, float* b, int ldb,
                          double beta, double* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_z_c_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, void* a, int lda, void* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_c_z_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, void* a, int lda, void* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_c_c_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, void* a, int lda, void* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_cgemm_c_s_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, void* a, int lda, float* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_cgemm_s_c_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, float* a, int lda, void* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_cgemm_s_s_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, float* a, int lda, float* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_z_d_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, void* a, int lda, double* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_d_z_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, double* a, int lda, void* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);
extern void BLAS_zgemm_d_d_x(enum blas_order_type order,
                          enum blas_trans_type transa,
                          enum blas_trans_type transb, int m, int n, int k,
                          void* alpha, double* a, int lda, double* b, int ldb,
                          void* beta, void* c, int ldc,
                          enum blas_prec_type prec);

extern void BLAS_dsymm_d_s(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        double alpha, double* a, int lda,
                        float* b, int ldb, double beta, double* c, int ldc);
extern void BLAS_dsymm_s_d(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        double alpha, float* a, int lda,
                        double* b, int ldb, double beta, double* c, int ldc );
extern void BLAS_dsymm_s_s(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        double alpha, float* a, int lda,
                        float* b, int ldb, double beta, double* c, int ldc );
extern void BLAS_zsymm_z_c(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, void* a, int lda,
                        void* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_zsymm_c_z(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, void* a, int lda,
                        void* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_zsymm_c_c(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, void* a, int lda,
                        void* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_csymm_c_s(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, void* a, int lda,
                        float* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_csymm_s_c(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, float* a, int lda,
                        void* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_csymm_s_s(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, float* a, int lda,
                        float* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_zsymm_z_d(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, void* a, int lda,
                        double* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_zsymm_d_z(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, double* a, int lda,
                        void* b, int ldb, void* beta, void* c, int ldc );
extern void BLAS_zsymm_d_d(enum blas_order_type order, enum blas_side_type side,
                        enum blas_uplo_type uplo, int m, int n,
                        void* alpha, double* a, int lda,
                        double* b, int ldb, void* beta, void* c, int ldc );

extern void BLAS_ssymm_x(enum blas_order_type order, enum blas_side_type side,
                      enum blas_uplo_type uplo, int m, int n,
                      float alpha, float* a, int lda,
                      float* b, int ldb, float beta,
                      float* c, int ldc , enum blas_prec_type prec);
extern void BLAS_dsymm_x(enum blas_order_type order, enum blas_side_type side,
                      enum blas_uplo_type uplo, int m, int n,
                      double alpha, double* a, int lda,
                      double* b, int ldb, double beta,
                      double* c, int ldc , enum blas_prec_type prec);
extern void BLAS_csymm_x(enum blas_order_type order, enum blas_side_type side,
                      enum blas_uplo_type uplo, int m, int n,
                      void* alpha, void* a, int lda,
                      void* b, int ldb, void* beta,
                      void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_x(enum blas_order_type order, enum blas_side_type side,
                      enum blas_uplo_type uplo, int m, int n,
                      void* alpha, void* a, int lda,
                      void* b, int ldb, void* beta,
                      void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_dsymm_d_s_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          double alpha, double* a, int lda,
                          float* b, int ldb, double beta,
                          double* c, int ldc, enum blas_prec_type prec);
extern void BLAS_dsymm_s_d_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          double alpha, float* a, int lda,
                          double* b, int ldb, double beta,
                          double* c, int ldc, enum blas_prec_type prec);
extern void BLAS_dsymm_s_s_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          double alpha, float* a, int lda,
                          float* b, int ldb, double beta,
                          double* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_z_c_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, void* a, int lda,
                          void* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_c_z_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, void* a, int lda,
                          void* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_c_c_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, void* a, int lda,
                          void* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_csymm_c_s_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, void* a, int lda,
                          float* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_csymm_s_c_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, float* a, int lda,
                          void* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_csymm_s_s_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, float* a, int lda,
                          float* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_z_d_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, void* a, int lda,
                          double* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_d_z_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, double* a, int lda,
                          void* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);
extern void BLAS_zsymm_d_d_x(enum blas_order_type order,
                          enum blas_side_type side,
                          enum blas_uplo_type uplo, int m, int n,
                          void* alpha, double* a, int lda,
                          double* b, int ldb, void* beta,
                          void* c, int ldc, enum blas_prec_type prec);

void BLAS_zhemm_z_c(enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, void* alpha,
                 void* a, int lda, void* b, int ldb, void* beta,
                 void* c, int ldc);
void BLAS_zhemm_c_z(enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, void* alpha,
                 void* a, int lda, void* b, int ldb, void* beta,
                 void* c, int ldc);
void BLAS_zhemm_c_c(enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, void* alpha,
                 void* a, int lda, void* b, int ldb, void* beta,
                 void* c, int ldc);
void BLAS_chemm_c_s(enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, void* alpha,
                 void* a, int lda, float* b, int ldb, void* beta,
                 void* c, int ldc);
void BLAS_zhemm_z_d(enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, void* alpha,
                 void* a, int lda, double* b, int ldb, void* beta,
                 void* c, int ldc);

void BLAS_chemm_x(enum blas_order_type order, enum blas_side_type side,
               enum blas_uplo_type uplo, int m, int n, void* alpha,
               void* a, int lda, void* b, int ldb, void* beta,
               void* c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_x(enum blas_order_type order, enum blas_side_type side,
               enum blas_uplo_type uplo, int m, int n, void* alpha,
               void* a, int lda, void* b, int ldb, void* beta,
               void* c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_z_c_x(enum blas_order_type order, enum blas_side_type side,
                   enum blas_uplo_type uplo, int m, int n, void* alpha,
                   void* a, int lda, void* b, int ldb, void* beta,
                   void* c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_c_z_x(enum blas_order_type order, enum blas_side_type side,
                   enum blas_uplo_type uplo, int m, int n, void* alpha,
                   void* a, int lda, void* b, int ldb, void* beta,
                   void* c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_c_c_x(enum blas_order_type order, enum blas_side_type side,
                   enum blas_uplo_type uplo, int m, int n, void* alpha,
                   void* a, int lda, void* b, int ldb, void* beta,
                   void* c, int ldc, enum blas_prec_type prec);
void BLAS_chemm_c_s_x(enum blas_order_type order, enum blas_side_type side,
                   enum blas_uplo_type uplo, int m, int n, void* alpha,
                   void* a, int lda, float* b, int ldb, void* beta,
                   void* c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_z_d_x(enum blas_order_type order, enum blas_side_type side,
                   enum blas_uplo_type uplo, int m, int n, void* alpha,
                   void* a, int lda, double* b, int ldb, void* beta,
                   void* c, int ldc, enum blas_prec_type prec);

extern void CBLAS_ERROR(const char* routine_name);

#endif /* __CBLAS_H_ */

