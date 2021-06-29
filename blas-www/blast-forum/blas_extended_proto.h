#ifndef BLAS_EXTENDED_PROTO_H
#define BLAS_EXTENDED_PROTO_H

/* prototypes */

/****************************
 *   Level 1 routines       *
 ****************************/

extern void BLAS_sdot(enum blas_conj_type conj, int n, float alpha, 
			const float* x, int incx, float beta,
			const float* y, int incy, float* r);
extern void BLAS_ddot(enum blas_conj_type conj, int n, double alpha, 
			const double* x, int incx, double beta,
			const double* y, int incy, double* r);
extern void BLAS_cdot(enum blas_conj_type conj, int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_zdot(enum blas_conj_type conj, int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_ddot_s_s(enum blas_conj_type conj,int n, double alpha, 
			const float* x, int incx, double beta,
			const float* y, int incy, double* r);
extern void BLAS_ddot_s_d(enum blas_conj_type conj,int n, double alpha, 
			const float* x, int incx, double beta,
			const double* y, int incy, double* r);
extern void BLAS_ddot_d_s(enum blas_conj_type conj,int n, double alpha, 
			const double* x, int incx, double beta,
			const float* y, int incy, double* r);
extern void BLAS_zdot_c_c(enum blas_conj_type conj,int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_zdot_c_z(enum blas_conj_type conj,int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_zdot_z_c(enum blas_conj_type conj,int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_cdot_s_s(enum blas_conj_type conj,int n, const void* alpha, 
			const float* x, int incx, const void* beta,
			const float* y, int incy, void* r);
extern void BLAS_cdot_s_c(enum blas_conj_type conj,int n, const void* alpha, 
			const float* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_cdot_c_s(enum blas_conj_type conj,int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const float* y, int incy, void* r);
extern void BLAS_zdot_d_d(enum blas_conj_type conj,int n, const void* alpha, 
			const double* x, int incx, const void* beta,
			const double* y, int incy, void* r);
extern void BLAS_zdot_z_d(enum blas_conj_type conj,int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const double* y, int incy, void* r);
extern void BLAS_zdot_d_z(enum blas_conj_type conj,int n, const void* alpha, 
			const double* x, int incx, const void* beta,
			const void* y, int incy, void* r);
extern void BLAS_sdot_x(enum blas_conj_type conj, int n, float alpha, 
			const float* x, int incx, float beta,
			const float* y, int incy, float* r,
			enum blas_prec_type prec);
extern void BLAS_ddot_x(enum blas_conj_type conj, int n, double alpha, 
			const double* x, int incx, double beta,
			const double* y, int incy, double* r,
			enum blas_prec_type prec);
extern void BLAS_cdot_x(enum blas_conj_type conj, int n, const void* alpha, 
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_x(enum blas_conj_type conj, int n, const void* alpha, 
			const void* x, int incx, const void* beta,
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
extern void BLAS_zdot_c_c_x(enum blas_conj_type conj, int n, const void* alpha,
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_c_z_x(enum blas_conj_type conj, int n, const void* alpha,
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_z_c_x(enum blas_conj_type conj, int n, const void* alpha,
			const void* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_cdot_s_s_x(enum blas_conj_type conj, int n, const void* alpha,
			const float* x, int incx, const void* beta,
			const float* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_cdot_s_c_x(enum blas_conj_type conj, int n, const void* alpha,
			const float* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_cdot_c_s_x(enum blas_conj_type conj, int n, const void* alpha,
			const void* x, int incx, const void* beta,
			const float* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_d_d_x(enum blas_conj_type conj, int n, const void* alpha,
			const double* x, int incx, const void* beta,
			const double* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_d_z_x(enum blas_conj_type conj, int n, const void* alpha,
			const double* x, int incx, const void* beta,
			const void* y, int incy, void* r,
			enum blas_prec_type prec);
extern void BLAS_zdot_z_d_x(enum blas_conj_type conj, int n, const void* alpha,
			const void* x, int incx, const void* beta,
			const double* y, int incy, void* r,
			enum blas_prec_type prec);


extern void BLAS_csum( int n, const void* x, int incx, void* sum );
extern void BLAS_csum_x( int n, const void* x, int incx, void* sum,
		      enum blas_prec_type prec );
extern void BLAS_dsum( int n, const double* x, int incx, double* sum );
extern void BLAS_dsum_x( int n, const double* x, int incx, double* sum,
		      enum blas_prec_type prec );
extern void BLAS_ssum( int n, const float* x, int incx, float* sum );
extern void BLAS_ssum_x( int n, const float* x, int incx, float* sum,
		      enum blas_prec_type prec );
extern void BLAS_zsum( int n, const void* x, int incx, void* sum );
extern void BLAS_zsum_x( int n, const void* x, int incx, void* sum,
		      enum blas_prec_type prec );


extern void BLAS_saxpby_x(int n, float alpha, const float *x, int incx,
		       float beta, float *y, int incy,
		       enum blas_prec_type prec);
extern void BLAS_daxpby_x(int n, double alpha, const double *x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
extern void BLAS_caxpby_x(int n, const void *alpha, const void *x, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
extern void BLAS_zaxpby_x(int n, const void *alpha, const void *x, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
extern void BLAS_daxpby_s_x(int n, double alpha, const float *x, int incx,
			 double beta, double *y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zaxpby_c_x(int n, const void *alpha, const void *x, int incx,
			 const void *beta, void *y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_caxpby_s_x(int n, const void *alpha, const float *x, int incx,
			 const void *beta, void *y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zaxpby_d_x(int n, const void *alpha, const double *x, int incx,
			const void *beta, void *y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_saxpby(int n, float alpha, const float *x, int incx,
		     float beta, float *y, int incy);
extern void BLAS_daxpby(int n, double alpha, const double *x, int incx,
		     double beta, double *y, int incy);
extern void BLAS_caxpby(int n, const void *alpha, const void *x, int incx,
		     const void *beta, void *y, int incy);
extern void BLAS_zaxpby(int n, const void *alpha, const void *x, int incx,
		     const void *beta, void *y, int incy);
extern void BLAS_daxpby_s(int n, double alpha, const float *x, int incx,
		       double beta, double *y, int incy);
extern void BLAS_caxpby_s(int n, const void *alpha, const float *x, int incx,
		       const void *beta, void *y, int incy);
extern void BLAS_zaxpby_c(int n, const void *alpha, const void *x, int incx,
		       const void *beta, void *y, int incy);
extern void BLAS_zaxpby_d(int n, const void *alpha, const double *x, int incx,
		       const void *beta, void *y, int incy);


extern void BLAS_swaxpby_x(int n, float alpha, const float *x, int incx,
			float beta, const float *y, int incy, float *w,
			int incw, enum blas_prec_type prec);
extern void BLAS_dwaxpby_x(int n, double alpha, const double *x, int incx,
			double beta, const double *y, int incy, double *w,
			int incw, enum blas_prec_type prec);
extern void BLAS_cwaxpby_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
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
extern void BLAS_cwaxpby_s_s_x(int n, const void *alpha, const float *x, int incx,
			    const void *beta, const float *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_cwaxpby_c_s_x(int n, const void *alpha, const void *x, int incx,
			    const void *beta, const float *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_cwaxpby_s_c_x(int n, const void *alpha, const float *x, int incx,
			    const void *beta, const void *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_c_c_x(int n, const void *alpha, const void *x, int incx,
			    const void *beta, const void *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_z_c_x(int n, const void *alpha, const void *x, int incx,
			    const void *beta, const void *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_c_z_x(int n, const void *alpha, const void *x, int incx,
			    const void *beta, const void *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_d_d_x(int n, const void *alpha, const double *x, int incx,
			    const void *beta, const double *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_z_d_x(int n, const void *alpha, const void *x, int incx,
			    const void *beta, const double *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_zwaxpby_d_z_x(int n, const void *alpha, const double *x, int incx,
			    const void *beta, const void *y, int incy, void *w,
			    int incw, enum blas_prec_type prec);
extern void BLAS_swaxpby(int n, float alpha, const float *x, int incx,
		      float beta, const float *y, int incy, float *w,
		      int incw);
extern void BLAS_dwaxpby(int n, double alpha, const double *x, int incx,
		      double beta, const double *y, int incy, double *w,
		      int incw);
extern void BLAS_cwaxpby(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w, int incw);
extern void BLAS_zwaxpby(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w, int incw);
extern void BLAS_dwaxpby_s_s(int n, double alpha, const float *x, int incx,
			  double beta, const float *y, int incy, double *w,
			  int incw);
extern void BLAS_dwaxpby_s_d(int n, double alpha, const float *x, int incx,
			  double beta, const double *y, int incy, double *w,
			  int incw);
extern void BLAS_dwaxpby_d_s(int n, double alpha, const double *x, int incx,
			  double beta, const float *y, int incy, double *w,
			  int incw);
extern void BLAS_cwaxpby_s_s(int n, const void *alpha, const float *x, int incx,
			  const void *beta, const float *y, int incy, void *w,
			  int incw);
extern void BLAS_cwaxpby_c_s(int n, const void *alpha, const void *x, int incx,
			  const void *beta, const float *y, int incy, void *w,
			  int incw);
extern void BLAS_cwaxpby_s_c(int n, const void *alpha, const float *x, int incx,
			  const void *beta, const void *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_c_c(int n, const void *alpha, const void *x, int incx,
			  const void *beta, const void *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_z_c(int n, const void *alpha, const void *x, int incx,
			  const void *beta, const void *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_c_z(int n, const void *alpha, const void *x, int incx,
			  const void *beta, const void *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_d_d(int n, const void *alpha, const double *x, int incx,
			  const void *beta, const double *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_z_d(int n, const void *alpha, const void *x, int incx,
			  const void *beta, const double *y, int incy, void *w,
			  int incw);
extern void BLAS_zwaxpby_d_z(int n, const void *alpha, const double *x, int incx,
			  const void *beta, const void *y, int incy, void *w,
			  int incw);



/****************************
 *   Level 2 routines       *
 ****************************/

extern void BLAS_sgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        float alpha,
                        const float* a, int lda,
                        const float* x, int incx,
                        float beta,
                        float* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_dgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        double alpha,
                        const double* a, int lda,
                        const double* x, int incx,
                        double beta,
                        double* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_zgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        const void* alpha,
                        const void* a, int lda,
                        const void* x, int incx,
                        const void* beta,
                        void* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_cgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        const void* alpha,
                        const void* a, int lda,
                        const void* x, int incx,
                        const void* beta,
                        void* y, int incy,
                        enum blas_prec_type prec);
void BLAS_dgemv_s_s_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              double alpha,
                              const float* a, int lda,
                              const float* x, int incx,
                              double beta,
                              double* y, int incy,
                              enum blas_prec_type prec);
void BLAS_dgemv_s_d_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              double alpha,
                              const float* a, int lda,
                              const double* x, int incx,
                              double beta,
                              double* y, int incy,
                              enum blas_prec_type prec);
void BLAS_dgemv_d_s_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              double alpha,
                              const double* a, int lda,
                              const float* x, int incx,
                              double beta,
                              double* y, int incy,
                              enum blas_prec_type prec);

void BLAS_zgemv_c_c_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const void* a, int lda,
                              const void* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);           
void BLAS_zgemv_c_z_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const void* a, int lda,
                              const void* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_zgemv_z_c_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const void* a, int lda,
                              const void* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_cgemv_s_s_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const float* a, int lda,
                              const float* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_cgemv_s_c_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const float* a, int lda,
                              const void* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_cgemv_c_s_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const void* a, int lda,
                              const float* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_zgemv_d_d_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const double* a, int lda,
                              const double* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec); 
void BLAS_zgemv_d_z_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const double* a, int lda,
                              const void* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);
void BLAS_zgemv_z_d_x(enum blas_order_type order,
                              enum blas_trans_type trans,
                              int m, int n,
                              const void* alpha,
                              const void* a, int lda,
                              const double* x, int incx,
                              const void* beta,
                              void* y, int incy,
                              enum blas_prec_type prec);

extern void BLAS_sgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        float alpha,
                        const float* a, int lda,
                        const float* x, int incx,
                        float beta,
                        float* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_dgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        double alpha,
                        const double* a, int lda,
                        const double* x, int incx,
                        double beta,
                        double* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_zgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        const void* alpha,
                        const void* a, int lda,
                        const void* x, int incx,
                        const void* beta,
                        void* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_cgemv_x(enum blas_order_type order,
                        enum blas_trans_type trans,
                        int m, int n,
                        const void* alpha,
                        const void* a, int lda,
                        const void* x, int incx,
                        const void* beta,
                        void* y, int incy,
                        enum blas_prec_type prec);
extern void BLAS_sgemv(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            float alpha,
                            const float* a, int lda,
                            const float* x, int incx,
                            float beta,
                            float* y, int incy);
extern void BLAS_dgemv(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            double alpha,
                            const double* a, int lda,
                            const double* x, int incx,
                            double beta,
                            double* y, int incy);
extern void BLAS_zgemv(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);
extern void BLAS_cgemv(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);
extern void BLAS_dgemv_s_s(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            double alpha,
                            const float* a, int lda,
                            const float* x, int incx,
                            double beta,
                            double* y, int incy);
extern void BLAS_dgemv_s_d(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            double alpha,
                            const float* a, int lda,
                            const double* x, int incx,
                            double beta,
                            double* y, int incy);
extern void BLAS_dgemv_d_s(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            double alpha,
                            const double* a, int lda,
                            const float* x, int incx,
                            double beta,
                            double* y, int incy);

extern void BLAS_zgemv_c_c(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);           
extern void BLAS_zgemv_c_z(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);              
extern void BLAS_zgemv_z_c(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);	 
extern void BLAS_cgemv_s_s(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const float* a, int lda,
                            const float* x, int incx,
                            const void* beta,
                            void* y, int incy);
extern void BLAS_cgemv_s_c(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const float* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);
extern void BLAS_cgemv_c_s(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const float* x, int incx,
                            const void* beta,
                            void* y, int incy);

extern void BLAS_zgemv_d_d(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const double* a, int lda,
                            const double* x, int incx,
                            const void* beta,
                            void* y, int incy); 
extern void BLAS_zgemv_d_z(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const double* a, int lda,
                            const void* x, int incx,
                            const void* beta,
                            void* y, int incy);
extern void BLAS_zgemv_z_d(enum blas_order_type order,
                            enum blas_trans_type trans,
                            int m, int n,
                            const void* alpha,
                            const void* a, int lda,
                            const double* x, int incx,
                            const void* beta,
                            void* y, int incy);


extern void BLAS_sgbmv(enum blas_order_type order,
		  	        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                float alpha,
                                const float* a,
                                int lda, 
                                const float* x, int incx,
                                float beta,
                                float* y, int incy);
extern void BLAS_dgbmv(enum blas_order_type order,
		  	        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                double alpha,
                                const double* a, int lda, 
                                const double* x, int incx,
                                double beta,
                                double* y, int incy);
extern void BLAS_cgbmv(enum blas_order_type order,
		  	        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv(enum blas_order_type order,
		  	        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_dgbmv_s_s(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                double alpha,
                                const float* a, int lda, 
                                const float* x, int incx,
                                double beta,
                                double* y, int incy);
extern void BLAS_dgbmv_s_d(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                double alpha,
                                const float* a, int lda, 
                                const double* x, int incx,
                                double beta,
                                double* y, int incy);
extern void BLAS_dgbmv_d_s(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                double alpha,
                                const double* a, int lda, 
                                const float* x, int incx,
                                double beta,
                                double* y, int incy);
extern void BLAS_zgbmv_c_c(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv_c_z(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv_z_c(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_cgbmv_s_s(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const float* a, int lda, 
                                const float* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_cgbmv_s_c(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const float* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_cgbmv_c_s(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const float* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv_d_d(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const double* a, int lda, 
                                const double* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv_z_d(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const void* a, int lda, 
                                const double* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_zgbmv_d_z(enum blas_order_type order,
			        enum blas_trans_type trans,
                                int m, int n, int kl, int ku,
                                const void* alpha,
                                const double* a, int lda, 
                                const void* x, int incx,
                                const void* beta,
                                void* y, int incy);
extern void BLAS_sgbmv_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         float alpha,
                         const float* a, int lda, 
                         const float* x, int incx,
                         float beta,
                         float* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_dgbmv_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha,
                         const double* a, int lda, 
                         const double* x, int incx,
                         double beta,
                         double* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_cgbmv_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_dgbmv_s_s_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha,
                         const float* a, int lda, 
                         const float* x, int incx,
                         double beta,
                         double* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_dgbmv_s_d_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha,
                         const float* a, int lda, 
                         const double* x, int incx,
                         double beta,
                         double* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_dgbmv_d_s_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         double alpha,
                         const double* a, int lda, 
                         const float* x, int incx,
                         double beta,
                         double* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_c_c_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_c_z_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_z_c_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_cgbmv_s_s_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const float* a, int lda, 
                         const float* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_cgbmv_s_c_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const float* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_cgbmv_c_s_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const float* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_d_d_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const double* a, int lda, 
                         const double* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_d_z_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const double* a, int lda, 
                         const void* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);
extern void BLAS_zgbmv_z_d_x(enum blas_order_type order,
			 enum blas_trans_type trans,
                         int m, int n, int kl, int ku,
                         const void* alpha,
                         const void* a, int lda, 
                         const double* x, int incx,
                         const void* beta,
                         void* y, int incy,
			 enum blas_prec_type prec);

   
extern void BLAS_ssymv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   float alpha, const float* a, int lda,
   const float* x, int incx,
   float beta,
   float* y, int incy    );
extern void BLAS_dsymv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* a, int lda,
   const double* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_csymv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zsymv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_dsymv_d_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* a, int lda,
   const float* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_dsymv_s_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* a, int lda,
   const double* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_dsymv_s_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* a, int lda,
   const float* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_zsymv_z_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zsymv_c_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zsymv_c_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_csymv_c_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_csymv_s_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_csymv_s_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zsymv_z_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zsymv_d_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zsymv_d_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_ssymv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   float alpha, const float* a, int lda,
   const float* x, int incx,
   float beta,
   float* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dsymv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* a, int lda,
   const double* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_csymv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dsymv_d_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* a, int lda,
   const float* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dsymv_s_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* a, int lda,
   const double* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dsymv_s_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* a, int lda,
   const float* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec); 
extern void BLAS_zsymv_z_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_c_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_c_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_csymv_c_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_csymv_s_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_csymv_s_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_z_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_d_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zsymv_d_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);


extern void BLAS_sspmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   float alpha, const float* ap,
   const float* x, int incx,
   float beta,
   float* y, int incy    );
extern void BLAS_dspmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* ap,
   const double* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_cspmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zspmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_dspmv_d_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* ap,
   const float* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_dspmv_s_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* ap,
   const double* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_dspmv_s_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* ap,
   const float* x, int incx,
   double beta,
   double* y, int incy    );
extern void BLAS_zspmv_z_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zspmv_c_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zspmv_c_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_cspmv_c_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_cspmv_s_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_cspmv_s_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zspmv_z_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zspmv_d_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zspmv_d_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_sspmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   float alpha, const float* ap,
   const float* x, int incx,
   float beta,
   float* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dspmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* ap,
   const double* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_cspmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dspmv_d_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const double* ap,
   const float* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dspmv_s_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* ap,
   const double* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec);
extern void BLAS_dspmv_s_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   double alpha, const float* ap,
   const float* x, int incx,
   double beta,
   double* y, int incy    , enum blas_prec_type prec); 
extern void BLAS_zspmv_z_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_c_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_c_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_cspmv_c_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_cspmv_s_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_cspmv_s_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const float* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_z_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_d_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zspmv_d_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const double* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);


extern void BLAS_chemv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhemv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zhemv_z_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhemv_c_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zhemv_c_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_chemv_c_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhemv_z_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_chemv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhemv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhemv_z_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhemv_c_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhemv_c_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_chemv_c_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhemv_z_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* a, int lda,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);


extern void BLAS_chpmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhpmv(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zhpmv_z_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhpmv_c_z(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_zhpmv_c_c(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_chpmv_c_s(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    );
extern void BLAS_zhpmv_z_d(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    ); 
extern void BLAS_chpmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhpmv_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhpmv_z_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhpmv_c_z_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhpmv_c_c_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const void* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_chpmv_c_s_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const float* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);
extern void BLAS_zhpmv_z_d_x(enum blas_order_type order,
   enum blas_uplo_type uplo, int n,
   const void* alpha, const void* ap,
   const double* x, int incx,
   const void* beta,
   void* y, int incy    , enum blas_prec_type prec);


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

extern void BLAS_strsv(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     float alpha, const float* T, int ldt, float* x, int incx);

extern void BLAS_dtrsv(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     double alpha, const double* T, int ldt, double* x, int incx);

extern void BLAS_dtrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     double alpha, const float* T, int ldt, double* x, int incx);

extern void BLAS_ctrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx, 
	     enum blas_prec_type prec);

extern void BLAS_ztrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx, 
	     enum blas_prec_type prec);

extern void BLAS_ztrsv_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx,
	     enum blas_prec_type prec);

extern void BLAS_ctrsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const float* T, int ldt, void* x, int incx,
	     enum blas_prec_type prec);

extern void BLAS_ztrsv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const double* T, int ldt, void* x, int incx,
	     enum blas_prec_type prec);                                                                                
extern void BLAS_ctrsv(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx);

extern void BLAS_ztrsv(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx);

extern void BLAS_ztrsv_c(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const void* T, int ldt, void* x, int incx);

extern void BLAS_ctrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const float* T, int ldt, void* x, int incx);

extern void BLAS_ztrsv_d(enum blas_order_type order, enum blas_uplo_type uplo,
             enum blas_trans_type trans, enum blas_diag_type diag, int n,
	     const void* alpha, const double* T, int ldt, void* x, int incx);


/****************************
 *   Level 3 routines       *
 ****************************/

extern void BLAS_sgemm(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   float alpha, const float* a, int lda, const float* b, int ldb,
   float beta, float* c, int ldc    );
extern void BLAS_dgemm(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const double* a, int lda, const double* b, int ldb,
   double beta, double* c, int ldc    );
extern void BLAS_cgemm(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_dgemm_d_s(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const double* a, int lda, const float* b, int ldb,
   double beta, double* c, int ldc    );
extern void BLAS_dgemm_s_d(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const float* a, int lda, const double* b, int ldb,
   double beta, double* c, int ldc    );
extern void BLAS_dgemm_s_s(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const float* a, int lda, const float* b, int ldb,
   double beta, double* c, int ldc    );
extern void BLAS_zgemm_z_c(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm_c_z(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm_c_c(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_cgemm_c_s(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const float* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_cgemm_s_c(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const float* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_cgemm_s_s(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const float* a, int lda, const float* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm_z_d(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const double* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm_d_z(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const double* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_zgemm_d_d(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const double* a, int lda, const double* b, int ldb,
   const void* beta, void* c, int ldc    );
extern void BLAS_sgemm_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   float alpha, const float* a, int lda, const float* b, int ldb,
   float beta, float* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_dgemm_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const double* a, int lda, const double* b, int ldb,
   double beta, double* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_cgemm_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_dgemm_d_s_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const double* a, int lda, const float* b, int ldb,
   double beta, double* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_dgemm_s_d_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const float* a, int lda, const double* b, int ldb,
   double beta, double* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_dgemm_s_s_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const float* a, int lda, const float* b, int ldb,
   double beta, double* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_z_c_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_c_z_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_c_c_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_cgemm_c_s_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const float* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_cgemm_s_c_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const float* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_cgemm_s_s_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const float* a, int lda, const float* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_z_d_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const void* a, int lda, const double* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_d_z_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const double* a, int lda, const void* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);
extern void BLAS_zgemm_d_d_x(enum blas_order_type order,
   enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   const void* alpha, const double* a, int lda, const double* b, int ldb,
   const void* beta, void* c, int ldc    , enum blas_prec_type prec);

   
extern void BLAS_ssymm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   float alpha, const float* a, int lda, 
   const float* b, int ldb, float beta, 
   float* c, int ldc );
extern void BLAS_dsymm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const double* a, int lda, 
   const double* b, int ldb, double beta, 
   double* c, int ldc );
extern void BLAS_csymm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zsymm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_dsymm_d_s(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const double* a, int lda, 
   const float* b, int ldb, double beta, 
   double* c, int ldc );
extern void BLAS_dsymm_s_d(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const float* a, int lda, 
   const double* b, int ldb, double beta, 
   double* c, int ldc );
extern void BLAS_dsymm_s_s(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const float* a, int lda, 
   const float* b, int ldb, double beta, 
   double* c, int ldc );
extern void BLAS_zsymm_z_c(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zsymm_c_z(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_zsymm_c_c(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_csymm_c_s(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_csymm_s_c(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const float* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_csymm_s_s(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const float* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zsymm_z_d(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_zsymm_d_z(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const double* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zsymm_d_d(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const double* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_ssymm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   float alpha, const float* a, int lda, 
   const float* b, int ldb, float beta, 
   float* c, int ldc , enum blas_prec_type prec);
extern void BLAS_dsymm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const double* a, int lda, 
   const double* b, int ldb, double beta, 
   double* c, int ldc , enum blas_prec_type prec);
extern void BLAS_csymm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_dsymm_d_s_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const double* a, int lda, 
   const float* b, int ldb, double beta, 
   double* c, int ldc , enum blas_prec_type prec);
extern void BLAS_dsymm_s_d_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const float* a, int lda, 
   const double* b, int ldb, double beta, 
   double* c, int ldc , enum blas_prec_type prec);
extern void BLAS_dsymm_s_s_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   double alpha, const float* a, int lda, 
   const float* b, int ldb, double beta, 
   double* c, int ldc , enum blas_prec_type prec); 
extern void BLAS_zsymm_z_c_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_c_z_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_c_c_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_csymm_c_s_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_csymm_s_c_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const float* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_csymm_s_s_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const float* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_z_d_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_d_z_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const double* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zsymm_d_d_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const double* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);

   
extern void BLAS_chemm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zhemm(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_zhemm_z_c(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zhemm_c_z(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_zhemm_c_c(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_chemm_c_s(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc );
extern void BLAS_zhemm_z_d(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc ); 
extern void BLAS_chemm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zhemm_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zhemm_z_c_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zhemm_c_z_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zhemm_c_c_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const void* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_chemm_c_s_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const float* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);
extern void BLAS_zhemm_z_d_x(enum blas_order_type order,
   enum blas_side_type side, 
   enum blas_uplo_type uplo, int m, int n, 
   const void* alpha, const void* a, int lda, 
   const double* b, int ldb, const void* beta, 
   void* c, int ldc , enum blas_prec_type prec);


extern int c_fpinfo_x(enum blas_cmach_type cmach, enum blas_prec_type prec);

#endif
 /* BLAS_EXTENDED_PROTO_H */
