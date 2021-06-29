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

  /* Prototypes for Chapter 3 */

  /* Level 1 Computational Routines */

float *c_susdot( int conj, int nz, const float *x, const int
                 *indx, const float *y, int incy, float *r,
                 enum blas_base_type index_base );
double *c_dusdot( int conj, int nz, const double *x, const int
                 *indx, const double *y, int incy, double *r,
                 enum blas_base_type index_base );
void *c_cusdot( int conj, int nz, const void *x, const int
                *indx, const void *y, int incy, void *r,
                enum blas_base_type index_base );
void *c_zusdot( int conj, int nz, const void *x, const int
                *indx, const void *y, int incy, void *r,
                enum blas_base_type index_base );

void BLAS_susaxpy( int nz, float alpha, const float *x, const int *indx,
                 float *y, int incy, enum blas_base_type index_base );
void BLAS_dusaxpy( int nz, double alpha, const double *x, const int *indx,
                 double *y, int incy, enum blas_base_type index_base );
void BLAS_cusaxpy( int nz, const void *alpha, const void *x, const int *indx,
                 void *y, int incy, enum blas_base_type index_base );
void BLAS_zusaxpy( int nz, const void *alpha, const void *x, const int *indx,
                 void *y, int incy, enum blas_base_type index_base );

void BLAS_susga( int nz, const float *y, int incy, float *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_dusga( int nz, const double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_cusga( int nz, const void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_zusga( int nz, const void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );

void BLAS_susgz( int nz, float *y, int incy, float *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_dusgz( int nz, double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_cusgz( int nz, void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_zusgz( int nz, void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );

void BLAS_sussc( int nz, const float *x, float *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_dussc( int nz, const double *x, double *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_cussc( int nz, const void *x, void *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_zussc( int nz, const void *x, void *y, int incy, const int *indx,
              enum blas_base_type index_base );

               /* Level 2 Computational Routines */

void BLAS_susmv( enum blas_trans_type transa, float alpha, int A,
              const float *x, int incx, float *y, int incy );
void BLAS_dusmv( enum blas_trans_type transa, double alpha, int A,
              const double *x, int incx, double *y, int incy );
void BLAS_cusmv( enum blas_trans_type transa, const void *alpha, int A,
              const void *x, int incx, void *y, int incy );
void BLAS_zusmv( enum blas_trans_type transa, const void *alpha, int A,
              const void *x, int incx, void *y, int incy );

void BLAS_sussv( enum blas_trans_type transt, float alpha, int T,
              float *x, int incx );
void BLAS_dussv( enum blas_trans_type transt, double alpha, int T,
              double *x, int incx );
void BLAS_cussv( enum blas_trans_type transt, const void *alpha, int T,
              void *x, int incx );
void BLAS_zussv( enum blas_trans_type transt, const void *alpha, int T,
              void *x, int incx );

               /* Level 3 Computational Routines */

void BLAS_susmm( enum blas_order_type order, enum blas_trans_type transa,
              int k, float alpha, int a, const float *b, int ldb,
              float *c, int ldc );
void BLAS_dusmm( enum blas_order_type order, enum blas_trans_type transa,
              int k, double alpha, int a, const double *b, int ldb,
              double *c, int ldc );
void BLAS_cusmm( enum blas_order_type order, enum blas_trans_type transa,
              int k, const void *alpha, int a, const void *b, int ldb,
              void *c, int ldc );
void BLAS_zusmm( enum blas_order_type order, enum blas_trans_type transa,
              int k, const void *alpha, int a, const void *b, int ldb,
              void *c, int ldc );

void BLAS_sussm( enum blas_order_type order, enum blas_trans_type transt,
              int k, float alpha, int t, float *b, int ldb );
void BLAS_dussm( enum blas_order_type order, enum blas_trans_type transt,
              int k, double alpha, int t, double *b, int ldb );
void BLAS_cussm( enum blas_order_type order, enum blas_trans_type transt,
              int k, const void *alpha, int t, void *b, int ldb );
void BLAS_zussm( enum blas_order_type order, enum blas_trans_type transt,
              int k, const void *alpha, int t, void *b, int ldb );

               /* Handle Management Routines */

               /* Creation Routines */

int c_suscr_begin( int m, int n );
int c_duscr_begin( int m, int n );
int c_cuscr_begin( int m, int n );
int c_zuscr_begin( int m, int n );

int c_suscr_block_begin( int Mb, int Nb, int k, int l );
int c_duscr_block_begin( int Mb, int Nb, int k, int l );
int c_cuscr_block_begin( int Mb, int Nb, int k, int l );
int c_zuscr_block_begin( int Mb, int Nb, int k, int l );

int c_suscr_variable_block_begin( int Mb, int Nb, const int *k,
                                  const int *l );
int c_duscr_variable_block_begin( int Mb, int Nb, const int *k,
                                  const int *l );
int c_cuscr_variable_block_begin( int Mb, int Nb, const int *k,
                                  const int *l );
int c_zuscr_variable_block_begin( int Mb, int Nb, const int *k,
                                  const int *l );

               /* Insertion Routines */

int c_suscr_insert( int a, float val, int i, int j );
int c_duscr_insert( int a, double val, int i, int j );
int c_cuscr_insert( int a, void *val, int i, int j );
int c_zuscr_insert( int a, void *val, int i, int j );

int c_suscr_insert_entries( int a, int nz, const float *val,
                            const int *indx, const int *jndx );
int c_duscr_insert_entries( int a, int nz, const double *val,
                            const int *indx, const int *jndx );
int c_cuscr_insert_entries( int a, int nz, const void *val,
                            const int *indx, const int *jndx );
int c_zuscr_insert_entries( int a, int nz, const void *val,
                            const int *indx, const int *jndx );

int c_suscr_insert_col( int a, int nz, const float *val, const int *indx );
int c_duscr_insert_col( int a, int nz, const double *val, const int *indx );
int c_cuscr_insert_col( int a, int nz, const void *val, const int *indx );
int c_zuscr_insert_col( int a, int nz, const void *val, const int *indx );

int c_suscr_insert_row( int a, int nz, const float *val, const int *indx );
int c_duscr_insert_row( int a, int nz, const double *val, const int *indx );
int c_cuscr_insert_row( int a, int nz, const void *val, const int *indx );
int c_zuscr_insert_row( int a, int nz, const void *val, const int *indx );

int c_suscr_insert_clique( int a, const int k, const int l, const float *val,
                           const int row_stride, const int col_stride,
                           const int *indx, const int *jndx );
int c_duscr_insert_clique( int a, const int k, const int l, const double *val,
                           const int row_stride, const int col_stride,
                           const int *indx, const int *jndx );
int c_cuscr_insert_clique( int a, const int k, const int l, const void *val,
                           const int row_stride, const int col_stride,
                           const int *indx, const int *jndx );
int c_zuscr_insert_clique( int a, const int k, const int l, const void *val,
                           const int row_stride, const int col_stride,
                           const int *indx, const int *jndx );

int c_suscr_insert_block( int a, const float *val, int row_stride,
                        int col_stride, int bi, int bj );
int c_duscr_insert_block( int a, const double *val, int row_stride,
                        int col_stride, int bi, int bj );
int c_cuscr_insert_block( int a, const void *val, int row_stride,
                        int col_stride, int bi, int bj );
int c_zuscr_insert_block( int a, const void *val, int row_stride,
                        int col_stride, int bi, int bj );

               /* Completion of Construction Routines */

int c_suscr_end( int a );
int c_duscr_end( int a );
int c_cuscr_end( int a );
int c_zuscr_end( int a );

               /* Matrix Property Routines */

int c_usgp( int a, int pname );

int c_ussp( int a, int pname );

               /* Destruction Routine */

void BLAS_usds( int a );
