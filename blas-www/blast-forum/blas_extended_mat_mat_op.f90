    MODULE blas_extended_mat_mat_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 4.?.? -- Matrix-Matrix Operations
      ! .  using external procedures.
      ! Written on 18 December 2000
      ! .  Susan Blackford, UT-ICL

      PRIVATE
      PUBLIC :: gemm, symm, hemm, trmm, trsm, syrk, herk, syr2k, her2k

      INTERFACE gemm
        ! gemm  a, b and c rank 2
        ! gemv  a rank 2; b and c rank 1
        ! ger,geru,gerc   a and b rank 1; c rank 2 (called ger)
        ! For c, rank 2, c has shape (m,n)
        ! .   a has shape (m,k) if transa = blas_no_trans (the default)
        ! .               (k,m) if transa /= blas_no_trans
        ! .               (m) if rank 1
        ! .   b has shape (k,n) if transb = blas_no_trans (the default)
        ! .               (n,k) if transb /= blas_no_trans
        ! .               (n) if rank 1
        ! For c, rank 1, c has shape (m)
        ! .   a has shape (m,n) if transa = blas_no_trans (the default)
        ! .               (n,m) if transa /= blas_no_trans
        ! .   b has shape (n)
        SUBROUTINE gemm_d(a,b,c,transa,transb,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa, transb
        END SUBROUTINE gemm_d
        SUBROUTINE gemv_d(a,b,c,transa,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa
        END SUBROUTINE gemv_d
        SUBROUTINE ger_d(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ger_d
        SUBROUTINE gemm_z(a,b,c,transa,transb,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa, transb
        END SUBROUTINE gemm_z
        SUBROUTINE gemv_z(a,b,c,transa,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa
        END SUBROUTINE gemv_z
        SUBROUTINE ger_z(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ger_z
        SUBROUTINE gemm_s(a,b,c,transa,transb,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa, transb
        END SUBROUTINE gemm_s
        SUBROUTINE gemv_s(a,b,c,transa,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa
        END SUBROUTINE gemv_s
        SUBROUTINE ger_s(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ger_s
        SUBROUTINE gemm_c(a,b,c,transa,transb,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa, transb
        END SUBROUTINE gemm_c
        SUBROUTINE gemv_c(a,b,c,transa,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: transa
        END SUBROUTINE gemv_c
        SUBROUTINE ger_c(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ger_c
      END INTERFACE

      INTERFACE symm
        ! symm  a, b and c rank 2
        ! symv  a rank 2; b and c rank 1
        ! For c, rank 2, c has shape (m,n)
        ! .   b has same shape as c
        ! .   a has shape (m,m) if  side = blas_left_side (the default)
        ! .     has shape (n,n) if  side /= blas_left_side
        ! For c, rank 1, c has shape (m)
        ! .   b has same shape as c
        ! .   a has shape (m,m)
        SUBROUTINE symm_d(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symm_d
        SUBROUTINE symv_d(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symv_d
        SUBROUTINE symm_z(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symm_z
        SUBROUTINE symv_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symv_z
        SUBROUTINE symm_s(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symm_s
        SUBROUTINE symv_s(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symv_s
        SUBROUTINE symm_c(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symm_c
        SUBROUTINE symv_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE symv_c
      END INTERFACE

      INTERFACE hemm
        ! hemm  a, b and c rank 2
        ! hemv  a rank 2; b and c rank 1
        ! For c, rank 2, c has shape (m,n)
        ! .   b has same shape as c
        ! .   a has shape (m,m) if  side = blas_left_side (the default)
        ! .     has shape (n,n) if  side /= blas_left_side
        ! For c, rank 1, c has shape (m)
        ! .   b has same shape as c
        ! .   a has shape (m,m)
        SUBROUTINE hemm_z(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hemm_z
        SUBROUTINE hemv_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hemv_z
        SUBROUTINE hemm_c(a,b,c,side,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hemm_c
        SUBROUTINE hemv_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hemv_c
      END INTERFACE

      INTERFACE trmm
        ! trsm  t and b rank 2
        ! trmv  t rank 2; b rank 1
        ! For b, rank 2, b has shape (m,n)
        ! .   t has shape (m,m) if  side = blas_left_side (the default)
        ! .     has shape (n,n) if  side /= blas_left_side
        ! For b, rank 1, b has shape (m)
        ! .   t has shape (m,m)
        SUBROUTINE trmm_d(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmm_d
        SUBROUTINE trmv_d(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmv_d
        SUBROUTINE trmm_z(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmm_z
        SUBROUTINE trmv_z(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmv_z
        SUBROUTINE trmm_s(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmm_s
        SUBROUTINE trmv_s(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmv_s
        SUBROUTINE trmm_c(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmm_c
        SUBROUTINE trmv_c(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trmv_c
      END INTERFACE

      INTERFACE trsm
        ! trsm  t and b rank 2
        ! trsv  t rank 2; b rank 1
        ! For b, rank 2, b has shape (m,n)
        ! .   t has shape (m,m) if  side = blas_left_side (the default)
        ! .     has shape (n,n) if  side /= blas_left_side
        ! For b, rank 1, b has shape (m)
        ! .   t has shape (m,m)
        SUBROUTINE trsm_d(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsm_d
        SUBROUTINE trsv_d(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsv_d
        SUBROUTINE trsm_z(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsm_z
        SUBROUTINE trsv_z(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsv_z
        SUBROUTINE trsm_s(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsm_s
        SUBROUTINE trsv_s(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: b(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsv_s
        SUBROUTINE trsm_c(t,b,side,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_side_type, &
           blas_trans_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsm_c
        SUBROUTINE trsv_c(t,b,uplo,trans,diag,alpha)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE trsv_c
      END INTERFACE

      INTERFACE syrk
        ! syrk  a and c rank 2
        ! syr  c rank 2; a rank 1
        ! c has shape (n,n)
        ! . For a, rank 2
        ! .   a has shape (n,k) if  trans = blas_no_trans (the default)
        ! .     has shape (k,n) if  trans /= blas_no_trans
        ! . For a, rank 1
        ! .   a has shape (n)
        SUBROUTINE syrk_d(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syrk_d
        SUBROUTINE syr_d(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr_d
        SUBROUTINE syrk_z(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syrk_z
        SUBROUTINE syr_z(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr_z
        SUBROUTINE syrk_s(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syrk_s
        SUBROUTINE syr_s(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr_s
        SUBROUTINE syrk_c(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syrk_c
        SUBROUTINE syr_c(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr_c
      END INTERFACE

      INTERFACE herk
        ! herk  a and c rank 2
        ! her  c rank 2; a rank 1
        ! c has shape (n,n)
        ! . For a, rank 2
        ! .   a has shape (n,k) if  trans = blas_no_trans (the default)
        ! .     has shape (k,n) if  trans /= blas_no_trans
        ! . For a, rank 1
        ! .   a has shape (n)
        SUBROUTINE herk_z(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE herk_z
        SUBROUTINE her_z(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her_z
        SUBROUTINE herk_c(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE herk_c
        SUBROUTINE her_c(a,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her_c
      END INTERFACE

      INTERFACE syr2k
        ! syr2k  a, b and c rank 2
        ! syr2  c rank 2; a and b rank 1
        ! c has shape (n,n)
        ! . For a, rank 2
        ! .   a has shape (n,k) if  trans = blas_no_trans (the default)
        ! .     has shape (k,n) if  trans /= blas_no_trans
        ! .   b has same shape as a
        ! . For a, rank 1
        ! .   a has shape (n)
        ! .   b has same shape as a
        SUBROUTINE syr2k_d(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2k_d
        SUBROUTINE syr2_d(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2_d
        SUBROUTINE syr2k_z(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2k_z
        SUBROUTINE syr2_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2_z
        SUBROUTINE syr2k_s(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2k_s
        SUBROUTINE syr2_s(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:), b(:)
          REAL (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2_s
        SUBROUTINE syr2k_c(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2k_c
        SUBROUTINE syr2_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE syr2_c
      END INTERFACE

      INTERFACE her2k
        ! her2k  a, b and c rank 2
        ! her2  c rank 2; a and b rank 1
        ! c has shape (n,n)
        ! . For a, rank 2
        ! .   a has shape (n,k) if  trans = blas_no_trans (the default)
        ! .     has shape (k,n) if  trans /= blas_no_trans
        ! .   b has same shape as a
        ! . For a, rank 1
        ! .   a has shape (n)
        ! .   b has same shape as a
        SUBROUTINE her2k_z(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: beta
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her2k_z
        SUBROUTINE her2_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: beta
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her2_z
        SUBROUTINE her2k_c(a,b,c,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: beta
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her2k_c
        SUBROUTINE her2_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:), b(:)
          COMPLEX (wp), INTENT (INOUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: beta
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE her2_c
      END INTERFACE

    END MODULE blas_extended_mat_mat_op
