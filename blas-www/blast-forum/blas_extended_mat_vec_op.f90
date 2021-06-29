    MODULE blas_dense_mat_vec_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 4.?.? -- Matrix-Vector Operations
      ! .  using external procedures.
      ! Written on 18 December 2000
      ! .  Susan Blackford, UT-ICL

      PRIVATE
      PUBLIC :: gbmv, sbmv, spmv, hbmv, hpmv, tbmv, tpmv, &
        ge_sum_mv, tbsv, tpsv

      INTERFACE gbmv
        ! m is not needed as an argument
        ! if trans = blas_no_trans then
        ! . x has shape (n)
        ! . y has shape (m)
        ! . a has shape (kl+ku+1,n), SIZE(a,1) > kl
        ! if trans =/ blas_no_trans then
        ! . x has shape (m)
        ! . y has shape (n)
        ! . a has shape (kl+ku+1,m), SIZE(a,1) > kl
        SUBROUTINE gbmv_d(a,m,kl,x,y,trans,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_trans_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), x(:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE gbmv_d
        SUBROUTINE gbmv_z(a,m,kl,x,y,trans,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_trans_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE gbmv_z
        SUBROUTINE gbmv_s(a,m,kl,x,y,trans,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_trans_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), x(:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE gbmv_s
        SUBROUTINE gbmv_c(a,m,kl,x,y,trans,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_trans_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE gbmv_c
      END INTERFACE

      INTERFACE sbmv
        ! x and y have shape (n)
        ! a has shape (k+1,n), (k=band width = SIZE(a,1)-1)
        SUBROUTINE sbmv_d(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE sbmv_d
        SUBROUTINE sbmv_z(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE sbmv_z
        SUBROUTINE sbmv_s(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE sbmv_s
        SUBROUTINE sbmv_c(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE sbmv_c
      END INTERFACE

      INTERFACE spmv
        ! x and y have shape (n)
        ! ap has shape ((n+1)*n/2)
        SUBROUTINE spmv_d(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:), x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE spmv_d
        SUBROUTINE spmv_z(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE spmv_z
        SUBROUTINE spmv_s(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:), x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE spmv_s
        SUBROUTINE spmv_c(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE spmv_c
      END INTERFACE

      INTERFACE hbmv_a_x
        ! ???????
        ! x and y have shape (n)
        ! a has shape (k+1,n), (k=band width = SIZE(a,1)-1)
        SUBROUTINE hbmv_z(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE hbmv_z
        SUBROUTINE hbmv_c(a,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE hbmv_c
      END INTERFACE

      INTERFACE hpmv_a_x
        ! ???????
        ! x and y have shape (n)
        ! ap has shape ((n+1)*n/2)
        SUBROUTINE hpmv_z(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE hpmv_z
        SUBROUTINE hpmv_c(ap,x,y,uplo,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:), x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE hpmv_c
      END INTERFACE

      INTERFACE tbmv
        ! x has shape (n)
        ! t has shape (k+1,n), (k=band width = SIZE(t,1)-1)
        SUBROUTINE tbmv_d(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbmv_d
        SUBROUTINE tbmv_z(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbmv_z
        SUBROUTINE tbmv_s(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbmv_s
        SUBROUTINE tbmv_c(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbmv_c
      END INTERFACE

      INTERFACE tpmv
        ! x has shape (n)
        ! tp has shape ((n+1)*n/2)
        SUBROUTINE tpmv_d(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: tp(:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpmv_d
        SUBROUTINE tpmv_z(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: tp(:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpmv_z
        SUBROUTINE tpmv_s(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: tp(:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpmv_s
        SUBROUTINE tpmv_c(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: tp(:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpmv_c
      END INTERFACE

      INTERFACE ge_sum_mv
        ! x has shape (n)
        ! y has shape (m)
        ! a and b have shape (m,n)
        SUBROUTINE ge_sum_mv_d(a,x,b,y,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:), x(:)
          REAL (wp), INTENT (OUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE ge_sum_mv_d
        SUBROUTINE ge_sum_mv_z(a,x,b,y,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:), x(:)
          COMPLEX (wp), INTENT (OUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE ge_sum_mv_z
        SUBROUTINE ge_sum_mv_s(a,x,b,y,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:), x(:)
          REAL (wp), INTENT (OUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE ge_sum_mv_s
        SUBROUTINE ge_sum_mv_c(a,x,b,y,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:), x(:)
          COMPLEX (wp), INTENT (OUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE ge_sum_mv_c
      END INTERFACE

      INTERFACE tbsv
        ! x has shape (n)
        ! t has shape (k+1,n), (k=band width = SIZE(t,1)-1)
        SUBROUTINE tbsv_d(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbsv_d
        SUBROUTINE tbsv_z(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbsv_z
        SUBROUTINE tbsv_s(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: t(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbsv_s
        SUBROUTINE tbsv_c(t,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: t(:,:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tbsv_c
      END INTERFACE

      INTERFACE tpsv
        ! x has shape (n)
        ! tp has shape ((n+1)*n/2)
        SUBROUTINE tpsv_d(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: tp(:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpsv_d
        SUBROUTINE tpsv_z(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: tp(:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpsv_z
        SUBROUTINE tpsv_s(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: tp(:)
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpsv_s
        SUBROUTINE tpsv_c(tp,x,uplo,trans,diag,alpha,prec)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type, &
           blas_diag_type, blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: tp(:)
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
          TYPE (blas_prec_type), INTENT(IN), OPTIONAL :: prec
        END SUBROUTINE tpsv_c
      END INTERFACE

    END MODULE blas_extended_mat_vec_op
