    MODULE blas_dense_mat_mov

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.9 -- Data Movement with Matrices
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: ge_copy, gb_copy, sy_copy, sb_copy, sp_copy, tr_copy, tb_copy, &
       tp_copy, he_copy, hb_copy, hp_copy, ge_trans, ge_permute

      INTERFACE ge_copy
        ! a has shape (m,n)
        ! b has same shape as a if trans = blas_no_trans
        ! . has shape (n,m)     if trans /= blas_no_trans
        SUBROUTINE ge_copy_d(a,b,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_copy_d
        SUBROUTINE ge_copy_z(a,b,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_copy_z
        SUBROUTINE ge_copy_s(a,b,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_copy_s
        SUBROUTINE ge_copy_c(a,b,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_copy_c
      END INTERFACE

      INTERFACE gb_copy
        ! a has shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        ! b has same shape as a     if trans = blas_no_trans
        ! . has shape ((kl+ku+1,m)  if trans /= blas_no_trans
        SUBROUTINE gb_copy_d(a,b,m,kl,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE gb_copy_d
        SUBROUTINE gb_copy_z(a,b,m,kl,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE gb_copy_z
        SUBROUTINE gb_copy_s(a,b,m,kl,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE gb_copy_s
        SUBROUTINE gb_copy_c(a,b,m,kl,trans)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE gb_copy_c
      END INTERFACE

      INTERFACE sy_copy
        ! a and b have shape (n,n)
        SUBROUTINE sy_copy_d(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_copy_d
        SUBROUTINE sy_copy_z(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_copy_z
        SUBROUTINE sy_copy_s(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_copy_s
        SUBROUTINE sy_copy_c(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_copy_c
      END INTERFACE

      INTERFACE sb_copy
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE sb_copy_d(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_copy_d
        SUBROUTINE sb_copy_z(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_copy_z
        SUBROUTINE sb_copy_s(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_copy_s
        SUBROUTINE sb_copy_c(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_copy_c
      END INTERFACE

      INTERFACE sp_copy
        ! ap and bp have shape ((n+1)*n/2)
        SUBROUTINE sp_copy_d(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_copy_d
        SUBROUTINE sp_copy_z(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_copy_z
        SUBROUTINE sp_copy_s(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_copy_s
        SUBROUTINE sp_copy_c(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_copy_c
      END INTERFACE

      INTERFACE tr_copy
        ! a and b have shape (n,n)
        SUBROUTINE tr_copy_d(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_copy_d
        SUBROUTINE tr_copy_z(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_copy_z
        SUBROUTINE tr_copy_s(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_copy_s
        SUBROUTINE tr_copy_c(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_copy_c
      END INTERFACE

      INTERFACE tb_copy
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE tb_copy_d(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_copy_d
        SUBROUTINE tb_copy_z(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_copy_z
        SUBROUTINE tb_copy_s(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_copy_s
        SUBROUTINE tb_copy_c(a,b,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_copy_c
      END INTERFACE

      INTERFACE tp_copy
        ! ap and bp have shape ((n+1)*n/2)
        SUBROUTINE tp_copy_d(ap,bp,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_copy_d
        SUBROUTINE tp_copy_z(ap,bp,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_copy_z
        SUBROUTINE tp_copy_s(ap,bp,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_copy_s
        SUBROUTINE tp_copy_c(ap,bp,uplo,trans,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type, &
           blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_copy_c
      END INTERFACE

      INTERFACE he_copy
        ! a and b have shape (n,n)
        SUBROUTINE he_copy_z(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE he_copy_z
        SUBROUTINE he_copy_c(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE he_copy_c
      END INTERFACE

      INTERFACE hb_copy
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE hb_copy_z(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hb_copy_z
        SUBROUTINE hb_copy_c(a,b,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (OUT) :: b(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hb_copy_c
      END INTERFACE

      INTERFACE hp_copy
        ! ap and bp have shape ((n+1)*n/2)
        SUBROUTINE hp_copy_z(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hp_copy_z
        SUBROUTINE hp_copy_c(ap,bp,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (OUT) :: bp(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hp_copy_c
      END INTERFACE

      INTERFACE ge_trans
        ! a has shape (n,n)
        SUBROUTINE ge_trans_d(a,conj)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
        END SUBROUTINE ge_trans_d
        SUBROUTINE ge_trans_z(a,conj)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
        END SUBROUTINE ge_trans_z
        SUBROUTINE ge_trans_s(a,conj)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
        END SUBROUTINE ge_trans_s
        SUBROUTINE ge_trans_c(a,conj)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
        END SUBROUTINE ge_trans_c
      END INTERFACE

      INTERFACE ge_permute
        ! a has shape (m,n)
        ! p has shape (m) if side = blas_left_side
        ! .           (n) if side = blas_right_side
        SUBROUTINE ge_permute_d(p,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          INTEGER, INTENT (IN) :: p(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_permute_d
        SUBROUTINE ge_permute_z(p,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          INTEGER, INTENT (IN) :: p(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_permute_z
        SUBROUTINE ge_permute_s(p,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          INTEGER, INTENT (IN) :: p(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_permute_s
        SUBROUTINE ge_permute_c(p,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          INTEGER, INTENT (IN) :: p(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_permute_c
      END INTERFACE
    END MODULE blas_dense_mat_mov
