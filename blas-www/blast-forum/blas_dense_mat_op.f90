    MODULE blas_dense_mat_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.7 -- Matrix Operations
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: ge_norm, gb_norm, sy_norm, he_norm, sb_norm, hb_norm, sp_norm, &
       hp_norm, tr_norm, tb_norm, tp_norm, ge_diag_scale, gb_diag_scale, &
       ge_lrscale, gb_lrscale, sy_lrscale, sb_lrscale, sp_lrscale, he_lrscale, &
       hb_lrscale, hp_lrscale, ge_diag_scale_acc, gb_diag_scale_acc, ge_acc, &
       gb_acc, sy_acc, sb_acc, sp_acc, tr_acc, tb_acc, tp_acc, ge_add, gb_add, &
       sy_add, sb_add, sp_add, tr_add, tb_add, tp_add

      INTERFACE ge_norm
        ! a has shape (m,n)
        FUNCTION ge_norm_d(a,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: ge_norm_d
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION ge_norm_d
        FUNCTION ge_norm_z(a,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: ge_norm_z
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION ge_norm_z
        FUNCTION ge_norm_s(a,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: ge_norm_s
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION ge_norm_s
        FUNCTION ge_norm_c(a,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: ge_norm_c
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION ge_norm_c
      END INTERFACE

      INTERFACE gb_norm
        ! a has shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        FUNCTION gb_norm_d(a,m,kl,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp) :: gb_norm_d
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION gb_norm_d
        FUNCTION gb_norm_z(a,m,kl,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp) :: gb_norm_z
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION gb_norm_z
        FUNCTION gb_norm_s(a,m,kl,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp) :: gb_norm_s
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION gb_norm_s
        FUNCTION gb_norm_c(a,m,kl,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp) :: gb_norm_c
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION gb_norm_c
      END INTERFACE

      INTERFACE sy_norm
        ! a has shape (n,n)
        FUNCTION sy_norm_d(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sy_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sy_norm_d
        FUNCTION sy_norm_z(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sy_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sy_norm_z
        FUNCTION sy_norm_s(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sy_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sy_norm_s
        FUNCTION sy_norm_c(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sy_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sy_norm_c
      END INTERFACE

      INTERFACE he_norm
        ! a has shape (n,n)
        FUNCTION he_norm_z(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: he_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION he_norm_z
        FUNCTION he_norm_c(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: he_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION he_norm_c
      END INTERFACE

      INTERFACE sb_norm
        ! a has shape (k+1,n), k=band width=size(a,1)-1
        FUNCTION sb_norm_d(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sb_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sb_norm_d
        FUNCTION sb_norm_z(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sb_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sb_norm_z
        FUNCTION sb_norm_s(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sb_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sb_norm_s
        FUNCTION sb_norm_c(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: sb_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sb_norm_c
      END INTERFACE

      INTERFACE hb_norm
        ! a has shape (k+1,n), k=band width=size(a,1)-1
        FUNCTION hb_norm_z(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: hb_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION hb_norm_z
        FUNCTION hb_norm_c(a,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: hb_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION hb_norm_c
      END INTERFACE

      INTERFACE sp_norm
        ! ap has shape ((n*1)*n/2)
        FUNCTION sp_norm_d(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: sp_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sp_norm_d
        FUNCTION sp_norm_z(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: sp_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sp_norm_z
        FUNCTION sp_norm_s(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: sp_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sp_norm_s
        FUNCTION sp_norm_c(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: sp_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION sp_norm_c
      END INTERFACE

      INTERFACE hp_norm
        ! ap has shape ((n*1)*n/2)
        FUNCTION hp_norm_z(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: hp_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION hp_norm_z
        FUNCTION hp_norm_c(ap,norm,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: hp_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
        END FUNCTION hp_norm_c
      END INTERFACE

      INTERFACE tr_norm
        ! a has shape (n,n)
        FUNCTION tr_norm_d(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tr_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tr_norm_d
        FUNCTION tr_norm_z(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tr_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tr_norm_z
        FUNCTION tr_norm_s(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tr_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tr_norm_s
        FUNCTION tr_norm_c(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tr_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tr_norm_c
      END INTERFACE

      INTERFACE tb_norm
        ! a has shape (k+1,n), k=band width=size(a,1)-1
        FUNCTION tb_norm_d(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tb_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tb_norm_d
        FUNCTION tb_norm_z(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tb_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tb_norm_z
        FUNCTION tb_norm_s(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tb_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tb_norm_s
        FUNCTION tb_norm_c(a,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          REAL (wp) :: tb_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tb_norm_c
      END INTERFACE

      INTERFACE tp_norm
        ! ap has shape ((n+1)*n/2)
        FUNCTION tp_norm_d(ap,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: tp_norm_d
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tp_norm_d
        FUNCTION tp_norm_z(ap,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: tp_norm_z
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tp_norm_z
        FUNCTION tp_norm_s(ap,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: tp_norm_s
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tp_norm_s
        FUNCTION tp_norm_c(ap,norm,uplo,diag)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_norm_type, &
           blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          REAL (wp) :: tp_norm_c
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END FUNCTION tp_norm_c
      END INTERFACE

      INTERFACE ge_diag_scale
        ! a has shape (m,n)
        ! d has shape (m), if side = blas_left_side
        ! .     shape (n), if side = blas_right_side
        SUBROUTINE ge_diag_scale_d(d,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_diag_scale_d
        SUBROUTINE ge_diag_scale_z(d,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_diag_scale_z
        SUBROUTINE ge_diag_scale_s(d,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_diag_scale_s
        SUBROUTINE ge_diag_scale_c(d,a,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE ge_diag_scale_c
      END INTERFACE

      INTERFACE gb_diag_scale
        ! a has shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        ! d has shape (m), if side = blas_left_side
        ! .     shape (n), if side = blas_right_side
        SUBROUTINE gb_diag_scale_d(d,a,m,kl,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE gb_diag_scale_d
        SUBROUTINE gb_diag_scale_z(d,a,m,kl,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE gb_diag_scale_z
        SUBROUTINE gb_diag_scale_s(d,a,m,kl,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE gb_diag_scale_s
        SUBROUTINE gb_diag_scale_c(d,a,m,kl,side)
          USE blas_operator_arguments, ONLY : blas_side_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          TYPE (blas_side_type), INTENT (IN), OPTIONAL :: side
        END SUBROUTINE gb_diag_scale_c
      END INTERFACE

      INTERFACE ge_lrscale
        ! a has shape (m,n)
        ! dl has shape (m)
        ! dr has shape (n)
        SUBROUTINE ge_lrscale_d(dl,dr,a)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: dl(:), dr(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_lrscale_d
        SUBROUTINE ge_lrscale_z(dl,dr,a)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: dl(:), dr(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_lrscale_z
        SUBROUTINE ge_lrscale_s(dl,dr,a)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: dl(:), dr(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_lrscale_s
        SUBROUTINE ge_lrscale_c(dl,dr,a)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: dl(:), dr(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_lrscale_c
      END INTERFACE

      INTERFACE gb_lrscale
        ! a has shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        ! dl has shape (m)
        ! dr has shape (n)
        SUBROUTINE gb_lrscale_d(dl,dr,a,m,kl)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: dl(:), dr(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_lrscale_d
        SUBROUTINE gb_lrscale_z(dl,dr,a,m,kl)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: dl(:), dr(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_lrscale_z
        SUBROUTINE gb_lrscale_s(dl,dr,a,m,kl)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: dl(:), dr(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_lrscale_s
        SUBROUTINE gb_lrscale_c(dl,dr,a,m,kl)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: dl(:), dr(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_lrscale_c
      END INTERFACE

      INTERFACE sy_lrscale
        ! a has shape (n,n)
        ! d has shape (n)
        SUBROUTINE sy_lrscale_d(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_lrscale_d
        SUBROUTINE sy_lrscale_z(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_lrscale_z
        SUBROUTINE sy_lrscale_s(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_lrscale_s
        SUBROUTINE sy_lrscale_c(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_lrscale_c
      END INTERFACE

      INTERFACE sb_lrscale
        ! a has shape (k+1,n), k=band width=size(a,1)-1
        ! d has shape (n)
        SUBROUTINE sb_lrscale_d(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_lrscale_d
        SUBROUTINE sb_lrscale_z(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_lrscale_z
        SUBROUTINE sb_lrscale_s(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_lrscale_s
        SUBROUTINE sb_lrscale_c(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_lrscale_c
      END INTERFACE

      INTERFACE sp_lrscale
        ! ap has shape ((n+1)*n/2)
        ! d has shape (n)
        SUBROUTINE sp_lrscale_d(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_lrscale_d
        SUBROUTINE sp_lrscale_z(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_lrscale_z
        SUBROUTINE sp_lrscale_s(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: d(:)
          REAL (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_lrscale_s
        SUBROUTINE sp_lrscale_c(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_lrscale_c
      END INTERFACE

      INTERFACE he_lrscale
        ! a has shape (n,n)
        ! d has shape (n)
        SUBROUTINE he_lrscale_z(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE he_lrscale_z
        SUBROUTINE he_lrscale_c(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE he_lrscale_c
      END INTERFACE

      INTERFACE hb_lrscale
        ! a has shape (k+1,n), k=band width=size(a,1)-1
        ! d has shape (n)
        SUBROUTINE hb_lrscale_z(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hb_lrscale_z
        SUBROUTINE hb_lrscale_c(d,a,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hb_lrscale_c
      END INTERFACE

      INTERFACE hp_lrscale
        ! ap has shape ((n+1)*n/2)
        ! d has shape (n)
        SUBROUTINE hp_lrscale_z(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hp_lrscale_z
        SUBROUTINE hp_lrscale_c(d,ap,uplo)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: d(:)
          COMPLEX (wp), INTENT (INOUT) :: ap(:)
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE hp_lrscale_c
      END INTERFACE

      INTERFACE ge_diag_scale_acc
        ! a and b have shape (m,n)
        ! d has shape (n)
        SUBROUTINE ge_diag_scale_acc_d(b,d,a)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: b(:,:), d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_diag_scale_acc_d
        SUBROUTINE ge_diag_scale_acc_z(b,d,a)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: b(:,:), d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_diag_scale_acc_z
        SUBROUTINE ge_diag_scale_acc_s(b,d,a)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: b(:,:), d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_diag_scale_acc_s
        SUBROUTINE ge_diag_scale_acc_c(b,d,a)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: b(:,:), d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
        END SUBROUTINE ge_diag_scale_acc_c
      END INTERFACE

      INTERFACE gb_diag_scale_acc
        ! a has shape (m,n)
        ! b has shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        ! d has shape (n)
        SUBROUTINE gb_diag_scale_acc_d(b,m,kl,d,a)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: b(:,:), d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_diag_scale_acc_d
        SUBROUTINE gb_diag_scale_acc_z(b,m,kl,d,a)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: b(:,:), d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_diag_scale_acc_z
        SUBROUTINE gb_diag_scale_acc_s(b,m,kl,d,a)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: b(:,:), d(:)
          REAL (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_diag_scale_acc_s
        SUBROUTINE gb_diag_scale_acc_c(b,m,kl,d,a)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: b(:,:), d(:)
          COMPLEX (wp), INTENT (INOUT) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
        END SUBROUTINE gb_diag_scale_acc_c
      END INTERFACE

      INTERFACE ge_acc
        ! a and b have shape (m,n)
        SUBROUTINE ge_acc_d(a,b,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_acc_d
        SUBROUTINE ge_acc_z(a,b,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_acc_z
        SUBROUTINE ge_acc_s(a,b,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_acc_s
        SUBROUTINE ge_acc_c(a,b,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE ge_acc_c
      END INTERFACE

      INTERFACE gb_acc
        ! a and b have shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        SUBROUTINE gb_acc_d(a,m,kl,b,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_acc_d
        SUBROUTINE gb_acc_z(a,m,kl,b,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_acc_z
        SUBROUTINE gb_acc_s(a,m,kl,b,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_acc_s
        SUBROUTINE gb_acc_c(a,m,kl,b,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_acc_c
      END INTERFACE

      INTERFACE sy_acc
        ! a and b have shape (n,n)
        SUBROUTINE sy_acc_d(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sy_acc_d
        SUBROUTINE sy_acc_z(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sy_acc_z
        SUBROUTINE sy_acc_s(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sy_acc_s
        SUBROUTINE sy_acc_c(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sy_acc_c
      END INTERFACE

      INTERFACE sb_acc
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE sb_acc_d(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sb_acc_d
        SUBROUTINE sb_acc_z(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sb_acc_z
        SUBROUTINE sb_acc_s(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sb_acc_s
        SUBROUTINE sb_acc_c(a,b,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sb_acc_c
      END INTERFACE

      INTERFACE sp_acc
        ! ap and bp have shape ((n+1)*n/2)
        SUBROUTINE sp_acc_d(ap,bp,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (INOUT) :: bp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sp_acc_d
        SUBROUTINE sp_acc_z(ap,bp,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (INOUT) :: bp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sp_acc_z
        SUBROUTINE sp_acc_s(ap,bp,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (INOUT) :: bp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sp_acc_s
        SUBROUTINE sp_acc_c(ap,bp,uplo,trans,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_trans_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (INOUT) :: bp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_trans_type), INTENT (IN), OPTIONAL :: trans
        END SUBROUTINE sp_acc_c
      END INTERFACE

      INTERFACE tr_acc
        ! a and b have shape (n,n)
        SUBROUTINE tr_acc_d(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_acc_d
        SUBROUTINE tr_acc_z(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_acc_z
        SUBROUTINE tr_acc_s(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_acc_s
        SUBROUTINE tr_acc_c(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_acc_c
      END INTERFACE

      INTERFACE tb_acc
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE tb_acc_d(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_acc_d
        SUBROUTINE tb_acc_z(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_acc_z
        SUBROUTINE tb_acc_s(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: b(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_acc_s
        SUBROUTINE tb_acc_c(a,b,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:)
          COMPLEX (wp), INTENT (INOUT) :: b(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_acc_c
      END INTERFACE

      INTERFACE tp_acc
        ! a and b have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE tp_acc_d(ap,bp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (INOUT) :: bp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_acc_d
        SUBROUTINE tp_acc_z(ap,bp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (INOUT) :: bp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_acc_z
        SUBROUTINE tp_acc_s(ap,bp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:)
          REAL (wp), INTENT (INOUT) :: bp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_acc_s
        SUBROUTINE tp_acc_c(ap,bp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:)
          COMPLEX (wp), INTENT (INOUT) :: bp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_acc_c
      END INTERFACE

      INTERFACE ge_add
        ! a, b and c have shape (m,n)
        SUBROUTINE ge_add_d(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ge_add_d
        SUBROUTINE ge_add_z(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ge_add_z
        SUBROUTINE ge_add_s(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ge_add_s
        SUBROUTINE ge_add_c(a,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE ge_add_c
      END INTERFACE

      INTERFACE gb_add
        ! a, b and c have shape (kl+ku+1,n), ku=SIZE(a,1)-1-kl
        SUBROUTINE gb_add_d(a,m,kl,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_add_d
        SUBROUTINE gb_add_z(a,m,kl,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_add_z
        SUBROUTINE gb_add_s(a,m,kl,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          INTEGER, INTENT (IN) :: m, kl
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_add_s
        SUBROUTINE gb_add_c(a,m,kl,b,c,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          INTEGER, INTENT (IN) :: m, kl
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE gb_add_c
      END INTERFACE

      INTERFACE sy_add
        ! a, b and c have shape (n,n)
        SUBROUTINE sy_add_d(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_add_d
        SUBROUTINE sy_add_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_add_z
        SUBROUTINE sy_add_s(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_add_s
        SUBROUTINE sy_add_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sy_add_c
      END INTERFACE

      INTERFACE sb_add
        ! a, b and c have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE sb_add_d(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_add_d
        SUBROUTINE sb_add_z(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_add_z
        SUBROUTINE sb_add_s(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_add_s
        SUBROUTINE sb_add_c(a,b,c,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sb_add_c
      END INTERFACE

      INTERFACE sp_add
        ! ap, bp and cp have shape ((n+1)*n/2)
        SUBROUTINE sp_add_d(ap,bp,cp,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:), bp(:)
          REAL (wp), INTENT (OUT) :: cp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_add_d
        SUBROUTINE sp_add_z(ap,bp,cp,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:), bp(:)
          COMPLEX (wp), INTENT (OUT) :: cp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_add_z
        SUBROUTINE sp_add_s(ap,bp,cp,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:), bp(:)
          REAL (wp), INTENT (OUT) :: cp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_add_s
        SUBROUTINE sp_add_c(ap,bp,cp,uplo,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:), bp(:)
          COMPLEX (wp), INTENT (OUT) :: cp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
        END SUBROUTINE sp_add_c
      END INTERFACE

      INTERFACE tr_add
        ! a, b and c have shape (n,n)
        SUBROUTINE tr_add_d(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_add_d
        SUBROUTINE tr_add_z(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_add_z
        SUBROUTINE tr_add_s(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_add_s
        SUBROUTINE tr_add_c(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tr_add_c
      END INTERFACE

      INTERFACE tb_add
        ! a, b and c have shape (k+1,n), k=band width=size(a,1)-1
        SUBROUTINE tb_add_d(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_add_d
        SUBROUTINE tb_add_z(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_add_z
        SUBROUTINE tb_add_s(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a(:,:), b(:,:)
          REAL (wp), INTENT (OUT) :: c(:,:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_add_s
        SUBROUTINE tb_add_c(a,b,c,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a(:,:), b(:,:)
          COMPLEX (wp), INTENT (OUT) :: c(:,:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tb_add_c
      END INTERFACE

      INTERFACE tp_add
        ! ap, bp and cp have shape ((n+1)*n/2)
        SUBROUTINE tp_add_d(ap,bp,cp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: ap(:), bp(:)
          REAL (wp), INTENT (OUT) :: cp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_add_d
        SUBROUTINE tp_add_z(ap,bp,cp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: ap(:), bp(:)
          COMPLEX (wp), INTENT (OUT) :: cp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_add_z
        SUBROUTINE tp_add_s(ap,bp,cp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: ap(:), bp(:)
          REAL (wp), INTENT (OUT) :: cp(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_add_s
        SUBROUTINE tp_add_c(ap,bp,cp,uplo,diag,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_uplo_type, blas_diag_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: ap(:), bp(:)
          COMPLEX (wp), INTENT (OUT) :: cp(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_uplo_type), INTENT (IN), OPTIONAL :: uplo
          TYPE (blas_diag_type), INTENT (IN), OPTIONAL :: diag
        END SUBROUTINE tp_add_c
      END INTERFACE
    END MODULE blas_dense_mat_op
