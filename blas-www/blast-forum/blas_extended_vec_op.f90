    MODULE blas_extended_vec_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 4.?.? -- Vector Operations
      ! .  using external procedures.
      ! Written on 18 December 2000
      ! .  Susan Blackford, UT-ICL

      PRIVATE
      PUBLIC :: axpby, waxpby

      INTERFACE axpby
        ! x and y have shape (n)
        SUBROUTINE axpby_d(x,y,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE axpby_d
        SUBROUTINE axpby_z(x,y,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE axpby_z
        SUBROUTINE axpby_s(x,y,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE axpby_s
        SUBROUTINE axpby_c(x,y,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE axpby_c
      END INTERFACE

      INTERFACE waxpby
        ! x, y and w have shape (n)
        SUBROUTINE waxpby_d(x,y,w,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (OUT) :: w(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE waxpby_d
        SUBROUTINE waxpby_z(x,y,w,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (OUT) :: w(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE waxpby_z
        SUBROUTINE waxpby_s(x,y,w,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (OUT) :: w(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE waxpby_s
        SUBROUTINE waxpby_c(x,y,w,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (OUT) :: w(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec
        END SUBROUTINE waxpby_c
      END INTERFACE

    END MODULE blas_extended_vec_op
