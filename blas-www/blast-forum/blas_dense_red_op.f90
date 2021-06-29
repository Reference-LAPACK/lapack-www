    MODULE blas_dense_red_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.2 -- Reduction Operations
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: dot, norm, sum, min_val, max_val, amin_val, amax_val, sumsq

      INTERFACE dot
        ! x and y must have same shape
        SUBROUTINE dot_d(x,y,r,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (INOUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE dot_d
        SUBROUTINE dot_z(x,y,r,conj,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (INOUT) :: r
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE dot_z
        SUBROUTINE dot_s(x,y,r,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (INOUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE dot_s
        SUBROUTINE dot_c(x,y,r,conj,alpha,beta)
          USE blas_operator_arguments, ONLY : blas_conj_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (INOUT) :: r
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE dot_c
      END INTERFACE

      INTERFACE norm
        FUNCTION norm_d(x,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          REAL (wp) :: norm_d
        END FUNCTION norm_d
        FUNCTION norm_z(x,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          REAL (wp) :: norm_z
        END FUNCTION norm_z
        FUNCTION norm_s(x,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          REAL (wp) :: norm_s
        END FUNCTION norm_s
        FUNCTION norm_c(x,norm)
          USE blas_operator_arguments, ONLY : blas_norm_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          TYPE (blas_norm_type), INTENT (IN), OPTIONAL :: norm
          REAL (wp) :: norm_c
        END FUNCTION norm_c
      END INTERFACE

      INTERFACE sum
        FUNCTION sum_d(x)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp) :: sum_d
        END FUNCTION sum_d
        FUNCTION sum_z(x)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp) :: sum_z
        END FUNCTION sum_z
        FUNCTION sum_s(x)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp) :: sum_s
        END FUNCTION sum_s
        FUNCTION sum_c(x)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp) :: sum_c
        END FUNCTION sum_c
      END INTERFACE

      INTERFACE min_val
        SUBROUTINE min_val_d(x,k,r)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE min_val_d
        SUBROUTINE min_val_s(x,k,r)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE min_val_s
      END INTERFACE

      INTERFACE max_val
        SUBROUTINE max_val_d(x,k,r)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE max_val_d
        SUBROUTINE max_val_s(x,k,r)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE max_val_s
      END INTERFACE

      INTERFACE amin_val
        SUBROUTINE amin_val_d(x,k,r)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amin_val_d
        SUBROUTINE amin_val_z(x,k,r)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amin_val_z
        SUBROUTINE amin_val_s(x,k,r)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amin_val_s
        SUBROUTINE amin_val_c(x,k,r)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amin_val_c
      END INTERFACE

      INTERFACE amax_val
        SUBROUTINE amax_val_d(x,k,r)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amax_val_d
        SUBROUTINE amax_val_z(x,k,r)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amax_val_z
        SUBROUTINE amax_val_s(x,k,r)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amax_val_s
        SUBROUTINE amax_val_c(x,k,r)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: r
          INTEGER, INTENT (OUT) :: k
        END SUBROUTINE amax_val_c
      END INTERFACE

      INTERFACE sumsq
        SUBROUTINE sumsq_d(x,ssq,scl)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: ssq, scl
        END SUBROUTINE sumsq_d
        SUBROUTINE sumsq_z(x,ssq,scl)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: ssq, scl
        END SUBROUTINE sumsq_z
        SUBROUTINE sumsq_s(x,ssq,scl)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: ssq, scl
        END SUBROUTINE sumsq_s
        SUBROUTINE sumsq_c(x,ssq,scl)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: ssq, scl
        END SUBROUTINE sumsq_c
      END INTERFACE
    END MODULE blas_dense_red_op
