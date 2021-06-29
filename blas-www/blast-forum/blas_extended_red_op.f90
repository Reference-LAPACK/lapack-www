    MODULE blas_extended_red_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 4.?.? -- Reduction Operations
      ! .  using external procedures.
      ! Written on 18 December 2000
      ! .  Susan Blackford, UT-ICL

      PRIVATE
      PUBLIC :: dot, sum

      INTERFACE dot
        ! x and y must have same shape
        SUBROUTINE dot_d(x,y,r,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (INOUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
        END SUBROUTINE dot_d
        SUBROUTINE dot_z(x,y,r,conj,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_conj_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (INOUT) :: r
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
        END SUBROUTINE dot_z
        SUBROUTINE dot_s(x,y,r,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (INOUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
        END SUBROUTINE dot_s
        SUBROUTINE dot_c(x,y,r,conj,alpha,beta,prec)
          USE blas_operator_arguments, ONLY : blas_conj_type, &
                                              blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (INOUT) :: r
          TYPE (blas_conj_type), INTENT (IN), OPTIONAL :: conj
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
        END SUBROUTINE dot_c
      END INTERFACE

      INTERFACE sum
        FUNCTION sum_d(x,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
          REAL (wp) :: sum_d
        END FUNCTION sum_d
        FUNCTION sum_z(x,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
          COMPLEX (wp) :: sum_z
        END FUNCTION sum_z
        FUNCTION sum_s(x,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
          REAL (wp) :: sum_s
        END FUNCTION sum_s
        FUNCTION sum_c(x,prec)
          USE blas_operator_arguments, ONLY : blas_prec_type
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          TYPE (blas_prec_type), INTENT (IN), OPTIONAL :: prec        
          COMPLEX (wp) :: sum_c
        END FUNCTION sum_c
      END INTERFACE

    END MODULE blas_extended_red_op
