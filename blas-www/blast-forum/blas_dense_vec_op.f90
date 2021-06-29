    MODULE blas_dense_vec_op

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.4 -- Vector Operations
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: rscale, axpby, waxpby, axpy_dot, apply_grot

      INTERFACE rscale
        ! x has shape (n)
        SUBROUTINE rscale_d(alpha,x)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: alpha
          REAL (wp), INTENT (INOUT) :: x(:)
        END SUBROUTINE rscale_d
        SUBROUTINE rscale_z(alpha,x)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: alpha
          COMPLEX (wp), INTENT (INOUT) :: x(:)
        END SUBROUTINE rscale_z
        SUBROUTINE rscale_s(alpha,x)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: alpha
          REAL (wp), INTENT (INOUT) :: x(:)
        END SUBROUTINE rscale_s
        SUBROUTINE rscale_c(alpha,x)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: alpha
          COMPLEX (wp), INTENT (INOUT) :: x(:)
        END SUBROUTINE rscale_c
      END INTERFACE

      INTERFACE axpby
        ! x and y have shape (n)
        SUBROUTINE axpby_d(x,y,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE axpby_d
        SUBROUTINE axpby_z(x,y,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE axpby_z
        SUBROUTINE axpby_s(x,y,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (INOUT) :: y(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE axpby_s
        SUBROUTINE axpby_c(x,y,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (INOUT) :: y(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE axpby_c
      END INTERFACE

      INTERFACE waxpby
        ! x, y and w have shape (n)
        SUBROUTINE waxpby_d(x,y,w,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (OUT) :: w(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE waxpby_d
        SUBROUTINE waxpby_z(x,y,w,alpha,beta)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (OUT) :: w(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE waxpby_z
        SUBROUTINE waxpby_s(x,y,w,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:), y(:)
          REAL (wp), INTENT (OUT) :: w(:)
          REAL (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE waxpby_s
        SUBROUTINE waxpby_c(x,y,w,alpha,beta)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:), y(:)
          COMPLEX (wp), INTENT (OUT) :: w(:)
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha, beta
        END SUBROUTINE waxpby_c
      END INTERFACE

      INTERFACE axpy_dot
        ! u, v and w have shape (n)
        SUBROUTINE axpy_dot_d(w,v,u,r,alpha)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: u(:), v(:)
          REAL (wp), INTENT (INOUT) :: w(:)
          REAL (wp), INTENT (OUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
        END SUBROUTINE axpy_dot_d
        SUBROUTINE axpy_dot_z(w,v,u,r,alpha)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: u(:), v(:)
          COMPLEX (wp), INTENT (INOUT) :: w(:)
          COMPLEX (wp), INTENT (OUT) :: r
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
        END SUBROUTINE axpy_dot_z
        SUBROUTINE axpy_dot_s(w,v,u,r,alpha)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: u(:), v(:)
          REAL (wp), INTENT (INOUT) :: w(:)
          REAL (wp), INTENT (OUT) :: r
          REAL (wp), INTENT (IN), OPTIONAL :: alpha
        END SUBROUTINE axpy_dot_s
        SUBROUTINE axpy_dot_c(w,v,u,r,alpha)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: u(:), v(:)
          COMPLEX (wp), INTENT (INOUT) :: w(:)
          COMPLEX (wp), INTENT (OUT) :: r
          COMPLEX (wp), INTENT (IN), OPTIONAL :: alpha
        END SUBROUTINE axpy_dot_c
      END INTERFACE

      INTERFACE apply_grot
        ! x and y have shape (n)
        SUBROUTINE apply_grot_d(c,s,x,y)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: c, s
          REAL (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE apply_grot_d
        SUBROUTINE apply_grot_z(c,s,x,y)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: c
          COMPLEX (wp), INTENT (IN) :: s
          COMPLEX (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE apply_grot_z
        SUBROUTINE apply_grot_s(c,s,x,y)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: c, s
          REAL (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE apply_grot_s
        SUBROUTINE apply_grot_c(c,s,x,y)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: c
          COMPLEX (wp), INTENT (IN) :: s
          COMPLEX (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE apply_grot_c
      END INTERFACE
    END MODULE blas_dense_vec_op
