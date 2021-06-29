    MODULE blas_dense_vec_mov

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.5 -- Data Movement with Vectors
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: copy, swap, sortv, permute

      INTERFACE copy
        ! x and y have shape (n)
        SUBROUTINE copy_d(x,y)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: y(:)
        END SUBROUTINE copy_d
        SUBROUTINE copy_z(x,y)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (OUT) :: y(:)
        END SUBROUTINE copy_z
        SUBROUTINE copy_s(x,y)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: x(:)
          REAL (wp), INTENT (OUT) :: y(:)
        END SUBROUTINE copy_s
        SUBROUTINE copy_c(x,y)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: x(:)
          COMPLEX (wp), INTENT (OUT) :: y(:)
        END SUBROUTINE copy_c
      END INTERFACE

      INTERFACE swap
        ! x and y have shape (n)
        SUBROUTINE swap_d(x,y)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE swap_d
        SUBROUTINE swap_z(x,y)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE swap_z
        SUBROUTINE swap_s(x,y)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE swap_s
        SUBROUTINE swap_c(x,y)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (INOUT) :: x(:), y(:)
        END SUBROUTINE swap_c
      END INTERFACE

      INTERFACE sortv
        ! x and p have shape (n)
        SUBROUTINE sortv_d(x,sort,p)
          USE blas_operator_arguments, ONLY : blas_sort_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: x(:)
          TYPE (blas_sort_type), INTENT (IN), OPTIONAL :: sort
          INTEGER, INTENT (OUT), OPTIONAL :: p(:)
        END SUBROUTINE sortv_d
        SUBROUTINE sortv_s(x,sort,p)
          USE blas_operator_arguments, ONLY : blas_sort_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: x(:)
          TYPE (blas_sort_type), INTENT (IN), OPTIONAL :: sort
          INTEGER, INTENT (OUT), OPTIONAL :: p(:)
        END SUBROUTINE sortv_s
      END INTERFACE

      INTERFACE permute
        ! x and p have shape (n)
        SUBROUTINE permute_d(x,p)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: x(:)
          INTEGER, INTENT (IN) :: p(:)
        END SUBROUTINE permute_d
        SUBROUTINE permute_z(x,p)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          INTEGER, INTENT (IN) :: p(:)
        END SUBROUTINE permute_z
        SUBROUTINE permute_s(x,p)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: x(:)
          INTEGER, INTENT (IN) :: p(:)
        END SUBROUTINE permute_s
        SUBROUTINE permute_c(x,p)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          INTEGER, INTENT (IN) :: p(:)
        END SUBROUTINE permute_c
      END INTERFACE
    END MODULE blas_dense_vec_mov
