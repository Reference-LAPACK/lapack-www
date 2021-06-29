    MODULE blas_dense_gen_trans

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.3 -- Generate Transformations
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: gen_grot, gen_jrot, gen_house

      INTERFACE gen_grot
        SUBROUTINE gen_grot_d(a,b,c,s,r)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: a, b
          REAL (wp), INTENT (OUT) :: c, s, r
        END SUBROUTINE gen_grot_d
        SUBROUTINE gen_grot_z(a,b,c,s,r)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (IN) :: a, b
          REAL (wp), INTENT (OUT) :: c
          COMPLEX (wp), INTENT (OUT) :: s, r
        END SUBROUTINE gen_grot_z
        SUBROUTINE gen_grot_s(a,b,c,s,r)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: a, b
          REAL (wp), INTENT (OUT) :: c, s, r
        END SUBROUTINE gen_grot_s
        SUBROUTINE gen_grot_c(a,b,c,s,r)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (IN) :: a, b
          REAL (wp), INTENT (OUT) :: c
          COMPLEX (wp), INTENT (OUT) :: s, r
        END SUBROUTINE gen_grot_c
      END INTERFACE

      INTERFACE gen_jrot
        SUBROUTINE gen_jrot_d(x,y,z,c,s,jrot)
          USE blas_operator_arguments, ONLY : blas_jrot_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: x, z
          REAL (wp), INTENT (IN) :: y
          REAL (wp), INTENT (OUT) :: c, s
          TYPE (blas_jrot_type), INTENT (IN), OPTIONAL :: jrot
        END SUBROUTINE gen_jrot_d
        SUBROUTINE gen_jrot_z(x,y,z,c,s,jrot)
          USE blas_operator_arguments, ONLY : blas_jrot_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: x, z
          COMPLEX (wp), INTENT (IN) :: y
          REAL (wp), INTENT (OUT) :: c
          COMPLEX (wp), INTENT (OUT) :: s
          TYPE (blas_jrot_type), INTENT (IN), OPTIONAL :: jrot
        END SUBROUTINE gen_jrot_z
        SUBROUTINE gen_jrot_s(x,y,z,c,s,jrot)
          USE blas_operator_arguments, ONLY : blas_jrot_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: x, z
          REAL (wp), INTENT (IN) :: y
          REAL (wp), INTENT (OUT) :: c, s
          TYPE (blas_jrot_type), INTENT (IN), OPTIONAL :: jrot
        END SUBROUTINE gen_jrot_s
        SUBROUTINE gen_jrot_c(x,y,z,c,s,jrot)
          USE blas_operator_arguments, ONLY : blas_jrot_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: x, z
          COMPLEX (wp), INTENT (IN) :: y
          REAL (wp), INTENT (OUT) :: c
          COMPLEX (wp), INTENT (OUT) :: s
          TYPE (blas_jrot_type), INTENT (IN), OPTIONAL :: jrot
        END SUBROUTINE gen_jrot_c
      END INTERFACE

      INTERFACE gen_house
        ! x has shape (n)
        SUBROUTINE gen_house_d(alpha,x,tau)
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (INOUT) :: alpha
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (OUT) :: tau
        END SUBROUTINE gen_house_d
        SUBROUTINE gen_house_z(alpha,x,tau)
          USE blas_precision, ONLY : wp => dp
          COMPLEX (wp), INTENT (INOUT) :: alpha
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (OUT) :: tau
        END SUBROUTINE gen_house_z
        SUBROUTINE gen_house_s(alpha,x,tau)
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (INOUT) :: alpha
          REAL (wp), INTENT (INOUT) :: x(:)
          REAL (wp), INTENT (OUT) :: tau
        END SUBROUTINE gen_house_s
        SUBROUTINE gen_house_c(alpha,x,tau)
          USE blas_precision, ONLY : wp => sp
          COMPLEX (wp), INTENT (INOUT) :: alpha
          COMPLEX (wp), INTENT (INOUT) :: x(:)
          COMPLEX (wp), INTENT (OUT) :: tau
        END SUBROUTINE gen_house_c
      END INTERFACE
    END MODULE blas_dense_gen_trans
