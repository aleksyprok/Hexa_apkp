MODULE tridiag

  USE shared_data

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE define_tridiag(u)
  ! Define the elements of the tridiagonal matrix, A, which represents the system of 
  ! equations which need to be solved i.e. A * v = d

  ! Input: u = temperature(:,n-1)

  ! Output: a, b, c = vectors defining the tridiagonal matrix

    REAL(num), DIMENSION(0:), INTENT(INOUT) :: u
    REAL(num) :: um, up

    DO ix = 0, nx

      ixm = ix - 1
      ixp = ix + 1

      um = (0.5_num * (u(ix ) + u(ixm))) ** 2.5_num
      up = (0.5_num * (u(ixp) + u(ix ))) ** 2.5_num

      ! a(ix) = -kappa / dx ** 2.0_num
      ! b(ix) = 1.0_num / dt + 2.0_num * kappa / dx ** 2.0_num
      ! c(ix) = -kappa / dx ** 2.0_num
      a(ix) = - um / dx ** 2.0_num
      b(ix) = 1.0_num / dt + (up + um) / dx ** 2.0_num
      c(ix) = - up / dx ** 2.0_num

    END DO

  END SUBROUTINE

  SUBROUTINE tridiagLU
  ! Obtain the the LU factorisation of a tridiagonal matrix
  
  ! Input: a, b, c = vectors defining the tridiagonal matrix. a is 
  !          the subdiagonal, b is the main diagonal and c is the superdiagonal
  
  ! Output: e, f = vectors defining the L and U factors of the tridiagonal matrix

    e(0) = b(0)
    f(0) = c(0) / b(0)
    DO ix = 1, nx
      ixm = ix - 1
      e(ix) = b(ix) - a(ix) * f(ixm)
      f(ix) = c(ix) / e(ix)
    END DO

  END SUBROUTINE tridiagLU

  SUBROUTINE tridiagLUsolve(v)
  ! Solve L * (U * v) = d where L and U are LU factors of a tridiagonal matrix

  ! Input: d    = right hand side of a system of equations
  !        e, f = vectors defining L and U factors of the tridiagonal matrix
  
  ! Output: v = solution vector at current time step

    REAL(num), DIMENSION(0:), INTENT(INOUT) :: v

    ! Forward substitution to solve L * w = d
    v(0) = d(0) / e(0)
    DO ix = 1, nx
      ixm = ix - 1
      v(ix) = (d(ix) - a(ix) * v(ixm)) / e(ix)
    END DO

    ! Backward substitution to solve U * v = w
    DO ix = nx - 1, 0, -1
      ixp = ix + 1
      v(ix) = v(ix) - f(ix) * v(ixp)
    END DO

  END SUBROUTINE tridiagLUsolve

END MODULE tridiag