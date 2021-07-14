MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

    length_cm = 2.95436E9_num
    time_s = 5760.0_num

    nx = 512
    nz = 512

    delx = 6.0_num / REAL(nx, num)
    delz = 6.0_num / REAL(nz + 1, num)

    ALLOCATE(xc(0:nx+1))
    ALLOCATE(zc(0:nz+1))
    ALLOCATE(xb(nx+1))
    ALLOCATE(zb(nz+1))

    DO ix = 0, nx + 1
      xc(ix) = ix * delx - 0.5_num * delx
    END DO
    DO ix = 1, nx + 1
      xb(ix) = (ix - 1) * delx
    END DO

    DO iz = 0, nz + 1
      zc(iz) = iz * delz + 0.5_num * delz
    END DO
    DO iz = 1, nz + 1
      zb(iz) = iz * delz
    END DO

  END SUBROUTINE grid_setup

  FUNCTION bx_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: bx_nlff

    bx_nlff = -(l / kx) * SIN(kx * x) * COSH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION bx_nlff

  FUNCTION by_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: by_nlff

    by_nlff = SQRT(1 - l * l / (kx * kx)) * SIN(kx * x) * SINH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION by_nlff

  FUNCTION bz_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: bz_nlff

    bz_nlff = COS(kx * x) * SINH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION bz_nlff

  FUNCTION ay_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: ay_nlff

    ay_nlff = SIN(kx * x) * SINH(l * (z - Lz)) / kx / SINH(-l * Lz)

  END FUNCTION ay_nlff

  FUNCTION az_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: az_nlff

    az_nlff = SQRT(1.0_num - l * l / (kx * kx)) * &
              COS(kx * x) * SINH(l * (z - Lz)) / kx / SINH(-l * Lz)

  END FUNCTION az_nlff

  SUBROUTINE calc_initial_field

    REAL(num) :: alpha

    aax = 0

    alpha = 0.0_num
    l = SQRT(kx * kx - alpha * alpha)
    DO iz = 1, nz + 1
      DO ix = 1, nx + 1
        aay(ix, iz) = ay_nlff(xb(ix), zb(iz), kx)
      END DO
    END DO

    aaz = 0

    bbx(:, 1:nz) =-(aay(:, 2:nz+1) - aay(:, 1:nz)) / delz

    bby(1:nx, 1:nz) = (aax(:     , 2:nz+1) - aax(:   , 1:nz)) / delz  &
                    - (aaz(2:nx+1, :     ) - aaz(1:nx, :   )) / delx

    bbz(1:nx, :) = (aay(2:nx+1, :) - aay(1:nx, :)) / delx

    CALL boundary_conditions(0.0_num)

    alpha = 0.5_num * kx
    DO iz = 1, nz
      DO ix = 1, nx + 1
        bx_a(ix, iz) = bx_nlff(xb(ix), zc(iz), SQRT(kx * kx - alpha * alpha))
      END DO
    END DO
    DO iz = 1, nz
      DO ix = 1, nx
        by_a(ix, iz) = by_nlff(xc(ix), zc(iz), SQRT(kx * kx - alpha * alpha))
      END DO
    END DO
    DO iz = 1, nz + 1
      DO ix = 1, nx
        bz_a(ix, iz) = bz_nlff(xc(ix), zb(iz), SQRT(kx * kx - alpha * alpha))
      END DO
    END DO

  END SUBROUTINE calc_initial_field

  SUBROUTINE boundary_conditions(t)

    REAL(num), INTENT(IN) :: t
    REAL(num) :: alpha

    ! Boundary conditions

    ! Do z = z_min BCs first
    ! alpha = 0.5_num * kx * ramp_up(0.1_num * t)
    alpha = 0.5_num * kx
    l = SQRT(kx * kx - alpha * alpha)
    aax0 = 0.0_num
    DO ix = 1, nx + 1
      aay0(ix) = ay_nlff(xb(ix), 0.0_num       , l)
      aaz0(ix) = az_nlff(xb(ix), 0.5_num * delz, l)
    END DO
    aax1 = 0.0_num
    CALL calc_aay1
    bbx(:   , 0) =-(aay1 - aay0) / delz
    bby(1:nx, 0) = (aax1 - aax0) / delz  &
                 - (aaz0(2:nx+1) - aaz0(1:nx)) / delx
    ! DO ix = 1, nx + 1
    !  bbx(ix, 0) = bx_nlff(xb(ix), 0.5_num * delz, l)
    ! END DO
    ! DO ix = 0, nx + 1
    !  bby(ix, 0) = by_nlff(xc(ix), 0.5_num * delz, l)
    ! END DO

    bby(0, :) = bby(nx, :)
    bbz(0, :) = bbz(nx, :)

    bby(nx+1, :) = bby(1, :)
    bbz(nx+1, :) = bbz(1, :)

    bbx(:, nz+1) = bbx(:, nz)
    bby(:, nz+1) = bby(:, nz)

  END SUBROUTINE boundary_conditions

  FUNCTION ramp_up(t)

    REAL(num), INTENT(IN) :: t
    REAL(num) :: ramp_up

    IF (t .LE. 1.0_num) THEN
      ramp_up = SIN(0.5_num * pi * t) ** 2
    ELSE
      ramp_up = 1.0_num
    END IF

  END FUNCTION ramp_up

  SUBROUTINE calc_aay1

    ! Need to invert tridiagonal matrix with non-zero corner points.
    ! To do this we use:
    ! https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)

    REAL(num), DIMENSION(nx) :: a, b, u, v, y, q, p
    REAL(num), DIMENSION(nx, 3) :: ab
    REAL(num) :: b1, vTy, vTq

    a = 1.0_num / (delx * delx)
    b = -2.0_num / (delx * delx)
    b1 = b(1)
    b(1) = b(1) + b1
    b(nx) = b(nx) + a(1) * a(1) / b1

    ab(:, 1) = a
    ab(:, 2) = b
    ab(:, 3) = a

    u = 0.0_num
    v = 0.0_num
    u(1) = -b1
    u(nx) = a(1)
    v(1) = 1
    v(nx) = -a(1) / b1

    y = solve_tridiag(ab, -bbz(1:nx, 1))
    q = solve_tridiag(ab, u)

    vTy = v(1) * y(1) + v(nx) * y(nx)
    vTq = v(1) * q(1) + v(nx) * q(nx)
    p = y - vTy / (1.0_num + vTq) * q

    aay1(2:nx) = -(p(2:nx) - p(1:nx-1)) / delx
    aay1(1) = (p(1) - p(nx)) / delx
    aay1(nx+1) = (p(1) - p(nx)) / delx

  END SUBROUTINE calc_aay1

  FUNCTION solve_tridiag(ab, d)

    ! We use the algorithm descirbed in:
    ! https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Thomas.f90
    REAL(num), DIMENSION(:, :), INTENT(IN) :: ab
    REAL(num), DIMENSION(:), INTENT(IN) :: d
    REAL(num), DIMENSION(SIZE(d)) :: solve_tridiag
    REAL(num), DIMENSION(SIZE(ab, 1), SIZE(ab, 2)) :: c
    REAL(num), DIMENSION(SIZE(d)) :: b, x
    REAL(num) :: coeff
    INTEGER :: i, n

    n = SIZE(d)
    c = ab
    b = d

    ! Step 1: Forward elimination
    DO i = 2,n
      coeff = c(i, 1) / c(i-1, 2)
      c(i, 2) = c(i, 2) - coeff * c(i-1, 3)
      b(i) = b(i) - coeff * b(i-1)
    END DO

    ! Step 2: Back substitution
    x(n) = b(n) / c(n, 2)
    DO i = n-1, 1, -1
       x(i) = (b(i)- c(i, 3) * x(i+1)) / c(i, 2)
    END DO

    solve_tridiag = x

  END FUNCTION solve_tridiag

  SUBROUTINE arrayaloc

    ALLOCATE(aax(1:nx, 1:nz+1))
    ALLOCATE(aay(1:nx+1, 1:nz+1))
    ALLOCATE(aaz(1:nx+1, 1:nz))

    ALLOCATE(bbx(1:nx+1, 0:nz+1))
    ALLOCATE(bby(0:nx+1, 0:nz+1))
    ALLOCATE(bbz(0:nx+1, 1:nz+1))

    ALLOCATE(ccx(0:nx+1, 1:nz+1))
    ALLOCATE(ccy(1:nx+1, 1:nz+1))
    ALLOCATE(ccz(1:nx+1, 0:nz+1))

    ALLOCATE(bx(1:nx+1, 1:nz+1))
    ALLOCATE(by(1:nx+1, 1:nz+1))
    ALLOCATE(bz(1:nx+1, 1:nz+1))

    ALLOCATE( bb(1:nx+1, 1:nz+1))
    ALLOCATE(bbm(1:nx+1, 1:nz+1))

    ALLOCATE(bx_a(1:nx+1, 1:nz))
    ALLOCATE(by_a(1:nx, 1:nz))
    ALLOCATE(bz_a(1:nx, 1:nz+1))

    ALLOCATE(cx(1:nx+1, 1:nz+1))
    ALLOCATE(cy(1:nx+1, 1:nz+1))
    ALLOCATE(cz(1:nx+1, 1:nz+1))

    ALLOCATE(vx(1:nx+1, 1:nz+1))
    ALLOCATE(vy(1:nx+1, 1:nz+1))
    ALLOCATE(vz(1:nx+1, 1:nz+1))

    ALLOCATE(ex(1:nx+1, 1:nz+1))
    ALLOCATE(ey(1:nx+1, 1:nz+1))
    ALLOCATE(ez(1:nx+1, 1:nz+1))

    ALLOCATE(eex(1:nx  , 1:nz+1))
    ALLOCATE(eey(1:nx+1, 1:nz+1))
    ALLOCATE(eez(1:nx+1, 1:nz  ))

    ALLOCATE(aax0(1:nx))
    ALLOCATE(aay0(1:nx+1))
    ALLOCATE(aaz0(1:nx+1))
    ALLOCATE(aax1(1:nx))
    ALLOCATE(aay1(1:nx+1))

  END SUBROUTINE arrayaloc

  SUBROUTINE arraydealoc

    DEALLOCATE(aax, aay, aaz)
    DEALLOCATE(bbx, bby, bbz)
    DEALLOCATE(ccx, ccy, ccz)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(bb, bbm)
    DEALLOCATE(bx_a, by_a, bz_a)
    DEALLOCATE(cx, cy, cz)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(ex, ey, ez)
    DEALLOCATE(eex, eey, eez)

    DEALLOCATE(aax0, aay0, aaz0)
    DEALLOCATE(aax1, aay1)

    DEALLOCATE(xc, xb, zc, zb)

  END SUBROUTINE arraydealoc

END MODULE grid
