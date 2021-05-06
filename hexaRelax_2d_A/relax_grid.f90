MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

    length_cm = 2.95436E9_num
    time_s = 5760.0_num

    nx = 256
    nz = 256

    delx = 6.0_num / REAL(nx, num)
    delz = 6.0_num / REAL(nz, num)

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
      zc(iz) = iz * delz - 0.5_num * delz
    END DO
    DO iz = 1, nz + 1
      zb(iz) = (iz - 1) * delz
    END DO

  END SUBROUTINE grid_setup

  FUNCTION bx_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: bx_nlff

    bx_nlff = (l / kx) * SIN(kx * x) * COSH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION bx_nlff

  FUNCTION by_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: by_nlff

    by_nlff = -SQRT(1 - l * l / (kx * kx)) * SIN(kx * x) * SINH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION by_nlff

  FUNCTION bz_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: bz_nlff

    bz_nlff = -COS(kx * x) * SINH(l * (z - Lz)) / SINH(-l * Lz)

  END FUNCTION bz_nlff

  FUNCTION ay_poten(x, z)

    REAL(num), INTENT(IN) :: x, z
    REAL(num) :: ay_poten

    ay_poten = -SIN(kx * x) * SINH(kx * (z - Lz)) / kx / SINH(-kx * Lz)

  END FUNCTION ay_poten

  SUBROUTINE calc_initial_field

    aax = 0

    DO iz = 1, nz + 1
      DO ix = 1, nx + 1
        aay(ix, iz) = ay_poten(xb(ix), zb(iz))
      END DO
    END DO

    aaz = 0

    bbx(:, 1:nz) =-(aay(:, 2:nz+1) - aay(:, 1:nz)) / delz

    bby(1:nx, 1:nz) = (aax(:     , 2:nz+1) - aax(:   , 1:nz)) / delz  &
                    - (aaz(2:nx+1, :     ) - aaz(1:nx, :   )) / delx

    bbz(1:nx, :) = (aay(2:nx+1, :) - aay(1:nx, :)) / delx

    CALL boundary_conditions(0.0_num)

  END SUBROUTINE calc_initial_field

  SUBROUTINE boundary_conditions(t)

    REAL(num), INTENT(IN) :: t
    REAL(num) :: alpha

    ! Bounary conditions

    bby(0, :) = bby(nx, :)
    bbz(0, :) = bbz(nx, :)

    bby(nx+1, :) = bby(1, :)
    bbz(nx+1, :) = bbz(1, :)

    alpha = 0.5_num * kx * ramp_up(0.1_num * t)
    l = SQRT(kx * kx - alpha * alpha)
    DO ix = 1, nx + 1
      bbx(ix, 0   ) = bx_nlff(xb(ix), 0.0_num, l)
      ! bbx(ix, nz+1) = bx_nlff(xb(ix), 6.0_num, l)
    END DO
    DO ix = 0, nx + 1
      bby(ix, 0   ) = by_nlff(xc(ix), 0.0_num, l)
      ! bby(ix, nz+1) = by_nlff(xc(ix), 6.0_num, l)
    END DO

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

  SUBROUTINE arrayaloc

    ALLOCATE(aax(1:nx, 1:nz+1))
    ALLOCATE(aay(1:nx+1, 1:nz+1))
    ALLOCATE(aaz(1:nx+1, 1:nz))
    ALLOCATE(diva(2:nx, 2:nz))

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

  END SUBROUTINE arrayaloc

  SUBROUTINE arraydealoc

    DEALLOCATE(aax, aay, aaz, diva)
    DEALLOCATE(bbx, bby, bbz)
    DEALLOCATE(ccx, ccy, ccz)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(bb, bbm)
    DEALLOCATE(cx, cy, cz)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(ex, ey, ez)
    DEALLOCATE(eex, eey, eez)

    DEALLOCATE(xc, xb, zc, zb)

  END SUBROUTINE arraydealoc

END MODULE grid
