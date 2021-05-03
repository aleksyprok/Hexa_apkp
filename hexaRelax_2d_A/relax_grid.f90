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

    bx_nlff = (l / kx) * SIN(kx * x) * EXP(-l * z)

  END FUNCTION bx_nlff

  FUNCTION by_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: by_nlff

    by_nlff = SQRT(1 - l * l / (kx * kx)) * SIN(kx * x) * EXP(-l * z)

  END FUNCTION by_nlff

  FUNCTION bz_nlff(x, z, l)

    REAL(num), INTENT(IN) :: x, z, l
    REAL(num) :: bz_nlff

    bz_nlff = COS(kx * x) * EXP(-l * z)

  END FUNCTION bz_nlff

  SUBROUTINE calc_initial_field

    DO iz = 0, nz + 1
      DO ix = 1, nx + 1
        bbx(ix, iz) = bx_nlff(xb(ix), zc(iz), l)
      END DO
    END DO

    DO iz = 0, nz + 1
      DO ix = 0, nx + 1
        bby(ix, iz) = by_nlff(xc(ix), zc(iz), l)
      END DO
    END DO

    DO iz = 1, nz + 1
      DO ix = 0, nx + 1
        bbz(ix, iz) = bz_nlff(xc(ix), zb(iz), l)
      END DO
    END DO

  END SUBROUTINE calc_initial_field

  SUBROUTINE arrayaloc

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
