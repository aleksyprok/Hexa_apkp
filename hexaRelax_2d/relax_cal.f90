
MODULE cal

  USE var_global
  USE io
  USE grid

  IMPLICIT none

CONTAINS

  SUBROUTINE relax_routine

    REAL(num) :: bbmax
    REAL(num) :: t_eneg_shot, dt_eneg_shot, t_snapshot, dt_snapshots
    INTEGER :: n, snapshot_num

    ! Define frictional coefficient in dimensionless units
    frc = frc_coef * 1.E10 / (length_cm * length_cm) * time_s

    t = 0.0_num ! Set local time variable to zero
    t_snapshot = 0.0_num
    dt_snapshots = 10.0_num
    dt_eneg_shot = 1.0_num
    t_snapshot = t_snapshot + dt_snapshots
    snapshot_num = 0
    n = 0

    DO WHILE (t .LT. 1.E4_num)

      n = n + 1
      t = t + dt

      ccx =-(bby(:, 1:nz+1) - bby(:, 0:nz)) / delz
      ccy = (bbx(:, 1:nz+1) - bbx(:, 0:nz)) / delz  &
          - (bbz(1:nx+1, :) - bbz(0:nx, :)) / delx
      ccz = (bby(1:nx+1, :) - bby(0:nx, :)) / delx

      bx = 0.5_num * (bbx(:, 0:nz) + bbx(:, 1:nz+1))
      by = 0.25_num * (bby(0:nx, 0:nz  ) + bby(1:nx+1, 0:nz  ) + &
                       bby(0:nx, 1:nz+1) + bby(1:nx+1, 1:nz+1))
      bz = 0.5_num * (bbz(0:nx, :) + bbz(1:nx+1, :))

      bb = bx * bx + by * by + bz * bz

      bbmax = 1.e-8_num * MAXVAL(bb)
      DO iz = 1, nz + 1
        DO ix = 1, nx + 1
           bbm(ix, iz) = MAX(bb(ix, iz), bbmax)
        ENDDO
      ENDDO

      cx = 0.5 * (ccx(0:nx, :) + ccx(1:nx+1, :))
      cy = ccy
      cz = 0.5 * (ccz(:, 0:nz) + ccz(:, 1:nz+1))

      vx = frc * (cy * bz - cz * by) / bbm
      vy = frc * (cz * bx - cx * bz) / bbm
      vz = frc * (cx * by - cy * bx) / bbm

      ex = vz * by - vy * bz
      ey = vx * bz - vz * bx
      ez = vy * bx - vx * by

      eex = 0.5 * (ex(1:nx, :) + ex(2:nx+1, :))
      eey = ey
      eez = 0.5 * (ez(:, 1:nz) + ez(:, 2:nz+1))

      bbx(:, 1:nz) = bbx(:, 1:nz) + dt * (eey(:, 2:nz+1) - eey(:, 1:nz)) / delz
      bby(1:nx, 1:nz) = bby(1:nx, 1:nz) + dt * ((eez(2:nx+1, :) - eez(1:nx, :)) / delx - &
                                                (eex(:, 2:nz+1) - eex(:, 1:nz)) / delz )
      bbz(1:nx, :) = bbz(1:nx, :) - dt * (eey(2:nx+1, :) - eey(1:nx, :)) / delx

      CALL boundary_conditions(t)

      IF (MOD(n, 100) .EQ. 0) PRINT*, n, t
      IF (t .GE. t_eneg_shot) THEN
        t_eneg_shot = t_eneg_shot + dt_eneg_shot
        CALL output_diag
      END IF
      IF (t .GE. t_snapshot) THEN
        t_snapshot = t_snapshot + dt_snapshots
        snapshot_num = snapshot_num + 1
        CALL writedata(snapshot_num)
      END IF

      ! Set timestep
      dt = 0.1_num *  MINVAL( (/ delx / MAXVAL( (/ ABS(vx), TINY(0.0_num) /) ), &
                                 delz / MAXVAL( (/ ABS(vz), TINY(0.0_num) /) ), &
                                 delx * delz / etad /) )

    END DO

  END SUBROUTINE relax_routine

  SUBROUTINE boundary_conditions(t)

    REAL(num), INTENT(IN) :: t

    ! Bounary conditions

    bby(0, :) = bby(1, :)
    bbz(0, :) = bbz(1, :)

    bby(nx+1, :) = bby(nx, :)
    bbz(nx+1, :) = bbz(nx, :)

    l = kx * (1.0_num - 0.5_num * ramp_up(0.1_num * t))
    DO ix = 1, nx + 1
      bbx(ix, 0   ) = bx_nlff(xb(ix), 0.0_num, l)
      bbx(ix, nz+1) = bx_nlff(xb(ix), 6.0_num, l)
    END DO
    DO ix = 0, nx + 1
      bby(ix, 0   ) = by_nlff(xc(ix), 0.0_num, l)
      bby(ix, nz+1) = by_nlff(xc(ix), 6.0_num, l)
    END DO

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

END MODULE cal
