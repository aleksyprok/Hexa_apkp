
MODULE cal

  USE var_global
  USE io

  IMPLICIT none

CONTAINS

  SUBROUTINE relax_routine

    REAL(num) :: frc, bbmax, t, local_max
    REAL(num) :: localmin, globalmin
    REAL(num), DIMENSION(nx, ny, nz) :: b_c, bm_c ! bb and bbm at cell centre
    REAL(num) :: b2local, b2global
    REAL(num), PARAMETER :: c = 3.e10_num ! cm/s
    REAL(num), DIMENSION(nx, ny) :: exbase, eybase, bxbase, bybase, poynting
    REAL(num) :: szlocal, szglobal
    REAL(num) :: qlocal, qglobal
    REAL(num), DIMENSION(nx+1, ny+1, nz+1) :: vv ! Velocity squared
    REAL(num), DIMENSION(nx, ny, nz) :: diss
    REAL(num) :: min_diva_local, min_diva_global, max_diva_local, max_diva_global
    REAL(num) :: mean_diva_local, mean_diva_global
    REAL(num) :: t_snapshot, dt_snapshots
    INTEGER :: i, j, k, n
    INTEGER :: snapshot_num

    ! Define frictional coefficient in dimensionless units
    frc = frc_coef * 1.E10 / (length_cm * length_cm) * time_s

    t = 0. ! Set local time variable to zero
    t_snapshot = 0.
    dt_snapshots = 1.0
    t_snapshot = t_snapshot + dt_snapshots
    snapshot_num = 0
    n = 0

    DO WHILE (t .LT. 100.0)
    ! DO WHILE (n .LE. 100)

      n = n + 1
      t = t + dt

      ccx = (bbz(:, 1:ny+1, :) - bbz(:, 0:ny, :)) / dely  &
          - (bby(:, :, 1:nz+1) - bby(:, :, 0:nz)) / delz
      ccy = (bbx(:, :, 1:nz+1) - bbx(:, :, 0:nz)) / delz  &
          - (bbz(1:nx+1, :, :) - bbz(0:nx, :, :)) / delx
      ccz = (bby(1:nx+1, :, :) - bby(0:nx, :, :)) / delx  &
          - (bbx(:, 1:ny+1, :) - bbx(:, 0:ny, :)) / dely

      bx = 0.25 * (bbx(:, 0:ny, 0:nz  ) + bbx(:, 1:ny+1, 0:nz  ) + &
                   bbx(:, 0:ny, 1:nz+1) + bbx(:, 1:ny+1, 1:nz+1))
      by = 0.25 * (bby(0:nx, :, 0:nz  ) + bby(1:nx+1, :, 0:nz  ) + &
                   bby(0:nx, :, 1:nz+1) + bby(1:nx+1, :, 1:nz+1))
      bz = 0.25 * (bbz(0:nx, 0:ny  , :) + bbz(1:nx+1, 0:ny  , :) + &
                   bbz(0:nx, 1:ny+1, :) + bbz(1:nx+1, 1:ny+1, :))

      bb = bx * bx + by * by + bz * bz

      ! Field strength with minimum value for magneto-friction
      DO k = 1, nz + 1
        local_max = 0.0
        DO j = 1, ny + 1
          DO i = 1, nx + 1
             local_max = MAX(bb(i, j, k), local_max)
          ENDDO
        ENDDO
        CALL MPI_ALLREDUCE(local_max, bbmax, 1, mpi_num, MPI_MAX, comm, ierr)
        bbmax = 1.e-4 * bbmax
        DO j = 1, ny + 1
          DO i = 1, nx + 1
             bbm(i, j, k) = MAX(bb(i, j, k), bbmax)
          ENDDO
        ENDDO
      ENDDO

      cx = 0.5 * (ccx(0:nx, :, :) + ccx(1:nx+1, :, :))
      cy = 0.5 * (ccy(:, 0:ny, :) + ccy(:, 1:ny+1, :))
      cz = 0.5 * (ccz(:, :, 0:nz) + ccz(:, :, 1:nz+1))

      vx = frc * (cy * bz - cz * by) / bbm
      vy = frc * (cz * bx - cx * bz) / bbm
      vz = frc * (cx * by - cy * bx) / bbm

      ex = vz * by - vy * bz
      ey = vx * bz - vz * bx
      ez = vy * bx - vx * by

      eex = 0.5 * (ex(1:nx, :, :) + ex(2:nx+1, :, :))
      eey = 0.5 * (ey(:, 1:ny, :) + ey(:, 2:ny+1, :))
      eez = 0.5 * (ez(:, :, 1:nz) + ez(:, :, 2:nz+1))

      ! Overwrite electric Field
      eex(:, :, 1) = 0.0_num
      eey(:, :, 1) = 0.0_num
      ! IF (t .GE. 11.0_num) THEN
      !   eex(:, :, 1) = 0.0_num
      !   eey(:, :, 1) = 0.0_num
      ! END IF

      aax = aax - dt * eex
      aay = aay - dt * eey
      aaz = aaz - dt * eez

      CALL diva_diffusion

      bbx(:, 1:ny, 1:nz) = (aaz(:, 2:ny+1, :) - aaz(:, 1:ny, :)) / dely &
                         - (aay(:, :, 2:nz+1) - aay(:, :, 1:nz)) / delz
      bby(1:nx, :, 1:nz) = (aax(:, :, 2:nz+1) - aax(:, :, 1:nz)) / delz &
                         - (aaz(2:nx+1, :, :) - aaz(1:nx, :, :)) / delx
      bbz(1:nx, 1:ny, :) = (aay(2:nx+1, :, :) - aay(1:nx, :, :)) / delx &
                         - (aax(:, 2:ny+1, :) - aax(:, 1:ny, :)) / dely
      CALL horizontal_transfer
      CALL vertical_transfer
      CALL boundary_conditions(t)

      !####################################################
      !                  Output values
      !####################################################

      ! Calculate magnetic energy
      b_c = 0.125_num * (bb(1:nx  , 1:ny  , 2:nz  ) + bb(2:nx+1, 1:ny  , 2:nz  ) + &
                         bb(1:nx  , 2:ny+1, 2:nz  ) + bb(2:nx+1, 2:ny+1, 2:nz  ) + &
                         bb(1:nx  , 1:ny  , 3:nz+1) + bb(2:nx+1, 1:ny  , 3:nz+1) + &
                         bb(1:nx  , 2:ny+1, 3:nz+1) + bb(2:nx+1, 2:ny+1, 3:nz+1))
      bm_c = 0.125_num * (bbm(1:nx  , 1:ny  , 2:nz  ) + bbm(2:nx+1, 1:ny  , 2:nz  ) + &
                          bbm(1:nx  , 2:ny+1, 2:nz  ) + bbm(2:nx+1, 2:ny+1, 2:nz  ) + &
                          bbm(1:nx  , 1:ny  , 3:nz+1) + bbm(2:nx+1, 1:ny  , 3:nz+1) + &
                          bbm(1:nx  , 2:ny+1, 3:nz+1) + bbm(2:nx+1, 2:ny+1, 3:nz+1))
      b2local = 0.5_num * SUM(b_c) * delx * dely * delz
      CALL MPI_ALLREDUCE(b2local, b2global, 1, mpi_num, MPI_SUM, MPI_COMM_WORLD, ierr)

      !poynting flux (ExBy-EyBx)

      bxbase = 0.25_num * (bbx(1:nx, 1:ny, 1) + bbx(2:nx+1, 1:ny, 1) + &
                           bbx(1:nx, 1:ny, 2) + bbx(2:nx+1, 1:ny, 2))
      bybase = 0.25_num * (bby(1:nx, 1:ny, 1) + bby(1:nx, 2:ny+1, 1) + &
                           bby(1:nx, 1:ny, 2) + bby(1:nx, 2:ny+1, 2))

      exbase = 0.5_num * (eex(1:nx, :, 2) + eex(2:nx+1, :, 2))
      eybase = 0.5_num * (eey(:, 1:ny, 2) + eey(:, 2:ny+1, 2))

      poynting = (exbase * bybase - eybase * bxbase) * delx * dely
      szlocal = SUM(poynting)

      CALL MPI_ALLREDUCE(szlocal, szglobal, 1, mpi_num, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Dissipation Q = B ^ 2 * nu * v ^ 2
      vv = vx * vx + vy * vy + vz * vz
      diss = 0.125_num * (vv(1:nx  , 1:ny  , 2:nz  ) + vv(2:nx+1, 1:ny  , 2:nz  ) + &
                          vv(1:nx  , 2:ny+1, 2:nz  ) + vv(2:nx+1, 2:ny+1, 2:nz  ) + &
                          vv(1:nx  , 1:ny  , 3:nz+1) + vv(2:nx+1, 1:ny  , 3:nz+1) + &
                          vv(1:nx  , 2:ny+1, 3:nz+1) + vv(2:nx+1, 2:ny+1, 3:nz+1))
      diss = diss / frc
      diss = diss * bm_c
      qlocal = SUM(diss) * delx * dely * delz
      CALL MPI_ALLREDUCE(qlocal, qglobal, 1, mpi_num, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! div A
      diva = ABS(diva)
      min_diva_local = MINVAL(diva)
      max_diva_local = MAXVAL(diva)
      mean_diva_local = SUM(diva) / (nxglobal - 1) / (nyglobal - 1) / (nzglobal - 1)
      CALL MPI_ALLREDUCE(min_diva_local, min_diva_global, 1, mpi_num, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(max_diva_local, max_diva_global, 1, mpi_num, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(mean_diva_local, mean_diva_global, 1, mpi_num, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (rank .EQ. rankstart) THEN
        IF (MOD(n, 100) .EQ. 0) PRINT*, n, t
        WRITE(50, *) t, &
                     dt, &
                     b2global, &
                     szglobal, &
                     qglobal, &
                     min_diva_global, &
                     max_diva_global, &
                     mean_diva_global
      ENDIF

      IF (t .GE. t_snapshot) THEN
        t_snapshot = t_snapshot + dt_snapshots
        snapshot_num = snapshot_num + 1
        CALL writedata(snapshot_num)
      END IF
      ! snapshot_num = snapshot_num + 1
      ! CALL writedata(snapshot_num)

      !####################################################
      !                 Determine timestep
      !####################################################

      ! Determine minimum cell crossing time for advection terms

      localmin = MINVAL( (/ delx / MAXVAL(ABS(vx)), &
                            dely / MAXVAL(ABS(vy)), &
                            delz / MAXVAL(ABS(vz)), &
                            delx * delx / etad /) )
      CALL MPI_ALLREDUCE(localmin, globalmin, 1, mpi_num, MPI_MIN, MPI_COMM_WORLD, ierr)

      ! Set timestep
      dt = 0.1_num * globalmin

    END DO

  END SUBROUTINE relax_routine

  SUBROUTINE diva_diffusion

    REAL(num), DIMENSION(0:nx+1, 1:ny+1, 1:nz+1) :: aax_dummy
    REAL(num), DIMENSION(1:nx+1, 0:ny+1, 1:nz+1) :: aay_dummy
    REAL(num), DIMENSION(1:nx+1, 1:ny+1, 0:nz+1) :: aaz_dummy

    aax_dummy = 0.0_num
    aay_dummy = 0.0_num
    aaz_dummy = 0.0_num

    aax_dummy(1:nx, :, :) = aax
    aay_dummy(:, 1:ny, :) = aay
    aaz_dummy(:, :, 1:nz) = aaz

    ! Horizontal transfer
    CALL MPI_SENDRECV(aax_dummy(nx, :, :), (ny + 1) * (nz + 1), mpi_num, right, tag, &
                      aax_dummy(0,  :, :), (ny + 1) * (nz + 1), mpi_num, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(aax_dummy(1   , :, :), (ny + 1) * (nz + 1), mpi_num, left,  tag, &
                      aax_dummy(nx+1, :, :), (ny + 1) * (nz + 1), mpi_num, right, tag, &
                      comm, stat, ierr)
    ! Vertical transfer
    CALL MPI_SENDRECV(aay_dummy(:, ny, :), (nx + 1) * (nz + 1), mpi_num, up,   tag, &
                      aay_dummy(:, 0,  :), (nx + 1) * (nz + 1), mpi_num, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(aay_dummy(:, 1,    :), (nx + 1) * (nz + 1), mpi_num, down, tag, &
                      aay_dummy(:, ny+1, :), (nx + 1) * (nz + 1), mpi_num, up,   tag, &
                      comm, stat, ierr)

    diva = (aax_dummy(1:nx+1, :, :) - aax_dummy(0:nx, :, :)) / delx &
         + (aay_dummy(:, 1:ny+1, :) - aay_dummy(:, 0:ny, :)) / dely &
         + (aaz_dummy(:, :, 1:nz+1) - aaz_dummy(:, :, 0:nz)) / delz

    ! Apply boundary conditions
    IF (left .EQ. MPI_PROC_NULL) diva(1, :, :) = diva(2, :, :)
    IF (right .EQ. MPI_PROC_NULL) diva(nx+1, :, :) = diva(nx, :, :)
    IF (down .EQ. MPI_PROC_NULL) diva(:, 1, :) = diva(:, 2, :)
    IF (up .EQ. MPI_PROC_NULL) diva(:, ny+1, :) = diva(:, ny, :)
    diva(:, :, 1) = diva(:, :, 2)
    diva(:, :, nz+1) = diva(:, :, nz)

    aax = aax + etad * dt * (diva(2:nx+1, :, :) - diva(1:nx, :, :)) / delx
    aay = aay + etad * dt * (diva(:, 2:ny+1, :) - diva(:, 1:ny, :)) / dely
    aaz = aaz + etad * dt * (diva(:, :, 2:nz+1) - diva(:, :, 1:nz)) / delz

  END SUBROUTINE diva_diffusion

  SUBROUTINE calc_initial_field

    CALL readdata(potential_field_file)

    bbx(:, 1:ny, 1:nz) =  &
        (aaz(:, 2:ny+1, :     ) - aaz(:, 1:ny, :   )) / dely  &
      - (aay(:, :     , 2:nz+1) - aay(:, :   , 1:nz)) / delz

    bby(1:nx, :, 1:nz) =  &
        (aax(:     , :, 2:nz+1) - aax(:   , :, 1:nz)) / delz  &
      - (aaz(2:nx+1, :, :     ) - aaz(1:nx, :, :   )) / delx

    bbz(1:nx, 1:ny, :) =  &
        (aay(2:nx+1, :     , :) - aay(1:nx, :   , :)) / delx  &
      - (aax(:     , 2:ny+1, :) - aax(:   , 1:ny, :)) / dely

    CALL horizontal_transfer
    CALL vertical_transfer

    CALL boundary_conditions(0.0_num)

  END SUBROUTINE calc_initial_field

  SUBROUTINE boundary_conditions(t)

    REAL(num), INTENT(IN) :: t
    REAL(num), DIMENSION(1:nx+1, 0:ny+1) :: bbx_p
    REAL(num), DIMENSION(0:nx+1, 1:ny+1) :: bby_p
    INTEGER :: ix, iy

    ! Bounary conditions

    IF (left .EQ. MPI_PROC_NULL) THEN
      bby(0, :, :) = bby(1, :, :)
      bbz(0, :, :) = bbz(1, :, :)
    END IF

    IF (right .EQ. MPI_PROC_NULL) THEN
      bby(nx+1, :, :) = bby(nx, :, :)
      bbz(nx+1, :, :) = bbz(nx, :, :)
    END IF

    IF (down .EQ. MPI_PROC_NULL) THEN
      bbx(:, 0, :) = bbx(:, 1, :)
      bbz(:, 0, :) = bbz(:, 1, :)
    END IF

    IF (up .EQ. MPI_PROC_NULL) THEN
      bbx(:, ny+1, :) = bbx(:, ny, :)
      bbz(:, ny+1, :) = bbz(:, ny, :)
    END IF

    bbx_p = bbx(:, :, 1) - delz &
          * (bbz(1:nx+1, :, 1) - bbz(0:nx, :, 1)) / delx
    bby_p = bby(:, :, 1) - delz &
          * (bbz(:, 1:ny+1, 1) - bbz(:, 0:ny, 1)) / dely
    bbx(:, :, 0) = bbx_p + (2. * bbx_fix - bbx(:, :, 1) - bbx_p) * ramp_up(t)
    bby(:, :, 0) = bby_p + (2. * bby_fix - bby(:, :, 1) - bby_p) * ramp_up(t)

    bbx(:, :, nz+1) = bbx(:, :, nz)
    bby(:, :, nz+1) = bby(:, :, nz)

    ! IF (left .EQ. MPI_PROC_NULL) THEN
    !   bby(0, :, :) = bby(1, :, :)
    !   bbz(0, :, :) = bbz(1, :, :)
    ! END IF
    !
    ! IF (right .EQ. MPI_PROC_NULL) THEN
    !   bby(nx+1, :, :) = bby(nx, :, :)
    !   bbz(nx+1, :, :) = bbz(nx, :, :)
    ! END IF
    !
    ! IF (down .EQ. MPI_PROC_NULL) THEN
    !   bbx(:, 0, :) = bbx(:, 1, :)
    !   bbz(:, 0, :) = bbz(:, 1, :)
    ! END IF
    !
    ! IF (up .EQ. MPI_PROC_NULL) THEN
    !   bbx(:, ny+1, :) = bbx(:, ny, :)
    !   bbz(:, ny+1, :) = bbz(:, ny, :)
    ! END IF
    !
    ! bbx(:, :, nz+1) = bbx(:, :, nz)
    ! bby(:, :, nz+1) = bby(:, :, nz)
    !
    ! IF (t .LE. 10.0_num) THEN
    !   bbx(:, :, 0) = 0.0_num
    !   DO iy = 1, ny + 1
    !     DO ix = 0, nx + 1
    !       bby(ix, iy, 0) = driver(xc(ix), yb(iy), t)
    !     END DO
    !   END DO
    ! END IF
    !
    ! bbx(:, :, nz+1) = bbx(:, :, nz) + delz &
    !                 * (bbz(1:nx+1, :, nz+1) - bbz(0:nx, :, nz+1)) / delx
    !
    ! bby(:, :, nz+1) = bby(:, :, nz) + delz &
    !                 * (bbz(:, 1:ny+1, nz+1) - bbz(:, 0:ny, nz+1)) / dely

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

  FUNCTION ramp_down(t)

    REAL(num), INTENT(IN) :: t
    REAL(num) :: ramp_down

    IF (t .LE. 10.0_num) THEN
      ramp_down = 1.0_num
    ELSE IF ((t .GT. 10.0_num) .AND. (t .LE. 11.0_num)) THEN
      ramp_down = 1.0_num - SIN(0.5_num * pi * t) ** 2
    ELSE
      ramp_down = 0.0_num
    END IF

  END FUNCTION ramp_down

  FUNCTION driver(x, y, t)

    REAL(num), INTENT(in) :: x, y, t
    REAL(num) :: driver
    REAL(num) :: r

    r = SQRT((x - 3.0_num) ** 2.0_num + (y - 3.0_num) ** 2.0_num)
    IF (r .LE. 3.0_num) THEN
      driver = COS(pi * r / 6.0_num) ** 2.0_num * SIN(0.1_num * pi * t)
    ELSE
      driver = 0.0_num
    END IF
    ! driver = SIN(pi * t) * ramp_up(t)

  END FUNCTION driver

  SUBROUTINE calc_boundary_field

    CALL readdata(evolution_field_file)

    ! Calculate bbx1

    bbx1(:, 1:ny) =  &
        (aaz(:, 2:ny+1, 1) - aaz(:, 1:ny, 1)) / dely  &
      - (aay(:, :     , 2) - aay(:, :   , 1)) / delz

    ! Vertical transfer
    CALL MPI_SENDRECV(bbx1(:, ny), nx + 1, mpi_num, up,   tag, &
                      bbx1(:, 0 ), nx + 1, mpi_num, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx1(:, 1   ), nx + 1, mpi_num, down, tag, &
                      bbx1(:, ny+1), nx + 1, mpi_num, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (down .EQ. MPI_PROC_NULL) bbx1(:, 0   ) = bbx1(:, 1 )
    IF (up   .EQ. MPI_PROC_NULL) bbx1(:, ny+1) = bbx1(:, ny)

    ! Calculate bby1

    bby1(1:nx, :) =  &
        (aax(:,      :, 2) - aax(:,    :, 1)) / delz  &
      - (aaz(2:nx+1, :, 1) - aaz(1:nx, :, 1)) / delx
    ! Horizontal transfer
    CALL MPI_SENDRECV(bby1(nx, :), ny + 1, mpi_num, right, tag, &
                      bby1(0,  :), ny + 1, mpi_num, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby1(1   , :), ny + 1, mpi_num, left,  tag, &
                      bby1(nx+1, :), ny + 1, mpi_num, right, tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bby1(0,    :) = bby1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bby1(nx+1, :) = bby1(nx, :)

    ! Calculate bbz1

    bbz1(1:nx, 1:ny) =  &
        (aay(2:nx+1,      :, 1) - aay(1:nx,    :, 1)) / delx  &
      - (aax(     :, 2:ny+1, 1) - aax(   :, 1:ny, 1)) / dely
    ! Horizontal transfer
    CALL MPI_SENDRECV(bbz1(nx, :), ny + 2, mpi_num, right, tag, &
                      bbz1(0,  :), ny + 2, mpi_num, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(1   , :), ny + 2, mpi_num, left,  tag, &
                      bbz1(nx+1, :), ny + 2, mpi_num, right, tag, &
                      comm, stat, ierr)
    ! Vertical transfer
    CALL MPI_SENDRECV(bbz1(:, ny), nx + 2, mpi_num, up,   tag, &
                      bbz1(:, 0 ), nx + 2, mpi_num, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(:, 1   ), nx + 2, mpi_num, down, tag, &
                      bbz1(:, ny+1), nx + 2, mpi_num, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bbz1(0,    :) = bbz1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bbz1(nx+1, :) = bbz1(nx, :)
    IF (down  .EQ. MPI_PROC_NULL) bbz1(:, 0   ) = bbz1(:, 1 )
    IF (up    .EQ. MPI_PROC_NULL) bbz1(:, ny+1) = bbz1(:, ny)

    ! Calculate bbx0, bby0
    bbx0 = bbx1 - delz / delx * (bbz1(1:nx+1, :) - bbz1(0:nx, :))
    bby0 = bby1 - delz / dely * (bbz1(:, 1:ny+1) - bbz1(:, 0:ny))

    bbx_fix = 0.5 * (bbx0 + bbx1)
    bby_fix = 0.5 * (bby0 + bby1)

  END SUBROUTINE calc_boundary_field

  SUBROUTINE horizontal_transfer

    CALL MPI_SENDRECV(bby(nx, :, :), (ny + 1) * (nz + 2), mpi_num, right, tag, &
                      bby(0,  :, :), (ny + 1) * (nz + 2), mpi_num, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby(1   , :, :), (ny + 1) * (nz + 2), mpi_num, left,  tag, &
                      bby(nx+1, :, :), (ny + 1) * (nz + 2), mpi_num, right, tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(nx, :, :), (ny + 2) * (nz + 1), mpi_num, right, tag, &
                      bbz(0,  :, :), (ny + 2) * (nz + 1), mpi_num, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(1   , :, :), (ny + 2) * (nz + 1), mpi_num, left,  tag, &
                      bbz(nx+1, :, :), (ny + 2) * (nz + 1), mpi_num, right, tag, &
                      comm, stat, ierr)

  END SUBROUTINE horizontal_transfer

  SUBROUTINE vertical_transfer

    CALL MPI_SENDRECV(bbx(:, ny, :), (nx + 1) * (nz + 2), mpi_num, up,   tag, &
                      bbx(:, 0 , :), (nx + 1) * (nz + 2), mpi_num, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx(:, 1   , :), (nx + 1) * (nz + 2), mpi_num, down, tag, &
                      bbx(:, ny+1, :), (nx + 1) * (nz + 2), mpi_num, up,   tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(:, ny, :), (nx + 2) * (nz + 1), mpi_num, up,   tag, &
                      bbz(:, 0 , :), (nx + 2) * (nz + 1), mpi_num, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(:, 1   , :), (nx + 2) * (nz + 1), mpi_num, down, tag, &
                      bbz(:, ny+1, :), (nx + 2) * (nz + 1), mpi_num, up,   tag, &
                      comm, stat, ierr)

  END SUBROUTINE vertical_transfer


END MODULE cal
