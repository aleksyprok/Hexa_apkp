
MODULE cal

  USE var_global
  USE io

  IMPLICIT none

CONTAINS

  SUBROUTINE relax_routine

    REAL :: frc, bbmax, t, local_max
    REAL :: localmin, globalmin
    REAL :: b2local, b2global
    REAL :: qlocal, qglobal
    REAL, DIMENSION(nx+1, ny+1, nz+1) :: diss
    REAL :: min_divb_local, min_divb_global, max_divb_local, max_divb_global
    REAL :: mean_divb_local, mean_divb_global
    REAL, DIMENSION(nx, ny, nz) :: divb
    REAL :: t_snapshot, dt_snapshots
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
        CALL MPI_ALLREDUCE(local_max, bbmax, 1, MPI_REAL8, MPI_MAX, comm, ierr)
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
      eex(:, :, 1) = 0.
      eey(:, :, 1) = 0.

      bbx(:, 1:ny, 1:nz) = bbx(:, 1:ny, 1:nz) + dt * ( &
                           (eey(:, :, 2:nz+1) - eey(:, :, 1:nz)) / delz - &
                           (eez(:, 2:ny+1, :) - eez(:, 1:ny, :)) / dely )
      bby(1:nx, :, 1:nz) = bby(1:nx, :, 1:nz) + dt * ( &
                           (eez(2:nx+1, :, :) - eez(1:nx, :, :)) / delx - &
                           (eex(:, :, 2:nz+1) - eex(:, :, 1:nz)) / delz )
      bbz(1:nx, 1:ny, :) = bbz(1:nx, 1:ny, :) + dt * ( &
                           (eex(:, 2:ny+1, :) - eex(:, 1:ny, :)) / dely - &
                           (eey(2:nx+1, :, :) - eey(1:nx, :, :)) / delx )

      CALL horizontal_transfer
      CALL vertical_transfer
      CALL boundary_conditions(t)

      CALL divb_diffusion

      !####################################################
      !                 Determine timestep
      !####################################################

      ! Determine minimum cell crossing time for advection terms

      localmin = MINVAL( (/ delx / MAXVAL(ABS(vx)), &
                            dely / MAXVAL(ABS(vy)), &
                            delz / MAXVAL(ABS(vz)), &
                            delx * delx / etad /) )
      CALL MPI_ALLREDUCE(localmin, globalmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      ! Set timestep
      dt = 0.2 * globalmin
      ! IF (dt .LT. basedt * 1.e-3) dt = 1.e-3 * basedt ! If timestep is too small, set to 0.001*basedt

      !####################################################
      !                  Output values
      !####################################################

      ! Calculate magnetic energy
      b2local = SUM(bb(1:nx, 1:ny, 1:nz+1)) / 8. / pi * delx ** 3
      CALL MPI_ALLREDUCE(b2local, b2global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Dissipation Q = B ^ 2 / (4 * pi) * nu * v ^ 2
      diss = vx * vx + vy * vy + vz * vz
      diss = diss / frc / 4. / pi
      diss = diss * bb
      qlocal = SUM(diss(1:nx,1:ny,1:nz+1)) * delx ** 3
      CALL MPI_ALLREDUCE(qlocal, qglobal, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! div B
      divb = (bbx(2:nx+1, 1:ny, 1:nz) - bbx(1:nx, 1:ny, 1:nz)) / delx &
           + (bby(1:nx, 2:ny+1, 1:nz) - bby(1:nx, 1:ny, 1:nz)) / dely &
           + (bbz(1:nx, 1:ny, 2:nz+1) - bbz(1:nx, 1:ny, 1:nz)) / delz
      divb = ABS(divb)
      min_divb_local = MINVAL(divb)
      max_divb_local = MAXVAL(divb)
      mean_divb_local = SUM(divb) / (nxglobal * nyglobal * nzglobal)
      CALL MPI_ALLREDUCE(min_divb_local, min_divb_global, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(max_divb_local, max_divb_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(mean_divb_local, mean_divb_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (rank .EQ. rankstart) THEN
        IF (MOD(n, 100) .EQ. 0) PRINT*, n, t
        WRITE(50, *) t, &
                     b2global * length_cm ** 3., &        ! Magnetic energy (erg)
                     qglobal / time_s * length_cm ** 3, & ! Frictional dissipation (erg/s)
                     min_divb_global, &
                     max_divb_global, &
                     mean_divb_global
      ENDIF

      ! IF (MOD(n+1, 1000) .EQ. 0) THEN
      !   CALL writedata(n)
      !   CALL write_electric(n-1)
      ! END IF
      IF (t .GE. t_snapshot) THEN
        t_snapshot = t_snapshot + dt_snapshots
        snapshot_num = snapshot_num + 1
        CALL writedata(snapshot_num)
      END IF

    END DO

  END SUBROUTINE relax_routine

  SUBROUTINE divb_diffusion

    REAL, DIMENSION(0:nx+1, 0:ny+1, 0:nz+1) :: divb

    divb(1:nx, 1:ny, 1:nz) = &
        (bbx(2:nx+1, 1:ny, 1:nz) - bbx(1:nx, 1:ny, 1:nz)) / delx &
      + (bby(1:nx, 2:ny+1, 1:nz) - bby(1:nx, 1:ny, 1:nz)) / dely &
      + (bbz(1:nx, 1:ny, 2:nz+1) - bbz(1:nx, 1:ny, 1:nz)) / delz

    ! Horizontal transfer
    CALL MPI_SENDRECV(divb(nx, :, :), (ny + 2) * (nz + 2), MPI_REAL8, right, tag, &
                      divb(0,  :, :), (ny + 2) * (nz + 2), MPI_REAL8, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(divb(1   , :, :), (ny + 2) * (nz + 2), MPI_REAL8, left,  tag, &
                      divb(nx+1, :, :), (ny + 2) * (nz + 2), MPI_REAL8, right, tag, &
                      comm, stat, ierr)
    ! Vertical transfer
    CALL MPI_SENDRECV(divb(:, ny, :), (nx + 2) * (nz + 2), MPI_REAL8, up,   tag, &
                      divb(:, 0,  :), (nx + 2) * (nz + 2), MPI_REAL8, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(divb(:, 1,    :), (nx + 2) * (nz + 2), MPI_REAL8, down, tag, &
                      divb(:, ny+1, :), (nx + 2) * (nz + 2), MPI_REAL8, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left .EQ. MPI_PROC_NULL) divb(0, :, :) = divb(1, :, :)
    IF (right .EQ. MPI_PROC_NULL) divb(nx+1, :, :) = divb(nx, :, :)
    IF (down .EQ. MPI_PROC_NULL) divb(:, 0, :) = divb(:, 1, :)
    IF (up .EQ. MPI_PROC_NULL) divb(:, ny+1, :) = divb(:, ny, :)
    divb(:, :, 0) = divb(:, :, 1)
    divb(:, :, nz+1) = divb(:, :, nz)

    bbx = bbx + etad * dt * (divb(1:nx+1, :, :) - divb(0:nx, :, :)) / delx
    bby = bby + etad * dt * (divb(:, 1:ny+1, :) - divb(:, 0:ny, :)) / dely
    bbz = bbz + etad * dt * (divb(:, :, 1:nz+1) - divb(:, :, 0:nz)) / delz

  END SUBROUTINE divb_diffusion

  SUBROUTINE boundary_conditions(t)

    REAL, INTENT(IN) :: t
    REAL, DIMENSION(1:nx+1, 0:ny+1) :: bbx_p
    REAL, DIMENSION(0:nx+1, 1:ny+1) :: bby_p

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

  END SUBROUTINE boundary_conditions

  FUNCTION ramp_up(t)

    REAL, INTENT(IN) :: t
    REAL :: ramp_up

    IF (t .LE. 1.) THEN
      ramp_up = SIN(0.5 * pi * t) ** 2
    ELSE
      ramp_up = 1.
    END IF

  END FUNCTION ramp_up

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

    CALL boundary_conditions(0.)

  END SUBROUTINE calc_initial_field

  SUBROUTINE calc_boundary_field

    CALL readdata(evolution_field_file)

    ! Calculate bbx1

    bbx1(:, 1:ny) =  &
        (aaz(:, 2:ny+1, 1) - aaz(:, 1:ny, 1)) / dely  &
      - (aay(:, :     , 2) - aay(:, :   , 1)) / delz

    ! Vertical transfer
    CALL MPI_SENDRECV(bbx1(:, ny), nx + 1, MPI_REAL8, up,   tag, &
                      bbx1(:, 0 ), nx + 1, MPI_REAL8, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx1(:, 1   ), nx + 1, MPI_REAL8, down, tag, &
                      bbx1(:, ny+1), nx + 1, MPI_REAL8, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (down .EQ. MPI_PROC_NULL) bbx1(:, 0   ) = bbx1(:, 1 )
    IF (up   .EQ. MPI_PROC_NULL) bbx1(:, ny+1) = bbx1(:, ny)

    ! Calculate bby1

    bby1(1:nx, :) =  &
        (aax(:,      :, 2) - aax(:,    :, 1)) / delz  &
      - (aaz(2:nx+1, :, 1) - aaz(1:nx, :, 1)) / delx
    ! Horizontal transfer
    CALL MPI_SENDRECV(bby1(nx, :), ny + 1, MPI_REAL8, right, tag, &
                      bby1(0,  :), ny + 1, MPI_REAL8, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby1(1   , :), ny + 1, MPI_REAL8, left,  tag, &
                      bby1(nx+1, :), ny + 1, MPI_REAL8, right, tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bby1(0,    :) = bby1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bby1(nx+1, :) = bby1(nx, :)

    ! Calculate bbz1

    bbz1(1:nx, 1:ny) =  &
        (aay(2:nx+1,      :, 1) - aay(1:nx,    :, 1)) / delx  &
      - (aax(     :, 2:ny+1, 1) - aax(   :, 1:ny, 1)) / dely
    ! Horizontal transfer
    CALL MPI_SENDRECV(bbz1(nx, :), ny + 2, MPI_REAL8, right, tag, &
                      bbz1(0,  :), ny + 2, MPI_REAL8, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(1   , :), ny + 2, MPI_REAL8, left,  tag, &
                      bbz1(nx+1, :), ny + 2, MPI_REAL8, right, tag, &
                      comm, stat, ierr)
    ! Vertical transfer
    CALL MPI_SENDRECV(bbz1(:, ny), nx + 2, MPI_REAL8, up,   tag, &
                      bbz1(:, 0 ), nx + 2, MPI_REAL8, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(:, 1   ), nx + 2, MPI_REAL8, down, tag, &
                      bbz1(:, ny+1), nx + 2, MPI_REAL8, up,   tag, &
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

    CALL MPI_SENDRECV(bby(nx, :, :), (ny + 1) * (nz + 2), MPI_REAL8, right, tag, &
                      bby(0,  :, :), (ny + 1) * (nz + 2), MPI_REAL8, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby(1   , :, :), (ny + 1) * (nz + 2), MPI_REAL8, left,  tag, &
                      bby(nx+1, :, :), (ny + 1) * (nz + 2), MPI_REAL8, right, tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(nx, :, :), (ny + 2) * (nz + 1), MPI_REAL8, right, tag, &
                      bbz(0,  :, :), (ny + 2) * (nz + 1), MPI_REAL8, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(1   , :, :), (ny + 2) * (nz + 1), MPI_REAL8, left,  tag, &
                      bbz(nx+1, :, :), (ny + 2) * (nz + 1), MPI_REAL8, right, tag, &
                      comm, stat, ierr)

  END SUBROUTINE horizontal_transfer

  SUBROUTINE vertical_transfer

    CALL MPI_SENDRECV(bbx(:, ny, :), (nx + 1) * (nz + 2), MPI_REAL8, up,   tag, &
                      bbx(:, 0 , :), (nx + 1) * (nz + 2), MPI_REAL8, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx(:, 1   , :), (nx + 1) * (nz + 2), MPI_REAL8, down, tag, &
                      bbx(:, ny+1, :), (nx + 1) * (nz + 2), MPI_REAL8, up,   tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(:, ny, :), (nx + 2) * (nz + 1), MPI_REAL8, up,   tag, &
                      bbz(:, 0 , :), (nx + 2) * (nz + 1), MPI_REAL8, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(:, 1   , :), (nx + 2) * (nz + 1), MPI_REAL8, down, tag, &
                      bbz(:, ny+1, :), (nx + 2) * (nz + 1), MPI_REAL8, up,   tag, &
                      comm, stat, ierr)

  END SUBROUTINE vertical_transfer

END MODULE cal
