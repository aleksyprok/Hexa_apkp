MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

    INTEGER :: namelen, num_hex_cells
    CHARACTER (LEN = MPI_MAX_PROCESSOR_NAME) :: procname
    INTEGER :: ix, iy, iz

    CALL MPI_DIMS_CREATE(mpisize, mpidir, nproc, ierr) !nproc = dims

    IF (periodic .EQ. 1) THEN
      periods(1) = .TRUE.
      periods(2) = .TRUE.
    ELSE
      periods(1) = .FALSE.
      periods(2) = .FALSE.
    END IF

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpidir, nproc, periods, .TRUE., comm, ierr)

    CALL MPI_COMM_RANK(comm, rank, ierr)

    CALL MPI_CART_COORDS(comm, rank, mpidir, coords, ierr)

    CALL MPI_GET_PROCESSOR_NAME(procname, namelen, ierr)

    CALL MPI_BARRIER(comm, ierr)

    CALL MPI_CART_SHIFT(comm, 0, 1, left, right, ierr)
    CALL MPI_CART_SHIFT(comm, 1, 1, down, up, ierr)

    ! Get rank of process at (0,0)
    CALL MPI_CART_RANK(comm, (/0, 0/), rankstart, ierr)

    IF (rank .EQ. rankstart) THEN
      OPEN(UNIT = 42, &
           FILE = parameters_file, &
           FORM = 'FORMATTED', &
           STATUS = 'OLD')
      READ(42, *) num_hex_cells
      READ(42, *) nxglobal, nyglobal, nzglobal
      READ(42, *) length_cm
      READ(42, *) time_s
      CLOSE(42)
    END IF
    CALL MPI_BCAST(num_hex_cells, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nxglobal, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nyglobal , 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nzglobal, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(length_cm, 1, mpi_num, rankstart, comm, ierr)
    CALL MPI_BCAST(time_s, 1, mpi_num, rankstart, comm, ierr)
                                                      !
    CALL MPI_BARRIER(comm,ierr)

    nz = nzglobal            ! no mpi in z direction
    nx = nxglobal / nproc(1) ! no of cells in x associated with each process
    ny = nyglobal / nproc(2) ! no of cells in y associated with each process

    IF ( ((nx * nproc(1)) .NE. nxglobal) .OR.  &
         ((ny * nproc(2)) .NE. nyglobal) ) THEN
        IF ( (nx * nproc(1)) .NE. nxglobal) THEN
           PRINT*,'Unable to subdivide equally in x. Fix grid'
        ENDIF
        IF ( (ny * nproc(2)) .NE. nyglobal) THEN
           PRINT*,'Unable to subdivide equally in y. Fix grid'
        ENDIF
        CALL MPI_FINALIZE(ierr)
        STOP
    ENDIF

    ! Use width of 2
    delx = 2.0_num / nxglobal
    dely = 2.0_num / nyglobal
    delz = 2.0_num / nzglobal

    ALLOCATE(xc(0:nx+1))
    ALLOCATE(yc(0:ny+1))
    ALLOCATE(zc(0:nz+1))
    ALLOCATE(xb(nx+1))
    ALLOCATE(yb(ny+1))
    ALLOCATE(zb(nz+1))

    ! xb ranges from -1 to 1
    DO ix = 0, nx + 1
      xc(ix) = (coords(1) * nx + ix) * delx - 0.5_num * delx - 1.0_num
    END DO
    DO ix = 1, nx + 1
      xb(ix) = (coords(1) * nx + ix - 1) * delx - 1.0_num
    END DO

    ! yb ranges from -1 to 1
    DO iy = 0, ny + 1
      yc(iy) = (coords(2) * ny + iy) * dely - 0.5_num * dely - 1.0_num
    END DO
    DO iy = 1, ny + 1
      yb(iy) = (coords(2) * ny + iy - 1) * dely - 1.0_num
    END DO

    DO iz = 0, nz + 1
      zc(iz) = iz * delz - 0.5_num * delz
    END DO
    DO iz = 1, nz + 1
      zb(iz) = (iz - 1) * delz
    END DO

  END SUBROUTINE grid_setup

  SUBROUTINE setup_param

    ! This subroutine reads the setup file to check if periodic / open
    ! boundary conditions are used.

    ! Dummy strings to read files
    CHARACTER (LEN = 10) :: s1
    CHARACTER (LEN = 6 ) :: s2
    INTEGER :: i

    OPEN(UNIT   = 42, &
         FILE   = setup_file, &
         FORM   = 'FORMATTED', &
         STATUS = 'OLD', &
         ACTION = 'READ')

    ! Skip first 6 lines of file
    DO i = 1, 6
      READ(42, *)
    END DO

    ! Read periodic and open values
    READ(42, *) s1
    READ(42, *) s2

    ! Get integer values from the strings
    READ(s1(10:10), *) periodic
    READ(s2(6:6)  , *) open

    CLOSE(42)

  END SUBROUTINE setup_param

  SUBROUTINE arrayaloc

    ALLOCATE(aax(nx,   ny+1, nz+1))
    ALLOCATE(aay(nx+1, ny,   nz+1))
    ALLOCATE(aaz(nx+1, ny+1, nz  ))

    ALLOCATE(bbx(1:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(bby(0:nx+1, 1:ny+1, 0:nz+1))
    ALLOCATE(bbz(0:nx+1, 0:ny+1, 1:nz+1))

    ALLOCATE(bbx_read(1:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(bby_read(0:nx+1, 1:ny+1, 0:nz+1))
    ALLOCATE(bbz_read(0:nx+1, 0:ny+1, 1:nz+1))

    ALLOCATE(bbx_po(1:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(bby_po(0:nx+1, 1:ny+1, 0:nz+1))
    ALLOCATE(bbz_po(0:nx+1, 0:ny+1, 1:nz+1))

    ALLOCATE(bbx_ff(1:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(bby_ff(0:nx+1, 1:ny+1, 0:nz+1))
    ALLOCATE(bbz_ff(0:nx+1, 0:ny+1, 1:nz+1))

    ALLOCATE(ccx(0:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(ccy(1:nx+1, 0:ny+1, 1:nz+1))
    ALLOCATE(ccz(1:nx+1, 1:ny+1, 0:nz+1))

    ALLOCATE(bx(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(by(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(bz(1:nx+1, 1:ny+1, 1:nz+1))

    ALLOCATE( bb(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(bbm(1:nx+1, 1:ny+1, 1:nz+1))

    ALLOCATE(cx(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(cy(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(cz(1:nx+1, 1:ny+1, 1:nz+1))

    ALLOCATE(vx(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(vy(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(vz(1:nx+1, 1:ny+1, 1:nz+1))

    ALLOCATE(ex(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(ey(1:nx+1, 1:ny+1, 1:nz+1))
    ALLOCATE(ez(1:nx+1, 1:ny+1, 1:nz+1))

    ALLOCATE(eex(1:nx  , 1:ny+1, 1:nz+1))
    ALLOCATE(eey(1:nx+1, 1:ny  , 1:nz+1))
    ALLOCATE(eez(1:nx+1, 1:ny+1, 1:nz  ))

  END SUBROUTINE arrayaloc

  SUBROUTINE arraydealoc

    DEALLOCATE(aax, aay, aaz)
    DEALLOCATE(bbx, bby, bbz)
    DEALLOCATE(ccx, ccy, ccz)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(bb, bbm)
    DEALLOCATE(cx, cy, cz)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(ex, ey, ez)
    DEALLOCATE(eex, eey, eez)

    DEALLOCATE(xc, xb, yc, yb, zc, zb)
    DEALLOCATE(bbx_read, bby_read, bbz_read)
    DEALLOCATE(bbx_po, bby_po, bbz_po)
    DEALLOCATE(bbx_ff, bby_ff, bbz_ff)

  END SUBROUTINE arraydealoc

END MODULE grid
