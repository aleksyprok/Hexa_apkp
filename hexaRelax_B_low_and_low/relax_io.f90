MODULE io

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE readdata(filename)

    CHARACTER (LEN = *), INTENT(IN) :: filename
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bbx_global, bby_global, bbz_global
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bbx_dummy, bby_dummy, bbz_dummy
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j

    ALLOCATE(bbx_dummy(1:nx+1, 0:ny+1, 0:nz+1))
    ALLOCATE(bby_dummy(0:nx+1, 1:ny+1, 0:nz+1))
    ALLOCATE(bbz_dummy(0:nx+1, 0:ny+1, 1:nz+1))

    ! Get bb_global and calculate bb for processor at (0, 0) cartesian coordinate
    IF (rank .EQ. rankstart) THEN

      ALLOCATE(bbx_global(1:nxglobal+1, 0:nyglobal+1, 0:nzglobal+1))
      ALLOCATE(bby_global(0:nxglobal+1, 1:nyglobal+1, 0:nzglobal+1))
      ALLOCATE(bbz_global(0:nxglobal+1, 0:nyglobal+1, 1:nzglobal+1))

      PRINT*, 'Reading 3D model from ' // filename

      OPEN(UNIT = 42, &
           FILE = filename, &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD')

      READ(42) bbx_global
      READ(42) bby_global
      READ(42) bbz_global

      CLOSE(42)

      bbx_dummy = bbx_global(1:nx+1, 0:ny+1, 0:nz+1)
      bby_dummy = bby_global(0:nx+1, 1:ny+1, 0:nz+1)
      bbz_dummy = bbz_global(0:nx+1, 0:ny+1, 1:nz+1)

    END IF

    ! Get bb for the remaining processors
    IF(rank .EQ. rankstart) THEN

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated bb for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_SEND(bbx_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny)     : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 1) * (ny + 2) * (nz + 2), &
                        mpi_num, nextrank, tag, comm, ierr)
          CALL MPI_SEND(bby_global((i * nx)     : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 1) * (nz + 2), &
                        mpi_num, nextrank, tag, comm, ierr)
          CALL MPI_SEND(bbz_global((i * nx)     : (i + 1) * nx + 1, &
                                   (j * ny)     : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 2) * (nz + 1), &
                        mpi_num, nextrank, tag, comm, ierr)
        END DO
      END DO

    ELSE

      CALL MPI_RECV(bbx_dummy, (nx + 1) * (ny + 2) * (nz + 2), mpi_num, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(bby_dummy, (nx + 2) * (ny + 1) * (nz + 2), mpi_num, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(bbz_dummy, (nx + 2) * (ny + 2) * (nz + 1), mpi_num, rankstart, tag, comm, stat, ierr)

    END IF

    bbx_read = bbx_dummy
    bby_read = bby_dummy
    bbz_read = bbz_dummy

    CALL MPI_BARRIER(comm, ierr)
    IF (rank .EQ. rankstart) DEALLOCATE(bbx_global, bby_global, bbz_global)
    DEALLOCATE(bbx_dummy, bby_dummy, bbz_dummy)

  END SUBROUTINE readdata
  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bbx_global, bby_global, bbz_global
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j

    IF (rank .EQ. rankstart) THEN

       ALLOCATE(bbx_global(1:nxglobal+1, 0:nyglobal+1, 0:nzglobal+1))
       ALLOCATE(bby_global(0:nxglobal+1, 1:nyglobal+1, 0:nzglobal+1))
       ALLOCATE(bbz_global(0:nxglobal+1, 0:nyglobal+1, 1:nzglobal+1))

       bbx_global(1:nx+1, 0:ny+1, 0:nz+1) = bbx
       bby_global(0:nx+1, 1:ny+1, 0:nz+1) = bby
       bbz_global(0:nx+1, 0:ny+1, 1:nz+1) = bbz

    END IF

    IF (rank .NE. rankstart) THEN

      CALL MPI_SEND(bbx, (nx + 1) * (ny + 2) * (nz + 2), mpi_num, rankstart, tag, comm, ierr)
      CALL MPI_SEND(bby, (nx + 2) * (ny + 1) * (nz + 2), mpi_num, rankstart, tag, comm, ierr)
      CALL MPI_SEND(bbz, (nx + 2) * (ny + 2) * (nz + 1), mpi_num, rankstart, tag, comm, ierr)

    ELSE

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated bb for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_RECV(bbx_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny)     : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 1) * (ny + 2) * (nz + 2), &
                        mpi_num, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(bby_global((i * nx)     : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 1) * (nz + 2), &
                        mpi_num, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(bbz_global((i * nx) : (i + 1) * nx + 1, &
                                   (j * ny) : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 2) * (nz + 1), &
                        mpi_num, nextrank, tag, comm, stat, ierr)
        END DO
      END DO

    END IF

    IF (rank .EQ. rankstart) THEN

      WRITE (filename, FMT = '(a, i5.5)')  output_file, n

      OPEN(UNIT = 42, &
           FILE = FILENAME, &
           FORM = 'UNFORMATTED', &
           STATUS = 'UNKNOWN')
      WRITE (42) bbx_global
      WRITE (42) bby_global
      WRITE (42) bbz_global
      CLOSE (42)

      DEALLOCATE(bbx_global, bby_global, bbz_global)

    ENDIF

  END SUBROUTINE writedata

END MODULE io
