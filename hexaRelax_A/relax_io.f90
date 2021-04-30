MODULE io

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE readdata(filename)

    CHARACTER (LEN = *), INTENT(IN) :: filename
    REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: aax_global, aay_global, aaz_global
    REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: aax_dummy, aay_dummy, aaz_dummy
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j, k, opt

    ALLOCATE(aax_dummy(nx,     ny + 1, nz + 1))
    ALLOCATE(aay_dummy(nx + 1, ny,     nz + 1))
    ALLOCATE(aaz_dummy(nx + 1, ny + 1, nz    ))

    ! Get aa_global and calculate aa for processor at (0, 0) cartesian coordinate
    IF (rank .EQ. rankstart) THEN

      ALLOCATE(aax_global(nxglobal,     nyglobal + 1, nzglobal + 1))
      ALLOCATE(aay_global(nxglobal + 1, nyglobal,     nzglobal + 1))
      ALLOCATE(aaz_global(nxglobal + 1, nyglobal + 1, nzglobal    ))

      PRINT*, 'Reading 3D model from ' // filename

      OPEN(UNIT = 42, &
           FILE = filename, &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD')

      READ(42) opt

      IF (.NOT. (opt .EQ. 1)) STOP 'Invalid option'

      READ(42) (((aax_global(i, j, k), i = 1, nxglobal    ), j = 1, nyglobal + 1), k = 1, nzglobal + 1)
      READ(42) (((aay_global(i, j, k), i = 1, nxglobal + 1), j = 1, nyglobal    ), k = 1, nzglobal + 1)
      READ(42) (((aaz_global(i, j, k), i = 1, nxglobal + 1), j = 1, nyglobal + 1), k = 1, nzglobal    )

      CLOSE(42)

      aax_dummy = aax_global(1:nx  , 1:ny+1, 1:nz+1)
      aay_dummy = aay_global(1:nx+1, 1:ny,   1:nz+1)
      aaz_dummy = aaz_global(1:nx+1, 1:ny+1, 1:nz  )

    END IF


    ! Get aa for the remaining processors
    IF(rank .EQ. rankstart) THEN

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated aa for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_SEND(aax_global((i * nx) + 1 : (i + 1) * nx, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   1 : nz + 1), &
                        nx * (ny + 1) * (nz + 1), &
                        MPI_REAL4, nextrank, tag, comm, ierr)
          CALL MPI_SEND(aay_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny, &
                                   1 : nz + 1), &
                        (nx + 1) * ny * (nz + 1), &
                        MPI_REAL4, nextrank, tag, comm, ierr)
          CALL MPI_SEND(aaz_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   1 : nz), &
                        (nx + 1) * (ny + 1) * nz, &
                        MPI_REAL4, nextrank, tag, comm, ierr)
        END DO
      END DO

    ELSE

      CALL MPI_RECV(aax_dummy, nx * (ny + 1) * (nz + 1), MPI_REAL4, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aay_dummy, (nx + 1) * ny * (nz + 1), MPI_REAL4, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aaz_dummy, (nx + 1) * (ny + 1) * nz, MPI_REAL4, rankstart, tag, comm, stat, ierr)

    END IF

    aax = aax_dummy
    aay = aay_dummy
    aaz = aaz_dummy

    CALL MPI_BARRIER(comm, ierr)
    IF (rank .EQ. rankstart) DEALLOCATE(aax_global, aay_global, aaz_global)
    DEALLOCATE(aax_dummy, aay_dummy, aaz_dummy)

  END SUBROUTINE readdata

  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: aax_global, aay_global, aaz_global
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j

    IF (rank .EQ. rankstart) THEN

       ALLOCATE(aax_global(nxglobal, nyglobal+1, nzglobal+1))
       ALLOCATE(aay_global(nxglobal+1, nyglobal, nzglobal+1))
       ALLOCATE(aaz_global(nxglobal+1, nyglobal+1, nzglobal))

       aax_global(1:nx, 1:ny+1, 1:nz+1) = aax
       aay_global(1:nx+1, 1:ny, 1:nz+1) = aay
       aaz_global(1:nx+1, 1:ny+1, 1:nz) = aaz

    END IF

    IF (rank .NE. rankstart) THEN

      CALL MPI_SEND(aax, nx * (ny + 1) * (nz + 1), mpi_num, rankstart, tag, comm, ierr)
      CALL MPI_SEND(aay, (nx + 1) * ny * (nz + 1), mpi_num, rankstart, tag, comm, ierr)
      CALL MPI_SEND(aaz, (nx + 1) * (ny + 1) * nz, mpi_num, rankstart, tag, comm, ierr)

    ELSE

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated aa for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_RECV(aax_global((i * nx) + 1 : (i + 1) * nx, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        nx * (ny + 1) * (nz + 1), &
                        mpi_num, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(aay_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny, &
                                   :), &
                        (nx + 1) * ny * (nz + 1), &
                        mpi_num, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(aaz_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 1) * (ny + 1) * nz, &
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
      WRITE (42) aax_global
      WRITE (42) aay_global
      WRITE (42) aaz_global
      CLOSE (42)

      DEALLOCATE(aax_global, aay_global, aaz_global)

    ENDIF

  END SUBROUTINE writedata

END MODULE io
