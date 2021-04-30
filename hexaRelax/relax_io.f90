MODULE io

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE readdata(filename)

    CHARACTER (LEN = *), INTENT(IN) :: filename
    REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: aax_global, aay_global, aaz_global
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j, k, opt

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

      aax = aax_global(1:nx  , 1:ny+1, 1:nz+1)
      aay = aay_global(1:nx+1, 1:ny,   1:nz+1)
      aaz = aaz_global(1:nx+1, 1:ny+1, 1:nz  )

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

      CALL MPI_RECV(aax, nx * (ny + 1) * (nz + 1), MPI_REAL4, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aay, (nx + 1) * ny * (nz + 1), MPI_REAL4, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aaz, (nx + 1) * (ny + 1) * nz, MPI_REAL4, rankstart, tag, comm, stat, ierr)

    END IF

    IF (rank .EQ. rankstart) THEN
      DEALLOCATE(aax_global, aay_global, aaz_global)
    ENDIF

  END SUBROUTINE readdata

  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: bbx_global, bby_global, bbz_global
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

      CALL MPI_SEND(bbx, (nx + 1) * (ny + 2) * (nz + 2), MPI_REAL8, rankstart, tag, comm, ierr)
      CALL MPI_SEND(bby, (nx + 2) * (ny + 1) * (nz + 2), MPI_REAL8, rankstart, tag, comm, ierr)
      CALL MPI_SEND(bbz, (nx + 2) * (ny + 2) * (nz + 1), MPI_REAL8, rankstart, tag, comm, ierr)

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
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(bby_global((i * nx)     : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 1) * (nz + 2), &
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(bbz_global((i * nx) : (i + 1) * nx + 1, &
                                   (j * ny) : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 2) * (ny + 2) * (nz + 1), &
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
        END DO
      END DO

    END IF

    IF (rank .EQ. rankstart) THEN

      WRITE (filename, FMT = '(a, i5.5)')  output_file, n

      OPEN(UNIT = 42, &
           FILE = FILENAME, &
           FORM = 'UNFORMATTED', &
           STATUS = 'UNKNOWN')
      ! WRITE (42) (((bbx_global(i,j,k), i = 1, nxglobal + 1), j = 0, nyglobal + 1), k = 0, nzglobal + 1)
      ! WRITE (42) (((bby_global(i,j,k), i = 0, nxglobal + 1), j = 1, nyglobal + 1), k = 0, nzglobal + 1)
      ! WRITE (42) (((bbz_global(i,j,k), i = 0, nxglobal + 1), j = 0, nyglobal + 1), k = 1, nzglobal + 1)
      WRITE (42) bbx_global
      WRITE (42) bby_global
      WRITE (42) bbz_global
      CLOSE (42)

      DEALLOCATE(bbx_global, bby_global, bbz_global)

    ENDIF

  END SUBROUTINE writedata

  SUBROUTINE write_electric(n)

    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: eex_global, eey_global, eez_global
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j

    IF (rank .EQ. rankstart) THEN

       ALLOCATE(eex_global(1:nxglobal  , 1:nyglobal+1, 1:nzglobal+1))
       ALLOCATE(eey_global(1:nxglobal+1, 1:nyglobal  , 1:nzglobal+1))
       ALLOCATE(eez_global(1:nxglobal+1, 1:nyglobal+1, 1:nzglobal  ))

       eex_global(1:nx  , 1:ny+1, 1:nz+1) = eex
       eey_global(1:nx+1, 1:ny  , 1:nz+1) = eey
       eez_global(1:nx+1, 1:ny+1, 1:nz  ) = eez

    END IF

    IF(rank .NE. rankstart) THEN

      CALL MPI_SEND(eex, nx * (ny + 1) * (nz + 1), MPI_REAL8, rankstart, tag, comm, ierr)
      CALL MPI_SEND(eey, (nx + 1) * ny * (nz + 1), MPI_REAL8, rankstart, tag, comm, ierr)
      CALL MPI_SEND(eez, (nx + 1) * (ny + 1) * nz, MPI_REAL8, rankstart, tag, comm, ierr)

    ELSE

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated ee for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_RECV(eex_global((i * nx) + 1 : (i + 1) * nx, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        nx * (ny + 1) * (nz + 1), &
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(eey_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny, &
                                   :), &
                        (nx + 1) * ny * (nz + 1), &
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
          CALL MPI_RECV(eez_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   :), &
                        (nx + 1) * (ny + 1) * nz, &
                        MPI_REAL8, nextrank, tag, comm, stat, ierr)
        END DO
      END DO

    END IF

    IF (rank .EQ. rankstart) THEN

      WRITE (filename, FMT = '(a, i5.5)')  electric_file, n

      OPEN(UNIT = 42, &
           FILE = FILENAME, &
           FORM = 'UNFORMATTED', &
           STATUS = 'UNKNOWN')
      WRITE (42) eex_global
      WRITE (42) eey_global
      WRITE (42) eez_global
      CLOSE (42)

      DEALLOCATE(eex_global, eey_global, eez_global)

    ENDIF

  END SUBROUTINE write_electric

END MODULE io
