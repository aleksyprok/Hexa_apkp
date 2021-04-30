PROGRAM test

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: aax
  INTEGER, PARAMETER :: nx=2, ny=2, nz=3
  INTEGER :: ix, iy, iz

  ALLOCATE(aax(1:nx,0:ny,0:nz))

  DO iz = 0, nz
    DO iy = 0, ny
      DO ix = 1, nx
        aax(ix,iy,iz) = ix + iy + iz
      END DO
    END DO
  END DO

  DO iz = 0, nz
    DO iy = 0, ny
      PRINT*, aax(:,iy,iz)
    END DO
    PRINT*, ' '
    PRINT*, ' '
  END DO

  OPEN(UNIT = 1, &
       FILE = 'test.dat', &
       FORM = 'unformatted', &
       STATUS = 'unknown')
  WRITE(1) (((aax(ix,iy,iz), ix = 1, nx), &
                             iy = 0, ny), &
                             iz = 0, nz)
  CLOSE(1)

  DEALLOCATE(aax)

END PROGRAM
