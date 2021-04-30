MODULE var_global

  IMPLICIT NONE

  INCLUDE "mpif.h"

  INTEGER, PARAMETER :: num = KIND(0.0D0)
  ! INTEGER, PARAMETER :: num = 16
  REAL(num), PARAMETER :: pi = 4.0_num * ATAN(1.0_num)
  REAL(num), PARAMETER :: frc_coef = 300.0_num ! Frictional coefficient (km^2/s)

  ! Grid parameters:
  REAL(num) :: delx, dely, delz

  ! Number of cells (not corners):
  INTEGER :: nxglobal, nyglobal, nzglobal
  INTEGER :: nx, ny, nz

  INTEGER :: periodic, open

  ! File and directory names
  CHARACTER (LEN = *), PARAMETER :: potential_field_file = 'setup_files/smothed_128x128/run1_00020p'
  CHARACTER (LEN = *), PARAMETER :: evolution_field_file = 'setup_files/smothed_128x128/run1_00020p'
  CHARACTER (LEN = *), PARAMETER :: parameters_file      = 'setup_files/smothed_128x128/param1'       ! Used to get nx, ny, nz
  CHARACTER (LEN = *), PARAMETER :: setup_file           = 'setup_files/smothed_128x128/run1_setup'   ! Used to check for periodic/open boundaries
  CHARACTER (LEN = *), PARAMETER :: output_file          = 'run1/relax_'
  CHARACTER (LEN = 50) :: filename

  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: aax, aay, aaz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bbx, bby, bbz, bx, by, bz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: ccx, ccy, ccz, cx, cy, cz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eex, eey, eez, ex, ey, ez
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: vx, vy, vz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bb, bbm

  REAL(num), DIMENSION(:, :), ALLOCATABLE :: bbx0, bby0
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: bbx1, bby1, bbz1
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: bbx_fix, bby_fix

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb, yc, yb, zc, zb

  ! REAL(num), PARAMETER :: kx = pi / 3.0_num, ky = pi / 3.0_num
  ! REAL(num), PARAMETER :: kz = SQRT(2.0_num) * pi / 3.0_num
  REAL(num), PARAMETER :: kx = pi / 3.0_num
  REAL(num) :: l = 0.5_num * kx

  !  Model parameters:
  REAL(num) :: etad

  ! Variables associated with MPI
  INTEGER :: mpisize, ierr, rank, comm
  INTEGER :: left, right, up, down, nextrank, rankstart, rankend
  INTEGER, PARAMETER :: mpidir = 2, tag = 100
  INTEGER, DIMENSION(mpidir) :: nproc, coords
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
  LOGICAL, DIMENSION(mpidir) :: periods
  INTEGER :: mpi_num

  ! Variables associated with length and time
  REAL(num) :: length_cm  ! Length of one Hexa length unit in cm
  REAL(num) :: time_s     ! Time of one hexa time unit in seconds
  REAL(num) :: dt, basedt ! Timestep and default timestep (hexa units)
  REAL(num) :: timestep_s ! Default timestep in seconds

END MODULE var_global
