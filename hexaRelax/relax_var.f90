MODULE var_global

  IMPLICIT NONE

  INCLUDE "mpif.h"

  REAL, PARAMETER :: pi = 3.141592654
  REAL, PARAMETER :: frc_coef = 3000 ! Frictional coefficient (km^2/s)

  ! Grid parameters:
  REAL :: delx, dely, delz

  ! Number of cells (not corners):
  INTEGER :: nxglobal, nyglobal, nzglobal
  INTEGER :: nx, ny, nz

  INTEGER :: periodic, open

  ! File and directory names
  CHARACTER (LEN = *), PARAMETER :: potential_field_file = 'setup_files/poten_00003p'
  CHARACTER (LEN = *), PARAMETER :: evolution_field_file = 'setup_files/run1_00003p'
  CHARACTER (LEN = *), PARAMETER :: parameters_file      = 'setup_files/param1'       ! Used to get nx, ny, nz
  CHARACTER (LEN = *), PARAMETER :: setup_file           = 'setup_files/run1_setup'   ! Used to check for periodic/open boundaries
  CHARACTER (LEN = *), PARAMETER :: output_file          = 'run1/relax_'
  CHARACTER (LEN = *), PARAMETER :: electric_file        = 'run1/electric_'
  CHARACTER (LEN = 50) :: filename

  REAL(4), DIMENSION(:, :, :), ALLOCATABLE :: aax, aay, aaz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: bbx, bby, bbz, bx, by, bz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: ccx, ccy, ccz, cx, cy, cz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: eex, eey, eez, ex, ey, ez
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vx, vy, vz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: bb, bbm

  REAL, DIMENSION(:, :), ALLOCATABLE :: bbx0, bby0
  REAL, DIMENSION(:, :), ALLOCATABLE :: bbx1, bby1, bbz1
  REAL, DIMENSION(:, :), ALLOCATABLE :: bbx_fix, bby_fix

  !  Model parameters:
  REAL :: etad

  ! Variables associated with MPI
  INTEGER :: mpisize, ierr, rank, comm
  INTEGER :: left, right, up, down, nextrank, rankstart, rankend
  INTEGER, PARAMETER :: mpidir = 2, tag = 100
  INTEGER, DIMENSION(mpidir) :: nproc, coords
  integer, DIMENSION(MPI_STATUS_SIZE) :: stat
  LOGICAL, DIMENSION(mpidir) :: periods

  ! Variables associated with length and time
  REAL(4) :: length_cm ! Length of one Hexa length unit in cm
  REAL(4) :: time_s ! Time of one hexa time unit in seconds
  REAL :: dt, basedt ! Timestep and default timestep (hexa units)
  REAL :: timestep_s ! Default timestep in seconds

END MODULE var_global
