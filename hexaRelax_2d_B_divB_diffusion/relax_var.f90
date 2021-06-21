MODULE var_global

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(0.0D0)
  REAL(num), PARAMETER :: pi = 4.0_num * ATAN(1.0_num)
  REAL(num), PARAMETER :: frc_coef = 300.0_num ! Frictional coefficient (km^2/s)

  ! Grid parameters:
  REAL(num) :: delx, delz
  REAL(num) :: t, frc

  ! Number of cells (not corners):
  INTEGER :: nx, nz
  INTEGER :: ix, iz

  CHARACTER (LEN = *), PARAMETER :: output_file          = 'run1/relax_'
  CHARACTER (LEN = *), PARAMETER :: output_hexa          = 'run1/run1_'
  CHARACTER (LEN = 50) :: filename

  REAL(num), DIMENSION(:, :), ALLOCATABLE :: aax, aay, aaz
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: bbx, bby, bbz, bx, by, bz
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: ccx, ccy, ccz, cx, cy, cz
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: eex, eey, eez, ex, ey, ez
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: vx, vy, vz
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: bb, bbm, bx_a, by_a, bz_a

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb, yc, yb, zc, zb

  !  Model parameters:
  REAL(num) :: etad, iter_no = 0.0_num
  REAL(num), PARAMETER :: kx = pi / 3.0_num, Lz = 6.0_num
  REAL(num) :: l = kx

  ! Variables associated with length and time
  REAL(num) :: length_cm  ! Length of one Hexa length unit in cm
  REAL(num) :: time_s     ! Time of one hexa time unit in seconds
  REAL(num) :: dt, basedt ! Timestep and default timestep (hexa units)
  REAL(num) :: timestep_s ! Default timestep in seconds

END MODULE var_global
