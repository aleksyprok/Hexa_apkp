PROGRAM relax

  USE var_global
  USE grid
  USE io
  USE cal

  IMPLICIT NONE

  timestep_s = 60. ! Timestep in seconds.

  CALL MPI_INIT(ierr) ! Initiate MPI
  IF (num .EQ. 4) mpi_num = MPI_REAL4
  IF (num .EQ. 8) mpi_num = MPI_REAL8
  IF (num .EQ. 16) mpi_num = MPI_REAL16

  CALL MPI_COMM_SIZE(comm, mpisize, ierr) ! Get # of processes

  CALL setup_param ! Check if periodic or open boundary conditions are used
  CALL grid_setup  ! Sets up computational grid and MPI (and reads in param1)
  CALL arrayaloc ! Allocates arrays

  basedt = timestep_s / time_s
  etad = 0.05 * delx * delx / basedt
  dt = 0.2 * delx * delx / etad

  CALL calc_boundary_field
  CALL calc_initial_field

  CALL writedata(0)
  IF (rank .EQ. rankstart) THEN
    OPEN (UNIT = 50, FILE = 'run1/diagnostic', STATUS = 'unknown')
  ENDIF

  CALL relax_routine

  CALL MPI_BARRIER(comm, ierr)
  IF (rank .EQ. rankstart) THEN
    CLOSE(50)
  ENDIF
  CALL MPI_BARRIER(comm, ierr)
  CALL arraydealoc
  CALL MPI_FINALIZE(ierr)

END PROGRAM relax
