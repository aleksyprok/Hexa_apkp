PROGRAM relax

  USE var_global
  USE grid
  USE io
  USE cal

  IMPLICIT NONE

  timestep_s = 60. ! Timestep in seconds.

  CALL grid_setup
  CALL arrayaloc

  basedt = timestep_s / time_s
  etad = 0.05 * delx * delx / basedt
  dt = 0.2 * delx * delx / etad

  CALL calc_initial_field

  CALL writedata(0)
  OPEN (UNIT = 50, FILE = 'run1/diagnostic', STATUS = 'unknown')

  CALL relax_routine

  CALL arraydealoc
  CLOSE(50)

END PROGRAM relax
