MODULE io

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n

    WRITE (filename, FMT = '(a, i5.5)')  output_file, n

    OPEN(UNIT = 42, &
         FILE = FILENAME, &
         FORM = 'UNFORMATTED', &
         STATUS = 'UNKNOWN')
    WRITE (42) t
    WRITE (42) bbx
    WRITE (42) bby
    WRITE (42) bbz
    CLOSE (42)

  END SUBROUTINE writedata

  SUBROUTINE output_diag

    REAL(num), DIMENSION(nx, 2:nz) :: b_c, bm_c ! bb and bbm at cell centre
    REAL(num), DIMENSION(nx) :: exbase, eybase, bxbase, bybase, poynting
    REAL(num), DIMENSION(nx+1, nz+1) :: vv ! Velocity squared
    REAL(num), DIMENSION(nx, 2:nz) :: diss
    REAL(num), DIMENSION(nx, 2:nz) :: divb
    REAL(num) :: mag_eng, sz, q, min_divb, max_divb, mean_divb

    ! Calculate magnetic energy
    b_c = 0.25_num * (bb(1:nx  , 2:nz  ) + bb(2:nx+1, 2:nz  ) + &
                      bb(1:nx  , 3:nz+1) + bb(2:nx+1, 3:nz+1))
    bm_c = 0.25_num * (bbm(1:nx  , 2:nz  ) + bbm(2:nx+1, 2:nz  ) + &
                       bbm(1:nx  , 3:nz+1) + bbm(2:nx+1, 3:nz+1))
    mag_eng = 0.5_num * SUM(b_c) * delx * delz

    !poynting flux (ExBy-EyBx)
    bxbase = 0.25_num * (bbx(1:nx, 2) + bbx(2:nx+1, 2) + &
                         bbx(1:nx, 3) + bbx(2:nx+1, 3))
    bybase = 0.25_num * (bby(1:nx, 2) + bby(1:nx, 2) + &
                         bby(1:nx, 3) + bby(1:nx, 3))
    exbase = eex(:, 2) + eex(:, 2)
    eybase = 0.5_num * (eey(1:nx, 2) + eey(2:nx+1, 2))
    poynting = (exbase * bybase - eybase * bxbase) * delx
    sz = SUM(poynting)

    ! Dissipation Q = B ^ 2 * nu * v ^ 2
    vv = vx * vx + vy * vy + vz * vz
    diss = 0.25_num * (vv(1:nx  , 2:nz  ) + vv(2:nx+1, 2:nz  ) + &
                       vv(1:nx  , 3:nz+1) + vv(2:nx+1, 3:nz+1))
    diss = diss / frc
    diss = diss * bm_c
    q = SUM(diss) * delx * delz

    ! div B
    divb = (bbx(2:nx+1, 2:nz) - bbx(1:nx, 2:nz)) / delx &
         + (bbz(1:nx, 3:nz+1) - bbz(1:nx, 2:nz)) / delz
    divb = ABS(divb)
    min_divb = MINVAL(divb)
    max_divb = MAXVAL(divb)
    mean_divb = SUM(divb) / REAL(nx * nz, num)

    WRITE(50, *) t, &
                 dt, &
                 mag_eng, &
                 sz, &
                 q, &
                 min_divb, &
                 max_divb, &
                 mean_divb

  END SUBROUTINE output_diag

END MODULE io
