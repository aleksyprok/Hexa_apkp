MODULE io

  USE var_global
  USE grid

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

  SUBROUTINE write_hexa(n)

    INTEGER, INTENT(IN) :: n
    INTEGER :: opt = 1
    REAL(4), DIMENSION(nx, 2, nz+1) :: aax_dummy
    REAL(4), DIMENSION(nx+1, 1, nz+1) :: aay_dummy
    REAL(4), DIMENSION(nx+1, 2, nz) :: aaz_dummy

    aax_dummy(:, 1, :) = REAL(aax, 4)
    aax_dummy(:, 2, :) = REAL(aax, 4)
    aay_dummy(:, 1, :) = REAL(aay, 4)
    aaz_dummy(:, 1, :) = REAL(aaz, 4)
    aaz_dummy(:, 2, :) = REAL(aaz, 4)

    WRITE (filename, FMT = '(a, i5.5, "p")')  output_hexa, n

    OPEN(UNIT = 42, &
         FILE = FILENAME, &
         FORM = 'UNFORMATTED', &
         STATUS = 'UNKNOWN')
    WRITE (42) opt
    ! WRITE (42) aax_dummy
    ! WRITE (42) aay_dummy
    ! WRITE (42) aaz_dummy
    WRITE (42) REAL(aax, 4)
    WRITE (42) REAL(aay, 4)
    WRITE (42) REAL(aaz, 4)
    CLOSE (42)

  END SUBROUTINE write_hexa

  SUBROUTINE output_diag

    REAL(num), DIMENSION(nx, 2:nz) :: b_c, bm_c ! bb and bbm at cell centre
    REAL(num), DIMENSION(nx) :: exbase, eybase, bxbase, bybase, poynting
    REAL(num), DIMENSION(nx+1, nz+1) :: vv ! Velocity squared
    REAL(num), DIMENSION(nx, 2:nz) :: diss
    REAL(num), DIMENSION(nx+1, nz) :: bx_err
    REAL(num), DIMENSION(nx, nz) :: by_err
    REAL(num), DIMENSION(nx, nz+1) :: bz_err
    REAL(num), DIMENSION(nx, 2:nz) :: err_array
    REAL(num) :: mag_eng, sz, q, min_diva, max_diva, mean_diva, err

    ! Calculate magnetic energy
    b_c = 0.25_num * (bb(1:nx  , 2:nz  ) + bb(2:nx+1, 2:nz  ) + &
                      bb(1:nx  , 3:nz+1) + bb(2:nx+1, 3:nz+1))
    bm_c = 0.25_num * (bbm(1:nx  , 2:nz  ) + bbm(2:nx+1, 2:nz  ) + &
                       bbm(1:nx  , 3:nz+1) + bbm(2:nx+1, 3:nz+1))
    mag_eng = 0.5_num * SUM(b_c) * delx * delz

    !poynting flux (ExBy-EyBx)
    bxbase = 0.25_num * (bbx(1:nx, 1) + bbx(2:nx+1, 1) + &
                         bbx(1:nx, 2) + bbx(2:nx+1, 2))
    bybase = 0.25_num * (bby(1:nx, 1) + bby(1:nx, 1) + &
                         bby(1:nx, 2) + bby(1:nx, 2))
    exbase = eex(:, 2)
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

    ! div A
    diva = (aax(2:nx, 2:nz) - aax(1:nx-1, 2:nz)) / delx &
         + (aaz(2:nx, 2:nz) - aaz(2:nx, 1:nz-1)) / delz
    diva = ABS(diva)
    min_diva = MINVAL(diva)
    max_diva = MAXVAL(diva)
    mean_diva = SUM(diva) / REAL(nx * nz, num)

    ! Analytic solutions vs. numerical solution

    bx_err = (bbx(:, 1:nz) - bx_a) ** 2.0_num
    by_err = (bby(1:nx, 1:nz) - by_a) ** 2.0_num
    bz_err = (bbz(1:nx, :) - bz_a) ** 2.0_num
    ! print*, MAXVAL(bx_err), MAXVAL(by_err), MAXVAL(bz_err)
    ! print*, ' -------------------------------------'
    ! print*, bx_a(1:5, 1)
    ! print*, ' -------------------------------------'
    ! print*, bbx(1:5, 1)



    err_array = 0.5_num * (bx_err(1:nx, 2:nz) + bx_err(2:nx+1, 2:nz)) &
              + by_err(:, 2:nz) &
              + 0.5_num * (bz_err(:, 2:nz) + bz_err(:, 3:nz+1))
    err = SUM(err_array) * delx * delz

    WRITE(50, *) t, &
                 dt, &
                 mag_eng, &
                 sz, &
                 q, &
                 min_diva, &
                 max_diva, &
                 mean_diva, &
                 err

  END SUBROUTINE output_diag

END MODULE io
